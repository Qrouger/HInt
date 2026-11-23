"""
HPC (Slurm) entry point of HInt

Author: Quentin Rouger

The single-node pipeline of HInt.py is split into Slurm jobs chained with dependencies :

    fetch    (login node)   UniProt sequences, downloaded by HInt-slurm submit (the login node has internet access)
    prepare  (1 job)        sequences, DeepLoc, SignalP, filters, job list of the first round
    setup    (1 job)        Jean Zay : prepare, afdb, msa and features in one job on prepost (1 GPU, internet access)
    afdb     (1 job)        MSAs of the AlphaFold database, only when compute nodes are offline (run on a partition with internet)
    msa      (1 job, GPU)   MSAs (AlphaFold database + MMseqs2-GPU colabfold_search)
    features (job array)    AlphaFold feature pickle files (one CPU job with the local species database)
    predict  (job array)    3D models of one round, jobs are packed on the GPUs of each task by VRAM
    score    (1 job)        scoring of one round, job list of the next round, final results after the last round

A round is one bait of Interact_with (in the given order, the next bait is only tested against the preys validated
with the previous one, as in HInt.py) or the homo-oligomer modeling.
The state shared between jobs is saved in log_file/slurm/state.pkl.

Usage (from the directory containing HInt.txt) :
    HInt-slurm submit --partition <partition> --gpu_vram <GiB> [options]
    HInt-slurm submit --cluster jean-zay --jz_project <project> --jz_gpu a100 [options]
    HInt-slurm status
"""

import os
import sys
import json
import time
import glob
import math
import shlex
import pickle
import shutil
import logging
import argparse
import subprocess
import multiprocessing

from concurrent.futures import ThreadPoolExecutor

SLURM_DIR = "log_file/slurm"
SCRIPT_DIR = "slurm_scripts"
STATE_FILE = f"{SLURM_DIR}/state.pkl"
CONFIG_FILE = f"{SLURM_DIR}/config.json"
NEED_PKL_FILE = f"{SLURM_DIR}/need_pkl.txt"
JOB_IDS_FILE = f"{SLURM_DIR}/job_ids.json"
LOCAL_SPECIES_DB = "uniref30_2302_species.sqlite" # in Path_Database, see set_species_env
LOCAL_SPECIES_FEATURES_CPUS = 6
HINT_LOG = "./log_file/HInt.log"
STAGES = ["fetch", "prepare", "afdb", "msa", "features", "predict", "score"]

# Utils_HInt opens ./log_file/HInt.log when imported
os.makedirs(SLURM_DIR, exist_ok=True)

from .Utils_HInt import *
from .File_proteins import File_proteins
from .Scoring_HInt import Score_interaction, Resume_file, Create_figures

logger = logging.getLogger()


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def setup_logging (name) :
    """
    Log every Slurm job in log_file/HInt.log, each line tagged with the job name.
    The file is appended (several jobs write in it) and reset by HInt-slurm submit.

    Parameters :
    ----------
    name : str
    """
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(f"%(asctime)s - %(levelname)s - [{name}] %(message)s")
    file_handler = logging.FileHandler(HINT_LOG, mode="a")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)


def format_time (seconds) :
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{seconds:06.3f}"


def allocated_cpus () :
    """
    Number of CPUs allocated by Slurm (multiprocessing.cpu_count() returns all CPUs of the node).
    """
    return int(os.environ.get("SLURM_CPUS_PER_TASK", multiprocessing.cpu_count()))


def allocated_gpus () :
    """
    GPUs allocated by Slurm, as seen by CUDA.
    """
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if visible == "" or visible == "NoDevFiles" :
        return ["0"]
    return [gpu.strip() for gpu in visible.split(",") if gpu.strip() != ""]


def local_gpu_vram () :
    """
    Smallest VRAM (GiB) of the GPUs visible on the node.
    """
    pynvml.nvmlInit()
    vram = min((pynvml.nvmlDeviceGetMemoryInfo(pynvml.nvmlDeviceGetHandleByIndex(i)).total / 1024**2) * 0.001 for i in range(pynvml.nvmlDeviceGetCount()))
    pynvml.nvmlShutdown()
    return vram


def clean_previous_run () :
    """
    Remove the logs and the state of a previous run before a new run starts from the beginning
    (old status files of the predict tasks would be read as finished tasks by the score stage).
    """
    for path in glob.glob(f"{SLURM_DIR}/*.log") + glob.glob(f"{SLURM_DIR}/*.out") + [STATE_FILE, f"{STATE_FILE}.tmp", NEED_PKL_FILE, f"{SLURM_DIR}/job_ids.txt"] :
        if os.path.isfile(path) :
            os.remove(path)
    for path in glob.glob(f"{SLURM_DIR}/round_*") :
        shutil.rmtree(path)


def load_config () :
    with open(CONFIG_FILE, "r") as f :
        return json.load(f)


def load_state () :
    with open(STATE_FILE, "rb") as f :
        return pickle.load(f)


def save_state (state) :
    """
    Atomic write of the state shared between Slurm jobs.
    """
    tmp = f"{STATE_FILE}.tmp"
    with open(tmp, "wb") as f :
        pickle.dump(state, f)
    os.replace(tmp, STATE_FILE)


def round_dir (round_id) :
    path = f"{SLURM_DIR}/round_{round_id}"
    os.makedirs(path, exist_ok=True)
    return path


def define_rounds (Informations_dict, compounds) :
    """
    List of modeling rounds : one per bait (in the order of Interact_with), then homo-oligomers, then custom interactions.

    Returns :
    ----------
    rounds : list of tuple (interaction_type, bait)
    """
    rounds = []
    if Informations_dict["Interact_with"] != [''] :
        for bait in Informations_dict["Multimer_bait"] :
            rounds.append(("Compounds" if compounds else "PPI_int", bait))
    if int(Informations_dict["Homo-oligomer"]) > 1 :
        rounds.append(("homo_int", ""))
    if Informations_dict["Custom_interactions"] != "" :
        rounds.append(("custom", ""))
    return rounds


def write_round_jobs (state, round_id, config) :
    """
    Generate and save the job list of a round from the current prey list.
    """
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    interaction_type, bait = state["rounds"][round_id]
    job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, interaction_type, bait, vram=state["gpu_vram"])
    with open(f"{round_dir(round_id)}/jobs.tsv", "w") as f :
        for job, vram in job_with_vram_length :
            f.write(f"{job.strip()}\t{vram}\n")
    n_tasks = predict_task_count(len(job_with_vram_length), config)
    logger.info(f"Round {round_id} ({interaction_type} {bait}) : {len(job_with_vram_length)} jobs, spread over {n_tasks} Slurm tasks of {predict_gpus(len(job_with_vram_length), config)} GPUs")
    HInt_object.Make_save_dict()


def predict_task_count (n_jobs, config) :
    """
    Number of predict tasks of a round : one task per jobs_per_task jobs (no maximum), the jobs are spread equally.
    """
    return math.ceil(n_jobs / config["jobs_per_task"])


def predict_gpus (n_jobs, config) :
    """
    GPUs of the predict tasks of a round, from the number of jobs per task : with 100 jobs and 4 GPUs per task,
    less than 25 jobs -> 1 GPU, less than 50 -> 2, less than 75 -> 3, otherwise 4.
    """
    n_tasks = predict_task_count(n_jobs, config)
    if n_tasks == 0 :
        return 0
    jobs_per_task = math.ceil(n_jobs / n_tasks)
    jobs_per_gpu = config["jobs_per_task"] / config["gpus_per_task"]
    return min(config["gpus_per_task"], int(jobs_per_task // jobs_per_gpu) + 1)


def record_job_id (name, job_id) :
    """
    Save a submitted job (log_file/slurm/job_ids.json and job_ids.txt, used by status and to cancel the pipeline).
    """
    named_ids = dict()
    if os.path.isfile(JOB_IDS_FILE) :
        with open(JOB_IDS_FILE, "r") as f :
            named_ids = json.load(f)
    named_ids[name] = job_id
    with open(JOB_IDS_FILE, "w") as f :
        json.dump(named_ids, f)
    with open(f"{SLURM_DIR}/job_ids.txt", "a") as f :
        f.write(f"{job_id}\n")


def wait_predict_dependency (timeout=300) :
    """
    Dependency of the predict jobs of round 0 on the last preparation job (msa or features), written by HInt-slurm submit
    once all its jobs are submitted. Empty when no preparation job follows the prepare stage.
    """
    start = time.time()
    while time.time() - start < timeout :
        if os.path.isfile(JOB_IDS_FILE) :
            with open(JOB_IDS_FILE, "r") as f :
                named_ids = json.load(f)
            if "predict_dependency" in named_ids :
                return named_ids["predict_dependency"] or None
        time.sleep(5)
    logger.warning(f"No dependency found in {JOB_IDS_FILE} (scripts submitted by hand ?) : the predict jobs of round 0 are submitted without dependency")
    return None


def submit_round (round_id, config, dependency=None, dry_run=False) :
    """
    Submit the predict job array of one round, with the number of tasks and GPUs set by its job list, then its score job.
    Called by prepare (round 0), score (next round) and HInt-slurm submit (resumed run).

    Parameters :
    ----------
    round_id : int
    config : dict
    dependency : str or None
        Dependency of the predict jobs, e.g. afterany:123.
    dry_run : bool
    """
    args = argparse.Namespace(**config["submit_args"])
    args.dry_run = dry_run
    n_jobs = len(read_round_jobs(round_id))
    n_tasks = predict_task_count(n_jobs, config)
    hint = "HInt-slurm"
    score_dependency = dependency
    if n_tasks :
        gpus = predict_gpus(n_jobs, config)
        throttle = f"%{args.max_parallel}" if args.max_parallel else ""
        script = sbatch_script(args, f"predict_r{round_id}", f"{hint} predict --round {round_id} --task ${{SLURM_ARRAY_TASK_ID}}", stage_settings(args, "predict", gpus, False), array=f"0-{n_tasks - 1}{throttle}")
        predict_id = sbatch(args, script, dependency)
        score_dependency = f"afterany:{predict_id}" # score even if some predictions failed
        logger.info(f"Round {round_id} : {n_jobs} jobs submitted in {n_tasks} predict tasks of {gpus} GPUs (job {predict_id})")
        if not dry_run :
            record_job_id(f"predict_r{round_id}", predict_id)
    else :
        logger.info(f"Round {round_id} : all models already done, no predict job")
    # internet for --send_to_database
    script = sbatch_script(args, f"score_r{round_id}", f"{hint} score --round {round_id}", stage_settings(args, "score", 0, True))
    score_id = sbatch(args, script, score_dependency)
    if not dry_run :
        record_job_id(f"score_r{round_id}", score_id)


def write_need_pkl (state, config) :
    """
    List the proteins without pkl file (read by the features stage).
    """
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]
    targets = list(HInt_object.get_possible_prey())
    for bait in Informations_dict["Interact_with"] :
        if bait != "" and bait not in targets :
            targets.append(bait)
    need_pkl = [prot for prot in targets if not os.path.isfile(f"{Path_Pickle_Feature}/{prot}.pkl")]
    without_msa = [prot for prot in need_pkl if not os.path.isfile(f"{Path_Pickle_Feature}/{prot}.a3m")]
    if without_msa :
        logger.warning(f"{len(without_msa)} proteins have no MSA after the MSA stage : {', '.join(without_msa[:20])}{' ...' if len(without_msa) > 20 else ''}")
    with open(NEED_PKL_FILE, "w") as f :
        f.write("\n".join(need_pkl) + ("\n" if need_pkl else ""))
    logger.info(f"{len(need_pkl)} proteins need pkl files, spread over {config['features_tasks']} Slurm tasks")


def missing_work (Informations_dict) :
    """
    Proteins that still need each step, from log_file/save_dict.pkl and Path_Pickle_Feature.
    The filters of the prepare stage are applied in the same order with the saved results : DeepLoc only runs on
    the baits and the preys of the right length, SignalP on those that also pass the DeepLoc filter, MSAs and pkl
    files on those that also pass the SignalP filter. A prey whose result is not known yet is kept.

    Returns :
    ----------
    missing : dict of list (keys : fetch, deeploc, signalp, msa, pkl)
    """
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"], Informations_dict["Interact_with"], Informations_dict["AlphaFold"], Informations_dict["Custom_interactions"])
    proteins = HInt_object.get_proteins()
    fasta_sequences = HInt_object.get_proteins_sequence_SP()
    saved = dict()
    if os.path.isfile("log_file/save_dict.pkl") :
        with open("log_file/save_dict.pkl", "rb") as f :
            saved = pickle.load(f)
    saved_sequences = saved.get("sequence_SP", {})
    saved_deeploc = saved.get("deeploc", {})
    saved_SP = saved.get("Signal_peptide", {})
    saved_no_SP = saved.get("sequence_no_SP", {})
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]

    baits = [prot for prot in Informations_dict["Interact_with"] if prot != ""]
    preys = [prot for prot in proteins if prot not in baits]

    # Length filter (filter_length)
    if Informations_dict["Min_protein_length"] != "0" or Informations_dict["Max_protein_length"] != "" :
        max_length = int(Informations_dict["Max_protein_length"]) if Informations_dict["Max_protein_length"] != "" else 100000
        min_length = int(Informations_dict["Min_protein_length"])
        kept = []
        for prot in preys :
            sequence = fasta_sequences.get(prot) or saved_sequences.get(prot)
            if not sequence or min_length < len(sequence) < max_length :
                kept.append(prot)
        preys = kept
    deeploc_proteins = baits + preys

    # DeepLoc filter (filter_deeploc)
    localisation = [loc.strip() for loc in Informations_dict["DeepLoc"].split(",")]
    if localisation != ["None"] :
        preys = [prot for prot in preys if prot not in saved_deeploc or any(loc in localisation for loc in saved_deeploc[prot].split("|"))]
    signalp_proteins = baits + preys

    # Signal peptide filter (filter_signalP)
    if Informations_dict["Signal_peptide"] != "None" :
        preys = [prot for prot in preys if prot not in saved_SP or saved_SP[prot] == Informations_dict["Signal_peptide"]]
    feature_proteins = baits + preys

    run_deeploc = Informations_dict["Organism"] != "None" and Informations_dict["Custom_interactions"] == ""
    return {"fetch" : [prot for prot in proteins if not fasta_sequences.get(prot) and prot not in saved_sequences],
            "deeploc" : [prot for prot in deeploc_proteins if prot not in saved_deeploc] if run_deeploc else [],
            "signalp" : [prot for prot in signalp_proteins if prot not in saved_SP or prot not in saved_no_SP] if Informations_dict["Organism"] != "None" else [],
            "msa" : [prot for prot in feature_proteins if not os.path.isfile(f"{Path_Pickle_Feature}/{prot}.a3m")],
            "pkl" : [prot for prot in feature_proteins if not os.path.isfile(f"{Path_Pickle_Feature}/{prot}.pkl")]}


def read_round_jobs (round_id) :
    jobs = []
    path = f"{SLURM_DIR}/round_{round_id}/jobs.tsv"
    if os.path.isfile(path) :
        with open(path, "r") as f :
            for line in f :
                if line.strip() :
                    job, vram = line.rstrip("\n").split("\t")
                    jobs.append((job + "\n", float(vram)))
    return jobs


def job_proteins (job, interaction_type) :
    """
    Proteins (or compound) of a job line : "BAIT,20-200;PREY", "BAIT1;BAIT2;PREY", "PREY:3" or "PROT1,20-200;PROT2:2" (custom).
    """
    return [obj.split(",")[0].split(":")[0] for obj in job.strip().split(";") if obj.strip()]


def job_output_dir (job, interaction_type) :
    """
    Output directory of a job, named as in gpu_job_runner.
    """
    if interaction_type == "homo_int" :
        name = job.strip().replace(":", "_homo_") + "er"
    elif interaction_type == "custom" :
        name = job.strip().replace(";", "_and_").replace(":", "_homo_")
        if "_homo_" in name :
            nbr = name.split("_homo_")[1].split("_")[0]
            name = name.replace(f"_{nbr}", f"_{nbr}er")
    else :
        name = job.strip().replace(";", "_and_")
    return f"./result_{interaction_type}/{name.replace(',', '_')}"


def model_exists (job, interaction_type) :
    outdir = job_output_dir(job, interaction_type)
    return bool(glob.glob(f"{outdir}/ranked_0.pdb") or glob.glob(f"{outdir}/*ranked_0*.pdb") or glob.glob(f"{outdir}/*_model.cif"))


def jax_cache_directory (config, Path_Database, gpu_id) :
    """
    Persistent JAX compilation cache shared by all jobs (one per GPU architecture).
    """
    if config.get("jax_cache_dir") :
        base = config["jax_cache_dir"]
    elif os.access(os.path.dirname(os.path.abspath(Path_Database)), os.W_OK) :
        base = f"{os.path.dirname(os.path.abspath(Path_Database))}/jax_cache"
    else :
        base = os.path.abspath("./jax_cache")
    try :
        cc = subprocess.check_output(["nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader", "-i", str(gpu_id)], text=True).strip().replace(".", "")
    except (subprocess.CalledProcessError, FileNotFoundError) :
        cc = "unknown"
    path = f"{base}/sm_{cc}"
    os.makedirs(path, exist_ok=True)
    return path


# ------------------------------------------------------------------
# Stages
# ------------------------------------------------------------------

def stage_prepare (args, submit_first_round=True) :
    """
    Sequences, DeepLoc, SignalP and all filters (same steps as HInt.py before feature generation).
    """
    setup_logging("prepare")
    config = load_config()
    GPU = allocated_gpus()
    CPU = allocated_cpus()
    Informations_dict = Define_informations()
    link_result_dirs(Informations_dict) # models in Path_Structures if set
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"], Informations_dict["Interact_with"], Informations_dict["AlphaFold"], Informations_dict["Custom_interactions"])

    logger.info(f"Slurm job {os.environ.get('SLURM_JOB_ID', 'local')} | GPUs : {GPU} | CPUs : {CPU}")
    logger.info(f"AlphaFold version set to: {Informations_dict['AlphaFold']}")

    timings = dict()
    summary = [Informations_dict["Interact_with"], Informations_dict["DeepLoc"], Informations_dict["Signal_peptide"], Informations_dict["Min_protein_length"], Informations_dict["Max_protein_length"], Informations_dict["Homo-oligomer"], str(len(HInt_object.get_possible_prey())), Informations_dict["Organism"]]

    for bait in Informations_dict["Interact_with"] : # Check that all bait proteins exist in the protein list
        if bait not in HInt_object.get_proteins() and bait != "" :
            raise Exception(f"Bait {bait} not found in the protein list {Informations_dict['Path_Uniprot_ID']}")

    start = time.time()
    need_msa, need_pkl, need_DeepLoc = (list(dict.fromkeys(need)) for need in HInt_object.check_save_dict(Informations_dict["Path_Pickle_Feature"]))
    timings["Sequences"] = [len(HInt_object.get_proteins()), format_time(time.time() - start)]

    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]])
    if Informations_dict["Min_protein_length"] != "0" or Informations_dict["Max_protein_length"] != "" :
        need_msa, need_pkl, need_DeepLoc = filter_length(HInt_object, Informations_dict, need_msa, need_pkl, need_DeepLoc)

    start = time.time()
    if Informations_dict["Organism"] == "None" :
        HInt_object.set_proteins_sequence_no_SP(HInt_object.get_proteins_sequence_SP()) #don't remove signal peptide
        need_DeepLoc = []
    if need_DeepLoc and Informations_dict["Custom_interactions"] == "" :
        if os.environ.get("CUDA_VISIBLE_DEVICES", "") in ["", "NoDevFiles"] : # prepare submitted without GPU (DeepLoc and SignalP found in log_file/save_dict.pkl)
            raise RuntimeError(f"DeepLoc is needed for {len(need_DeepLoc)} proteins ({', '.join(need_DeepLoc[:20])}) but this job has no GPU. "
                               "Submit the pipeline again, the prepare stage will get a GPU.")
        run_deeploc(HInt_object, Informations_dict["Organism"], need_DeepLoc, GPU)
    timings["Deeploc"] = [len(need_DeepLoc), format_time(time.time() - start)]

    if Informations_dict["DeepLoc"].split(",") != ["None"] :
        need_msa, need_pkl = filter_deeploc(HInt_object, Informations_dict, need_msa, need_pkl)
    HInt_object.Make_save_dict()

    start = time.time()
    need_SP = list()
    for protein in need_msa :
        if protein not in HInt_object.get_proteins_sequence_no_SP().keys() or protein not in HInt_object.get_prot_SP().keys() :
            need_SP.append(protein)
    if Informations_dict["Organism"] == "None" :
        need_SP = []
    if len(need_SP) > 0 :
        need_msa = run_SP(HInt_object, Informations_dict, need_SP, need_msa)
    timings["SignalP"] = [len(need_SP), format_time(time.time() - start)]

    need_msa = check_exist_MSA(HInt_object, Informations_dict, need_msa)
    HInt_object.Make_save_dict()

    if Informations_dict["Custom_interactions"] != "" : # Regions of interest of the custom interactions
        regions_dict = {}
        with open(Informations_dict["Custom_interactions"], "r") as custom_file :
            for line in custom_file :
                for prot in line.strip().split(";") :
                    if "(" in prot :
                        regions_dict[prot.split("(")[0].strip()] = prot.split("(")[1].split(")")[0].strip()
        Informations_dict["Regions"] = regions_dict

    for prot, region in Informations_dict["Regions"].items() : # Adjust protein lengths if specific regions are defined
        if region != "0-0" :
            region_start, region_end = map(int, region.split("-"))
            dict_length = HInt_object.get_length_prot()
            dict_length[prot] = region_end - region_start + 1
            HInt_object.set_length_prot(dict_length)

    need_msa, need_pkl = filter_signalP(HInt_object, Informations_dict, need_msa, need_pkl)
    HInt_object.Make_save_dict()

    if os.path.exists("log_file/shallow_MSA.txt") :
        os.remove("log_file/shallow_MSA.txt")

    gpu_vram = config["gpu_vram"]
    if gpu_vram is None :
        gpu_vram = local_gpu_vram()
        logger.warning(f"--gpu_vram not set, VRAM of this node used for the job lists : {gpu_vram:.1f} GiB")

    rounds = define_rounds(Informations_dict, HInt_object.get_compounds() != {})
    if len(rounds) != config["n_rounds"] :
        raise ValueError(f"{len(rounds)} rounds found but {config['n_rounds']} were submitted, submit the pipeline again.")

    state = {"HInt_object" : HInt_object,
             "Informations_dict" : Informations_dict,
             "need_msa" : need_msa,
             "rounds" : rounds,
             "gpu_vram" : gpu_vram,
             "timings" : timings,
             "summary" : summary}
    if rounds and not config.get("msa_stage", True) : # MSAs found by HInt-slurm submit
        missing_msa = [prot for prot in need_msa if not os.path.isfile(f"{Informations_dict['Path_Pickle_Feature']}/{prot}.a3m")]
        if missing_msa :
            raise RuntimeError(f"No MSA for {len(missing_msa)} proteins ({', '.join(missing_msa[:20])}), e.g. MSAs removed after signal peptide cleavage. "
                               "Submit the pipeline again, the MSA stage will be included.")
        write_need_pkl(state, config)
    if rounds :
        write_round_jobs(state, 0, config)
    save_state(state)

    if rounds :
        if submit_first_round :
            submit_round(0, config, wait_predict_dependency())
    else :
        finalize(state, CPU)
    logger.info("Prepare stage done")


def stage_setup (args) :
    """
    prepare, afdb and features in one job (Jean Zay : prepost partition, internet access, GPU of the node for DeepLoc),
    only the missing steps, then submission of the predict and score jobs of round 0.
    MMseqs2-GPU does not run on the GPU of the prepost nodes : when MSAs are missing, msa (and features) run in the
    msa job on the GPU nodes, and this job submits round 0 after it.
    """
    config = load_config()
    stage_prepare(args, submit_first_round=False)
    if not load_state()["rounds"] :
        return
    if config["afdb_stage"] and config["msa_stage"] :
        stage_afdb(args)
    if config["msa_stage"] :
        setup_logging("setup")
        submit_round(0, config, wait_predict_dependency())
        logger.info("Setup stage done, MSAs and pkl files computed by the msa job")
        return
    if config.get("features_stage", True) :
        args.task = 0
        try :
            stage_features(args)
        except SystemExit :
            logger.error("Some pkl files are missing, the predict jobs of their interactions will fail")
    setup_logging("setup")
    submit_round(0, config)
    logger.info("Setup stage done")


def stage_msa_features (args) :
    """
    msa then features in one job on the GPU nodes (MMseqs2-GPU), used with the local species database
    (the GPU nodes have no internet access).
    """
    config = load_config()
    stage_msa(args)
    if config.get("features_stage", True) :
        args.task = 0
        stage_features(args)


def fetch_sequences (Informations_dict) :
    """
    Download the UniProt sequences, they are saved in log_file/save_dict.pkl for the prepare stage.
    """
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"], Informations_dict["Interact_with"], Informations_dict["AlphaFold"], Informations_dict["Custom_interactions"])
    HInt_object.check_save_dict(Informations_dict["Path_Pickle_Feature"])
    logger.info(f"Sequences of {len(HInt_object.get_proteins())} proteins saved")


def stage_fetch (args) :
    """
    Download the UniProt sequences on a node with internet access (done by HInt-slurm submit on the login node).
    """
    setup_logging("fetch")
    fetch_sequences(Define_informations())


def stage_afdb (args) :
    """
    Download the MSAs available in the AlphaFold database on a node with internet access.
    """
    setup_logging("afdb")
    CPU = allocated_cpus()
    state = load_state()
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    start = time.time()
    need_msa = [prot for prot in state["need_msa"] if not os.path.isfile(f"{Informations_dict['Path_Pickle_Feature']}/{prot}.a3m")]
    if need_msa :
        create_feature(HInt_object, Informations_dict, [], CPU, need_msa, [], afdb_only=True)
    state["timings"]["AlphaFold database MSA"] = [len(state["need_msa"]) - len(need_msa), format_time(time.time() - start)]
    save_state(state)
    logger.info(f"{len(state['need_msa']) - len(need_msa)} MSAs found in the AlphaFold database")


def stage_msa (args) :
    """
    Generate all missing MSAs at once (MMseqs2-GPU is efficient on large batches), then list proteins without pkl file.
    """
    setup_logging("msa")
    config = load_config()
    GPU = allocated_gpus()
    CPU = allocated_cpus()
    state = load_state()
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]

    start = time.time()
    need_msa = [prot for prot in state["need_msa"] if not os.path.isfile(f"{Path_Pickle_Feature}/{prot}.a3m")] # MSAs downloaded by the afdb stage
    nb_msa = len(need_msa)
    if need_msa :
        use_afdb = config["use_afdb"] and not config.get("afdb_stage", False)
        create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, [], use_afdb=use_afdb, local_search_threshold=config["local_search_threshold"])
    state["timings"]["MSA"] = [nb_msa, format_time(time.time() - start)]
    write_need_pkl(state, config)
    save_state(state)


def stage_features (args) :
    """
    Generate the feature pickle files of one slice of the protein list (one Slurm array task).
    """
    task = args.task
    setup_logging(f"features_task_{task}")
    config = load_config()
    CPU = allocated_cpus()
    state = load_state()
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]
    re_search_templates_mmseqs2 = "False" # as in create_feature

    with open(NEED_PKL_FILE, "r") as f :
        need_pkl = [line.strip() for line in f if line.strip()]
    my_proteins = need_pkl[task::config["features_tasks"]]
    logger.info(f"Features task {task} : {len(my_proteins)} proteins with {CPU} CPUs")

    # One species cache per task : SQLite (WAL) cannot be shared by several nodes on a network file system
    cache_name = "uniprot_species_cache.sqlite" if config["features_tasks"] == 1 else f"uniprot_species_cache_task_{task}.sqlite"
    set_species_env(Path_Pickle_Feature, Informations_dict["Path_Database"], cache_name)
    if "ALPHAPULLDOWN_LOCAL_SPECIES_DB" not in os.environ :
        logger.warning(f"No local species database in {Informations_dict['Path_Database']}, species used to pair MSAs are fetched from UniProt (needs internet access)")

    prot_no_SP = HInt_object.get_proteins_sequence_no_SP()
    log_path = HINT_LOG
    with ThreadPoolExecutor(max_workers=max(1, CPU)) as executor :
        for protein in my_proteins :
            if os.path.isfile(f"{Path_Pickle_Feature}/{protein}.pkl") :
                continue
            pkl_file = f"./log_file/{protein}_pkl.fasta"
            cmd = feature_cmd(pkl_file, Informations_dict["Path_Database"], Path_Pickle_Feature, re_search_templates_mmseqs2)
            executor.submit(create_ind_feature, protein, pkl_file, prot_no_SP, cmd, log_path)

    missing = [prot for prot in my_proteins if not os.path.isfile(f"{Path_Pickle_Feature}/{prot}.pkl")]
    if missing :
        logger.error(f"No pkl file for {len(missing)} proteins : {', '.join(missing)}")
        sys.exit(1)
    logger.info(f"Features task {task} done")


def stage_predict (args) :
    """
    Run the AlphaFold jobs of one slice of a round on the GPUs allocated to this Slurm task.
    Jobs are packed on the GPUs according to their predicted VRAM, as in HInt.py.
    """
    round_id = args.round
    task = args.task
    setup_logging(f"predict_round_{round_id}_task_{task}")
    config = load_config()
    state = load_state()
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    Path_Database = Informations_dict["Path_Database"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]
    AF_version = Informations_dict["AlphaFold"]
    interaction_type, bait = state["rounds"][round_id]

    all_jobs = read_round_jobs(round_id)
    n_tasks = predict_task_count(len(all_jobs), config)
    jobs = all_jobs[task::n_tasks] if task < n_tasks else []
    status_file = f"{round_dir(round_id)}/task_{task}.json"
    if not jobs :
        logger.info(f"No job for task {task} of round {round_id}")
        with open(status_file, "w") as f :
            json.dump({"done" : [], "failed" : [], "skipped" : []}, f)
        return

    GPU = allocated_gpus()
    max_vram = config["gpu_vram"] if config["gpu_vram"] is not None else local_gpu_vram()
    multi_job_per_gpu = config["multi_job_per_gpu"] == "True"
    jax_cache_dir = jax_cache_directory(config, Path_Database, GPU[0])
    compounds_dict = HInt_object.get_compounds()
    seq_bait = HInt_object.get_proteins_sequence_no_SP()
    log_path = HINT_LOG
    logger.info(f"Round {round_id} ({interaction_type} {bait}) task {task} : {len(jobs)} jobs | GPUs : {GPU} | VRAM per GPU : {max_vram:.1f} GiB | multi job per GPU : {multi_job_per_gpu}")

    done, failed, skipped = [], [], []
    pending = []
    for job, vram in jobs :
        if model_exists(job, interaction_type) :
            skipped.append(job.strip())
        else :
            pending.append((job, vram, 0))

    result_queue = multiprocessing.Manager().Queue()
    gpu_vram_used = {gpu : 0.0 for gpu in GPU}
    jobs_running = {gpu : [] for gpu in GPU}
    start = time.time()

    while pending or any(jobs_running.values()) :
        while not result_queue.empty() :
            result_queue.get()

        # Collect finished jobs
        for gpu_id in GPU :
            for item in list(jobs_running[gpu_id]) :
                job, p, vram, attempt = item
                if p.exitcode is None :
                    continue
                p.join()
                jobs_running[gpu_id].remove(item)
                gpu_vram_used[gpu_id] = max(0.0, gpu_vram_used[gpu_id] - vram)
                if p.exitcode == 0 :
                    done.append(job.strip())
                elif attempt < config["max_retries"] :
                    logger.error(f"[GPU {gpu_id}] Job {job.strip()} crash (exit code {p.exitcode}), retry {attempt + 1}/{config['max_retries']}")
                    pending.append((job, vram, attempt + 1))
                else :
                    logger.error(f"[GPU {gpu_id}] Job {job.strip()} crash (exit code {p.exitcode}), abandoned")
                    failed.append(job.strip())

        # Launch new jobs
        launched = False
        for item in list(pending) :
            job, vram, attempt = item
            if interaction_type == "Compounds" :
                inputs = [f"{Path_Pickle_Feature}/{prot}.a3m" for prot in job_proteins(job, interaction_type) if prot in Informations_dict["Interact_with"]]
            else :
                inputs = [f"{Path_Pickle_Feature}/{prot}.pkl" for prot in job_proteins(job, interaction_type)]
            if not all(os.path.isfile(path) for path in inputs) :
                logger.error(f"Missing features for {job.strip()}, job abandoned")
                pending.remove(item)
                failed.append(job.strip())
                continue
            gpu_id = min(GPU, key=lambda g : gpu_vram_used[g])
            free = max_vram - gpu_vram_used[gpu_id]
            fits = vram <= free or not jobs_running[gpu_id] # a job larger than the GPU is run alone (unified memory)
            if fits and (multi_job_per_gpu or not jobs_running[gpu_id]) :
                compound = None
                if interaction_type == "Compounds" :
                    compound = compounds_dict[job.strip().split(";")[-1]]
                p = multiprocessing.Process(target=gpu_job_runner,
                    args=(gpu_id, job, vram, result_queue, Path_Database, Path_Pickle_Feature, interaction_type, AF_version, seq_bait, Informations_dict["Interact_with"], compound, jax_cache_dir, log_path))
                p.start()
                pending.remove(item)
                gpu_vram_used[gpu_id] += vram
                jobs_running[gpu_id].append((job, p, vram, attempt))
                logger.info(f"[GPU {gpu_id}] Job started: '{job.strip()}' | PID: {p.pid} | Predicted VRAM: {vram:.2f} GB | Predicted total VRAM usage: {gpu_vram_used[gpu_id]:.2f} GB | Active jobs: {len(jobs_running[gpu_id])}")
                launched = True
                break
        if not launched :
            time.sleep(2)

    with open(status_file, "w") as f :
        json.dump({"done" : done, "failed" : failed, "skipped" : skipped, "time" : format_time(time.time() - start)}, f)
    logger.info(f"Round {round_id} task {task} : {len(done)} done, {len(skipped)} already modeled, {len(failed)} failed")


def stage_score (args) :
    """
    Score one round, then prepare the next round or write the final results.
    """
    round_id = args.round
    setup_logging(f"score_round_{round_id}")
    config = load_config()
    CPU = allocated_cpus()
    state = load_state()
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    interaction_type, bait = state["rounds"][round_id]

    # Summary of the predict tasks
    done, failed = 0, []
    times = []
    unfinished = []
    for task in range(predict_task_count(len(read_round_jobs(round_id)), config)) :
        status_file = f"{SLURM_DIR}/round_{round_id}/task_{task}.json"
        if not os.path.isfile(status_file) :
            unfinished.append(str(task))
            continue
        with open(status_file, "r") as f :
            status = json.load(f)
        done += len(status["done"])
        failed.extend(status["failed"])
        if "time" in status :
            times.append(status["time"])
    if unfinished :
        # Unmodeled interactions would be scored as negative, stop the pipeline instead
        logger.error(f"Predict tasks {', '.join(unfinished)} of round {round_id} did not finish (timeout or crash), check their Slurm logs in {SLURM_DIR}/. "
                     f"Resume with : HInt-slurm submit <same options> --start_at predict --start_round {round_id} (models already done are skipped)")
        sys.exit(1)
    with open(f"{round_dir(round_id)}/failed_jobs.txt", "w") as f :
        f.write("\n".join(failed) + ("\n" if failed else ""))
    if failed :
        logger.warning(f"{len(failed)} jobs failed in round {round_id}, see {round_dir(round_id)}/failed_jobs.txt")
    state["timings"][f"PPI round {round_id} ({interaction_type} {bait})"] = [len(read_round_jobs(round_id)), f"max task time {max(times) if times else 'NA'}"]

    start = time.time()
    if os.path.exists("log_file/shallow_MSA.txt") :
        os.remove("log_file/shallow_MSA.txt")
    Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"], HInt_object.get_possible_prey())
    if interaction_type in ["homo_int", "custom"] :
        Score_interaction(HInt_object, Informations_dict, CPU, interaction_type, config["send_to_database"], "", config["multi_scoring"], "")
    else :
        Score_interaction(HInt_object, Informations_dict, CPU, interaction_type, config["send_to_database"], "", config["multi_scoring"], bait)
    HInt_object.Make_save_dict()
    state["timings"][f"Scoring round {round_id}"] = [len(HInt_object.get_possible_prey()), format_time(time.time() - start)]

    if round_id + 1 < len(state["rounds"]) :
        write_round_jobs(state, round_id + 1, config)
        save_state(state)
        submit_round(round_id + 1, config)
    else :
        save_state(state)
        finalize(state, CPU)
    logger.info(f"Score stage of round {round_id} done")


def finalize (state, CPU) :
    """
    Final results, figures and report (same outputs as HInt.py).
    """
    HInt_object = state["HInt_object"]
    Informations_dict = state["Informations_dict"]
    sorted_protein = Resume_file(HInt_object, Informations_dict)
    if Informations_dict["Interact_with"] != [''] and HInt_object.get_compounds() == {} :
        Create_figures(HInt_object, Informations_dict, Informations_dict["AlphaFold"], sorted_protein, CPU)
    summary = state["summary"]
    line_report = f"===== SUMMARY INFO =====\nBait : {summary[0]}\nCellular localization : {summary[1]}\nPeptide signal : {summary[2]}\nMinimum protein length : {summary[3]}\nMaximum protein length : {summary[4]}\nHomo-oligomer state : {summary[5]}\nNumber of preys : {summary[6]}\nOrganism : {summary[7]}\nExecution : Slurm\n\n===== TIMINGS =====\nStep\tNumber of proteins\tTime\n"
    for step, (number, step_time) in state["timings"].items() :
        line_report += f"{step}\t{number}\t{step_time}\n"
    with open("./log_file/HInt_report.txt", "w") as f :
        f.write(line_report)
    logger.info("All steps done")


# ------------------------------------------------------------------
# Submission
# ------------------------------------------------------------------

DEFAULT_RESOURCES = {
    "fetch" : {"cpus" : 2, "mem" : "8G", "time" : "00:30:00"},
    "prepare" : {"cpus" : 2, "mem" : "32G", "time" : "01:00:00"},
    "afdb" : {"cpus" : 8, "mem" : "16G", "time" : "12:00:00"},
    "msa" : {"cpus" : 32, "mem" : "128G", "time" : "24:00:00"},
    "features" : {"cpus" : 16, "mem" : "64G", "time" : "12:00:00"},
    "predict" : {"cpus" : 8, "mem" : "64G", "time" : "48:00:00"},
    "score" : {"cpus" : 32, "mem" : "64G", "time" : "12:00:00"},
}

# IDRIS Jean Zay : GPU nodes are selected with a constraint and an account <project>@<gpu>, the memory is proportional
# to the number of CPUs (no --mem), compute nodes have no internet access, the prepost partition has it (max 20h).
JEAN_ZAY_GPUS = {
    "a100" : {"constraint" : "a100", "cpus_per_gpu" : 8, "vram" : 80, "qos" : "qos_gpu_a100", "setup" : []},
    "h100" : {"constraint" : "h100", "cpus_per_gpu" : 24, "vram" : 80, "qos" : "qos_gpu_h100", "setup" : ["module load arch/h100"]},
    "v100" : {"constraint" : "v100-32g", "cpus_per_gpu" : 10, "vram" : 32, "qos" : "qos_gpu", "setup" : []},
}
JEAN_ZAY_PARTITIONS = {"afdb" : "archive"} # downloads only, the other CPU stages run on prepost (fetch runs on the login node)
JEAN_ZAY_RESOURCES = {
    "fetch" : {"cpus" : 2, "time" : "00:30:00"},
    "prepare" : {"cpus" : 2, "time" : "01:00:00"},
    "afdb" : {"cpus" : 8, "time" : "04:00:00"},
    "msa" : {"cpus" : None, "time" : "10:00:00"},
    "features" : {"cpus" : 16, "time" : "00:30:00"},
    "predict" : {"cpus" : None, "time" : "05:00:00"},
    "score" : {"cpus" : 16, "time" : "04:00:00"},
}


def hours (time_limit) :
    """
    Slurm time (D-HH:MM:SS, HH:MM:SS or MM) in hours.
    """
    days = 0
    if "-" in time_limit :
        days, time_limit = time_limit.split("-")
    parts = [int(p) for p in time_limit.split(":")]
    if len(parts) == 1 :
        parts = [0, parts[0], 0]
    while len(parts) < 3 :
        parts.append(0)
    return int(days) * 24 + parts[0] + parts[1] / 60 + parts[2] / 3600


def stage_settings (args, stage, gpus, internet) :
    """
    Slurm settings of one stage.

    Parameters :
    ----------
    args : argparse.Namespace
    stage : str
    gpus : int
    internet : bool
        The stage needs internet access (run on --internet_partition when compute nodes are offline).

    Returns :
    ----------
    settings : dict
    """
    jean_zay = args.cluster == "jean-zay"
    defaults = JEAN_ZAY_RESOURCES[stage] if jean_zay else DEFAULT_RESOURCES[stage]
    settings = {"cpus" : getattr(args, f"{stage}_cpus") or defaults["cpus"],
                "mem" : getattr(args, f"{stage}_mem") or (None if jean_zay else defaults["mem"]),
                "time" : getattr(args, f"{stage}_time") or defaults["time"],
                "gres" : None, "constraint" : None, "partition" : None, "account" : args.account, "qos" : args.qos,
                "extra" : list(args.sbatch_extra), "setup" : list(args.setup)}
    if settings["mem"] in ["none", "None"] :
        settings["mem"] = None
    on_internet = internet and args.offline_compute and not gpus

    if jean_zay :
        gpu = JEAN_ZAY_GPUS[args.jz_gpu]
        settings["extra"].insert(0, "--hint=nomultithread")
        if gpus :
            settings["constraint"] = gpu["constraint"]
            settings["gres"] = f"gpu:{gpus}"
            settings["account"] = args.account or f"{args.jz_project}@{args.jz_gpu}"
            settings["cpus"] = getattr(args, f"{stage}_cpus") or gpu["cpus_per_gpu"] * gpus
            if not settings["qos"] :
                settings["qos"] = f"{gpu['qos']}-t4" if hours(settings["time"]) > 20 else f"{gpu['qos']}-t3"
            settings["setup"] = gpu["setup"] + settings["setup"]
        else :
            if hours(settings["time"]) > 20 :
                raise ValueError(f"--{stage}_time : the prepost and archive partitions of Jean Zay are limited to 20h")
            on_internet = True # all CPU stages run on prepost (internet access, hours not deducted)
            settings["account"] = args.jz_internet_account or f"{args.jz_project}@cpu"
            settings["qos"] = None
        if on_internet :
            settings["partition"] = getattr(args, f"{stage}_partition") or JEAN_ZAY_PARTITIONS.get(stage) or args.internet_partition or "prepost"
        else :
            settings["partition"] = getattr(args, f"{stage}_partition")
        return settings

    if gpus :
        gpu_type = f"{args.gpu_type}:" if args.gpu_type else ""
        settings["gres"] = f"gpu:{gpu_type}{gpus}"
    if on_internet :
        settings["partition"] = getattr(args, f"{stage}_partition") or args.internet_partition
        settings["account"] = args.internet_account or args.account
        settings["qos"] = args.internet_qos
    else :
        settings["partition"] = getattr(args, f"{stage}_partition") or (args.gpu_partition if gpus else args.cpu_partition) or args.partition
    return settings


def sbatch_script (args, name, command, settings, array=None) :
    """
    Write a sbatch script and return its path.
    """
    workdir = os.getcwd()
    output = f"{SLURM_DIR}/%x_%A_%a.out" if array else f"{SLURM_DIR}/%x_%j.out"
    lines = ["#!/bin/bash",
             f"#SBATCH --job-name={args.job_prefix}_{name}",
             f"#SBATCH --chdir={workdir}",
             f"#SBATCH --output={output}",
             "#SBATCH --nodes=1",
             "#SBATCH --ntasks=1",
             f"#SBATCH --cpus-per-task={settings['cpus']}",
             f"#SBATCH --time={settings['time']}"]
    for key, option in [("mem", "--mem"), ("partition", "--partition"), ("account", "--account"), ("qos", "--qos"), ("constraint", "--constraint"), ("gres", "--gres")] :
        if settings[key] :
            lines.append(f"#SBATCH {option}={settings[key]}")
    if array :
        lines.append(f"#SBATCH --array={array}")
    for extra in settings["extra"] :
        lines.append(f"#SBATCH {extra}")
    lines += ["", "set -eo pipefail", ""]
    if args.cluster == "jean-zay" :
        lines.append("module purge")
    lines += settings["setup"]
    if not args.no_conda :
        conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda") or "conda"
        conda_env = args.conda_env or os.environ.get("CONDA_PREFIX") or "HInt"
        lines.append(f'source "$({shlex.quote(conda_exe)} info --base)/etc/profile.d/conda.sh"')
        lines.append(f"conda activate {shlex.quote(conda_env)}")
    if args.cluster == "jean-zay" : # models downloaded by DeepLoc (torch hub, Hugging Face) in $WORK, $HOME is limited to 3 GB
        lines += ['export TORCH_HOME="${TORCH_HOME:-$WORK/.cache/torch}"',
                  'export HF_HOME="${HF_HOME:-$WORK/.cache/huggingface}"']
    lines += ["", 'echo "HInt ${SLURM_JOB_ID} on $(hostname) | CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-none}"', command, ""]
    path = f"{SCRIPT_DIR}/{name}.sbatch"
    with open(path, "w") as f :
        f.write("\n".join(lines))
    return path


def sbatch (args, script, dependency=None) :
    """
    Submit a script and return the Slurm job id.
    """
    cmd = ["sbatch", "--parsable"]
    if dependency :
        cmd += [f"--dependency={dependency}", "--kill-on-invalid-dep=yes"]
    cmd.append(script)
    if args.dry_run :
        fake_id = f"<{os.path.basename(script).replace('.sbatch', '')}>"
        print("[dry run] " + " ".join(cmd))
        return fake_id
    # A job submitted from another job must not inherit its Slurm allocation (SLURM_*, SBATCH_*, GPUs)
    env = {key : value for key, value in os.environ.items() if not key.startswith(("SLURM_", "SBATCH_")) and key != "CUDA_VISIBLE_DEVICES"}
    job_id = subprocess.check_output(cmd, text=True, env=env).strip().split(";")[0]
    print(f"Submitted {os.path.basename(script):<28} job {job_id}" + (f"  (dependency {dependency})" if dependency else ""))
    return job_id


def stage_submit (args) :
    """
    Write all sbatch scripts and submit the whole pipeline with dependencies.
    """
    if not os.path.isfile("HInt.txt") :
        raise FileNotFoundError("HInt.txt not found, run HInt-slurm from the screen directory")
    if not args.dry_run and shutil.which("sbatch") is None :
        raise RuntimeError("sbatch not found, run HInt-slurm submit on a Slurm login node (or use --dry_run)")
    if args.gpus_per_task is None :
        args.gpus_per_task = 4 if args.cluster == "jean-zay" else 1
    if args.cluster == "jean-zay" :
        if not args.jz_project :
            raise ValueError("--jz_project is needed on Jean Zay (Unix group of the project, see idrproj)")
        args.offline_compute = True
        if args.gpu_vram is None :
            args.gpu_vram = JEAN_ZAY_GPUS[args.jz_gpu]["vram"]
    elif args.offline_compute and not args.internet_partition :
        raise ValueError("--offline_compute needs --internet_partition (a partition whose nodes have internet access)")
    Informations_dict = Define_informations()
    link_result_dirs(Informations_dict) # models in Path_Structures if set
    n_rounds = len(define_rounds(Informations_dict, False))
    # With the local species database, pkl files need no internet access and are fast : one CPU job
    local_species = os.path.isfile(os.path.join(Informations_dict["Path_Database"], LOCAL_SPECIES_DB))
    # Jean Zay : prepare, afdb, msa and features in one job on prepost (internet access and one GPU)
    single_setup = args.cluster == "jean-zay" and not args.split_setup and STAGES.index(args.start_at) <= STAGES.index("prepare")
    if single_setup and hours(args.setup_time) > 20 :
        raise ValueError("--setup_time : the prepost partition of Jean Zay is limited to 20h")
    if local_species or single_setup :
        args.features_tasks = 1
        args.features_cpus = args.features_cpus or LOCAL_SPECIES_FEATURES_CPUS
    use_afdb = not args.no_afdb
    if not args.dry_run :
        open(HINT_LOG, "w").close() # new log for each submission
        if args.start_at in ["fetch", "prepare"] : # new run from the beginning
            clean_previous_run()
    # Only the missing steps are submitted
    missing = missing_work(Informations_dict)
    # UniProt sequences downloaded here : the login node has internet access
    fetched_here = False
    if missing["fetch"] and STAGES.index(args.start_at) <= STAGES.index("fetch") :
        if args.dry_run :
            print(f"[dry run] {len(missing['fetch'])} UniProt sequences would be downloaded on this node")
        else :
            print(f"Downloading {len(missing['fetch'])} UniProt sequences on this node ...")
            setup_logging("fetch")
            fetch_sequences(Informations_dict)
            fetched_here = True
            missing = missing_work(Informations_dict)
    gpu_for_prepare = bool(missing["deeploc"] or missing["signalp"]) # prepare on CPU only when DeepLoc and SignalP are already done
    need_msa = bool(missing["msa"])
    need_features = bool(missing["pkl"])
    if args.gpu_vram is None and not gpu_for_prepare :
        raise ValueError("--gpu_vram is needed : the prepare stage runs without GPU (DeepLoc and SignalP already done)")

    for value, flag in [(args.multi_job_per_gpu, "--multi_job_per_gpu"), (args.multi_scoring, "--multi_scoring")] :
        if value not in ["True", "False"] :
            raise ValueError(f"Invalid value for {flag}. Need True or False.")
    if args.send_to_database not in ["All", "Negative", "None"] :
        raise ValueError("Invalid value for --send_to_database. Need All, Negative or None.")
    if args.start_at in ["predict", "score"] and not os.path.isfile(STATE_FILE) :
        raise FileNotFoundError(f"{STATE_FILE} not found, start the pipeline at prepare")

    config = {"gpu_vram" : args.gpu_vram,
              "multi_job_per_gpu" : args.multi_job_per_gpu,
              "multi_scoring" : args.multi_scoring,
              "send_to_database" : args.send_to_database,
              "max_retries" : args.max_retries,
              "use_afdb" : use_afdb,
              "afdb_stage" : use_afdb and args.offline_compute, # AlphaFold database downloaded on a node with internet access
              "local_search_threshold" : args.local_search_threshold,
              "jax_cache_dir" : args.jax_cache_dir,
              "features_tasks" : args.features_tasks,
              "jobs_per_task" : args.jobs_per_task,
              "gpus_per_task" : args.gpus_per_task,
              "msa_stage" : need_msa,
              "features_stage" : need_features,
              "n_rounds" : n_rounds,
              "submit_args" : vars(args)} # the predict and score jobs are submitted by the prepare and score jobs
    if not args.dry_run or not os.path.isfile(CONFIG_FILE) :
        with open(CONFIG_FILE, "w") as f :
            json.dump(config, f, indent=2)
    os.makedirs(SCRIPT_DIR, exist_ok=True)

    throttle = f"%{args.max_parallel}" if args.max_parallel else ""
    hint = "HInt-slurm"
    order = STAGES.index(args.start_at)
    first_round = args.start_round if args.start_at == "predict" else 0
    if order > STAGES.index("prepare") and n_rounds and not os.path.isfile(f"{SLURM_DIR}/round_{first_round}/jobs.tsv") :
        raise FileNotFoundError(f"No job list for round {first_round} : start the pipeline at prepare" + (" or at the previous round" if first_round else ""))
    last = None
    prepare_id = None
    if not args.dry_run :
        for path in [JOB_IDS_FILE, f"{SLURM_DIR}/job_ids.txt"] :
            if os.path.isfile(path) :
                os.remove(path)

    def submit (name, stage, command, gpus=0, internet=False, array=None, dependency="afterok", settings=None) :
        nonlocal last
        if last is None and array is None and not fetched_here : # the first job resets log_file/HInt.log, also when its script is submitted by hand
            command += " --new_log"
        script = sbatch_script(args, name, command, settings or stage_settings(args, stage, gpus, internet), array)
        last = sbatch(args, script, f"{dependency}:{last}" if last else None)
        if not args.dry_run :
            record_job_id(name, last)
        return last

    skipped = []
    if not need_msa :
        skipped.append("afdb/msa (all MSAs found)")
    if not need_features :
        skipped.append("features (all pkl files found)")
    if single_setup and n_rounds :
        setup_gpus = args.setup_gpus if gpu_for_prepare else 0 # DeepLoc (MMseqs2-GPU runs in the msa job on the GPU nodes)
        settings = stage_settings(args, "prepare", 0, True)
        settings.update({"cpus" : args.setup_cpus, "time" : args.setup_time, "partition" : args.setup_partition})
        if setup_gpus : # DeepLoc and MMseqs2-GPU use the GPU of the prepost node, reachable without --gres
            settings["setup"] = settings["setup"] + ['export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"']
        prepare_id = submit("setup", "prepare", f"{hint} setup", settings=settings)
    elif order <= STAGES.index("prepare") :
        prepare_id = submit("prepare", "prepare", f"{hint} prepare", gpus=3 if gpu_for_prepare else 0, internet=not args.offline_compute)
    if n_rounds == 0 :
        print("No bait, no homo-oligomer and no custom interaction : only the prepare stage is needed.")
    else :
        if config["afdb_stage"] and need_msa and order <= STAGES.index("afdb") and not single_setup :
            submit("afdb", "afdb", f"{hint} afdb", internet=True)
        # With the setup job, features runs in the msa job (GPU nodes, no internet access : local species database needed)
        msa_with_features = single_setup and need_msa and need_features and local_species
        if need_msa and order <= STAGES.index("msa") : # GPU nodes : MMseqs2-GPU
            submit("msa", "msa", f"{hint} msa_features" if msa_with_features else f"{hint} msa", gpus=args.msa_gpus)
        if need_features and order <= STAGES.index("features") and (not single_setup or need_msa) and not msa_with_features :
            if local_species :
                submit("features", "features", f"{hint} features --task 0")
            else : # species used to pair MSAs are fetched from UniProt
                submit("features", "features", f"{hint} features --task ${{SLURM_ARRAY_TASK_ID}}", internet=True, array=f"0-{args.features_tasks - 1}{throttle}")
            next_dependency = "afterany" # missing features are reported by the predict tasks
        elif msa_with_features :
            next_dependency = "afterany"
        else :
            next_dependency = "afterok"
        # The predict and score jobs of a round are submitted once its job list is known : by prepare for round 0,
        # by the score job of the previous round for the next ones, or here when prepare is not submitted (resumed run)
        predict_dependency = f"{next_dependency}:{last}" if last and last != prepare_id else ""
        if prepare_id is None :
            submit_round(first_round, config, predict_dependency or None, dry_run=args.dry_run)
        elif not args.dry_run : # read by the prepare job
            with open(JOB_IDS_FILE, "r") as f :
                named_ids = json.load(f)
            named_ids["predict_dependency"] = predict_dependency
            with open(JOB_IDS_FILE, "w") as f :
                json.dump(named_ids, f)

    if not args.dry_run :
        print(f"\nRounds : {n_rounds} | features tasks : {args.features_tasks}{' (local species database)' if local_species else ''} | predict : one task per {args.jobs_per_task} PPI, up to {args.gpus_per_task} GPUs per task, submitted once the job list of each round is known")
        if skipped :
            print(f"Skipped : {', '.join(skipped)} | prepare on {'GPU (DeepLoc/SignalP)' if gpu_for_prepare else 'CPU (DeepLoc and SignalP already done)'}")
        print(f"Follow the pipeline : squeue -u $USER   or   HInt-slurm status")
        print(f"Cancel the pipeline : scancel $(cat {SLURM_DIR}/job_ids.txt)")
    check_pi_score_env()


def check_pi_score_env () :
    """
    The scoring needs the pi_score conda environment, compute nodes may not have internet access to create it.
    """
    tool = next((t for t in ["micromamba", "mamba", "conda"] if shutil.which(t)), None)
    if tool is None :
        return
    try :
        envs = json.loads(subprocess.run([tool, "env", "list", "--json"], capture_output=True, text=True, check=True).stdout)["envs"]
    except Exception :
        return
    if not any("pi_score" in env for env in envs) :
        print(f"\nWARNING : conda environment pi_score not found. If compute nodes have no internet access, create it now :\n  {tool} create -y -n pi_score python=2.7 scikit-learn=0.20.4 biopython biopandas")


def stage_status (args) :
    """
    Print the progress of the pipeline.
    """
    if not os.path.isfile(CONFIG_FILE) :
        print("No Slurm run found in this directory.")
        return
    config = load_config()
    if os.path.isfile(f"{SLURM_DIR}/job_ids.txt") and shutil.which("squeue") :
        with open(f"{SLURM_DIR}/job_ids.txt") as f :
            ids = f.read().split()
        if ids :
            subprocess.run(["squeue", "--jobs", ",".join(ids), "--format=%.18i %.30j %.10T %.12M %.20R"])
    if not os.path.isfile(STATE_FILE) :
        print("prepare : not done")
        return
    state = load_state()
    print("prepare : done")
    if os.path.isfile(NEED_PKL_FILE) :
        with open(NEED_PKL_FILE) as f :
            need_pkl = [line.strip() for line in f if line.strip()]
        path = state["Informations_dict"]["Path_Pickle_Feature"]
        print(f"msa : done | features : {sum(os.path.isfile(f'{path}/{p}.pkl') for p in need_pkl)}/{len(need_pkl)} pkl files")
    else :
        print("msa : not done")
    for round_id, (interaction_type, bait) in enumerate(state["rounds"]) :
        jobs = read_round_jobs(round_id)
        if not os.path.isfile(f"{SLURM_DIR}/round_{round_id}/jobs.tsv") :
            print(f"round {round_id} ({interaction_type} {bait}) : waiting for the previous round")
            continue
        modeled = sum(model_exists(job, interaction_type) for job, _ in jobs)
        n_tasks = predict_task_count(len(jobs), config)
        tasks = sum(os.path.isfile(f"{SLURM_DIR}/round_{round_id}/task_{t}.json") for t in range(n_tasks))
        scored = os.path.isfile(f"{SLURM_DIR}/round_{round_id}/failed_jobs.txt")
        print(f"round {round_id} ({interaction_type} {bait}) : {modeled}/{len(jobs)} models | {tasks}/{n_tasks} predict tasks finished | scoring {'done' if scored else 'not done'}")
    if os.path.isfile("./log_file/HInt_report.txt") and os.path.getmtime("./log_file/HInt_report.txt") > os.path.getmtime(CONFIG_FILE) :
        print("final results : All_Final_result_HInt.csv")


# ------------------------------------------------------------------
# Command line
# ------------------------------------------------------------------

def add_submit_arguments (parser) :
    cluster = parser.add_argument_group("Cluster")
    cluster.add_argument("--cluster", help="Cluster profile: generic or jean-zay (IDRIS) (default: generic)", choices=["generic", "jean-zay"], default="generic")
    cluster.add_argument("--jz_project", help="Jean Zay: Unix group of the project (account <project>@<gpu>, see idrproj)", default=None)
    cluster.add_argument("--jz_gpu", help="Jean Zay: GPU partition, a100, h100 or v100 (default: a100)", choices=list(JEAN_ZAY_GPUS), default="a100")
    cluster.add_argument("--jz_internet_account", help="Jean Zay: account of the prepost jobs (default: <project>@cpu)", default=None)
    cluster.add_argument("--offline_compute", help="Compute nodes have no internet access: run sequence download, AlphaFold database MSAs, pkl files (UniProt species without local database) and scoring on --internet_partition", action="store_true")
    cluster.add_argument("--internet_partition", help="Partition with internet access (Jean Zay: prepost)", default=None)
    cluster.add_argument("--internet_account", help="Account of the jobs on --internet_partition (default: --account)", default=None)
    cluster.add_argument("--internet_qos", help="QOS of the jobs on --internet_partition", default=None)

    slurm = parser.add_argument_group("Slurm")
    slurm.add_argument("--partition", help="Partition used for all jobs", default=None)
    slurm.add_argument("--gpu_partition", help="Partition of GPU jobs (default: --partition)", default=None)
    slurm.add_argument("--cpu_partition", help="Partition of CPU jobs (default: --partition)", default=None)
    slurm.add_argument("--account", help="Slurm account (Jean Zay default: <project>@<gpu>)", default=None)
    slurm.add_argument("--qos", help="Slurm QOS of the GPU and CPU jobs (Jean Zay default: qos_gpu_<gpu>-t3, or -t4 above 20h)", default=None)
    slurm.add_argument("--gpu_type", help="GPU type in --gres, e.g. a100 (gres=gpu:a100:N)", default=None)
    slurm.add_argument("--sbatch_extra", help="Extra #SBATCH line, e.g. '--constraint=a100' (repeatable)", action="append", default=[])
    slurm.add_argument("--setup", help="Shell line run before HInt in each job, e.g. 'module load cuda/12.4' (repeatable)", action="append", default=[])
    slurm.add_argument("--conda_env", help="Conda environment activated in each job (default: the current one)", default=None)
    slurm.add_argument("--no_conda", help="Do not activate a conda environment (use --setup instead)", action="store_true")
    slurm.add_argument("--job_prefix", help="Prefix of the Slurm job names (default: HInt)", default="HInt")
    slurm.add_argument("--max_parallel", help="Maximum number of array tasks running at the same time", type=int, default=None)
    slurm.add_argument("--dry_run", help="Write the sbatch scripts without submitting them", action="store_true")

    work = parser.add_argument_group("Work distribution")
    work.add_argument("--features_tasks", help=f"Number of array tasks for pkl files (default: 10, one job with {LOCAL_SPECIES_FEATURES_CPUS} CPUs when Path_Database contains {LOCAL_SPECIES_DB})", type=int, default=10)
    work.add_argument("--jobs_per_task", help="Maximum models (PPI) per predict task, one task per jobs_per_task models (default: 100)", type=int, default=100)
    work.add_argument("--gpus_per_task", help="GPUs of a full predict task, fewer for fewer models: one GPU per jobs_per_task/gpus_per_task models (default: 1, Jean Zay: 4)", type=int, default=None)
    work.add_argument("--msa_gpus", help="GPUs of the MSA job (default: 3)", type=int, default=3)
    work.add_argument("--split_setup", help="Jean Zay: submit prepare, afdb, msa and features as separate jobs instead of one setup job on prepost", action="store_true")
    work.add_argument("--setup_partition", help="Jean Zay: partition of the setup job (default: prepost)", default="prepost")
    work.add_argument("--setup_gpus", help="Jean Zay: GPUs of the setup job, when DeepLoc, SignalP or MSAs have to run (default: 1)", type=int, default=1)
    work.add_argument("--setup_cpus", help="Jean Zay: CPUs of the setup job (default: 16)", type=int, default=16)
    work.add_argument("--setup_time", help="Jean Zay: time limit of the setup job (default: 10:00:00)", default="10:00:00")
    work.add_argument("--start_at", help="First stage to submit, to resume a run (default: fetch/prepare)", choices=STAGES[:6], default="fetch")
    work.add_argument("--start_round", help="First round to submit with --start_at predict (default: 0)", type=int, default=0)
    for stage in STAGES :
        work.add_argument(f"--{stage}_cpus", type=int, default=None, help=f"CPUs of the {stage} jobs (default: {DEFAULT_RESOURCES[stage]['cpus']}, Jean Zay: {JEAN_ZAY_RESOURCES[stage]['cpus'] or 'CPUs per GPU x GPUs'})")
        work.add_argument(f"--{stage}_mem", default=None, help=f"memory of the {stage} jobs, 'none' to omit (default: {DEFAULT_RESOURCES[stage]['mem']}, Jean Zay: omitted)")
        work.add_argument(f"--{stage}_time", default=None, help=f"time limit of the {stage} jobs (default: {DEFAULT_RESOURCES[stage]['time']}, Jean Zay: {JEAN_ZAY_RESOURCES[stage]['time']})")
        work.add_argument(f"--{stage}_partition", default=None, help=f"partition of the {stage} jobs" + (f" (Jean Zay default: {JEAN_ZAY_PARTITIONS[stage]})" if stage in JEAN_ZAY_PARTITIONS else ""))

    hint = parser.add_argument_group("HInt")
    hint.add_argument("--gpu_vram", help="VRAM (GiB) of the GPUs used for modeling, sets job lists and packing (default: read on the node, Jean Zay: 80 for a100/h100, 32 for v100)", type=float, default=None)
    hint.add_argument("--multi_job_per_gpu", help="Allow multiple jobs on the same GPU if VRAM allows it (default: True)", default="True")
    hint.add_argument("--multi_scoring", help="Score all models of each interaction (default: False)", default="False")
    hint.add_argument("--send_to_database", help="Publish models to QSProteome (https://qsproteome.org/). All, Negative or None", default="None")
    hint.add_argument("--max_retries", help="Retries of a crashed AlphaFold job (default: 1)", type=int, default=1)
    hint.add_argument("--no_afdb", help="Do not download MSAs from the AlphaFold database", action="store_true")
    hint.add_argument("--local_search_threshold", help="Use local MMseqs2-GPU (colabfold_search) from this number of MSAs, below use the ColabFold server (default: 1)", type=int, default=1)
    hint.add_argument("--jax_cache_dir", help="Shared JAX compilation cache (default: next to Path_Database)", default=None)


def main () :
    parser = argparse.ArgumentParser(prog="HInt-slurm", description="Run HInt on a Slurm cluster.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_submit_arguments(subparsers.add_parser("submit", help="Submit the whole pipeline"))
    subparsers.add_parser("status", help="Show the progress of the pipeline")
    for stage, help_text in [("fetch", "UniProt sequences download (run by submit on the login node)"), ("prepare", "(job) sequences, DeepLoc, SignalP and filters"),
                             ("setup", "(job) prepare, afdb and features in one job"),
                             ("msa_features", "(job) msa then features in one job"),
                             ("afdb", "(job) MSAs of the AlphaFold database"), ("msa", "(job) MSA generation")] :
        subparsers.add_parser(stage, help=help_text).add_argument("--new_log", help="Reset log_file/HInt.log (first job of a run)", action="store_true")
    features = subparsers.add_parser("features", help="(array job) pkl files")
    features.add_argument("--new_log", help="Reset log_file/HInt.log (first job of a run)", action="store_true")
    features.add_argument("--task", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    predict = subparsers.add_parser("predict", help="(array job) AlphaFold models of one round")
    predict.add_argument("--round", type=int, required=True)
    predict.add_argument("--task", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    score = subparsers.add_parser("score", help="(job) scoring of one round and final results")
    score.add_argument("--round", type=int, required=True)
    args = parser.parse_args()
    if getattr(args, "new_log", False) :
        open(HINT_LOG, "w").close()

    if args.command == "features" and not 0 <= args.task < load_config()["features_tasks"] :
        raise ValueError(f"--task must be between 0 and {load_config()['features_tasks'] - 1}")

    try :
        {"submit" : stage_submit,
         "status" : stage_status,
         "fetch" : stage_fetch,
         "prepare" : stage_prepare,
         "setup" : stage_setup,
         "afdb" : stage_afdb,
         "msa" : stage_msa,
         "msa_features" : stage_msa_features,
         "features" : stage_features,
         "predict" : stage_predict,
         "score" : stage_score}[args.command](args)
    except Exception :
        logger.exception(f"HInt-slurm {args.command} failed") # traceback in log_file/HInt.log
        raise


if __name__ == "__main__" :
    main()
