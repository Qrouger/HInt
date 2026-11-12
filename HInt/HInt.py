"""
Main entry point of HInt

Author: Quentin Rouger

This script orchestrates the full HInt pipeline, from input parsing to final results generation.
"""

import sys
import os
import time
import logging
import argparse
import threading
import multiprocessing

# ------------------------------------------------------------------
# Logging configuration
# ------------------------------------------------------------------

log_filename = "./log_file/HInt.log"
os.makedirs("log_file", exist_ok=True)  # must exist before importing HInt modules

from .Utils_HInt import *
from .File_proteins import *
from .Scoring_HInt import Score_interaction, Resume_file, Create_figures

# force=True removes existing handlers to avoid duplicated logs, mode "w" resets the log file
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(log_filename, mode="w"), logging.StreamHandler(sys.stdout)],
    force=True,
)
logger = logging.getLogger()


# ------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------

def add_arguments(parser) :
    """
    Define command-line arguments.

    Parameters
    ----------
    parser : argparse.ArgumentParser
    """
    parser.add_argument("--gpu", default="0", help="Comma-separated list of GPUs available for computation (default: 0)")
    parser.add_argument("--cpu", type=int, default=multiprocessing.cpu_count() // 2, help="Number of CPUs available for computation (default: half of the CPUs)")
    parser.add_argument("--multi_job_per_gpu", default="True", choices=["True", "False"], help="Allow multiple jobs to run on the same GPU if VRAM allows it (default: True)")
    parser.add_argument("--multi_scoring", default="False", choices=["True", "False"], help="Score all models of each interaction (default: False) and make a mean score.")
    parser.add_argument("--send_to_database", default="None", choices=["All", "Negative", "None"], help="Publish models to the public QSProteome database (https://qsproteome.org/). All, Negative or None")


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def format_time(seconds) :
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{seconds:06.3f}"


def get_priority_pkl(job_with_vram_length, First_batch, GPU, mode, need_pkl) :
    """One flat list in priority order: proteins that still need a pkl, without duplicates."""
    if First_batch is None :
        return []
    need_pkl_set = set(need_pkl)
    prio_list_MSA = prioritize_by_vram_fit(job_with_vram_length, First_batch, GPU, mode)
    return list(dict.fromkeys(prot for batch in prio_list_MSA for prot in batch if prot in need_pkl_set))


# ------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------

def main() :

    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()

    GPU = args.gpu.split(",")
    CPU = args.cpu
    if CPU > multiprocessing.cpu_count() :
        raise ValueError(f"Number of CPUs specified ({CPU}) exceeds the number of available CPUs ({multiprocessing.cpu_count()}).")
    multi_job_per_gpu = args.multi_job_per_gpu
    multi_scoring = args.multi_scoring
    send_to_database = args.send_to_database
    Informations_dict = Define_informations()
    baits = Informations_dict["Interact_with"]

    # --------------------------------------------------------------
    # Initialize protein container
    # --------------------------------------------------------------

    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"], baits, Informations_dict["AlphaFold"], Informations_dict["Custom_interactions"])

    logger.info(f"GPUs set to: {GPU}")
    logger.info(f"Multi-job per GPU set to: {multi_job_per_gpu}")
    logger.info(f"Number of CPUs set to: {CPU}")
    logger.info(f"AlphaFold version set to: {Informations_dict['AlphaFold']}")
    logger.info(f"Send to database option set to: {send_to_database}")

    time_dict = HInt_object.get_time_dict()
    time_dict["Summarize_info"] = [baits, Informations_dict["DeepLoc"], Informations_dict["Signal_peptide"],
                                   Informations_dict["Min_protein_length"], Informations_dict["Max_protein_length"],
                                   Informations_dict["Homo-oligomer"], str(len(HInt_object.get_possible_prey())),
                                   Informations_dict["Organism"]]

    for bait in baits : # Check that all bait proteins exist in the protein list
        if bait != "" and bait not in HInt_object.get_proteins() :
            raise Exception(f"Bait {bait} not found in the protein list {Informations_dict['Path_Uniprot_ID']}")

    def model_and_features(mode, job_with_vram_length, prio_list=None, preys=None) :
        """
        Run 3D modeling in a GPU thread. If prio_list is given, the remaining features (pkl)
        are computed at the same time, then the MSA depth figures of the preys are made.
        """
        gpu_thread = threading.Thread(target=Generate_3D_model, args=(HInt_object, CPU, multi_scoring, send_to_database, Informations_dict, mode, job_with_vram_length, GPU, multi_job_per_gpu))
        gpu_thread.start()
        if prio_list is not None :
            create_feature(HInt_object, Informations_dict, GPU, CPU, [], prio_list)
            Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], baits, preys)
        gpu_thread.join()

    # --------------------------------------------------------------
    # Checkpointing: determine which features need to be computed
    # --------------------------------------------------------------

    start = time.time()
    # need_msa also includes proteins that only require signal peptide information
    need_msa, need_pkl, need_DeepLoc = HInt_object.check_save_dict(Informations_dict["Path_Pickle_Feature"])
    time_dict["Sequences"] = [len(HInt_object.get_proteins()), format_time(time.time() - start)]

    # Remove bait proteins from the prey list
    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in baits])

    # --------------------------------------------------------------
    # Length-based filtering (default: remove proteins shorter than 20 AA)
    # --------------------------------------------------------------

    if Informations_dict["Min_protein_length"] != "0" or Informations_dict["Max_protein_length"] != "" :
        need_msa, need_pkl, need_DeepLoc = filter_length(HInt_object, Informations_dict, need_msa, need_pkl, need_DeepLoc)

    # --------------------------------------------------------------
    # DeepLoc filtering
    # --------------------------------------------------------------

    start = time.time()
    if Informations_dict["Organism"] == "None" :
        HInt_object.set_proteins_sequence_no_SP(HInt_object.get_proteins_sequence_SP()) # don't remove signal peptide
        need_DeepLoc = []
    if need_DeepLoc and Informations_dict["Custom_interactions"] == "" : # Run DeepLoc only for proteins without localization information
        run_deeploc(HInt_object, Informations_dict["Organism"], need_DeepLoc, GPU)
    time_dict["Deeploc"] = [len(need_DeepLoc), format_time(time.time() - start)]

    if Informations_dict["DeepLoc"].split(",") != ["None"] : # Apply DeepLoc-based filtering if enabled
        need_msa, need_pkl = filter_deeploc(HInt_object, Informations_dict, need_msa, need_pkl)
    HInt_object.Make_save_dict()  # Save DeepLoc results

    # --------------------------------------------------------------
    # Signal peptide processing
    # --------------------------------------------------------------

    start = time.time()
    need_SP = []
    if Informations_dict["Organism"] != "None" : # Run SignalP only if the organism is specified
        # protein needs MSA but can already have a sequence without SP
        need_SP = [protein for protein in need_msa if protein not in HInt_object.get_proteins_sequence_no_SP() or protein not in HInt_object.get_prot_SP()]
    if need_SP : # Run SignalP for proteins without signal peptide annotation
        need_msa = run_SP(HInt_object, Informations_dict, need_SP, need_msa)
    time_dict["SignalP"] = [len(need_SP), format_time(time.time() - start)]

    need_msa = check_exist_MSA(HInt_object, Informations_dict, need_msa) # Check for existing MSA files after SignalP processing
    HInt_object.Make_save_dict() # Save sequences without signal peptides

    # --------------------------------------------------------------
    # Regions of interest (custom interactions)
    # --------------------------------------------------------------

    if Informations_dict["Custom_interactions"] != "" :
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

    # Filter proteins based on signal peptide criteria
    need_msa, need_pkl = filter_signalP(HInt_object, Informations_dict, need_msa, need_pkl)
    HInt_object.Make_save_dict()
    if os.path.exists("log_file/shallow_MSA.txt") :
        os.remove("log_file/shallow_MSA.txt")

    # --------------------------------------------------------------
    # One vs all protein-protein interaction modeling
    # --------------------------------------------------------------

    if baits != [''] :
        for bait in Informations_dict["Multimer_bait"] :
            job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "PPI_int", bait)
            # first batch of proteins to process, based on the number of available GPUs and CPUs
            First_batch = Generate_first_batch(job_with_vram_length, GPU, multi_job_per_gpu)
            if First_batch is not None :
                First_batch[0].extend(bait.split(","))
            need_msa, need_pkl, first_need_pkl = First_need_sorted(First_batch, need_msa, need_pkl)

            start = time.time()
            if first_need_pkl or need_msa : # generate MSA for all and pkl only for first batch
                need_pkl.extend(create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, first_need_pkl))
            time_dict["MSA"] = [len(need_msa), format_time(time.time() - start)]

            logger.info("Generating MSA depth figures")
            Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], baits, []) # Just for bait prot
            HInt_object.Make_save_dict()
            prio_list = get_priority_pkl(job_with_vram_length, First_batch, GPU, "PPI_int", need_pkl)

            start = time.time()
            if HInt_object.get_compounds() != {} :
                model_and_features("Compounds", job_with_vram_length)
            else :
                model_and_features("PPI_int", job_with_vram_length, prio_list, HInt_object.get_possible_prey())
            time_dict["PPI"] = [len(job_with_vram_length), format_time(time.time() - start)]
            # if score is not done in Generate_3D_model, do it here
            Score_interaction(HInt_object, Informations_dict, CPU, "PPI_int", send_to_database, "", multi_scoring, bait)
    HInt_object.Make_save_dict()

    # --------------------------------------------------------------
    # Homo-oligomer modeling
    # --------------------------------------------------------------

    if int(Informations_dict["Homo-oligomer"]) > 1 :
        job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "homo_int", "")
        if baits == [''] :
            First_batch = Generate_first_batch(job_with_vram_length, GPU, multi_job_per_gpu)
            need_msa, need_pkl, first_need_pkl = First_need_sorted(First_batch, need_msa, need_pkl)
            if first_need_pkl and need_msa :
                need_pkl.extend(create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, first_need_pkl))
            logger.info("Generating MSA depth figures")
            Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], baits, []) # Just for bait prot
            prio_list = get_priority_pkl(job_with_vram_length, First_batch, GPU, "homo_int", need_pkl)
            model_and_features("homo_int", job_with_vram_length, prio_list, HInt_object.get_possible_prey())
        else :
            Generate_3D_model(HInt_object, CPU, multi_scoring, send_to_database, Informations_dict, "homo_int", job_with_vram_length, GPU, multi_job_per_gpu)
        Score_interaction(HInt_object, Informations_dict, CPU, "homo_int", send_to_database, "", multi_scoring)
        HInt_object.Make_save_dict()

    # --------------------------------------------------------------
    # Custom interactions modeling
    # --------------------------------------------------------------

    if Informations_dict["Custom_interactions"] != '' :
        job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "custom", "")
        First_batch = Generate_first_batch(job_with_vram_length, GPU, multi_job_per_gpu, True)
        need_msa, need_pkl, first_need_pkl = First_need_sorted(First_batch, need_msa, need_pkl)
        if first_need_pkl and need_msa :
            need_pkl.extend(create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, first_need_pkl))
        logger.info("Generating MSA depth figures")
        HInt_object.Make_save_dict()
        prio_list = get_priority_pkl(job_with_vram_length, First_batch, GPU, "custom", need_pkl)
        model_and_features("custom", job_with_vram_length, prio_list, "")
        Score_interaction(HInt_object, Informations_dict, CPU, "custom", send_to_database, "", multi_scoring)
        HInt_object.Make_save_dict()

    # --------------------------------------------------------------
    # Final summary and figures
    # --------------------------------------------------------------

    sorted_protein = Resume_file(HInt_object, Informations_dict)
    if baits != [''] :
        Create_figures(HInt_object, Informations_dict, Informations_dict["AlphaFold"], sorted_protein, CPU)

    labels = ["Bait", "Cellular localization", "Peptide signal", "Minimum protein length",
              "Maximum protein length", "Homo-oligomer state", "Number of preys", "Organism"]
    lines = ["===== SUMMARY INFO ====="]
    lines += [f"{label} : {value}" for label, value in zip(labels, time_dict["Summarize_info"])]
    lines += [f"GPU : {GPU}", f"CPU : {CPU}", "", "===== TIMINGS =====", "Step\tNumber of proteins\tTime"]
    lines += [f"{step}\t{values[0]}\t{values[1]}" for step, values in time_dict.items() if step != "Summarize_info"]
    with open("./log_file/HInt_report.txt", "w") as f :
        f.write("\n".join(lines) + "\n")
    logger.info("All steps done")
