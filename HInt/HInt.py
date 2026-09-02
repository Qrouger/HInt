"""
Main entry point of HInt

Author: Quentin Rouger

This script orchestrates the full HInt pipeline, from input parsing to final results generation.
"""

import sys
import logging
import argparse
import os
import threading

# ------------------------------------------------------------------
# Logging configuration
# ------------------------------------------------------------------

log_filename = "./log_file/HInt.log"

# Create log directory if it does not exist
if not os.path.exists("log_file") :
    os.system("mkdir log_file")

from .Utils_HInt import *
from .File_proteins import *
from .Scoring_HInt import Score_interaction, Resume_file, Create_figures

# Reset existing handlers to avoid duplicated logs
logging.getLogger().handlers.clear()
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler = logging.FileHandler(log_filename, mode='w')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(stream_handler)


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
    parser.add_argument("--gpu", help="Comma-separated list of GPUs available for computation (default: 0)", required=False, default="0")

    N_CPU = multiprocessing.cpu_count()
    default_cpu = N_CPU // 2  # by default, use half of the available CPUs
    parser.add_argument("--cpu", help="Number of CPUs available for computation (default: half of the CPUs)", required=False, default=default_cpu, type=int)
    parser.add_argument("--multi_job_per_gpu", help="Allow multiple jobs to run on the same GPU if VRAM allows it (default: True)", required=False, default="True")
    parser.add_argument("--multi_scoring", help="Score all models of each interaction (default: False) and make a mean score.", required=False, default="False")

# ------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------

def main() :

    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()

    GPU = [gpu for gpu in args.gpu.split(",")]
    CPU = args.cpu
    if CPU > multiprocessing.cpu_count() :
        raise ValueError(f"Number of CPUs specified ({CPU}) exceeds the number of available CPUs ({multiprocessing.cpu_count()}).")
    multi_job_per_gpu = args.multi_job_per_gpu
    if multi_job_per_gpu not in ["True", "False"] :
        raise ValueError("Invalid value for --multi_job_per_gpu. Need True or False.")
    multi_scoring = args.multi_scoring
    if multi_scoring not in ["True", "False"] :
        raise ValueError("Invalid value for --multi_scoring. Need True or False.")
    Informations_dict = Define_informations()
    
    # --------------------------------------------------------------
    # Initialize protein container
    # --------------------------------------------------------------
    
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"], Informations_dict["Interact_with"], Informations_dict["AlphaFold"])

    logger.info("GPUs set to: %s", GPU)
    logger.info("Number of CPUs set to: %s", CPU)

    time_dict = HInt_object.get_time_dict() 
    time_dict["Summarize_info"] = [Informations_dict["Interact_with"], Informations_dict["DeepLoc"], Informations_dict["Signal_peptide"], Informations_dict["Min_protein_length"], Informations_dict["Max_protein_length"], Informations_dict["Homo-oligomer"],str(len(HInt_object.get_possible_prey())),Informations_dict["Organism"]]


    
    for bait in Informations_dict["Interact_with"] : # Check that all bait proteins exist in the protein list
        if bait not in HInt_object.get_proteins() and bait != "" :
            raise Exception(f"Bait {bait} not found in the protein list {Informations_dict['Path_Uniprot_ID']}")

    # --------------------------------------------------------------
    # Checkpointing: determine which features need to be computed
    # --------------------------------------------------------------

    start_sequence = time.time()
    # need_msa also includes proteins that only require signal peptide information
    need_msa, need_pkl, need_DeepLoc = HInt_object.check_save_dict(Informations_dict["Path_Pickle_Feature"])
    time_sequence = format_time(time.time() - start_sequence)
    time_dict["Sequences"] = [len(HInt_object.get_proteins()),time_sequence]

    # Remove bait proteins from the prey list
    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]])
    # --------------------------------------------------------------
    # Length-based filtering
    # --------------------------------------------------------------

    # Filter proteins based on sequence length
    # (default: remove proteins shorter than 20 AA)
    need_msa, need_pkl, need_DeepLoc = filter_length(HInt_object, Informations_dict, need_msa, need_pkl, need_DeepLoc)

    # --------------------------------------------------------------
    # DeepLoc filtering
    # --------------------------------------------------------------

    start_DeepLoc = time.time()
    if Informations_dict["Organism"] == "None" :
        HInt_object.set_proteins_sequence_no_SP(HInt_object.get_proteins_sequence_SP()) #don't remove signal peptide
        need_DeepLoc = []


    if len(need_DeepLoc) > 0 : # Run DeepLoc only for proteins without localization information
        run_deeploc(HInt_object, Informations_dict["Organism"], need_DeepLoc, GPU)
    time_DeepLoc = format_time(time.time() - start_DeepLoc)
    time_dict["Deeploc"] = [len(need_DeepLoc), time_DeepLoc]


    if Informations_dict["DeepLoc"].split(",") != ["None"] : # Apply DeepLoc-based filtering if enabled
        need_msa, need_pkl = filter_deeploc(HInt_object, Informations_dict, need_msa, need_pkl)

    HInt_object.Make_save_dict()  # Save DeepLoc results

    # --------------------------------------------------------------
    # Signal peptide processing
    # --------------------------------------------------------------

    start_SP = time.time()
    need_SP = list()
    for protein in need_msa :
        if protein not in HInt_object.get_proteins_sequence_no_SP().keys() or protein not in HInt_object.get_prot_SP().keys() : #protein need MSA but can already have sequence without SP
            need_SP.append(protein)


    if Informations_dict["Organism"] == "None" : # Run SignalP only if the organism is specified
        need_SP = []
    if len(need_SP) > 0 : # Run SignalP for proteins without signal peptide annotation
        need_msa = run_SP(HInt_object, Informations_dict, need_SP, need_msa)

    time_SP = format_time(time.time() - start_SP)
    time_dict["SignalP"] = [len(need_SP), time_SP]

    need_msa = check_exist_MSA(HInt_object, Informations_dict, need_msa) # Check for existing MSA files after SignalP processing
    HInt_object.Make_save_dict() # Save sequences without signal peptides
    
    for bait in Informations_dict["Interact_with"] : # Adjust bait protein lengths if specific regions are defined
        if Informations_dict["Regions"][bait] != "0-0" :
            dict_length = HInt_object.get_length_prot()
            start = int(Informations_dict["Regions"][bait].split("-")[0])
            end = int(Informations_dict["Regions"][bait].split("-")[1])
            dict_length[bait] = end - start + 1
            HInt_object.set_length_prot(dict_length)

    # Filter proteins based on signal peptide criteria
    need_msa, need_pkl = filter_signalP(HInt_object, Informations_dict, need_msa, need_pkl)

    HInt_object.Make_save_dict()

    if os.path.exists("log_file/shallow_MSA.txt") :
        os.remove("log_file/shallow_MSA.txt")

    # --------------------------------------------------------------
    # Feature generation (MSA + pickle files)
    # --------------------------------------------------------------


    #Create batch just for loop on
    job_with_vram_length = []
    if Informations_dict["Interact_with"] != [''] or Informations_dict["Interact_with"] == ['']:
        for bait in Informations_dict["Multimer_bait"] :
            job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "PPI_int", bait)
            First_batch = Generate_first_batch(job_with_vram_length, GPU, multi_job_per_gpu) #correspond to the first batch of proteins to process, based on the number of available GPUs and CPUs
            if First_batch != None :
                First_batch[0].extend([b for b in bait.split(",")])
            need_msa, need_pkl, first_need_pkl = First_need_sorted(First_batch, need_msa, need_pkl)
            start_MSA = time.time()
            if first_need_pkl != [] or need_msa != [] :
                need_pkl.extend(create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, first_need_pkl))#generate MSA for all and pkl only for first batch
            time_MSA = format_time(time.time() - start_MSA)
            time_dict["MSA"] = [len(need_msa), time_MSA]
            logger.info("Generating MSA depth figures")
            Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"], []) # Just for bait prot
            HInt_object.Make_save_dict()
            prio_list_MSA = []
            if First_batch != None :
                prio_list_MSA = prioritize_by_vram_fit(job_with_vram_length, First_batch, GPU)
            batch_MSA = []
            if prio_list_MSA != [] :
                for batch in prio_list_MSA :
                    new_prio_list = []
                    for prot in batch :
                        if prot in need_pkl :
                            new_prio_list.append(prot)
                    if new_prio_list != [] :
                        batch_MSA.append(new_prio_list)
            start_PPI = time.time()
            if HInt_object.get_compounds() != {} :
                gpu_thread = threading.Thread(target=Generate_3D_model, args=(HInt_object, CPU, multi_scoring, Informations_dict, "Compounds", job_with_vram_length, GPU, multi_job_per_gpu))
                gpu_thread.start()
                gpu_thread.join()
            else :
                gpu_thread = threading.Thread(target=Generate_3D_model, args=(HInt_object, CPU, multi_scoring, Informations_dict, "PPI_int", job_with_vram_length, GPU, multi_job_per_gpu))
                gpu_thread.start()
                for batch in batch_MSA :
                    create_feature(HInt_object, Informations_dict, GPU, CPU, [], batch)
                Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"], HInt_object.get_possible_prey())
                gpu_thread.join()
            time_PPI = format_time(time.time() - start_PPI)
            time_dict["PPI"] = [len(job_with_vram_length), time_PPI]
            Score_interaction(HInt_object, Informations_dict, CPU, "PPI_int", "", multi_scoring, bait) #if score is not done in Generate_3D_model, do it here
    HInt_object.Make_save_dict()



    # --------------------------------------------------------------
    # Homo-oligomer modeling
    # --------------------------------------------------------------

    if int(Informations_dict["Homo-oligomer"]) > 1 :
        job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "homo_int", "")
        if Informations_dict["Interact_with"] == [''] :
            First_batch = Generate_first_batch(job_with_vram_length, GPU, multi_job_per_gpu)
            need_msa, need_pkl, first_need_pkl = First_need_sorted(First_batch, need_msa, need_pkl)
            if first_need_pkl != [] and need_msa != [] :
                need_pkl.extend(create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, first_need_pkl))
            logger.info("Generating MSA depth figures")
            Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"], []) # Just for bait prot
            prio_list_MSA = []
            if First_batch != None :
                prio_list_MSA = prioritize_by_vram_fit(job_with_vram_length, First_batch, GPU)
            batch_MSA = []
            if prio_list_MSA != [] :
                for batch in prio_list_MSA :
                    new_prio_list = []
                    for prot in batch :
                        if prot in need_pkl :
                            new_prio_list.append(prot)
                    if new_prio_list != [] :
                        batch_MSA.append(new_prio_list)
            gpu_thread = threading.Thread(target=Generate_3D_model, args=(HInt_object, CPU, multi_scoring, Informations_dict, "homo_int", job_with_vram_length, GPU, multi_job_per_gpu))
            gpu_thread.start()
            for batch in batch_MSA :
                create_feature(HInt_object, Informations_dict, GPU, CPU, [], batch)
            Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"], HInt_object.get_possible_prey())
            gpu_thread.join()
            Score_interaction(HInt_object, Informations_dict, CPU, "homo_int", "", multi_scoring) 
        else :
            Generate_3D_model(HInt_object, CPU, multi_scoring, Informations_dict, "homo_int", job_with_vram_length, GPU, multi_job_per_gpu)
            Score_interaction(HInt_object, Informations_dict, CPU, "homo_int", "", multi_scoring) 
            HInt_object.Make_save_dict()



    # --------------------------------------------------------------
    # Final summary and figures
    # --------------------------------------------------------------


    sorted_protein = Resume_file(HInt_object, Informations_dict)

    if Informations_dict["Interact_with"] != [''] :

        Create_figures(HInt_object, Informations_dict, Informations_dict["AlphaFold"], sorted_protein, CPU)


    line_report = f"===== SUMMARY INFO =====\nBait : {time_dict['Summarize_info'][0]}\nCellular localization : {time_dict['Summarize_info'][1]}\nPeptide signal : {time_dict['Summarize_info'][2]}\nMinimum protein length : {time_dict['Summarize_info'][3]}\nMaximum protein length : {time_dict['Summarize_info'][4]}\nHomo-oligomer state : {time_dict['Summarize_info'][5]}\nNumber of preys : {time_dict['Summarize_info'][6]}\nOrganism : {time_dict['Summarize_info'][7]}\nGPU : {GPU}\nCPU : {CPU}\n\n===== TIMINGS =====\nStep\tNumber of proteins\tTime\n"
    for step in time_dict.keys() :
        if step != "Summarize_info" :
            line_report += f"{step}\t{time_dict[step][0]}\t{time_dict[step][1]}\n"
    with open("./log_file/HInt_report.txt", "w") as f :
        f.write(line_report)
    logger.info("All steps done")


def format_time(seconds) :
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{seconds:06.3f}"
