"""
Main entry point of HInt

Author: Quentin Rouger

This script orchestrates the full HInt pipeline, from input parsing to final results generation.
"""

import sys
import logging
import argparse
import os
import time

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
    parser.add_argument("--multi_scoring", help="Score all models of each interactions (default: False)", required=False, default="False")

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

    time_dict = HInt_object.get_time_dict() 
    time_dict["Summarize_info"] = [Informations_dict["Interact_with"], Informations_dict["DeepLoc"], Informations_dict["Signal_peptide"], Informations_dict["Min_protein_length"], Informations_dict["Max_protein_length"], Informations_dict["Homo-oligomer"],str(len(HInt_object.get_possible_prey())),Informations_dict["Organism"]]

    logger.info("GPUs set to: %s", GPU)
    logger.info("Number of CPUs set to: %s", CPU)


    
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

    need_SP = list()
    for protein in need_msa :
        if protein not in HInt_object.get_proteins_sequence_no_SP().keys() or protein not in HInt_object.get_prot_SP().keys() : #protein need MSA but can already have sequence without SP
            need_SP.append(protein)

    start_SP = time.time()
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


    # --------------------------------------------------------------
    # Feature generation (MSA + pickle files)
    # --------------------------------------------------------------

    start_MSA = time.time()
    if (len(need_msa) > 0 or len(need_pkl) > 0) :#and Informations_dict["Interact_with"] != [''] :
        create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, need_pkl)
    time_MSA = format_time(time.time() - start_MSA)
    time_dict["MSA"] = [len(need_msa), time_MSA]

    HInt_object.Make_save_dict()  # Save SignalP and feature results


    # --------------------------------------------------------------
    # MSA depth analysis
    # --------------------------------------------------------------

    logger.info("Generating MSA depth figures")
    Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"])


    # --------------------------------------------------------------
    # Protein–protein interaction modeling
    # --------------------------------------------------------------

    if Informations_dict["Interact_with"] != [''] :
        for bait in Informations_dict["Multimer_bait"] :
            if HInt_object.get_compounds() != {} :
                job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "Compounds", bait)
                start_PPI = time.time()
                Generate_3D_model(Informations_dict, "Compounds", job_with_vram_length, GPU, multi_job_per_gpu, HInt_object.get_proteins_sequence_no_SP(), HInt_object.get_compounds())
                time_PPI = format_time(time.time() - start_PPI)
                time_dict["PPI"] = [len(job_with_vram_length), time_PPI]
                start_Scoring = time.time()
                nbr_new_score = Score_interaction(HInt_object, Informations_dict, CPU, "Compounds", bait, multi_scoring)
                time_Scoring = format_time(time.time() - start_Scoring)
                time_dict["Scoring"] = [nbr_new_score, time_Scoring]
            else :
                job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "PPI_int", bait)
                start_PPI = time.time()
                Generate_3D_model(Informations_dict, "PPI_int", job_with_vram_length, GPU, multi_job_per_gpu)
                time_PPI = format_time(time.time() - start_PPI)
                time_dict["PPI"] = [len(job_with_vram_length), time_PPI]
                start_Scoring = time.time()
                nbr_new_score = Score_interaction(HInt_object, Informations_dict, CPU, "PPI_int", bait, multi_scoring)
                time_Scoring = format_time(time.time() - start_Scoring)
                time_dict["Scoring_PPI"] = [nbr_new_score, time_Scoring]
            HInt_object.Make_save_dict()  # Save interaction scores


    # --------------------------------------------------------------
    # Homo-oligomer modeling
    # --------------------------------------------------------------


    if int(Informations_dict["Homo-oligomer"]) > 1 :
        job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "homo_int", "")
        start_homo = time.time()
        Generate_3D_model(Informations_dict, "homo_int", job_with_vram_length, GPU, multi_job_per_gpu)
        time_homo = format_time(time.time() - start_homo)
        time_dict["Homo"] = [len(job_with_vram_length), time_homo]
        start_Scoring = time.time()
        nbr_new_score = Score_interaction(HInt_object, Informations_dict, CPU, "homo_int")
        time_Scoring = format_time(time.time() - start_Scoring)
        time_dict["Scoring_Homo"] = [nbr_new_score, time_Scoring]
        HInt_object.Make_save_dict()


    # --------------------------------------------------------------
    # Final summary and figures
    # --------------------------------------------------------------


    sorted_protein = Resume_file(HInt_object, Informations_dict)

    if Informations_dict["Interact_with"] != [''] :
        start_Figure = time.time()
        nbr_fig = Create_figures(HInt_object, Informations_dict, Informations_dict["AlphaFold"], sorted_protein, CPU)
        time_Figure = format_time(time.time() - start_Figure)
        time_dict["Figures"] = [len(nbr_fig), time_Figure]

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
