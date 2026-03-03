"""
Main entry point of HInt

Author: Quentin Rouger

This script orchestrates the full HInt pipeline, from input parsing to final results generation.
"""

import sys
import logging
import argparse
import os

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
    parser.add_argument("--gpu", help="Comma-separated list of GPUs available for computation", required=False, default="0")

    N_CPU = multiprocessing.cpu_count()
    default_cpu = N_CPU // 2  # by default, use half of the available CPUs
    parser.add_argument("--cpu", help="Number of CPUs available for computation", required=False, default=default_cpu, type=int)
    parser.add_argument("--multi_job_per_gpu", help="Allow multiple jobs to run on the same GPU if VRAM allows it (default: True)", required=False, default="True")


# ------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------

if __name__ == "__main__" :

    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()

    GPU = [gpu for gpu in args.gpu.split(",")]
    CPU = args.cpu
    multi_job_per_gpu = args.multi_job_per_gpu
    if multi_job_per_gpu not in ["True", "False"] :
        raise ValueError("Invalid value for --multi_job_per_gpu.")
        
    Informations_dict = Define_informations()

    logger.info("GPUs set to: %s", GPU)
    logger.info("Number of CPUs set to: %s", CPU)

    # --------------------------------------------------------------
    # Initialize protein container
    # --------------------------------------------------------------

    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"])

    
    for bait in Informations_dict["Interact_with"] : # Check that all bait proteins exist in the protein list
        if bait not in HInt_object.get_proteins() and bait != "" :
            raise Exception(f"Bait {bait} not found in the protein list {Informations_dict['Path_Uniprot_ID']}")

    # --------------------------------------------------------------
    # Checkpointing: determine which features need to be computed
    # --------------------------------------------------------------


    # need_msa also includes proteins that only require signal peptide information
    need_msa, need_pkl, need_DeepLoc = HInt_object.check_save_dict(Informations_dict["Path_Pickle_Feature"])


    # Remove bait proteins from the prey list
    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]])

    # --------------------------------------------------------------
    # Length-based filtering
    # --------------------------------------------------------------

    # Filter proteins based on sequence length
    # (default: remove proteins shorter than 20 AA)
    need_msa, need_pkl, need_DeepLoc = filter_lenght(HInt_object, Informations_dict, need_msa, need_pkl, need_DeepLoc)

    # --------------------------------------------------------------
    # DeepLoc filtering
    # --------------------------------------------------------------

    if len(need_DeepLoc) > 0 : # Run DeepLoc only for proteins without localization information
        run_deeploc(HInt_object, Informations_dict["Organism"], need_DeepLoc, GPU)


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
    if len(need_SP) > 0 : # Run SignalP for proteins without signal peptide annotation
        need_msa = run_SP(HInt_object, Informations_dict, need_SP, need_msa)

    HInt_object.Make_save_dict() # Save sequences without signal peptides

    
    for bait in Informations_dict["Interact_with"] : # Adjust bait protein lengths if specific regions are defined
        if Informations_dict["Regions"][bait] != "0-0" :
            dict_lenght = HInt_object.get_lenght_prot()
            start = int(Informations_dict["Regions"][bait].split("-")[0])
            end = int(Informations_dict["Regions"][bait].split("-")[1])
            dict_lenght[bait] = end - start + 1
            HInt_object.set_lenght_prot(dict_lenght)


    # Filter proteins based on signal peptide criteria
    need_msa, need_pkl = filter_signalP(HInt_object, Informations_dict, need_msa, need_pkl)


    # --------------------------------------------------------------
    # Feature generation (MSA + pickle files)
    # --------------------------------------------------------------


    if len(need_msa) > 0 or len(need_pkl) > 0 :
        create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, need_pkl)

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
            job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "PPI_int", bait)
            Generate_3D_model(Informations_dict, "PPI_int", job_with_vram_length, GPU, multi_job_per_gpu)
            Score_interaction(HInt_object, Informations_dict, CPU, "PPI_int", bait)
            HInt_object.Make_save_dict()  # Save interaction scores


    # --------------------------------------------------------------
    # Homo-oligomer modeling
    # --------------------------------------------------------------


    if int(Informations_dict["Homo-oligomer"]) > 1 :
        job_with_vram_length = Generate_scripts(HInt_object, Informations_dict, "homo_int", "")
        Generate_3D_model(Informations_dict, "homo_int", job_with_vram_length, GPU, multi_job_per_gpu)
        Score_interaction(HInt_object, Informations_dict, CPU, "homo_int")


    # --------------------------------------------------------------
    # Final summary and figures
    # --------------------------------------------------------------


    sorted_protein = Resume_file(HInt_object, Informations_dict)
    if Informations_dict["Interact_with"] != [''] :
        Create_figures(HInt_object, Informations_dict, Informations_dict["AlphaFold"], sorted_protein, CPU)
