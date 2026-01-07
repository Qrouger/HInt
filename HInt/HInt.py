""" Main file of HInt

    Author: Quentin Rouger
"""

from Utils_HInt import *
from File_proteins import *
from Scoring_HInt import Score_interaction, Resume_file, Create_figures

import sys
import logging
from datetime import datetime
import argparse


log_filename = "./log_file/HInt.log"
if os.path.exists("log_file") == False :
    os.system("mkdir log_file")
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

def add_arguments(parser) :
    parser.add_argument("--gpu", help = "list of GPUs available for work", required = False, default = "0")
    N_CPU = multiprocessing.cpu_count()
    default_cpu = N_CPU//2 #by default use half of cpu available
    parser.add_argument("--cpu", help = "Number of CPUs availabe for work", required = False, default = default_cpu, type = int)



if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    GPU = [gpu for gpu in args.gpu.split(",")]
    Informations_dict = Define_informations()
    logger.info("GPU set to : "+str(GPU))
    CPU = args.cpu
    logger.info("Number of CPU set to : "+str(CPU))

    #Create objects
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"])

    #Error if bait not in protein list
    for bait in Informations_dict["Interact_with"] :
        if bait not in HInt_object.get_proteins() and bait != "" :
            raise Exception(f"Bait {bait} not in the protein list of {Informations_dict['Path_Uniprot_ID']}")


    #Look at the Checkpoint
    need_msa, need_pkl, need_DeepLoc = HInt_object.check_save_dict(Informations_dict["Path_Pickle_Feature"], Informations_dict["Regions"]) #need MSA represnt alos protein who just need SP informations

    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]]) #remove baits from prey list

    #Filter with lenght
    need_msa, need_pkl, need_DeepLoc = filter_lenght(HInt_object, Informations_dict, need_msa, need_pkl, need_DeepLoc) #filter proteins based on lenght, by default > 30 aa (remove junk prot in proteome data), update need lists, use length with SP

    #Filter with DeepLoc
    if len(need_DeepLoc) > 0 :
        run_deeploc(HInt_object, Informations_dict["Organism"], need_DeepLoc, GPU) #run DeepLoc for new proteins and set in dict
    if Informations_dict["DeepLoc"].split(",") != ["None"] :
        need_msa, need_pkl = filter_deeploc(HInt_object, Informations_dict, need_msa, need_pkl)
    HInt_object.Make_save_dict() #Save DeepLoc results

    #Filter with SignalP
    if len(need_msa) > 0 :
        need_msa = run_SP(HInt_object, Informations_dict, need_msa) #run SignalP for new proteins and set in dict, if protein have already msa but no SP informations don't return it
    HInt_object.Make_save_dict() #Save sequence without signal peptide

    for bait in Informations_dict["Interact_with"] :
        if Informations_dict["Regions"][bait] != "0-0" : #modifying real lenght of bait
            dict_lenght = HInt_object.get_lenght_prot()
            start = int(Informations_dict["Regions"][bait].split("-")[0])
            end = int(Informations_dict["Regions"][bait].split("-")[1])
            dict_lenght[bait] = end - start + 1
            HInt_object.set_lenght_prot(dict_lenght)

    #if Informations_dict["Signal_peptide"] != "None" :
    need_msa, need_pkl = filter_signalP(HInt_object,Informations_dict, need_msa, need_pkl) #filter proteins with SignalP, if None just describe in dict

    #Create MSA features
    if len(need_msa) > 0 or len(need_pkl) > 0 : 
        create_feature(HInt_object, Informations_dict, GPU, CPU, need_msa, need_pkl) #run MSA and create pkl files for new proteins

    HInt_object.Make_save_dict() #Save SignalP results

    logger.info(str(datetime.now())+" Generation of MSA depth figures")
    Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    
    #print(str(datetime.now()) + " Remove baits from prey list")
    #HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]]) #remove baits from prey list

    if Informations_dict["Interact_with"] != [""] :
        for bait in Informations_dict["Multimer_bait"] :
            Generate_scripts(HInt_object, Informations_dict, "PPI_int", bait, GPU) #generate bait_vs_prey with new preys
            Generate_3D_model(Informations_dict, "PPI_int", GPU)
            Score_interaction(HInt_object, Informations_dict, CPU, "PPI_int", bait)
            HInt_object.Make_save_dict() #Save scores #Maybe make it for each prot score
    if int(Informations_dict["Homo-oligomer"]) > 1 : #select proteins who can create an homo-oligomer
        Generate_scripts(HInt_object, Informations_dict, "homo_int","", GPU)
        Generate_3D_model(Informations_dict, "homo_int", GPU)
        Score_interaction(HInt_object, Informations_dict, "homo_int")
    Resume_file(HInt_object, Informations_dict)
    Create_figures(HInt_object, Informations_dict, Informations_dict["AlphaFold"])
    
