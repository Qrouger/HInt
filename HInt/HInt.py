""" Main file of HInt

    Author: Quentin Rouger
"""

from Utils_HInt import *
from File_proteins import *

import sys
import logging
import datetime
import argparse


log_filename = "./HInt.log"
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    GPU =[gpu for gpu in args.gpu.split(",")]
    Informations_dict = Define_informations()
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"])
    HInt_object.get_result_dict()
    need_msa, need_pkl = HInt_object.find_proteins_sequence(Informations_dict["Path_Pickle_Feature"]) #and real name of proteins
    HInt_object.create_fasta_file(need_msa, need_pkl)
    if len(need_msa) > 0 :
        remove_SP(HInt_object,Informations_dict, need_msa, "msa")
    if len(need_pkl) > 0 :
        remove_SP(HInt_object,Informations_dict, need_pkl, "pkl")
    filter_signalP(HInt_object,Informations_dict)
    if Informations_dict["DeepLoc"] != "None" :
        filter_deeploc(HInt_object, Informations_dict["DeepLoc"])
    create_feature(HInt_object,Informations_dict,GPU)
    print(str(datetime.datetime.now())+" Generation of MSA depth figures")
    Make_all_MSA_coverage(HInt_object,Informations_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    print(str(datetime.datetime.now()) + " Remove baits from prey list")
    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]]) #remove baits from prey list
    if int(Informations_dict["Homo-oligomer"]) > 1 : #select proteins who can create an homo-oligomer
        Use_RF2_PPI(HInt_object, Informations_dict, "RF2_homo_int", GPU, Informations_dict["Regions"]) #set better interactions all_vs_bait #Maybe remove
    if Informations_dict["Interact_with"] != "" :
        Generate_scripts(HInt_object, Informations_dict, "APD_PPI_int", GPU) #generate bait_vs_prey with new preys
        Generate_3D_model(Informations_dict, "APD_PPI_int", GPU)
        Score_interaction_APD(HInt_object, Informations_dict, "APD_PPI_int")
    if int(Informations_dict["Homo-oligomer"]) > 1 : #select proteins who can create an homo-oligomer
        Generate_scripts(HInt_object, Informations_dict, "APD_homo_int", GPU)
        Generate_3D_model(Informations_dict, "APD_homo_int", GPU)
        Score_interaction_APD(HInt_object, Informations_dict, "APD_homo_int")
    Resume_file(HInt_object)

    
