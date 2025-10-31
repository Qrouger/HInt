""" Main file of HInt

    Author: Quentin Rouger
"""

from Utils_HInt import *
from File_proteins import *

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



if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    GPU =[gpu for gpu in args.gpu.split(",")]
    Informations_dict = Define_informations()

    #Create objects and filtre proteins
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"])
    need_msa, need_pkl, need_DeepLoc = HInt_object.check_save_dict(Informations_dict["Path_Pickle_Feature"])
    HInt_object.Make_save_dict()
    if len(need_DeepLoc) > 0 :
        run_deeploc(HInt_object, Informations_dict["Organism"], need_DeepLoc, GPU) #run DeepLoc for new proteins and set in dict
    if Informations_dict["DeepLoc"].split(",") != ["None"] :
        need_msa, need_pkl = filter_deeploc(HInt_object, Informations_dict, need_msa, need_pkl)
    HInt_object.Make_save_dict()

    
    if len(need_msa) > 0 :
        run_SP(HInt_object,Informations_dict, need_msa) #run SignalP for new proteins and set in dict
    HInt_object.Make_save_dict()
    if Informations_dict["Signal_peptide"] != "None" :
        need_msa, need_pkl = filter_signalP(HInt_object,Informations_dict, need_msa, need_pkl) #filter proteins with SignalP
    if len(need_msa) > 0 or len(need_pkl) > 0 :
        create_feature(HInt_object, Informations_dict,GPU, need_msa, need_pkl) #run MSA and create pkl files for new proteins
    HInt_object.Make_save_dict()

    print(str(datetime.now())+" Generation of MSA depth figures")
    Make_all_MSA_coverage(HInt_object, Informations_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    
    print(str(datetime.now()) + " Remove baits from prey list")
    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]]) #remove baits from prey list

    #Predict PPI and score interactions
   # if int(Informations_dict["Homo-oligomer"]) > 1 : #select proteins who can create an homo-oligomer
        #Use_RF2_PPI(HInt_object, Informations_dict, "RF2_homo_int", GPU) #set better interactions all_vs_bait #Maybe remove
    if Informations_dict["Interact_with"] != [""] :
        for bait in Informations_dict["Interact_with"] :
            Generate_scripts(HInt_object, Informations_dict, "PPI_int", bait, GPU) #generate bait_vs_prey with new preys
            Generate_3D_model(Informations_dict, "PPI_int", GPU)
            Score_interaction_APD(HInt_object, Informations_dict, "PPI_int")
    if int(Informations_dict["Homo-oligomer"]) > 1 : #select proteins who can create an homo-oligomer
        Generate_scripts(HInt_object, Informations_dict, "homo_int","", GPU)
        Generate_3D_model(Informations_dict, "homo_int", GPU)
        Score_interaction_APD(HInt_object, Informations_dict, "homo_int")
    Resume_file(HInt_object, Informations_dict)

    #figures generations
    
