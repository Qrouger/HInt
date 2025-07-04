""" Main file of HInt

    Author: Quentin Rouger
"""

from Utils_HInt import *
from File_proteins import *

import sys
import logging
import datetime


log_filename = "./HInt.log"
logging.basicConfig(filename=log_filename, filemode="w", level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s",)
class Logger(object):
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, "a")  # Mode append pour conserver l'historique
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        self.terminal.flush()
        self.log.flush()
        
sys.stdout = Logger(log_filename)
sys.stderr = Logger(log_filename)


GPU = ["0","1"]
#enlever new_pickle ?
if __name__ == "__main__":
    Informations_dict = Define_informations()
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"])
    HInt_object.get_result_dict()
    need_msa = HInt_object.find_proteins_sequence(Informations_dict["Path_Pickle_Feature"]) #and real name of proteins
    HInt_object.create_fasta_file(need_msa)
    if len(need_msa) > 0 :
        remove_SP(HInt_object,Informations_dict, need_msa)
    filter_signalP(HInt_object,Informations_dict)
    create_feature(HInt_object,Informations_dict["Path_AlphaFold_Data"],Informations_dict["Path_Pickle_Feature"])
    print(str(datetime.datetime.now())+" Generation of MSA depth figures")
    Make_all_MSA_coverage(HInt_object,Informations_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    print(str(datetime.datetime.now()) + " Remove baits from prey list")
    HInt_object.set_possible_prey([protein for protein in HInt_object.get_possible_prey() if protein not in Informations_dict["Interact_with"]]) #remove baits from prey list
    if int(Informations_dict["Homo-oligomer"]) > 1 : #select proteins who can create an homo-oligomer
        Use_RF2_PPI(HInt_object, Informations_dict, "RF2_homo_int", GPU) #set better interactions all_vs_bait #Maybe remove
    if Informations_dict["Interact_with"] != "" :
        if len(HInt_object.get_possible_prey()) > 1 :
           if Informations_dict["Organism"] == "euk" :
              Use_RF2_PPI(HInt_object, Informations_dict, "RF2_PPI_int", GPU) #set better interactions all_vs_bait
           if Informations_dict["Organism"] in ["gram+","gram-"] :
              Use_RF2_Lite(HInt_object, Informations_dict, GPU)      
        Generate_scripts(HInt_object, 1500, Informations_dict, "APD_PPI_int", GPU) #generate bait_vs_prey with new preys
        Generate_3D_model(Informations_dict, "APD_PPI_int", GPU)
        Score_interaction_APD(HInt_object, Informations_dict, "APD_PPI_int")
    if int(Informations_dict["Homo-oligomer"]) > 1 : #select proteins who can create an homo-oligomer
        Generate_scripts(HInt_object, 1500, Informations_dict, "APD_homo_int", GPU)
        Generate_3D_model(Informations_dict, "APD_homo_int", GPU)
        Score_interaction_APD(HInt_object, Informations_dict, "APD_homo_int")
    print(HInt_object.get_result_dict())
    Resume_file(HInt_object)

    
