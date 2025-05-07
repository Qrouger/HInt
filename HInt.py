""" Main file of HInt

    Author: Quentin Rouger
"""

from Utils_HInt import *
from File_proteins import *

import sys
import logging



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

#enlever new_pickle ?
if __name__ == "__main__":
    informations_dict = define_informations()
    HInt_object = File_proteins(informations_dict["Path_Uniprot_ID"])
    HInt_object.find_proteins_sequence() #and real name of proteins
    HInt_object.create_fasta_file()
    remove_SP(HInt_object,informations_dict["Organism"])
    if informations_dict["Signal_peptide"] != "None" : #generated MSA just for protein of interest
        filtered_signalP(HInt_object,informations_dict)
    create_feature(HInt_object,informations_dict["Path_AlphaFold_Data"],informations_dict["Path_Pickle_Feature"])
    Make_all_MSA_coverage(HInt_object,informations_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    if informations_dict["Interact_with"] != "" :
        generate_bait_vs_prey(HInt_object,1500,informations_dict)
        while len(HInt_object.get_possible_prey()) > 15 :
            Rosetta_PPI #set better interactions all_vs_bait #need to set_possible_prey
        generate_bait_vs_prey(HInt_object,1500,informations_dict) #regenerate bait_vs_prey with new preys
        Make_bait_vs_prey(informations_dict)
 #   if informations_dict["Homo-oligomer"] >= 2 :
                    
 #   else :

    #generate_APD_script(PPI_object, args.max_aa)



    
    #add_iQ_score(informations_dict["Path_Singularity_Image"])
    #if informations_dict["homo-oligomerization"] == "" :
    #  Make_homo_oligo(informations_dict["Path_AlphaFold_Data"])
    #  add_hiQ_score(informations_dict["Path_Singularity_Image"])
    #  PPI_object.update_iQ_score_hiQ_score()
        
        
