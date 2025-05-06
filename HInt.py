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

if __name__ == "__main__":
    informations_dict = define_informations()
    print(informations_dict)
    HInt_object = File_proteins(informations_dict["Path_Uniprot_ID"])
    HInt_object.find_proteins_sequence() #and real name of proteins
    HInt_object.create_fasta_file()
    #remove_SP(HInt_object,informations_dict["Organism"])
    if len(HInt_object.already_pickle(informations_dict["Path_Pickle_Feature"])) > 0 : #if new feature pickle is need
        create_feature(HInt_object,informations_dict["Path_AlphaFold_Data"],informations_dict["Path_Pickle_Feature"])
    Make_all_MSA_coverage(HInt_object,informations_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    recover_prot_sequence(HInt_object,informations_dict["Path_Pickle_Feature"]) #set sequence dict without peptide signal #remplacer par une fonctions qui prends les séquences en fonction du SignalP
    HInt_object.find_prot_lenght()
    if informations_dict["Signal_peptide"] != "None" : #mettre avant la génération des MSA ?
        filtered_signalP(HInt_object,informations_dict["Signal_peptide"])
        print(HInt_object.get_possible_prey())
    if informations_dict["Interact_with"] != "" :
        if len(HInt_object.get_possible_prey()) > 15 :
            Rosetta_PPI #set better interactions all_vs_bait
        else :
            generate_bait_vs_prey(HInt_object,1500,informations_dict)

 #   if informations_dict["Homo-oligomer"] >= 2 :
                    
 #   else :

    #generate_APD_script(PPI_object, args.max_aa)



    
    #add_iQ_score(informations_dict["Path_Singularity_Image"])
    #if informations_dict["homo-oligomerization"] == "" :
    #  Make_homo_oligo(informations_dict["Path_AlphaFold_Data"])
    #  add_hiQ_score(informations_dict["Path_Singularity_Image"])
    #  PPI_object.update_iQ_score_hiQ_score()
        
        
