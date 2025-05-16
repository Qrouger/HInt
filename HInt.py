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
    Informations_dict = Define_informations()
    HInt_object = File_proteins(Informations_dict["Path_Uniprot_ID"])
    HInt_object.find_proteins_sequence() #and real name of proteins
    HInt_object.create_fasta_file()
    remove_SP(HInt_object,Informations_dict["Organism"])
    if Informations_dict["Signal_peptide"] != "None" : #select proteins with or without peptide signal
        filtered_signalP(HInt_object,Informations_dict)
    create_feature(HInt_object,Informations_dict["Path_AlphaFold_Data"],Informations_dict["Path_Pickle_Feature"])
    Make_all_MSA_coverage(HInt_object,Informations_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    if Informations_dict["Homo-oligomer"] != 1 : #select proteins who can create an homo-oligomer
        Use_Rosetta_PPI(HInt_object,Informations_dict,"RoseTTAFold_homo_int") #set better interactions all_vs_bait
    if Informations_dict["Interact_with"] != "" :
        if len(HInt_object.get_possible_prey()) > 15 :
            Use_Rosetta_PPI(HInt_object,Informations_dict,"RoseTTAFold_PPI_int") #set better interactions all_vs_bait
        generate_bait_vs_prey(HInt_object,1500,Informations_dict) #regenerate bait_vs_prey with new preys
        Generate_3D_model(Informations_dict,"bait_vs_prey")
        Score_PPI_interaction(Informations_dict)

 #   if int(informations_dict["Homo-oligomer"]) >= 2 :
                    
 #   else :

    #generate_APD_script(PPI_object, args.max_aa)



    
    #add_iQ_score(informations_dict["Path_Singularity_Image"])
    #if informations_dict["homo-oligomerization"] == "" :
    #  Make_homo_oligo(informations_dict["Path_AlphaFold_Data"])
    #  add_hiQ_score(informations_dict["Path_Singularity_Image"])
        
