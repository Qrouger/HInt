""" Main file of HInt

    Author: Quentin Rouger
"""
import argparse

from .Utils_HInt import *

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

def main() :
    path_dict = define_path()
    HInt_object = File_proteins(path_dict["Path_Uniprot_ID"])
    HInt_object.find_proteins_sequence() #and real name of proteins
    PPI_object.find_proteins_sequence() #and real name of proteins
    PPI_object.create_fasta_file()
    if args.use_signalP == True :
        remove_SP(PPI_object,args.org)
    if len(HInt_object.already_pickle(path_dict["Path_Pickle_Feature"])) > 0 : #if new feature pickle is need
        create_feature(PPI_object,path_dict["Path_AlphaFold_Data"],path_dict["Path_Pickle_Feature"],args.use_mmseq)
    Make_all_MSA_coverage(PPI_object,path_dict["Path_Pickle_Feature"]) #make MSA depth for new pickle and set shallow_MSA.txt
    recover_prot_sequence(PPI_object,path_dict["Path_Pickle_Feature"]) #set sequence dict without peptide signal
    HInt_object.find_prot_lenght()
    if bait_prot :
        #write_bait_file
        Rosetta_PPI #set better interactions
    
    #generate_APD_script(PPI_object, args.max_aa)
y


    
    #add_iQ_score(path_dict["Path_Singularity_Image"])
    #if path_dict["homo-oligomerization"] == "" :
    #  Make_homo_oligo(path_dict["Path_AlphaFold_Data"])
    #  add_hiQ_score(path_dict["Path_Singularity_Image"])
    #  PPI_object.update_iQ_score_hiQ_score()
        
        
 
