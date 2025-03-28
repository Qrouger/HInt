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


def add_arguments(parser) :
   # parser.add_argument("--use_mmseq", help = "Use MMseqs2 for feature generation (True/False)", required = False, default = True)
   # parser.add_argument("--make_multimers", help = "Enable or disable multimer model generation (True/False) and analyse", required = False, default = True)
    parser.add_argument("--max_aa" , help = "Maximum number of amino acids that can be generated per cluster", required = False, default = 2500, type = int)
    parser.add_argument("--org" , help = "Organism of interest: arch, gram+, gram-, or euk for SignalP", required = False, default = "gram-", type = str)

def main() :
    path_dict = define_path()
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    HInt_object = File_proteins(path_dict["Path_Uniprot_ID"])
    HInt_object.find_proteins_sequence() #and real name of proteins
    if len(HInt_object.already_pickle(path_dict["Path_Pickle_Feature"])) > 0 : #if new feature pickle is need
        HInt_object.create_fasta_file()
        remove_SP(HInt_object,args.org)
        if path_dict["Signal_P"] :
            for prot with signalP
        create_feature(HInt_object,path_dict["Path_AlphaFold_Data"])
    Make_all_MSA_coverage(PPI_object) #make MSA depth for new pickle and set shallow_MSA.txt
    recover_prot_sequence(PPI_object) #set sequence dict without peptide signal
    HInt_object.find_prot_lenght()
    generate_APD_script(PPI_object, args.max_aa)
    if nbr_all_int > 50 :
      Rosetta_PPI
    else :
      Make_all_int
    add_iQ_score(path_dict["Path_Singularity_Image"])
    if path_dict["homo-oligomerization"] == "" :
      Make_homo_oligo(path_dict["Path_AlphaFold_Data"])
      add_hiQ_score(path_dict["Path_Singularity_Image"])
      PPI_object.update_iQ_score_hiQ_score()
        
        
 
