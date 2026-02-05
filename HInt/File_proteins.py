""" Create File_proteins object

    Author: Quentin Rouger
"""
import urllib.request
import re
from Utils_HInt import *
import csv
import os
import copy

class File_proteins() :
    """
    Manipulate and save the file that contains all proteins.
    """
    def __init__ (self, path_txt_file) :
        """
        Constructor :
        Set attributes for a single entry file.

        Parameters:
    	-----------
        path_txt_file : string
        """
        self.set_all_att(path_txt_file)

    def set_proteins_sequence_SP (self, new_protein_sequence) :
        """
        Sets a dict of all sequences with Signal peptide.
        
        Parameters:
        ----------
        new_protein_sequence = dictionary
        
        Returns:
        ----------
        """
        self.protein_sequence_SP = new_protein_sequence

    def set_proteins_sequence_no_SP (self, new_protein_sequence) :
        """
        Sets a dict of all sequences without Signal peptide and set a lenght dict.
        
        Parameters:
        ----------
        new_protein_sequence = dictionary
        
        Returns:
        ----------
        """
        self.protein_sequence_no_SP = new_protein_sequence
        self.find_prot_lenght(new_protein_sequence)

    def set_proteins (self, new_protein) :
        """
        Sets a list of all proteins UniprotID.
        
        Parameters:
        ----------
        new_protein = list
        
        Returns:
        ----------
        """
        self.protein = new_protein

    def set_file_name (self, filename) :
        """
        Sets new filename for the txt file.
        
        Parameters:
        ----------
        filename = string
        
        Returns:
        ----------
        """
        self.file_name = filename

    def set_lenght_prot (self, lenght_prot) :
        """
        Sets lenght of all proteins.
        
        Parameters:
        ----------
        lenght_prot = dictionary
        
        Returns:
        ----------
        """
        self.lenght_prot = lenght_prot
    
    def set_result_dict (self, result_dict) :
        """
        Sets dict who will contains all results.
        
        Parameters:
        ----------
        result_dict = dict
        
        Returns:
        ----------
        """
        self.result_dict = result_dict
    
    def set_possible_prey (self, possible_prey) :
        """
        Sets list of all possible prey.
        
        Parameters:
        ----------
        possible_prey : list
        
        Returns:
        ----------
        """
        self.possible_prey = possible_prey

    def set_deeploc (self, deeploc) :
        """
        Sets a dict of all DeepLoc results.
        
        Parameters:
        ----------
        deeploc : dictionary
        
        Returns:
        ----------
        """
        self.deeploc = deeploc

    def set_int_score (self, int_score) :
        """
        Sets a dict of interaction score.
        
        Parameters:
        ----------
        int_score : dictionary
        
        Returns:
        ----------
        """
        self.int_score = int_score

    def set_prot_SP (self, dict_SP) :
        """
        Sets a dict indicating whether each protein has a signal peptide (SP) or not.
        
        Parameters:
        ----------
        int_score : dictionary
        
        Returns:
        ----------
        """
        self.prot_SP = dict_SP

    def get_proteins_sequence_SP (self) :
        """
        Return the new amino acid sequence dictionary with SP.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        proteins_sequence : dictionary
        """
        return self.protein_sequence_SP

    def get_proteins_sequence_no_SP (self) :
        """
        Return the new amino acid sequence dictionary without SP.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        proteins_sequence : dictionary
        """
        return self.protein_sequence_no_SP

    def get_proteins (self) :
        """
        Return the new proteins name list.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        protein : list
        """
        return self.protein
    
    def get_file_name (self) :
        """
        Return the name of the file.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        file_name : string
        """
        return self.file_name
    
    def get_lenght_prot (self) :
        """
        Return the lenght of proteins.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        lenght_prot : dictionary
        """
        return self.lenght_prot

    def get_result_dict (self) :
        """
        Return new result dict.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        result_dict : dict
        """
        return self.result_dict

    def get_possible_prey (self) :
        """
        Return list of all possible prey.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        possible_prey : list
        """
        return self.possible_prey
    
    def get_deeploc (self) :
        """
        Return the DeepLoc result dict.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        deeploc : dictionary
        """
        return self.deeploc
    
    def get_int_score (self) :
        """
        Return the interaction result dict.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        int_score : dictionary
        """
        return self.int_score

    def get_prot_SP (self) :
        """
        Return a dict indicating whether each protein has a signal peptide (SP) or not.
        
        Parameters:
        ----------
        int_score : dictionary
        
        Returns:
        ----------
        """
        return self.prot_SP

### Generating of features and pre-file to run multimer

    def set_all_att (self, path_txt) :
        """
        Initialize and populate all attributes from a protein input text file.

        This method parses a text file containing protein identifiers and/or FASTA sequences, cleans NCBI-formatted headers, validates protein uniqueness,
        checks sequence validity, and initializes all internal data structures required for downstream analyses.

        The input file may contain :
        - Comma-separated protein identifiers
        - FASTA-formatted protein sequences
        - NCBI-style headers, which are automatically cleaned and normalized

        Parameters :
        ----------
        path_txt : string

        Notes :
        ----------
        - Protein identifiers are converted to uppercase for consistency.
        - FASTA headers are simplified by extracting protein identifiers from NCBI annotations.
        - Sequences containing non-standard amino acids are rejected to ensure compatibility with MSA generation and peptid signal.
        """
        new_proteins = list()
        sequence_SP = dict()
        result_dict = dict()
        prot_SP = dict()
        already_fasta = dict()
        int_score = dict()
        save_prot = ""
        with open(path_txt,"r") as check_f : #clean ncbi file
            new_fasta = str()
            for line in check_f :
                if line[0] == ">" and  "[protein_id=" in line or "[locus_tag=" in line or "[gbkey=" in line : #clean ncbi file
                    new_fasta += ">" + line.split(" ")[1].split("=")[1][0:len(line.split(" ")[1].split("=")[1])-1] + "\n"
                else :
                    new_fasta += line.replace("*", "")
        with open(path_txt,"w") as w_file :
            w_file.write(new_fasta)
        with open(path_txt,"r") as in_file :
            for line in in_file :
                if "," in str(line) or (line[0] != ">" and save_prot == "" and line.strip() != "") :
                    save_prot = ""
                    new_line = (line.strip().split(","))
                    for prot in new_line :
                        if prot.upper().strip() in new_proteins :
                            raise ValueError(f"Protein {prot.upper().strip()} is duplicated in the input file.")
                        else :
                            new_proteins.append(prot.upper().strip())
                        result_dict[prot.upper().strip()] = dict()
                        int_score[prot.upper().strip()] = dict()
                elif line[0] == ">" :
                    save_prot = line[1:len(line)].strip("\n").strip(" ")
                    if save_prot in new_proteins :
                        raise ValueError(f"Protein {save_prot} is duplicated in the input file.")
                    else :
                        new_proteins.append(save_prot)
                    already_fasta[save_prot] = str()
                    sequence_SP[save_prot] = ""
                    result_dict[save_prot] = dict()
                    int_score[save_prot] = dict()
                elif len(line) > 1 and save_prot != "" :
                    sequence_SP[save_prot] = sequence_SP[save_prot] + line.strip("\n").strip("\t").replace(" ","")
                    for aa in ["O", "B", "Z", "J", "X"] :
                        if aa in line.strip("\n") :
                            raise ValueError(f"Sequence {save_prot} contains {aa}.")
        self.set_file_name(path_txt)
        self.set_proteins(new_proteins)
        self.set_possible_prey(new_proteins)
        self.set_proteins_sequence_SP(sequence_SP)
        self.set_int_score(int_score)
        self.set_result_dict(result_dict)
        self.set_prot_SP(prot_SP)

    def check_save_dict (self, Path_Pickle_Feature) :
        """
        Inspect previously saved computation states and determine missing features for each protein.

        This method compares the current protein set and sequences with the persistent save dictionary (save_dict.pkl) in order to identify which
        computational steps must be recomputed. It verifies the presence and consistency of :
        - Multiple sequence alignments (MSA, .a3m files)
        - Feature pickle files (.pkl)
        - DeepLoc subcellular localization predictions
        - Signal peptide annotations
        - Interaction scores

        Parameters :
        ----------
        Path_Pickle_Feature : str

        Returns :
        ----------
        need_msa : list
        need_pkl : list
        need_DeepLoc : list

        Notes :
        ----------
        - Sequence mismatches between MSAs and cached data trigger a full reset for the affected protein.
        - Interaction scores are invalidated if the corresponding structural models are missing.
        - UniProt sequences are retrieved using the API and cleaned to remove formatting artifacts.
        - The save dictionary is always updated at the end of the procedure.
        """
        need_msa = list()
        need_pkl = list()
        need_DeepLoc = list()
        protein_sequence_no_SP = dict()
        deeploc_prot = dict()
        int_score = self.get_int_score()
        result_dict = self.get_result_dict()
        pattern = r"SQ   SEQUENCE   .*  .*\n([\s\S]*)"
        del_car = ["\n"," ","//"]

        sequences_SP = self.get_proteins_sequence_SP()
        proteins = self.get_proteins()

        if os.path.isfile('log_file/save_dict.pkl') == True :
            with open('log_file/save_dict.pkl', 'rb') as save_dict :
                all_info = pickle.load(save_dict)
            deeploc_prot = copy.deepcopy(all_info["deeploc"])
            protein_sequence_no_SP = copy.deepcopy(all_info["sequence_no_SP"])
            int_score.update(copy.deepcopy(all_info["int_score"]))
            prot_SP = copy.deepcopy(all_info["Signal_peptide"])
            for protein in proteins :
                result_dict[protein] = dict()
                if protein in deeploc_prot.keys() : #if protein in deeploc of save dict
                    result_dict[protein]["DeepLoc"] = deeploc_prot[protein]
                if protein in prot_SP.keys() : #if protein in SP of save dict
                    result_dict[protein]["Signal_peptide"] = prot_SP[protein]
                if protein in all_info["sequence_SP"].keys() : #if protein in sequence_SP of save dict
                    if protein in sequences_SP.keys() : #check if two sequence match
                        if sequences_SP[protein] != all_info["sequence_SP"][protein] : #if not match, remove MSA files and start at zero
                            cmd = f"rm -rf {Path_Pickle_Feature}/*{protein}*"
                            cmd2 = f"rm -rf ./result_PPI_int/*{protein}*"
                            cmd3 = f"rm -rf ./result_homo_int/*{protein}*"
                            os.system(cmd)
                            os.system(cmd2)
                            os.system(cmd3)
                            need_msa.append(protein)
                            need_DeepLoc.append(protein)
                            int_score[protein] = dict() #remove int score
                        else : #sequences match, set all arguments
                            sequences_SP[protein] = all_info["sequence_SP"][protein]
                            if protein not in all_info["sequence_no_SP"].keys() :
                                need_msa.append(protein)
                            if protein not in all_info["deeploc"].keys() :
                                need_DeepLoc.append(protein)

                    else : #protein not in fasta format, so UniprotID is good
                        if protein in all_info["sequence_SP"].keys() :
                            sequences_SP[protein] = all_info["sequence_SP"][protein]
                        if protein not in all_info["sequence_no_SP"].keys() :
                            need_msa.append(protein)
                        if protein not in all_info["deeploc"].keys() :
                            need_DeepLoc.append(protein)
                if protein not in all_info["sequence_SP"].keys() :
                    if protein not in sequences_SP.keys() : #if protein need sequence in Uniprot
                        logger.info("Search sequence for " + protein)
                        try:
                            urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/"+protein+".txt","log_file/temp_file.txt")
                        except Exception as e :
                            raise Exception(f"{protein} is not a compliant UniprotID")
                            break
                        with open("log_file/temp_file.txt","r") as in_file:
                            for seq in re.finditer(pattern, in_file.read()) :
                                sequences_SP[protein] = seq.group(1)
                        #  with open("temp_file.txt","r") as in_file:
                        #      for name in re.finditer(pattern2, in_file.read()) :
                        #          names[protein] = name.group(1)
                        for car in del_car :
                            sequences_SP[protein] = sequences_SP[protein].replace(car,"")
                        os.remove("log_file/temp_file.txt")
                        need_DeepLoc.append(protein) #add of new protein in txt file
                        need_msa.append(protein)
                    else : #protein in fasta format, who don't have sequence without SP and DeepLoc
                        need_DeepLoc.append(protein)
                        need_msa.append(protein)
                if os.path.isfile(f"{Path_Pickle_Feature}/{protein}.a3m") == False and protein not in need_msa : #proteins without msa
                    need_msa.append(protein)
                    cmd = f"rm -rf {Path_Pickle_Feature}/*{protein}*" #remove residue files
                    os.system(cmd)
                    int_score[protein] = dict() #remove int score
                if os.path.isfile(f"{Path_Pickle_Feature}/{protein}.pkl") == False and os.path.isfile(f"{Path_Pickle_Feature}/{protein}.a3m") == True : #proteins with msa without pkl file
                    need_pkl.append(protein)
                if os.path.isfile(f"{Path_Pickle_Feature}/{protein}.a3m") == True and protein not in need_msa : #check sequence in msa file match with save dict
                    with open(f"{Path_Pickle_Feature}/{protein}.a3m","r") as msa_file :
                        index = 0
                        for line in msa_file :
                            if index == 1 :
                                msa_seq = line.strip("\n")
                                break
                            if line[0] == ">" :
                                index += 1
                    if msa_seq != all_info["sequence_no_SP"][protein] : #if not match, remove pkl file and start at zero
                        cmd = f"rm -rf {Path_Pickle_Feature}/*{protein}*"
                        os.system(cmd)
                        need_msa.append(protein)
                        need_DeepLoc.append(protein)
                        int_score[protein] = dict() #remove int score

            #Check is save score interaction are still valid or delete it
            for protein in all_info["int_score"].keys() :
                for key in all_info["int_score"][protein].keys() :
                    bait = key.split("iQ_score_vs_")[1]
                    if os.path.isfile(f"./result_PPI_int/{bait}_and_{protein}/ranked_0.pdb") == False :
                        if f"iQ_score_vs_{bait}" in int_score[protein].keys() :
                            del int_score[protein][f"iQ_score_vs_{bait}"]
                        else :
                            del int_score[protein]["iQ_score"]

        else : #no save dict, so check if protein have MSA or pkl
            for protein in proteins :
                if protein not in sequences_SP.keys() : #if protein need sequence in Uniprot
                    print("Search sequence for " + protein)
                    urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/"+protein+".txt","log_file/temp_file.txt")
                    if os.path.getsize("log_file/temp_file.txt") == 0 :
                        logger.error(f"{protein} is not a compliant UniprotID")
                    with open("log_file/temp_file.txt","r") as in_file :
                        for seq in re.finditer(pattern, in_file.read()):
                            sequences_SP[protein] = seq.group(1)
                    #  with open("temp_file.txt","r") as in_file :
                    #      for name in re.finditer(pattern2, in_file.read()) :
                    #          names[protein] = name.group(1)
                    for car in del_car :
                        sequences_SP[protein] = sequences_SP[protein].replace(car,"")

                    os.remove("log_file/temp_file.txt")
                if os.path.isfile(f"{Path_Pickle_Feature}/{protein}.a3m") == False : #proteins without msa
                    cmd = f"rm -rf {Path_Pickle_Feature}/*{protein}*" #remove residue files
                    os.system(cmd)
                if os.path.isfile(f"{Path_Pickle_Feature}/{protein}.pkl") == False and os.path.isfile(f"{Path_Pickle_Feature}/{protein}.a3m") == True : #proteins with msa without pkl file
                    need_pkl.append(protein)
                    
                need_msa.append(protein) #make SignalP and DeepLoc for all proteins, too create the save dict
                need_DeepLoc.append(protein)
        self.set_result_dict(result_dict)
        self.set_deeploc(deeploc_prot)
        self.set_proteins_sequence_no_SP(protein_sequence_no_SP)
        self.set_proteins_sequence_SP(sequences_SP)
        self.set_int_score(int_score)
        self.Make_save_dict() #Save modification of the save dict
        return need_msa, need_pkl, need_DeepLoc


    def Make_save_dict (self) :
        """
        Save all relevant intermediate results into a pickle file in order to avoid recomputing completed steps in future runs.
        """
        pkl_dict = dict()
        pkl_dict["sequence_SP"] = self.get_proteins_sequence_SP()
        pkl_dict["sequence_no_SP"] = self.get_proteins_sequence_no_SP()
        pkl_dict["deeploc"] = self.get_deeploc()
        pkl_dict["int_score"] = self.get_int_score()
        pkl_dict["Signal_peptide"] = self.get_prot_SP()
        with open('log_file/save_dict.pkl', 'wb') as out_file :
            pickle.dump(pkl_dict, out_file)

    def find_prot_lenght (self, prot_dict = None) :
        """
        Compute and store the length (number of amino acids) of each protein based on sequences without signal peptides.

        Parameters :
        ----------
        prot_dict = dict or None
        """
        if prot_dict == None :
            proteins = self.get_proteins()
        else :
            proteins = prot_dict
        sequences = self.get_proteins_sequence_no_SP()
        lenght_prot = dict()
        for protein in proteins :
            lenght_prot[protein] = len(sequences[protein])
        self.set_lenght_prot(lenght_prot)


    def create_fasta_file (self, with_SP, need_msa=[], need_pkl=[]) :
        """
        Generate FASTA files for proteins requiring MSA or pkl generation.

        Parameters :
        ----------
        with_SP : boolean
        need_msa : list
        need_pkl : list
        """
        line_msa = str()
        sequences_SP = self.get_proteins_sequence_SP()
        sequences_no_SP = self.get_proteins_sequence_no_SP()
        file_name = self.get_file_name().split("/")[-1]
        file_msa = file_name.replace(".txt","_msa.fasta")
        file_pkl = file_name.replace(".txt","_pkl.fasta")
        if os.path.isfile(f"log_file/{file_msa}") == True :
            os.remove(f"log_file/{file_msa}")
        if os.path.isfile(f"log_file/{file_pkl}") == True :
            os.remove(f"log_file/{file_pkl}")
        if len(need_msa) != 0 :
            if with_SP == True :
                sequences = sequences_SP
            if with_SP == False :
                sequences = sequences_no_SP
            for protein in need_msa :
                line_msa += ">" + protein + "\n" + sequences[protein] + "\n"
            with open(f"log_file/{file_msa}","w") as f_msa :
                f_msa.write(line_msa)

        line_pkl = str()
        if len(need_pkl) != 0 :
            for protein in need_pkl :
                line_pkl += ">" + protein + "\n" + sequences_no_SP[protein] + "\n"
            with open(f"log_file/{file_pkl}","w") as f_pkl :
                f_pkl.write(line_pkl)

