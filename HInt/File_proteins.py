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

    def set_proteins_sequence (self, new_protein_sequence) :
        """
        Sets a dict of all sequences and set a lenght dict.
        
        Parameters:
        ----------
        new_protein_sequence = dictionary
        
        Returns:
        ----------
        """
        self.protein_sequence = new_protein_sequence
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

    def set_fasta_file (self, filename) :
        """
        Sets new filename for the fasta file.
        
        Parameters:
        ----------
        filename = string
        
        Returns:
        ----------
        """
        self.fasta_file = filename

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

    def set_names (self, name) :
        """
        Sets names of all proteins.
        
        Parameters:
        ----------
        name = dictionary
        
        Returns:
        ----------
        """
        self.name = name
    
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

    def get_proteins_sequence (self) :
        """
        Return the new amino acid sequence list.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        proteins_sequence : dictionary
        """
        return self.protein_sequence
    
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
    
    def get_fasta_file (self) :
        """
        Return the name of the fasta file.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        fasta_file : string
        """
        return self.fasta_file
    
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
    
    def get_names (self) :
        """
        Return names of proteins.
        
        Parameters:
        ----------
        
        Returns:
        ----------
        name : dictionary
        """
        return self.name

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
    
### Generating of features and pre-file to run multimer

    def set_all_att (self, path_txt) :
        """
        Set all values for all attribut for one txt file.
        
        Parameters:
        ----------
        filename : string
        
        Returns:
        ----------
        """
        new_proteins = list()
        already_fasta = dict()
        result_dict = dict()
        save_prot = ""
        with open(path_txt,"r") as in_file :
            for line in in_file :
                if "," in str(line) or (line[0] != ">" and save_prot == "") :
                    save_prot = ""
                    new_line = (line.split(","))
                    for prot in new_line :
                        new_proteins.append(prot.upper().strip())
                        result_dict[prot.upper().strip()] = ""
                elif line[0] == ">" :
                    save_prot = line[1:len(line)].strip("\n").strip(" ")
                    new_proteins.append(save_prot)
                    already_fasta[save_prot] = str()
                    result_dict[save_prot] = ""
                elif len(line) > 2 and save_prot != "" :
                    already_fasta[save_prot] = already_fasta[save_prot] + line.strip("\n")
        self.set_file_name(path_txt)
        self.set_proteins(new_proteins)
        self.set_possible_prey(new_proteins)
        self.set_proteins_sequence(already_fasta)
        self.set_result_dict(result_dict)
    
    def find_proteins_sequence (self, Path_Pickle_Feature) :
        """
        Search for the amino acid sequence on the UniProt website and clean it. If the a3m file exist, take the sequence from it and return a list of proteins that do not have a MSA. Set if the sequence have Peptide signal or not.
        
        Parameters:
        ----------

        Returns:
        ----------
        need_msa : list
        """
        need_msa = list()
        sequences = self.get_proteins_sequence()
        names = dict()
        pattern = r"SQ   SEQUENCE   .*  .*\n([\s\S]*)"
        pattern2 = r"GN   Name=([\w]*)"
        del_car = ["\n"," ","//"]
        for proteins in self.get_proteins() :
            if proteins not in sequences.keys() :
                if os.path.isfile(f"{Path_Pickle_Feature}/{proteins}.a3m") == True :
                    with open(f"{Path_Pickle_Feature}/{proteins}.a3m","r") as file :
                        new_sequence = str()
                        for line in file :
                            if line[0] != "#" and line[0] != ">" :
                                new_sequence = line.strip("\n")
                            if new_sequence != "" :
                                break
                        sequences[proteins] = new_sequence
                        print(f"Sequence for {proteins} found in pickle feature")
                else :
                    print("Search sequence for " + proteins)
                    urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/"+proteins+".txt","temp_file.txt")
                    if os.path.getsize("temp_file.txt") == 0 :
                        print(f"{proteins} is not a compliant UniprotID")
                    with open("temp_file.txt","r") as in_file:
                        for seq in re.finditer(pattern, in_file.read()):
                            sequences[proteins] = seq.group(1)
                    with open("temp_file.txt","r") as in_file:
                        for name in re.finditer(pattern2, in_file.read()) :
                            names[proteins] = name.group(1)
                    for car in del_car :
                        sequences[proteins] = sequences[proteins].replace(car,"")
                    os.remove("temp_file.txt")
                    need_msa.append(proteins)
            else :
                if os.path.isfile(f"{Path_Pickle_Feature}/{proteins}.a3m") == False :
                    need_msa.append(proteins)
        self.set_proteins_sequence(sequences)
        self.set_names(names)
        return need_msa

    def find_prot_lenght (self, prot_dict = None) :
        """
        Set the lenght for all proteins from a dictionary of amino acid sequence.

        Parameters:
        ----------
        prot_dict = dictionary
        
        Returns:
        ----------
        """
        #lenght_prot = self.get_lenght_prot() #not already set
        if prot_dict == None :
            proteins = self.get_proteins()
        else :
            proteins = prot_dict
        sequences = self.get_proteins_sequence()
        lenght_prot = dict()
        for protein in proteins :
            lenght_prot[protein] = len(sequences[protein])
        self.set_lenght_prot(lenght_prot)

    def create_fasta_file (self, need_msa) :
        """
        Generate a FASTA file for protein who don't have MSA.

        Parameters:
        ----------
        need_msa : list
        
        Returns:
        ----------
        """
        line = str()
        sequences = self.get_proteins_sequence()
        for protein in need_msa :
            line = line + ">" + protein + "\n" + sequences[protein] + "\n"        
        file_name = self.get_file_name()
        file_out = file_name.replace("txt","fasta")
        with open(file_out,"w") as fh :
            fh.write(line)
        self.set_fasta_file(file_out)
