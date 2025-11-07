import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import csv
import logging
import subprocess
import multiprocessing
import pynvml
import copy
import sys
import get_good_inter_pae
import signal
import glob
import pandas as pd
from pathlib import Path
from datetime import datetime
from numpy import load
from Bio import SeqIO
from tqdm import tqdm
from Bio.PDB import MMCIFParser, PDBIO
from contextlib import closing

# Configure global logger
logging.basicConfig(
    filename="HInt.log",  # Log file name
    level=logging.INFO,  # Log level
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

# Use the logger configuration from HInt.py
logger = logging.getLogger()

def Define_informations() :
    """
    Extract all paths from HInt.txt and store them in a dictionary.
    
    Parameters:
    ----------

    Returns:
    ----------
    Informations_dict : dict
    """
    logger.info("Defining informations")
    Informations_dict = dict()
    list_inf = ["Signal_peptide", "Homo-oligomer", "Interact_with", "Organism", "DeepLoc", "Regions", "Path_Uniprot_ID", "Path_AlphaFold_Data", "Path_Pickle_Feature", "Path_Singularity_Image", "Path_MMseqs2_Data"]
    with open("HInt.txt", "r") as file :
        for lines in file :
            if ":" in lines :
                info = lines.split("#")[0]
                informations_name = info.split(":")[0].strip().strip("\n")
                informations = info.split(":")[1].strip().strip("\n")
                Informations_dict[informations_name] = informations
    for info in list_inf :
        if info not in Informations_dict.keys() : #if settings file is not authentic
            if info in ["Interact_with", "Organism","Path_Uniprot_ID", "Path_AlphaFold_Data", "Path_Pickle_Feature"] :
                raise ValueError(f"HInt.txt file is compromised, verify the file. {info} is missing")
            elif info in ["Signal_peptide","Homo-oligomer","Path_MMseqs2_Data","Regions","DeepLoc"] :
                Informations_dict[info] = ""
    for informations_key in Informations_dict.keys() : #verify all informations and set default value
        if type(Informations_dict[informations_key]) is str and Informations_dict[informations_key].endswith("/") : #avoid error in path
            Informations_dict[informations_key] = Informations_dict[informations_key][:-1]
        if len(Informations_dict[informations_key]) == 0 :
            logger.info(f'Informations : {informations_key} is empty')
            if informations_key == "Path_ccp4" :
                logger.info("Set ccp4 path by default on ./opt/xtal/ccp4-9")
                Informations_dict[informations_key] = "./opt/xtal/ccp4-9"
            elif informations_key == "Path_AlphaFold_Data" :
                logger.info("Set AlphaFold data by default on ./alphadata")
                Informations_dict[informations_key] = "./alphadata"
            elif informations_key == "Path_Pickle_Feature" :
                logger.info("Set pickle feature path by default on ./feature")
                Informations_dict[informations_key] = "./feature"
            elif informations_key == "Path_MMseqs2_Data" :
                logger.info("/!\ local MMseqs2 GPU will not be used")
            elif informations_key == "Signal_peptide" :
                Informations_dict[informations_key] = "None"
            elif informations_key == "Homo-oligomer" :
                Informations_dict[informations_key] = "1"
            elif informations_key == "DeepLoc" :
                Informations_dict[informations_key] = "None"
        if informations_key == "Interact_with" :
            regions_dict = dict()
            new_baits_list = list()
            list_baits = [prot for prot in Informations_dict["Interact_with"].split(",")]
            if len(list_baits) > 3 :
                raise ValueError("HInt don't support more than 3 protein baits")
            for prot in list_baits :
                if "-" in prot :
                    name_prot = prot.split("(")[0]
                    new_baits_list.append(name_prot)
                    regions_dict[name_prot] = prot.split("(")[1].strip(")")
                else :
                    regions_dict[prot] = "0-0"
                    new_baits_list.append(prot)
            Informations_dict["Regions"] = regions_dict
            Informations_dict["Interact_with"] = new_baits_list
        if informations_key == "Homo-oligomer" :
            if int(Informations_dict[informations_key]) == False :
                raise ValueError(f"Homo-oligomer is not an integer")
            if Informations_dict[informations_key] == "0" :
                Informations_dict[informations_key] = "1"
        if informations_key == "DeepLoc" :
            if Informations_dict["Organism"] == "euk" : #euk
                for value in Informations_dict[informations_key].split(",") :
                    if value not in ["Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion", "Plastid", "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgo apparatus", "Peroxisome","None"] :
                        raise ValueError(f"Incorrect DeepLoc value : {value}")
            else : #other
                for value in Informations_dict[informations_key].split(",") :
                    if value not in ["Cell wall & surface","Extracellular","Cytoplasmic","Cytoplasmic Membrane","Outer Membrane","Periplasmic","None"] :
                        raise ValueError(f"Incorrect DeepLocPro value : {value}")
    if len(Informations_dict["Signal_peptide"]) == 0 and len(Informations_dict["DeepLoc"]) == 0 and len(Informations_dict["Homo-oligomer"]) == 0 and len(Informations_dict["Interact_with"]) == 0 : #no info
        logger.info("need information to discriminate the potential homologue")
        exit()
    return(Informations_dict)

def run_deeploc(file, org, need_DeepLoc, GPU) :
    """
    Launch DeepLoc on sequence with Signal peptide and class protein in function.

    Parameters :
    ----------
    file : object of class File_proteins
    org : string
    need_DeepLoc : list
    GPU : list

    Returns:
    ----------
    """
    result_dict = file.get_result_dict()
    prot_seq = file.get_proteins_sequence_SP()
    deeploc = file.get_deeploc()
    GPU_str = ""
    for nbr_GPU in GPU :
        GPU_str += nbr_GPU + ","
    GPU_str = GPU_str.strip(",")
    file_name = file.get_file_name()
    fasta_file = file_name.replace(".txt","_msa.fasta")

    if org == "euk" :
        logger.info(str(datetime.now())+" Start DeepLoc eucaryote")
        software = "deeploc2"
    else :
        logger.info(str(datetime.now())+" Start DeepLocPro")
        software = "deeplocpro"

    dp_lines = str()
    for protein in need_DeepLoc :
        dp_lines += ">"+protein+"\n"+prot_seq[protein]+"\n"
    with open(f"log_file/{fasta_file}", "w") as dp_file :
        dp_file.write(dp_lines)

    if os.path.exists("log_file/result_deeploc") == True :
        os.system("rm -r log_file/result_deeploc")
    cmd = f"CUDA_VISIBLE_DEVICES={GPU_str} {software} -f log_file/{fasta_file} -o log_file/result_deeploc -d cuda"
    os.system(cmd)
    DeepLoc_file = os.listdir(f"log_file/result_deeploc")
    if software == "deeplocpro" :
        min_index = 3
        min_score = 0.1
        index_name = 1
    else :
        min_index = 4
        min_score = 0.4
        index_name = 0
    with open(f"log_file/result_deeploc/{DeepLoc_file[0]}", "r") as DL_file :
        reader = csv.reader(DL_file, delimiter=',')
        for index, line in enumerate(reader) :
            compartment = str()
            if index == 0 :
                first_line = line #save title name
            else :
                protein = line[index_name]
                for index, score in enumerate(line) :
                    if index >= min_index and float(score) > min_score :
                        compartment += first_line[index] + "|"
                compartment = compartment.strip("|")
                
                deeploc[protein] = compartment
                result_dict[protein]["DeepLoc"] = deeploc[protein]
    file.set_result_dict(result_dict)
    file.set_deeploc(deeploc)



def run_SP (file, Informations_dict, need_prot) :
    """
    Create a new FASTA file without the signal peptide using SignalP and filtred signal peptide.
    If protein have already msa but no sequence without SP don't return it.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dict
    need_prot : list

    Returns:
    ----------
    need_msa : list
    """
    final_file = str()
    SP_signal = 0
    prot_SP = dict()
    file_name = file.get_file_name()
    fasta_file = file_name.replace(".txt","_msa.fasta")
    output_file = fasta_file.replace(".fasta","")

    file.create_fasta_file(True, need_prot)
    cmd1 = f"signalp -fasta log_file/{fasta_file} -org {Informations_dict['Organism']} -prefix log_file/{output_file}"
    os.system(cmd1)
    file_signalp = fasta_file.replace(".fasta","_summary.signalp5")
    with open(f"log_file/{file_signalp}","r") as fh :
        for line in fh :
            new_line = line.split("\t")
            if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                prot_SP[new_line[0]] = new_line[len(new_line)-1].split("-")[1].split(".")[0]
    new_fasta_dict = file.get_proteins_sequence_no_SP()
    with open(f"log_file/{fasta_file}", "r") as fa_file :
        for line2 in fa_file :
            new_line2 = line2
            if SP_signal == 0 and line2[0] != ">" :
                new_line2 = line2
                new_fasta_dict[save_key] = line2.strip("\n")
            if int(SP_signal) > 0 :
                new_line2 = line2[int(SP_signal)-1:len(line2)]
                new_fasta_dict[save_key] = line2[int(SP_signal)-1:len(line2)].strip("\n")
                SP_signal = 0
            if line2[0] == ">" :
                save_key = line2[1:len(line2)-1]
                if str(line2[1:len(line2)-1]) in prot_SP.keys() :
                    SP_signal = prot_SP[line2[1:len(line2)-1]]
            final_file = final_file + new_line2
    need_msa = list()
    for prot in need_prot : 
        if os.path.isfile(f"{Informations_dict['Path_Pickle_Feature']}/{prot}.a3m") == False :
            need_msa.append(prot)
    file.set_proteins_sequence_no_SP(new_fasta_dict)
    return need_msa

def create_feature (file, Informations_dict, GPU, need_msa, need_pkl) :
    """
    Launch command to generate features.

    Parameters:
    ----------
    file : object of class File_proteins
    informations_dict : dict
    GPU : list
    need_msa : list
    need_pkl : list

    Returns:
    ----------
    """
    Path_AlphaFold_Data = Informations_dict["Path_AlphaFold_Data"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]
    baits = Informations_dict["Interact_with"]
    regions = Informations_dict["Regions"]
    Path_MMseqs2_Data = Informations_dict["Path_MMseqs2_Data"]
    prot_no_SP = file.get_proteins_sequence_no_SP()
    prot_SP = file.get_proteins_sequence_SP()
    file_name = file.get_file_name()
    msa_name = file_name.replace(".txt","_msa.fasta")
    pkl_name = file_name.replace(".txt","_pkl.fasta")
    GPU_str = ""
    for nbr_GPU in GPU :
        GPU_str += nbr_GPU + ","
    GPU_str = GPU_str.strip(",")
    logger.info(f"GPU use : {GPU_str}")
    
    logger.info(f"{len(need_msa)} proteins need msa")
    logger.info(f"{len(need_pkl)} proteins need pkl files")
    
    if os.path.exists(Path_Pickle_Feature) == False :
        os.system(f"mkdir {Path_Pickle_Feature}")

    #Look for MSA in AlphaFold database
    generated_msa = copy.deepcopy(need_msa)
    for protein in generated_msa : #check if prot have an MSA in alphafold database
        url = f"https://alphafold.ebi.ac.uk/files/msa/AF-{protein}-F1-msa_v6.a3m"
        name_file = Path_Pickle_Feature + "/" + protein +".a3m"
        outfile = os.path.basename(url)
        check = subprocess.run(["wget", "--spider", "-q", url])
        if check.returncode == 0 :
            subprocess.run(["wget", "-q", "-O",name_file , url], check=True)
            logger.info(f"MSA for {protein} found in AF database")
            need_msa.remove(protein) #msa found
            need_pkl.append(protein)
            #Cut SP for all MSA
            msa_in = f"{Path_Pickle_Feature}/{protein}.a3m"
            SP = len(prot_SP[protein])-len(prot_no_SP[protein])
            if SP > 0 : #if no SP don't modify the MSA
                logger.info(f"Remove SP from MSA")
                trimmed_records = []
                for rec in SeqIO.parse(msa_in, "fasta") :
                    new_seq = rec.seq[SP:]  #cut SP from MSA
                    if any(c.isupper() for c in str(new_seq)): #remove empty sequence
                        new_rec = rec[:]
                        new_rec.seq = new_seq
                        trimmed_records.append(new_rec)
                if trimmed_records : #if objects is not empty
                    with open(msa_in, "w") as msa_file : #overwrites the old MSA
                        for rec in trimmed_records:
                            msa_file.write(f">{rec.description}\n{rec.seq}\n")
                num_cpus_os = os.cpu_count()
                os.system(f"mafft --quiet --anysymbol --thread {num_cpus_os} --auto {Path_Pickle_Feature}/{protein}.a3m > {Path_Pickle_Feature}/{protein}.aln") #realign the new cut MSA
                os.system(f"reformat.pl fas a3m {Path_Pickle_Feature}/{protein}.aln {Path_Pickle_Feature}/{protein}.a3m")
                os.system(f"rm {Path_Pickle_Feature}/{protein}.aln")
                #Delete \n
                all_lines = str()
                with open(f"{Path_Pickle_Feature}/{protein}.a3m","r") as in_a3m :
                    seq = ""
                    header = None
                    for line in in_a3m:
                        if line.startswith(">") :
                            if header:
                                all_lines += f"{header}\n{seq}\n"
                            header = line.strip()
                            seq = ""
                        else:
                            seq += line.strip()
                    if header :
                        all_lines += f"{header}\n{seq}\n"
                with open(f"{Path_Pickle_Feature}/{protein}.a3m","w") as out_a3m :
                    out_a3m.write(all_lines)

    file.create_fasta_file(False, need_msa, need_pkl)

    #Create MSA files with ColabFold mmseq2 GPU accelerated for proteins without MSA
    if len(need_msa) > 100 : #if more than 50 sequences, use local colabfold_search GPU
        cmd = f"CUDA_VISIBLE_DEVICES={GPU_str} colabfold_search /log_file/{msa_name} {Path_MMseqs2_Data} {Path_Pickle_Feature} --db-load-mode 2 --gpu 1 "  #-e 0.1
        os.system(cmd)

    if len(need_msa) < 1 :
        logger.info("All MSAs have already been generated")

    #Create MSA files with ColabFold mmseq2 classic pipeline for proteins without MSA
    if os.path.isfile(f"log_file/{msa_name}") == True :
        cmd = ["create_individual_features.py",
        f"--fasta_paths=./log_file/{msa_name}",
        f"--data_dir={Path_AlphaFold_Data}",
        "--save_msa_files=False",
        f"--output_dir={Path_Pickle_Feature}",
        "--max_template_date=2024-05-02",
        "--skip_existing=True",
        "--use_mmseqs2=True",
        "--use_precomputed_msas=True"]
        process = subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
        stdout, stderr = process.communicate()
        
    #Create pkl files for proteins without pkl file
    if os.path.isfile(f"log_file/{pkl_name}") == True :
        cmd2 = ["create_individual_features.py",
        f"--fasta_paths=./log_file/{pkl_name}",
        f"--data_dir={Path_AlphaFold_Data}",
        "--save_msa_files=False",
        f"--output_dir={Path_Pickle_Feature}",
        "--max_template_date=2024-05-02",
        "--skip_existing=True",
        "--use_mmseqs2=True",
        "--use_precomputed_msas=True"]
        process = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
        for line in process.stdout:
            logger.info(line)
        stdout, stderr = process.communicate()

    ### if regions in bait, cut the MSA  No need if we made APD on a segment
 #   all_lines = str()
 #   for bait in baits :
 #       start = int(regions[bait].split("-")[0]) - 1
 #       end = int(regions[bait].split("-")[1]) - 1
 #       bait_name_file = f"{bait}_{str(start+1)}_{str(end+1)}"
 #       if regions[bait] != "0-0" and f"{bait_name_file}.a3m" not in os.listdir(Path_Pickle_Feature) : #if MSA regions of bait not exist
 #           with open(f"{Path_Pickle_Feature}/{bait}.a3m", "r") as a3m_file :
 #               for line in a3m_file:
 #                   if line[0] == ">" :
 #                       mem_name = line.strip().split("\t")[0] + "\n"
 #                   else :
 #                       if not all(cara == '-' for cara in line[start:end]) :
 #                           all_lines += mem_name + line[start:end] + "\n"
 #           with open(f"{Path_Pickle_Feature}/{bait_name_file}.a3m", "w") as new_file :
 #               new_file.write(all_lines)
 #           cmd1 = f"mafft --anysymbol {Path_Pickle_Feature}/{bait_name_file}.a3m > {Path_Pickle_Feature}/{bait_name_file}.aln"
 #           cmd2 = f"reformat.pl fas a3m {Path_Pickle_Feature}/{bait_name_file}.aln {Path_Pickle_Feature}/{bait_name_file}.a3m"
 #           os.system(cmd1)
 #           os.system(cmd2)
 #           with open(f"{Path_Pickle_Feature}/{bait_name_file}.a3m", "r") as a3m2_file :
 #               all_lines = str()
 #               all_seq = str()
 #               for line in a3m2_file:
 #                   if line[0] == ">" :
 #                       if all_seq != "" :
 #                           all_lines += mem_name + all_seq + "\n"
 #                       mem_name = line.strip().split("\t")[0] + "\n"
 #                       all_seq = str()
 #                   else :
 #                       all_seq += line.strip()
 #           with open(f"{Path_Pickle_Feature}/{bait_name_file}.a3m", "w") as a3m2_file :
 #               a3m2_file.write(all_lines)
 #           cmd3 = f"rm {Path_Pickle_Feature}/{bait_name_file}.aln"
 #           os.system(cmd3)
            
def filter_signalP(file, Informations_dict, need_msa, need_pkl) :
    """
    Filter proteins based on the presence of a signal peptide using SignalP results.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dict
    need_msa : list
    need_pkl : list

    Returns:
    ----------
    need_msa : list
    need_pkl : list
    """
    result_dict = file.get_result_dict()
    seq_dict = file.get_proteins_sequence_no_SP()
    possible_prey = file.get_possible_prey()
    new_possible_prey = list()
    SignalP = Informations_dict["Signal_peptide"]
    for protein in possible_prey :
        if seq_dict[protein][0] == "M" :
            result_dict[protein]["Signal_peptide"] = "No"
        else :
            result_dict[protein]["Signal_peptide"] = "Yes"
        if SignalP == "None" :
            new_possible_prey.append(protein)
        elif SignalP == result_dict[protein]["Signal_peptide"] :
            new_possible_prey.append(protein)
        elif SignalP != result_dict[protein]["Signal_peptide"] :
            result_dict[protein]["Reason_for_filtering"] = "Signal peptide : SignalP"
            if protein in need_msa :
                need_msa.remove(protein)
            if protein in need_pkl :
                need_pkl.remove(protein)
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)
    return need_msa, need_pkl

def Make_all_MSA_coverage (file, Path_Pickle_Feature) :
    """
    Generating MSA coverage for all proteins and write shallow_MSA text file.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_Pickle_Feature : string

    Returns:
    ----------
    """
    shallow_MSA = str()
    result_dict = file.get_result_dict()
    possible_prey = file.get_possible_prey()
    new_possible_prey = list()
    for prot in possible_prey :
        if Path(f'{Path_Pickle_Feature}/{prot}_coverage.pdf').exists() :
            pass
        else :
            pre_feature_dict = pickle.load(open(f'{Path_Pickle_Feature}/{prot}.pkl','rb'))
            feature_dict = pre_feature_dict.feature_dict
            msa = feature_dict['msa']
            seqid = (np.array(msa[0] == msa).mean(-1))
            seqid_sort = seqid.argsort()
            non_gaps = (msa != 21).astype(float)
            non_gaps[non_gaps == 0] = np.nan
            final = non_gaps[seqid_sort] * seqid[seqid_sort, None]
            plt.figure(figsize=(14, 4), dpi=100)
            plt.subplot(1, 2, 1)
            plt.title(f"Sequence coverage ({prot})")
            plt.imshow(final, interpolation='nearest', aspect='auto', cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
            plt.plot((msa != 21).sum(0), color='black')
            plt.xlim(-0.5, msa.shape[1] - 0.5)
            plt.ylim(-0.5, msa.shape[0] - 0.5)
            plt.colorbar(label="Sequence identity to query", )
            plt.xlabel("Positions")
            plt.ylabel("Sequences")
            plt.savefig(f"{Path_Pickle_Feature}/{prot+('_' if prot else '')}coverage.pdf")
            plt.close()
        a3m_file = open(f'{Path_Pickle_Feature}/{prot}.a3m', 'r')
        msa = subprocess.run(['wc', '-l', f'{Path_Pickle_Feature}/{prot}.a3m'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        line_msa = int(msa.stdout.split()[0])
        if line_msa/2 <= 100 :
            if line_msa/2 <= 2 :
               shallow_MSA += prot + " : " + str(int(line_msa/2)) + " sequences\n"
               result_dict[prot]["shallow_MSA"] = "No MSA"
               result_dict[prot]["Reason_for_filtering"] = "No MSA"
            else :
               shallow_MSA += prot + " : " + str(int(line_msa/2)) + " sequences\n"
               result_dict[prot]["shallow_MSA"] = "Shallow MSA"
               new_possible_prey.append(prot)
        else :
            new_possible_prey.append(prot)
    with open("log_file/shallow_MSA.txt", "w") as MSA_file :
        MSA_file.write(shallow_MSA)
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)


def filter_deeploc(file, Informations_dict, need_msa, need_pkl) :
    """
    Filter proteins based on cellular localisation and set new prey list and new need_msa/need_pkl list.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dict
    need_msa : list
    need_pkl : list

    Returns:
    ----------
    new_need_msa : list
    new_need_pkl : list
    """
    localisation = Informations_dict["DeepLoc"].split(",")
    deeploc = file.get_deeploc()
    possible_baits = Informations_dict["Interact_with"]
    possible_prey = file.get_possible_prey()
    result_dict = file.get_result_dict()
    new_possible_prey = list()
    new_need_msa = list()
    new_need_pkl = list()
    for protein in possible_prey :
        for loc in deeploc[protein].split("|") :
            if loc in localisation :
                new_possible_prey.append(protein)
                break #stop the loop if one localisation is correct
        if protein not in new_possible_prey :
            result_dict[protein]["Reason_for_filtering"] = "Cellular localisation : DeepLoc"
    for prot_msa in need_msa :
        if prot_msa in new_possible_prey or prot_msa in possible_baits :
            new_need_msa.append(prot_msa)
    for prot_pkl in need_pkl :
        if prot_pkl in new_possible_prey or prot_pkl in possible_baits :
            new_need_pkl.append(prot_pkl)
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)
    return(new_need_msa, new_need_pkl)


def Generate_scripts (file, Informations_dict, Interaction_file, bait, GPU) :
    """
    Write one local script to use AlphaPullDown. This script should be written based on the maximum number of amino acids.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionnary
    Interaction_file : string
    bait : string
    GPU : list

    Returns:
    ----------
    """
    all_lines = str()
    OOM_int = str()
    result_dict = file.get_result_dict()
    regions = Informations_dict["Regions"]
    possible_prey = file.get_possible_prey()
    lenght_prot = file.get_lenght_prot()
    save_lenght_line = dict()

    #found max amino acids for your GPU
    pynvml.nvmlInit()
    device_count = pynvml.nvmlDeviceGetCount()
    for i in range(len(GPU)):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        name = pynvml.nvmlDeviceGetName(handle)
        mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        vram = (mem_info.total / 1024**2) * 0.001 # in MiB
    pynvml.nvmlShutdown()
    max_aa = int(vram * 120) #120 aa per Go of vram
    start = 0
    end = 0
    if Interaction_file == "PPI_int" :
        if regions[bait] != "0-0" :
            start = int(regions[bait].split("-")[0])
            end = int(regions[bait].split("-")[1])
            bait_file = f"{bait}_{start}-{end}"
            lenght = end - start + 1
        else :
            lenght = lenght_prot[bait]
            bait_file = bait
        #nbr_prey = 0
        for prey in possible_prey :
            int_lenght = lenght + lenght_prot[prey]
            #if nbr_prey > 100 : #take 100 int if filtred is not sufficient # need to verify 100 first ?
            #    break
            path1 = glob.glob(f"./result_PPI_int/{bait_file}_and_{prey}/ranked_0*")
            path2 = glob.glob(f"./result_PPI_int/{prey}_and_{bait_file}/ranked_0*")
            if len(path1) == 0 and len(path2) == 0 : #if model don't exist
                if int_lenght <= max_aa: #make interaction if doesn't exist and is not too long
                    if start == 0 and end == 0 :
                        save_lenght_line[f"{bait}_and_{prey}"] = [int_lenght, f"{bait};{prey}\n"]
                    else :
                        save_lenght_line[f"{bait}_and_{prey}"] = [int_lenght, f"{bait},{start}-{end};{prey}\n"]
                    #nbr_prey += 1
                else : #if interaction is too large
                    OOM_int += OOM_int + bait + ";" + prey + "\n"
                    result_dict[prey][f"iQ_score_vs_{bait}"] = "Too big interactions : AF OOM"
            else :
                #nbr_prey += 1
                pass

    if Interaction_file == "homo_int" :
        nbr_oligo = Informations_dict["Homo-oligomer"]
        nbr_prey = 0
        for prey in possible_prey :
            int_lenght = lenght_prot[prey] * int(nbr_oligo)
            if nbr_prey > 100 : #take 100 int if filtred is not sufficient # need to verify 100 first ?
                break
            path = glob.glob(f"./result_homo_int/{prey}_homo_{nbr_oligo}er/ranked_0*")
            if len(path) == 0 : #if model don't exist
                if int_lenght <= max_aa: #make interaction if doesn't exist and is not too long
                    save_lenght_line[f"{prey}_homer_{nbr_oligo}er"] = [int_lenght, f"{prey}:{nbr_oligo}\n"]
                    nbr_prey += 1
                else : #if interaction is too large
                    OOM_int += prey + ":" + nbr_oligo + "\n"
                    result_dict[prey]["Reason_for_filtering"] = "Homo-oligomer too large for your GPU"
            else :
                nbr_prey += 1
                pass
    dict_split_GPU = Split_to_GPU(file, save_lenght_line, GPU)
    for GPU_i in dict_split_GPU.keys() :
        with open(f"log_file/{Interaction_file}_{GPU_i}.txt",'w') as final_file :
            final_file.write(dict_split_GPU[GPU_i])
    with open("log_file/OOM_interactions.txt", "w") as OOM_file :
        OOM_file.write(OOM_int)
    file.set_result_dict(result_dict)

def Generate_3D_model (Informations_dict, Interaction_file, GPU) :
    """
    Use Alphapulldown script to generate 3D interactions.

    Parameters:
    ----------
    Informations_dict : dictionnary
    Interaction_file : string
    GPU : list

    Returns:
    ----------
    """
    Path_AlphaFold_Data = Informations_dict["Path_AlphaFold_Data"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]
    start_time = datetime.now()
    processes = []
    for gpu_id in GPU:
        p = multiprocessing.Process(
            target=run_AF_on_gpu,
            args=(gpu_id, Interaction_file, Path_AlphaFold_Data, Path_Pickle_Feature)
        )
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
    end_time = datetime.now()
    logger.info("Time APD : %s\n", end_time - start_time)

def run_AF_on_gpu(gpu_id, Interaction_file, Path_AlphaFold_Data, Path_Pickle_Feature):
    """
    Run the AlphaFold script on a specific GPU.

    Parameters:
    ----------
    gpu_id : int
    Interaction_file : string
    Path_AlphaFold_Data : string
    Path_Pickle_Feature : string

    Returns:
    ----------
    """
    #for H100, cluster ?
    # JAX-specific optimizations
    #export JAX_ENABLE_X64=0  # Use float32 to save memory
    #export JAX_DEFAULT_MATMUL_PRECISION="high"
    #export JAX_TRACEBACK_FILTERING=off  # Better debugging
    #export TF_FORCE_UNIFIED_MEMORY=1
    #export XLA_PYTHON_CLIENT_PREALLOCATE=false
    #export XLA_CLIENT_MEM_FRACTION=4.0  # Allow oversubscription
    env = os.environ.copy()
    env['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    env['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
    env['TF_FORCE_UNIFIED_MEMORY'] = 'true'
    env['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '3.2'
    env['XLA_FLAGS'] = '--xla_gpu_enable_triton_gemm=false'
    cmd =f"run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=False \--output_path=./result_{Interaction_file} \--data_dir={Path_AlphaFold_Data} \--protein_lists=log_file/{Interaction_file}_GPU_{gpu_id}.txt \--monomer_objects_dir={Path_Pickle_Feature} \--remove_keys_from_pickles=False"
    log_file = f"log_file/log_GPU_{gpu_id}.txt"
    run_cmd(cmd, env=env)

def Score_interaction_APD (file, Informations_dict, Interaction, bait=None) :
    """
    Generate scores for all interactions and set a list of possible prey based on the scores.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary
    Interaction : string

    Returns:
    ----------
    """
    start_time = datetime.now()
    result_dict = file.get_result_dict()
    possible_prey = file.get_possible_prey()
    int_score = file.get_int_score()
    new_possible_prey = list()
    Path_ccp4 = Informations_dict["Path_ccp4"]
    regions = Informations_dict["Regions"]
    ppi_list = list()

    if bait is not None : #setup bait name
        if regions[bait] != "0-0" :
            start = int(regions[bait].split("-")[0])
            end = int(regions[bait].split("-")[1])
            bait = f"{bait}_{start}-{end}"

    #Multiprocessing to score all interactions
    N_CPU = multiprocessing.cpu_count()
    if os.path.isdir(f"./result_{Interaction}") == True :
        if bait is None : #for homo-oligomer
            for direc in os.listdir(f"./result_{Interaction}") :
                if "_homo_" in direc :
                    ppi_list.append(f"./result_{Interaction}/{direc}") #Found a solution for score only homo-oligomer without score
        else : #for one vs all
            for protein in possible_prey : #create pdb for AF3
                if os.path.isfile(f"./result_{Interaction}/{bait}_and_{protein}/ranked_0.pdb") == False :
                    parser = MMCIFParser(QUIET=True)
                    structure = parser.get_structure("model", f"./result_{Interaction}/{bait}_and_{protein}/ranked_0_model.cif")
                    io = PDBIO()
                    io.set_structure(structure)
                    io.save(f"./result_{Interaction}/{bait}_and_{protein}/ranked_0.pdb")
                if f"iQ_score_vs_{bait}" not in int_score[protein].keys() :
                    ppi_list.append(f"./result_{Interaction}/{bait}_and_{protein}")
                else :
                    result_dict[protein][f"iQ_score_vs_{bait}"] = int_score[protein][f"iQ_score_vs_{bait}"]
                    if int_score[protein][f"iQ_score_vs_{bait}"] > 0 :
                        new_possible_prey.append(protein)

        results = []
        with multiprocessing.Pool(N_CPU) as pool :
            tasks = [(ppi, "./", Path_ccp4) for ppi in ppi_list]
            results_iter = pool.imap_unordered(run_scoring, tasks)
            for df in tqdm(results_iter, total=len(ppi_list), desc="Scoring interactions") :
                if df is not None and not df.empty :
                    results.append(df)
            pool.close()
            pool.join()
        if results :
            merged_df = pd.concat(results, ignore_index=True)
            merged_df.to_csv(os.path.join(f"./result_{Interaction}", "predictions_with_good_interpae.csv"), index=False)

        #Resume all score and set new possible prey
        with open(f"result_{Interaction}/predictions_with_good_interpae.csv", "r") as result_file :
            reader = csv.DictReader(result_file)

            #For one vs all
            if Interaction == "PPI_int" :
                all_lines = "jobs,pi_score,iptm_ptm,pDockQ,iQ_score\n"
                for row in reader :
                    job = row['jobs']
                    if '_and_' in job and job.split("_and_")[1] in possible_prey and bait in job : #check if interaction is a PPI and if prey is in possible prey list
                        if row['pi_score'] == 'No interface detected' :
                            iQ_score = float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30 #pi_score don't detect interface so it's set on -2.63
                            line =f'{row["jobs"]},-2.63,{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                        else :
                            iQ_score = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                            line =f'{row["jobs"]},{row["pi_score"]},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'

                        int_score[job.split("_and_")[1]][f"iQ_score_vs_{bait}"] = iQ_score
                        result_dict[job.split("_and_")[1]][f"iQ_score_vs_{bait}"] = iQ_score
                        new_possible_prey.append(job.split("_and_")[1])
                        all_lines = all_lines + line
                        name_int = job.split("/")[-1]
                        os.system(f"cp result_{Interaction}/{job}/ranked_0.pdb result_{Interaction}/{job}/{name_int}_ranked_0.pdb") #rename pdb file with explicit name
                for protein in possible_prey :
                    if protein not in new_possible_prey :
                        int_score[protein][f"iQ_score_vs_{bait}"] = 0
                        result_dict[protein][f"iQ_score_vs_{bait}"] = 0 #if prey don't have interaction, set iQ_score to 0
                        result_dict[protein]["Reason_for_filtering"] = f"Bad interactions with {bait} : PAE < 10"

            #For homo-oligomer
            if Interaction == "homo_int" : #score homo-oligomer
                all_homo = dict()
                save_pi_score = dict()
                all_lines = "jobs,pi_score,iptm_ptm,hiQ_score\n"
                for row in reader :
                    job = row['jobs']
                    if row['pi_score'] != 'No interface detected' :
                        if job not in all_homo.keys() :
                            all_homo[job] = (row['pi_score'],1,row)
                            save_pi_score[job] = [float(row['pi_score'])]
                        else :
                            save_pi_score[job].append(float(row['pi_score']))
                            sum_pi_score = float(all_homo[job][0]) + float(row['pi_score'])
                            sum_int = all_homo[job][1] + 1
                            all_homo[job] = (sum_pi_score,sum_int,row)
                for key in all_homo.keys() :
                    row = all_homo[key][2]
                    number_oligo = row["jobs"].split("_")[2].replace("er","") #AFPD 2.0.4
                    prot_name = row["jobs"].split("_")[0]
                    if len(save_pi_score[key]) > int(number_oligo) : #if model have more interface than number of homo-oligomerization
                        new_sum_pi_score = 0
                        save_pi_score[key].sort(reverse=True)
                        for index in range(0,int(number_oligo)) :
                            new_sum_pi_score += save_pi_score[key][index]
                            hiQ_score = (((float(new_sum_pi_score)/int(number_oligo))+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm are always same for each interface
                        line =f'{key},{str(float(new_sum_pi_score)/int(number_oligo))},{row["iptm_ptm"]},{str(hiQ_score)}\n'
                        all_lines += line   
                    else :
                        hiQ_score = (((float(all_homo[key][0])/all_homo[key][1])+2.63)/5.26)*60+float(row['iptm_ptm'])*40
                        line =f'{key},{str(float(all_homo[key][0])/all_homo[key][1])},{row["iptm_ptm"]},{str(hiQ_score)}\n'
                        all_lines += line
                    new_possible_prey.append(prot_name)
                    result_dict[key.split("_homo_")[0]]["hiQ_score"] = hiQ_score
                    name_int = key.split("/")[-1]
                    os.system(f"cp {key}/ranked_0.pdb {key}/{name_int}_ranked_0.pdb") #rename pdb file
                for protein in possible_prey :
                    if protein not in new_possible_prey :
                        result_dict[protein]["hiQ_score"] = 0
                        result_dict[protein]["Reason_for_filtering"] = "Bad homo-oligomer PAE : AF"

            if len(all_lines.strip("\n")) > 1 : #if all_lines is not empty
                with open(f"result_{Interaction}/new_predictions_with_good_interpae.csv", "w") as file2 :
                    file2.write(all_lines)
        file.set_result_dict(result_dict)
        file.set_possible_prey(new_possible_prey)
    else :
        logger.info(f"result_{Interaction}/ don't exist")
    end_time = datetime.now()
    logger.info("Time scoring interactions : %s\n", end_time - start_time)
    file.set_int_score(int_score)

def run_scoring(args) :
    """
    Run get_good_inter_pae script.

    Parameters:
    ----------
    args : set

    Returns:
    ----------
    result : string
    """
    interaction, output_dir, Path_ccp4 = args
    result = get_good_inter_pae.main(interaction, output_dir, 10, 2, Path_ccp4)
    return result

def Split_to_GPU(file, save_lenght_line, GPU) :
    """
    Split interactions in function of their lenght on the least loaded GPU.

    Parameters:
    ----------
    file : object of class File_proteins
    save_lenght_line : dictionnary
    GPU : list

    Returns:
    ----------
    dict_split_GPU : dictionnary
    """
    sorted_list = list(sorted(save_lenght_line.items(), key=lambda item: item[1][0], reverse = True)) #sorted in function of interaction lenght
    dict_split_GPU = dict()
    gpu_loads = dict()
    for nbr_GPU in GPU :
        dict_split_GPU[f"GPU_{nbr_GPU}"] = ""
        gpu_loads[f"GPU_{nbr_GPU}"] = 0
    for interactions in sorted_list :
        target_gpu = min(gpu_loads, key=gpu_loads.get)
        dict_split_GPU[target_gpu] += interactions[1][1]
        gpu_loads[target_gpu] += interactions[1][0]
    return dict_split_GPU

def Resume_file(file, Informations_dict) :
    """
    Create a resume file with all interactions, their scores and their status.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary

    Returns:
    ----------
    """
    iQ_score_dict = dict()
    list_name_baits = list()
    result_dict = file.get_result_dict()
    possible_prey = file.get_possible_prey()
    possible_baits = Informations_dict["Interact_with"]
    regions = Informations_dict["Regions"]
    proteins = file.get_proteins()
    informations = ["DeepLoc","Signal_peptide"]
    big_csv_lines = "Name,DeepLoc,Signal_peptide\n"
    small_csv_lines = "Name, Reason_for_filtering\n"
    if Informations_dict["Interact_with"] != [""] : #sorted in function of all baits
        for bait in possible_baits : #remove bait from result dict
            del result_dict[bait]
            if regions[bait] != "0-0" :
                start = int(regions[bait].split("-")[0])
                end = int(regions[bait].split("-")[1])
                bait = f"{bait}_{start}-{end}"
            informations.append(f"iQ_score_vs_{bait}")
            list_name_baits.append(bait)
            if len(possible_baits) == 1 :
                big_csv_lines = big_csv_lines.strip("\n") + ",iQ_score\n"
            else :
                big_csv_lines = big_csv_lines.strip("\n") + f",iQ_score_vs_{bait}\n"
        sorted_proteins = sorted(result_dict.items(),key=lambda x: sum(x[1].get(f"iQ_score_vs_{bait}", 0.0) for bait in list_name_baits), reverse=True) #sorted in function of all baits

    if Informations_dict["Homo-oligomer"] != "1" :
        informations.append("hiQ_score")
        big_csv_lines = big_csv_lines.strip("\n") + "hiQ_score\n"
    if Informations_dict["Interact_with"] == [""] and Informations_dict["Homo-oligomer"] != "1" : #if no PPI interactions filter on hiQ_score
        sorted_proteins = sorted(result_dict.items(),key=lambda x: (x[1].get("hiQ_score", 0), len(x[1]),), reverse=True)

    sorted_dict = dict(sorted_proteins)
    logger.info(sorted_dict)
    for prot in sorted_dict.keys() :
        small_csv_lines += prot
        big_csv_lines += prot
        for info in informations :
            if info in result_dict[prot].keys() :
                big_csv_lines += "," + str(result_dict[prot][info])
            else :
                big_csv_lines += ","
        if "Reason_for_filtering" in result_dict[prot].keys() :
            small_csv_lines += ", " + str(result_dict[prot]["Reason_for_filtering"])
        else :
            small_csv_lines += ", " + str("possible hit")
        big_csv_lines += "\n"
        small_csv_lines += "\n"
    with open("All_Final_result.csv", "w") as All_result_file :
        All_result_file.write(big_csv_lines)
    with open("Summary_result.csv", "w") as summary :
        summary.write(small_csv_lines)

@staticmethod
def run_cmd(cmd, env=None):
    """
    Run a shell command and handle KeyboardInterrupt to terminate the process group.

    Parameters:
    ----------
    cmd : string
    env : dict

    Returns:
    ----------
    """
    p = subprocess.Popen(cmd, shell=True, env=env, preexec_fn=os.setsid)
    try:
        stdout, stderr = p.communicate()
    except KeyboardInterrupt:
        logger.info("Interrupt detected, subprocess stopped…")
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)
        stdout, stderr = p.communicate()
