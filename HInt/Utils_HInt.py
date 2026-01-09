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
import signal
import glob
import pandas as pd
from pathlib import Path
from datetime import datetime
from numpy import load
from Bio import SeqIO
import time
import queue
from concurrent.futures import ThreadPoolExecutor, as_completed


# Configure global logger
logging.basicConfig(
    filename="HInt.log",  # Log file name
    level=logging.INFO,  # Log level
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

# Use the logger configuration from HInt.py
logger = logging.getLogger()


#check if U, O, B, Z, J, X is not in a sequence
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
    list_inf = ["Signal_peptide", "Homo-oligomer", "Interact_with", "Organism", "DeepLoc", "Regions", "Multimer_bait", "AlphaFold", "Max_protein_lenght", "Min_protein_lenght", "Path_Uniprot_ID", "Path_AlphaFold_Data", "Path_Pickle_Feature", "Path_Singularity_Image", "Path_MMseqs2_Data"]
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
            elif info in ["Signal_peptide","Homo-oligomer","Path_MMseqs2_Data","Regions","Multimer_bait","DeepLoc","AlphaFold","Max_protein_lenght","Min_protein_lenght"] :
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
            elif informations_key == "AlphaFold" :
                Informations_dict[informations_key] = "2"
                logger.info("Set AlphaFold version by default on AlphaFold2")
            elif informations_key == "Min_protein_lenght" :
                Informations_dict[informations_key] = "20"
                logger.info("Minimum lenght for prey protein set default 20 AA")
            elif informations_key == "Max_protein_lenght" :
                Informations_dict[informations_key] = ""
        if informations_key == "Interact_with" :
            regions_dict = dict()
            new_baits_list = list()
            if "[" in Informations_dict["Interact_with"] : #create multimer for bait
                list_complex = [prot.strip("[").strip(",").strip() for prot in Informations_dict["Interact_with"].split("]")]
                if "" in list_complex :
                    list_complex.remove("")
                nbr_baits = len(list_complex)
                list_baits = [prot.strip("[").strip("]") for prot in Informations_dict["Interact_with"].split(",")]
                Informations_dict["Multimer_bait"] = list_complex
            else :
                list_baits = [prot.strip(",").strip() for prot in Informations_dict["Interact_with"].split(",")]
                nbr_baits = len(list_baits)
                Informations_dict["Multimer_bait"] = list_baits
            if nbr_baits > 3 :
                raise ValueError("HInt don't support more than 3 differents baits")
            for prot in list_baits :
                if "-" in prot :
                    name_prot = prot.split("(")[0]
                    new_baits_list.append(name_prot.strip())
                    regions_dict[name_prot.strip()] = prot.split("(")[1].strip(")")
                    for i,multimer in enumerate(Informations_dict["Multimer_bait"]) :
                        if prot in multimer :
                            new_multimer = multimer.replace(prot,name_prot.strip())
                            Informations_dict["Multimer_bait"][i] = new_multimer
                else :
                    if prot.strip() in regions_dict.keys() and regions_dict[prot.strip()] != "0-0" :
                        raise ValueError(f"HInt don't support multiple regions for protein bait : {prot.strip()}")
                    regions_dict[prot.strip()] = "0-0"
                    new_baits_list.append(prot.strip())
            Informations_dict["Regions"] = regions_dict
            Informations_dict["Interact_with"] = new_baits_list
        if informations_key == "Homo-oligomer" :
            if int(Informations_dict[informations_key]) == False :
                raise ValueError(f"Homo-oligomer is not an integer")
            if Informations_dict[informations_key] == "0" :
                Informations_dict[informations_key] = "1"
        if informations_key == "DeepLoc" :
            if Informations_dict["Organism"] == "euk" : #euk
                for value in Informations_dict[informations_key].split(","):
                    if value.strip()  not in ["Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion", "Plastid", "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgo apparatus", "Peroxisome","None"] :
                        raise ValueError(f"Incorrect DeepLoc value : {value}")
            else : #other
                for value in Informations_dict[informations_key].split(","):
                    if value.strip()  not in ["Cell wall & surface","Extracellular","Cytoplasmic","Cytoplasmic Membrane","Outer Membrane","Periplasmic","None"] :
                        raise ValueError(f"Incorrect DeepLocPro value : {value}")
    if len(Informations_dict["Signal_peptide"]) == 0 and len(Informations_dict["DeepLoc"]) == 0 and len(Informations_dict["Homo-oligomer"]) == 0 and len(Informations_dict["Interact_with"]) == 0 : #no info
        logger.info("Need information to discriminate the potential homolog/interolog")
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


    if os.path.exists("log_file/result_deeploc") == True :
        os.system("rm -r log_file/result_deeploc")
    dp_lines = str()
    for protein in need_DeepLoc :
        dp_lines += ">"+protein+"\n"+prot_seq[protein]+"\n"
    with open(f"log_file/{fasta_file}", "w") as dp_file :
        dp_file.write(dp_lines)

    if org == "euk" :
        logger.info(str(datetime.now())+" Start DeepLoc eucaryote")
        software = "deeploc2"
        cmd = f"CUDA_VISIBLE_DEVICES={GPU_str} {software} -f log_file/{fasta_file} -o log_file/result_deeploc -d cuda"
    else :
        logger.info(str(datetime.now())+" Start DeepLocPro")
        software = "deeplocpro"
        if org == "gram-" :
            group = "negative"
        if org == "gram+" :
            group = "positive"
        if org == "arch" :
            group = "archaea"
        cmd = f"CUDA_VISIBLE_DEVICES={GPU_str} {software} -f log_file/{fasta_file} -o log_file/result_deeploc -d cuda -g {group}" #deeplocpro don't use all GPU
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



def run_SP (file, Informations_dict, need_SP, need_msa) :
    """
    Create a new FASTA file without the signal peptide using SignalP.
    If protein have already msa but no sequence without SP don't return it.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dict
    need_SP : list #proteins needing SP removal and MSA
    need_msa : list #proteins needing MSA

    Returns:
    ----------
    new_need_msa : list
    """
    final_file = str()
    SP_signal = 0
    prot_SP = dict()
    file_name = file.get_file_name()
    fasta_file = file_name.replace(".txt","_msa.fasta")
    output_file = fasta_file.replace(".fasta","")
    new_fasta_dict = file.get_proteins_sequence_no_SP()

    file.create_fasta_file(True, need_SP)
    cmd1 = f"signalp -fasta log_file/{fasta_file} -org {Informations_dict['Organism']} -prefix log_file/{output_file}"
    os.system(cmd1)

    file_signalp = fasta_file.replace(".fasta","_summary.signalp5")
    with open(f"log_file/{file_signalp}","r") as fh :
        for line in fh :
            new_line = line.split("\t")
            if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                prot_SP[new_line[0]] = new_line[len(new_line)-1].split("-")[1].split(".")[0]
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
    new_need_msa = list()
    for prot in need_msa : 
        if os.path.isfile(f"{Informations_dict['Path_Pickle_Feature']}/{prot}.a3m") == False :
            new_need_msa.append(prot)
    file.set_proteins_sequence_no_SP(new_fasta_dict)
    return new_need_msa

def create_feature (file, Informations_dict, GPU, CPU, need_msa, need_pkl) :
    """
    Launch command to generate MSA and feature pickle.

    Parameters:
    ----------
    file : object of class File_proteins
    informations_dict : dict
    GPU : list
    CPU : intS
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
    #logger.info(f"GPU use : {GPU_str}")
    
    logger.info(f"{len(need_msa)} proteins need MSA")
    logger.info(f"{len(need_pkl)} proteins need pkl files")
    
    if os.path.exists(Path_Pickle_Feature) == False :
        os.system(f"mkdir {Path_Pickle_Feature}")

    #Look for MSA in AlphaFold database
    generated_msa = copy.deepcopy(need_msa)
    logger.info(f"Search MSA in AlphaFold database")

    # ThreadPool to parallelise
    cpu_per_mafft = 2
    max_workers = max(1, CPU // cpu_per_mafft)  # how many mafft in parallel


    start = time.time()

    futures_list = []
    found = []
    not_found = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor : #CPU parallelization
        for protein in generated_msa : #check if prot have an MSA in alphafold database
            l_p = len(protein)
            if l_p >= 5 and l_p <= 10 and "_" not in protein : #if not, is not an UniprotID
                futures_list.append(executor.submit(fetch_trim_mafft,protein,Path_Pickle_Feature,prot_SP,prot_no_SP,))

        for future in as_completed(futures_list) :
            protein, ok = future.result()
            if ok :
                logger.info(f"MSA for {protein} processed")
                need_msa.remove(protein) #msa found
                need_pkl.append(protein)
                    
    logger.info("MSA search in AlphaFold database complete")

    file.create_fasta_file(False, need_msa, need_pkl)

    #Create MSA files with ColabFold mmseq2 GPU accelerated for proteins without MSA
    if len(need_msa) > 10 :
        if len(GPU_str.split(",")) >= 4 : #Due to error by using mmseqGPU with more than 3 GPU 
            GPU_str = GPU_str[:-2]
        cmd = f"CUDA_VISIBLE_DEVICES={GPU_str} colabfold_search ./log_file/{msa_name} {Path_MMseqs2_Data} {Path_Pickle_Feature} --db-load-mode 1 --gpu 1"  #-e 0.1
        os.system(cmd)

    if len(need_msa) < 1 :
        logger.info("All MSAs have already been generated")

    #Create MSA files with ColabFold mmseq2 classic pipeline for proteins without MSA, and pickle files for proteins generated with MMseqs2 GPU
    if len(need_msa) > 0 and len(need_msa) <= 10 :
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

    end = time.time()
    elapsed = end - start
    logger.info("Create MSA take "+ str(elapsed/60)+" minutes")

    #Create pkl files for proteins without pkl file (from MSA found in AFdb or error during MSA generation with ColabFold mmseqs2)
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
        process = subprocess.Popen(cmd2, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
        stdout, stderr = process.communicate()
            
def fetch_trim_mafft(protein, Path_Pickle_Feature, prot_SP, prot_no_SP) :
    """
    Found MSA on AFdb ,if signal peptide cut it an re align with mafft.

    Parameters:
    ----------
    protein : str
    Path_Pickle_Feature : str
    prot_SP : dict
    prot_no_SP : dict

    Returns:
    ----------
    protein : str
    bool

    """
    url = f"https://alphafold.ebi.ac.uk/files/msa/AF-{protein}-F1-msa_v6.a3m"
    msa_in = f"{Path_Pickle_Feature}/{protein}.a3m"
    aln_out = f"{Path_Pickle_Feature}/{protein}.aln"


    check = subprocess.run(["wget", "--spider", "-q", url])
    if check.returncode != 0 : #normaly 0
        return protein, False  # MSA not found


    subprocess.run(["wget", "-q", "-O", msa_in, url], check=True)


    SP = len(prot_SP[protein]) - len(prot_no_SP[protein])
    if SP > 0 : #if no SP don't modify the MSA
        trimmed_records = []
        for rec in SeqIO.parse(msa_in, "fasta") :
            new_seq = rec.seq[SP:]
            len_aa_count = sum(1 for c in new_seq if c.isupper())
            if len_aa_count < 10 or not any(c.isupper() for c in str(new_seq)) : #skip short sequence after SP cut or sequence don't have aa
                continue
            new_rec = rec[:]
            new_rec.seq = new_seq
            trimmed_records.append(new_rec)

        if trimmed_records :
            with open(msa_in, "w") as msa_file :
                for rec in trimmed_records:
                    msa_file.write(f">{rec.description}\n{rec.seq}\n")

        #realign with mafft
        cmd_mafft = f"mafft --quiet --anysymbol --thread 1 --parttree --retree 1 --maxiterate 0 {msa_in} > {aln_out}" #monothread
        subprocess.run(cmd_mafft, shell=True, check=True)

        #reformat a3m
        cmd_reformat = f"reformat.pl fas a3m {aln_out} {msa_in}"
        subprocess.run(cmd_reformat, shell=True, check=True)
        cmd_rm_f = f"rm {aln_out}"
        subprocess.run(cmd_rm_f, shell=True, check=True)

        #delete \n
        all_lines = ""
        with open(msa_in, "r") as in_a3m:
            seq, header = "", None
            for line in in_a3m:
                if line.startswith(">"):
                    if header:
                        all_lines += f"{header}\n{seq}\n"
                    header = line.strip()
                    seq = ""
                else:
                    seq += line.strip()
            if header:
                all_lines += f"{header}\n{seq}\n"
        with open(msa_in, "w") as out_a3m:
            out_a3m.write(all_lines)

    return protein, True



def filter_signalP(file, Informations_dict, need_msa, need_pkl) :
    """
    Filter proteins based on the presence of a signal peptide using SignalP results, if not specified keep all proteins and just describe them in the result file.

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
    if SignalP != "None" :
        logger.info("Protein preys remaining after SignalP filtering : " + str(len(new_possible_prey)))
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)
    return need_msa, need_pkl

def Make_all_MSA_coverage(file, Path_Pickle_Feature) :
    """
    Generating MSA coverage for all proteins and write shallow_MSA text file.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_Pickle_Feature : string

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


def filter_lenght(file, Informations_dict, need_msa, need_pkl, need_DeepLoc) :
    """
    Filter proteins based on their lenght and set new prey list.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dict
    need_msa : list
    need_pkl : list
    need_DeepLoc : list

    Returns:
    ----------
    """
    result_dict = file.get_result_dict()
    sequence_dict = file.get_proteins_sequence_SP() #use sequence with SP for lenght filtering
    possible_prey = file.get_possible_prey()
    if Informations_dict["Max_protein_lenght"] == "" : #no set by default, depend of GPU memory
        max_lenght = 100000
    else :
        max_lenght = int(Informations_dict["Max_protein_lenght"])
    min_lenght = int(Informations_dict["Min_protein_lenght"]) #set by default to 30
    new_possible_prey = list()
    for protein in possible_prey :
        if len(sequence_dict[protein]) < max_lenght and len(sequence_dict[protein]) > min_lenght :
            new_possible_prey.append(protein)
        else :
            result_dict[protein]["Reason_for_filtering"] = "Lenght filtering"
            if protein in need_msa :
                need_msa.remove(protein)
            if protein in need_pkl :
                need_pkl.remove(protein)
            if protein in need_DeepLoc :
                need_DeepLoc.remove(protein)
    logger.info("Protein preys remaining after lenght filtering : " + str(len(new_possible_prey)))
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)
    return(need_msa, need_pkl, need_DeepLoc)


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
    logger.info("Protein prey remaining after DeepLoc filtering : " + str(len(new_possible_prey)))
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)
    return(new_need_msa, new_need_pkl)


def Generate_scripts(file, Informations_dict, Interaction_file, bait):
    """
    Prepare the list of interaction jobs for AlphaPullDown based on protein length
    and GPU VRAM constraints.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dict
    Interaction_file : str
    bait : str

    Returns:
    ----------
    job_list : list of str
    """
    job_list = []
    OOM_int = ""
    result_dict = file.get_result_dict()
    AF_version = Informations_dict["AlphaFold"]
    possible_prey = file.get_possible_prey()
    lenght_prot = file.get_lenght_prot()
    regions = Informations_dict["Regions"]

    # Estimate max amino acids based on GPU VRAM (choose first GPU as reference)
    pynvml.nvmlInit()
    handle = pynvml.nvmlDeviceGetHandleByIndex(0)
    mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
    vram = (mem_info.total / 1024**2) * 0.001  # GiB
    max_aa = int(vram * 120)  # 120 AA per GiB VRAM
    pynvml.nvmlShutdown()

    # Determine bait length and region
    complexe = "," in bait
    if complexe :
        save_multimer = bait
        bait_file = save_multimer.replace(",", "_and_")
        lenght = sum(lenght_prot[prot] for prot in save_multimer.split(","))
        for prot in save_multimer.split(",") :
            if regions[prot] != "0-0":
                start, end = int(regions[prot].split("-")[0]), int(regions[prot].split("-")[1])
                bait_file = bait_file.replace(prot, f"{prot}_{start}-{end}")
                bait_for_job = bait_for_job.replace(prot, f"{prot},{start}-{end}")
            else :
                bait_for_job = save_multimer.replace(",", ";")
    else :
        lenght = lenght_prot[bait]
        if regions[bait] != "0-0" :
            start, end = int(regions[bait].split("-")[0]), int(regions[bait].split("-")[1])
            bait_file = f"{bait}_{start}-{end}"
            bait_for_job = f"{bait},{start}-{end}"
        else: 
            bait_file = bait
            bait_for_job = bait

    # Build job list
    if Interaction_file == "PPI_int" :
        for prey in possible_prey :
            int_lenght = lenght + lenght_prot[prey]

            # Check if model already exists
            if AF_version == "3":
                path1 = glob.glob(f"./result_PPI_int/{bait_file}_and_{prey}/ranked_0_model.cif")
                path2 = glob.glob(f"./result_PPI_int/{prey}_and_{bait_file}/ranked_0_model.cif")
            else:  # AF_version == "2"
                path1 = glob.glob(f"./result_PPI_int/{bait_file}_and_{prey}/ranked_0.pdb")
                path2 = glob.glob(f"./result_PPI_int/{prey}_and_{bait_file}/ranked_0.pdb")

            if len(path1) == 0 and len(path2) == 0 :
                if int_lenght <= max_aa:
                    job_str = f"{bait_for_job};{prey}\n"
                    job_list.append(job_str)
                else :
                    OOM_int += f"{bait_for_job};{prey}\n"
                    result_dict[prey][f"iQ_score_vs_{bait}"] = "Too big interactions: AF OOM"

    elif Interaction_file == "homo_int" :
        nbr_oligo = Informations_dict.get("Homo-oligomer", 2)
        for prey in possible_prey :
            int_lenght = lenght_prot[prey] * int(nbr_oligo)
            path = glob.glob(f"./result_homo_int/{prey}_homo_{nbr_oligo}er/ranked_0*")
            if len(path) == 0 :
                if int_lenght <= max_aa :
                    job_str = f"{prey}:{nbr_oligo}\n"
                    job_list.append(job_str)
                else :
                    OOM_int += f"{prey}:{nbr_oligo}\n"
                    result_dict[prey]["Reason_for_filtering"] = "Homo-oligomer too large for your GPU"


    # Save OOM interactions
    with open("log_file/OOM_interactions.txt", "w") as OOM_file :
        OOM_file.write(OOM_int)

    file.set_result_dict(result_dict)

    return job_list


def Generate_3D_model(Informations_dict, interaction_type, job_list, GPU) :
    """
    Genrerate 3D models using multiple GPUs and multiprocessing.

    Parameters:
    ----------
    Informations_dict : dict
    interaction_files : str
    job_list : list
    GPU : list
    """
    AF_version = Informations_dict["AlphaFold"]
    Path_AlphaFold_Data = Informations_dict["Path_AlphaFold_Data"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]

    stop_flag = multiprocessing.Event()
    job_queue = multiprocessing.Queue()

    for interaction in job_list :
        job_queue.put(interaction)

    processes = []
    start_time = datetime.now()

    try:
        for gpu_id in GPU:
            p = multiprocessing.Process(target=gpu_worker, args=(gpu_id, job_queue,Path_AlphaFold_Data, Path_Pickle_Feature, interaction_type, AF_version, stop_flag))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

    except KeyboardInterrupt:
        logger.warning("Ctrl+C detected — stopping all GPU workers")
        stop_flag.set() 
        for p in processes:
            if p.is_alive():
                p.terminate()

        for p in processes:
            p.join(timeout=5)

        sys.exit(1)

    logger.info("Time modelisation : %s", datetime.now() - start_time)


def gpu_worker(gpu_id, job_queue, Path_AlphaFold_Data, Path_Pickle_Feature, interaction_type, AF_version, stop_flag) :
    """
    Run the AlphaFold script on a specific GPU, allows parallel processing of jobs on multiple GPUs.
    A new jobs is fetched from the job_queue when the current one is finished.

    Parameters:
    ----------
    gpu_id : str
    job_queue : list
    Path_AlphaFold_Data : str
    Path_Pickle_Feature : str
    interaction_type : str
    AF_version : str
    stop_flag : multiprocessing.Event
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


    stop_flag = multiprocessing.Event()

    while not stop_flag.is_set() :
        try :
            interaction_file = job_queue.get(timeout=5)
        except queue.Empty :
            break  #all jobs are done

        logger.info(f"[GPU {gpu_id}] Starting {interaction_file}")
        with open(f"log_file/{interaction_type}_GPU_{gpu_id}.txt", "w") as int_file :
            int_file.write(interaction_file)

        if AF_version == "2" :
            cmd = (f"run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_{interaction_type} \--data_dir={Path_AlphaFold_Data} \--protein_lists=log_file/{interaction_type}_GPU_{gpu_id}.txt \--monomer_objects_dir={Path_Pickle_Feature} \--remove_keys_from_pickles=False")


        elif AF_version == "3" :
            cmd = (
                "run_multimer_jobs.py --mode=custom "
                f"--output_path=./result_{interaction_type} "
                f"--data_dir={Path_AlphaFold_Data} "
                f"--protein_lists=log_file/{interaction_type}_GPU_{gpu_id}.txt "
                f"--monomer_objects_dir={Path_Pickle_Feature} "
                "--fold_backend=alphafold3"
            )


        current_process = subprocess.Popen(cmd,shell=True, env=env, preexec_fn=os.setsid)
        try :
            while True :
                if stop_flag.is_set() :
                    os.killpg(os.getpgid(current_process.pid), signal.SIGTERM)
                    logger.warning(f"[GPU {gpu_id}] Interrupted — killing AlphaFold")
                    return
                retcode = current_process.poll()
                if retcode is not None :
                    break
                time.sleep(0.5)


        except KeyboardInterrupt :
            os.killpg(os.getpgid(current_process.pid), signal.SIGTERM)
            logger.warning(f"[GPU {gpu_id}] Ctrl+C detected — killing AlphaFold")
            return

        cmd_rm = f"rm log_file/{interaction_type}_GPU_{gpu_id}.txt"
        os.system(cmd_rm)
        logger.info(f"[GPU {gpu_id}] Finished {interaction_file}")




@staticmethod
def run_cmd(cmd, env=None):
    """
    Run a shell command and handle KeyboardInterrupt to terminate the process group.

    Parameters:
    ----------
    cmd : str
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
