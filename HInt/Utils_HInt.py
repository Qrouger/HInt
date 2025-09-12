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
from pathlib import Path
from datetime import datetime
from numpy import load
from queue import Queue


def Define_informations() :
    """
    Extract all paths from HInt.txt and store them in a dictionary.
    
    Parameters:
    ----------

    Returns:
    ----------
    Informations_dict : dict
    """
    logging.basicConfig(level=logging.INFO,format="%(asctime)s - %(levelname)s - %(message)s",stream=sys.stdout)
    Informations_dict = dict()
    list_inf = ["Signal_peptide", "Homo-oligomer", "Interact_with", "Organism", "DeepLoc", "Regions", "Path_Uniprot_ID", "Path_AlphaFold_Data", "Path_Pickle_Feature", "Path_Singularity_Image", "Path_MMseqs2_Data"]
    with open("HInt.txt", "r") as file :
        for lines in file :
            if ":" in lines :
                informations_name = lines.split(":")[0].strip().strip("\n")
                informations = lines.split(":")[1].strip().strip("\n")
                Informations_dict[informations_name] = informations
    for info in list_inf :
        if info not in Informations_dict.keys() :
            if info in ["Interact_with", "Organism","Path_Uniprot_ID", "Path_AlphaFold_Data", "Path_Pickle_Feature", "Path_Singularity_Image"] :
                print(f"HInt.txt file is compromised, verify the file. {info} is missing")
                break
            elif info in ["Path_MMseqs2_Data","Regions","DeepLoc"] :
                Informations_dict[info] = ""
    for informations_key in Informations_dict.keys() :
        if len(Informations_dict[informations_key]) == 0 :
            logging.info(f'Informations : {informations_key} is empty')
            if informations_key == "Path_Uniprot_ID" or informations_key == "Path_Singularity_Image" :
                exit()
            elif informations_key == "Path_AlphaFold_Data" :
                logging.info("set by default on ./alphadata")
                Informations_dict[informations_key] = "./alphadata"
            elif informations_key == "Path_Pickle_Feature" :
                logging.info("set by default on ./feature")
                Informations_dict[informations_key] = "./feature"
            elif informations_key == "Signal_peptide" :
                Informations_dict[informations_key] = "None"
            elif informations_key == "Homo-oligomer" :
                Informations_dict[informations_key] = 1
            elif informations_key == "DeepLoc" :
                Informations_dict[informations_key] = "None"
            elif informations_key == "Organism" :
                exit()
            elif informations_key == "Path_MMseqs2_Data" :
                logging.info("/!\ local MMseqs2 GPU will not be used")
        if informations_key == "Interact_with" :
            regions_dict = dict()
            new_baits_list = list()
            list_baits = [prot for prot in Informations_dict["Interact_with"].split(",")]
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
    if Informations_dict["DeepLoc"] != "None" and Informations_dict["Organism"] != "euk" :
        logging.info("DeepLoc is only available for eukaryotes")
        Informations_dict["DeepLoc"] = "None"
    if len(Informations_dict["Signal_peptide"]) == 0 and len(Informations_dict["Homo-oligomer"]) == 0 and len(Informations_dict["Interact_with"]) == 0 : #no info
        logging.info("need information to discriminate the potential homologue")
        exit()
    return(Informations_dict)

def remove_SP (file, Informations_dict, need_msa) :
    """
    Create a new FASTA file without the signal peptide using SignalP and filtred signal peptide.

    Parameters:
    ----------
    file : object of class File_proteins
    org : string

    Returns:
    ----------
    """
    final_file = str()
    SP_signal = 0
    prot_SP = dict()
    Prot_Signal_string = str()
    fasta_file = file.get_fasta_file()
    cmd1 = "signalp -fasta " + fasta_file + " -org " + Informations_dict["Organism"]
    os.system(cmd1)
    SignalP = Informations_dict["Signal_peptide"]
    file_signalp = fasta_file.replace(".fasta","_summary.signalp5")
    with open(file_signalp,"r") as fh :
        for line in fh :
            new_line = line.split("\t")
            if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                prot_SP[new_line[0]] = new_line[len(new_line)-1].split("-")[1].split(".")[0]
    new_fasta_dict = file.get_proteins_sequence()
    with open(fasta_file, "r") as fa_file :
        for line2 in fa_file :
            new_line2 = line2
            if SP_signal == 0 and line2[0] != ">" :
                new_line2 = line2
                new_fasta_dict[save_key] = line2.strip("\n")
            if int(SP_signal) > 0 :
                new_line2 = line2[int(SP_signal)-1:len(line2)]
                new_fasta_dict[save_key] = line2[int(SP_signal)+1:len(line2)].strip("\n")
                SP_signal = 0
            if line2[0] == ">" :
                save_key = line2[1:len(line2)-1]
                if str(line2[1:len(line2)-1]) in prot_SP.keys() :
                    SP_signal = prot_SP[line2[1:len(line2)-1]]
            final_file = final_file + new_line2
    fasta_lines = str()
    for proteins in new_fasta_dict.keys() :
        if SignalP == "Yes" :
            if proteins in prot_SP.keys() and proteins in need_msa :
                fasta_lines += ">"+proteins+"\n"+new_fasta_dict[proteins]+"\n"
            else :
                pass
        if SignalP == "No" :
            if proteins not in prot_SP.keys() and proteins in need_msa :
                fasta_lines += ">"+proteins+"\n"+new_fasta_dict[proteins]+"\n"
            else :
                pass
        if SignalP == "None" :
            if proteins in need_msa :
                fasta_lines += ">"+proteins+"\n"+new_fasta_dict[proteins]+"\n"
    baits = Informations_dict["Interact_with"]
    for bait in baits :
        if SignalP == "Yes" and bait not in prot_SP.keys() and bait in need_msa :
            fasta_lines += ">"+bait+"\n"+new_fasta_dict[bait]+"\n" #add baits to the new fasta file
        if SignalP == "No" and bait in prot_SP.keys() and bait in need_msa :
            fasta_lines += ">"+bait+"\n"+new_fasta_dict[bait]+"\n" #add baits to the new fasta file
    with open(fasta_file, "w") as new_file2 :
        new_file2.write(fasta_lines)
    file.set_proteins_sequence(new_fasta_dict)

def create_feature (file, Informations_dict, GPU) :
    """
    Launch command to generate features.

    Parameters:
    ----------
    file : object of class File_proteins
    informations_dict : dict
    GPU : list

    Returns:
    ----------
    """
    Path_AlphaFold_Data = Informations_dict["Path_AlphaFold_Data"]
    Path_Pickle_Feature = Informations_dict["Path_Pickle_Feature"]
    baits = Informations_dict["Interact_with"]
    regions = Informations_dict["Regions"]
    Path_MMseqs2_Data = Informations_dict["Path_MMseqs2_Data"]
    fasta_file = file.get_fasta_file()
    a3m_file = open(fasta_file, 'r')
    nb_line = 0
    GPU_str = ""
    for nbr_GPU in GPU :
        GPU_str += nbr_GPU + ","
    GPU_str = GPU_str.strip(",")
    print(GPU_str)
    for line in a3m_file:
        nb_line += 1
    print(f"Number of sequences in {fasta_file} : {int(nb_line/2)}")
    if nb_line > 50 and Path_MMseqs2_Data!="" : #if more than 25 sequences, use local colabfold_search GPU
        cmd = f"CUDA_VISIBLE_DEVICES={GPU_str} colabfold_search {fasta_file} {Path_MMseqs2_Data} {Path_Pickle_Feature} --gpu 1 --db-load-mode 2" #-e 0.1
        os.system(cmd)
    if nb_line < 1 :
        logging.info("All MSAs have already been generated")
    cmd = ["create_individual_features.py",
    f"--fasta_paths=./{fasta_file}",
    f"--data_dir={Path_AlphaFold_Data}",
    "--save_msa_files=False",
    f"--output_dir={Path_Pickle_Feature}",
    "--max_template_date=2024-05-02",
    "--skip_existing=True",
    "--use_mmseqs2=True",
    "--use_precomputed_msas=True"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
    for line in process.stdout:
        print(line, end="")
    process.stdout.close()
    process.wait()
    all_lines = str()


    for bait in baits :
        start = int(regions[bait].split("-")[0]) - 1
        end = int(regions[bait].split("-")[1]) - 1
        bait_name_file = f"{bait}_{str(start+1)}_{str(end+1)}"
        if regions[bait] != "0-0" and f"{bait_name_file}.a3m" not in os.listdir(Path_Pickle_Feature) : #if MSA regions of bait not exist
            with open(f"{Path_Pickle_Feature}/{bait}.a3m", "r") as a3m_file :
                for line in a3m_file:
                    if line[0] == ">" :
                        mem_name = line.strip().split("\t")[0] + "\n"
                    else :
                        if not all(cara == '-' for cara in line[start:end]) :
                            all_lines += mem_name + line[start:end] + "\n"
            with open(f"{Path_Pickle_Feature}/{bait_name_file}.a3m", "w") as new_file :
                new_file.write(all_lines)
            cmd1 = f"mafft --anysymbol {Path_Pickle_Feature}/{bait_name_file}.a3m > {Path_Pickle_Feature}/{bait_name_file}.aln"
            cmd2 = f"reformat.pl fas a3m {Path_Pickle_Feature}/{bait_name_file}.aln {Path_Pickle_Feature}/{bait_name_file}.a3m"
            os.system(cmd1)
            os.system(cmd2)
            with open(f"{Path_Pickle_Feature}/{bait_name_file}.a3m", "r") as a3m2_file :
                all_lines = str()
                all_seq = str()
                for line in a3m2_file:
                    if line[0] == ">" :
                        if all_seq != "" :
                            all_lines += mem_name + all_seq + "\n"
                        mem_name = line.strip().split("\t")[0] + "\n"
                        all_seq = str()
                    else :
                        all_seq += line.strip()
            with open(f"{Path_Pickle_Feature}/{bait_name_file}.a3m", "w") as a3m2_file :
                a3m2_file.write(all_lines)
            cmd3 = f"rm {Path_Pickle_Feature}/{bait_name_file}.aln"
            os.system(cmd3)


def Make_all_MSA_coverage (file, Path_Pickle_Feature) : #relativement long pour parse toutes les prot #need to make this just with .a3m files
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
#    for prot in possible_prey : #just write shallow_MSA.txt
        a3m_file = open(f'{Path_Pickle_Feature}/{prot}.a3m', 'r')
        msa = 0
        for line in a3m_file:
            msa += 1
#        pre_feature_dict = pickle.load(open(f'{Path_Pickle_Feature}/{prot}.pkl','rb'))
#        feature_dict = pre_feature_dict.feature_dict
#        msa = feature_dict['msa']
        if msa/2 <= 100 : #name and sequence
            if msa <= 2 :
               shallow_MSA += prot + " : " + str(int(msa/2)) + " sequences\n"
               result_dict[prot] = result_dict[prot] + "No MSA"
            else :
               shallow_MSA += prot + " : " + str(int(msa/2)) + " sequences\n"
               result_dict[prot] = result_dict[prot] + "Shallow MSA | "
               new_possible_prey.append(prot)
        else :
            new_possible_prey.append(prot)
    with open("shallow_MSA.txt", "w") as MSA_file :
        MSA_file.write(shallow_MSA)
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)

def filter_signalP(file, Informations_dict) :
    """
    Filter proteins based on the presence of a signal peptide using SignalP results.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionnary

    Returns:
    ----------
    """
    result_dict = file.get_result_dict()
    seq_dict = file.get_proteins_sequence()
    possible_prey = file.get_possible_prey()
    new_possible_prey = list()
    SignalP = Informations_dict["Signal_peptide"]
    for protein in possible_prey :
        if seq_dict[protein][0] == "M" :
            if SignalP == "Yes" :
                result_dict[protein] = "Don't have a signal peptide"
            elif SignalP == "No" :
                new_possible_prey.append(protein)
            else :
                result_dict[protein] = "Don't have a signal peptide | "
                new_possible_prey.append(protein)
        else :
            if SignalP == "Yes" :
                new_possible_prey.append(protein)
            elif SignalP == "No" :
                result_dict[protein] = "Have a signal peptide"
            else :
                result_dict[protein] = "Have a signal peptide | "
                new_possible_prey.append(protein)
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)

def filter_deeploc(file, localisation) :
    """
    Filter proteins based on cellular localisation.

    Parameters:
    ----------
    file : object of class File_proteins
    localisation : string

    Returns:
    ----------
    """
    localisation = localisation.split(",")
    new_possible_prey = list()
    fasta_file = file.get_fasta_file()
    sequence_dict = file.get_proteins_sequence()
    possible_prey = file.get_possible_prey()
    all_seq = str()
    for prot in possible_prey :
        all_seq += ">" + prot + "\n" +sequence_dict[prot] + "\n"
    with open(f"{fasta_file}", "w") as fasta :
        fasta.write(all_seq)
    if os.path.exists("result_deeploc") == True :
        os.system("rm -r result_deeploc")
    cmd = f"deeploc2 -f {fasta_file} -o ./result_deeploc"
    os.system(cmd)
    file_result = os.listdir("./result_deeploc")[0]
    with open(f"result_deeploc/{file_result}", "r") as result :
        reader = csv.reader(result, delimiter=',')
        for line in reader :
            for cell_loc in localisation :
                if cell_loc.strip() in line[1] :
                    new_possible_prey.append(line[0])
    file.set_possible_prey(new_possible_prey)


def Generate_scripts (file, Informations_dict, Interaction_file, GPU) :
    """
    Write one local script to use AlphaPullDown. This script should be written based on the maximum number of amino acids.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionnary
    Interaction_file : string
    GPU : list

    Returns:
    ----------
    """
    all_lines = str()
    OOM_int = str()
    result_dict = file.get_result_dict()
    possible_baits = Informations_dict["Interact_with"]
    possible_prey = file.get_possible_prey()
    lenght_prot = file.get_lenght_prot()
    save_lenght_line = dict()
    max_aa = 2400 #need to bench
    if Interaction_file == "APD_PPI_int" :
        for bait in possible_baits :
            nbr_prey = 0
            lenght = lenght_prot[bait]
            for prey in possible_prey :
                int_lenght = lenght + lenght_prot[prey]
                if nbr_prey > 30 :
                    break
                if os.path.exists(f"./result_APD_PPI_int/{bait}_and_{prey}/ranked_0.pdb") == False and os.path.exists(f"./result_APD_PPI_int/{prey}_and_{bait}/ranked_0.pdb") == False : #if model don't exist
                    if int_lenght <= max_aa: #make interaction if doesn't exist and is not too long
                        save_lenght_line[f"{bait}_and_{prey}"] = [int_lenght, f"{bait};{prey}\n"]
                        nbr_prey += 1
                    else : #if interaction is too large
                        OOM_int = OOM_int + bait + ";" + prey + "\n"
                        result_dict[prey] = result_dict[prey] + " Too big interactions : AF OOM"
                else :
                    nbr_prey += 1
                    pass
    if Interaction_file == "APD_homo_int" :
        nbr_oligo = Informations_dict["Homo-oligomer"]
        nbr_prey = 0
        for prey in possible_prey :
            int_lenght = lenght_prot[prey] * int(nbr_oligo)
            if nbr_prey > 30 : #take 30 better interactions from RosseTTA
                break
            if os.path.exists(f"./result_APD_homo_int/{prey}_homo_{nbr_oligo}er/ranked_0.pdb") == False : #if model don't exist
                if int_lenght <= max_aa: #make interaction if doesn't exist and is not too long
                    save_lenght_line[f"{prey}_homer_{nbr_oligo}er"] = [int_lenght, f"{prey}:{nbr_oligo}\n"]
                    nbr_prey += 1
                else : #if interaction is too large
                    OOM_int = prey + ":" + nbr_oligo + "\n"
                    result_dict[prey] = result_dict[prey] + " | Homo-oligomer too large for your GPU"
            else :
                nbr_prey += 1
                pass
    dict_split_GPU = Split_to_GPU(file, save_lenght_line, GPU)
    for GPU_i in dict_split_GPU.keys() :
        with open(f"{Interaction_file}_{GPU_i}.txt",'w') as final_file :
            final_file.write(dict_split_GPU[GPU_i])
    with open("OOM_interactions.txt", "w") as OOM_file :
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
    print("Time APD:", end_time - start_time)

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
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    env['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
    env['TF_FORCE_UNIFIED_MEMORY'] = 'true'
    env['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '3.2'
    env['XLA_FLAGS'] = '--xla_gpu_enable_triton_gemm=false'
    cmd =f"run_multimer_jobs.py --mode=custom \--num_cycle=3 \--num_predictions_per_model=1 \--compress_result_pickles=True \--output_path=./result_{Interaction_file} \--data_dir={Path_AlphaFold_Data} \--protein_lists={Interaction_file}_GPU_{gpu_id}.txt \--monomer_objects_dir={Path_Pickle_Feature} \--remove_keys_from_pickles=False"
    log_file = f"log_GPU_{gpu_id}.txt"
    subprocess.run(cmd, shell=True, env=env)


def Prepare_RF2_PPI(file, Path_Pickle_Feature, possible_baits, Interaction, GPU, regions) :
    """
    Prepare script for RoseTTAFold2-PPI.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_Pickle_Feature : string
    possible_baits : list
    Interaction : string
    GPU : list
    regions : dict

    Returns:
    ----------
    total_int : integer
    """
    total_int = 0
    save_lenght_line = dict()
    possible_prey = file.get_possible_prey()
    new_possible_prey = list()
    lenght_prot = file.get_lenght_prot()
    result_dict = file.get_result_dict()
    index_file = 0
    pynvml.nvmlInit()
    device_count = pynvml.nvmlDeviceGetCount()
    for i in range(len(GPU)):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        name = pynvml.nvmlDeviceGetName(handle).decode("utf-8")
        mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        vram = (mem_info.total / 1024**2) * 0.001 # in MiB
    pynvml.nvmlShutdown()
    max_lenght_int = int((vram + 14.2) / 0.0223) #bench on data
    if os.path.exists(f"{Path_Pickle_Feature}/{Interaction}.txt") :
       os.remove(f"{Path_Pickle_Feature}/{Interaction}.txt")
    if Interaction == "RF2_homo_int" :
        for prey in possible_prey :
            lenght_first_prot = lenght_prot[prey]
            save_lenght_line[f"{prey}_and_{prey}.a3m"] = [lenght_prot[prey] + lenght_prot[prey],f"{Path_Pickle_Feature}/{prey}_and_{prey}.a3m {lenght_first_prot}\n"] #save lenght of interactions and line to write in the final
            total_int += 1
            if os.path.exists(f"{Path_Pickle_Feature}/{prey}_and_{prey}.a3m") == False :
                fusioned_MSA(prey+".a3m",prey+".a3m", Path_Pickle_Feature, "PPI", total_int)
    if Interaction == "RF2_PPI_int" :
        for bait in possible_baits :
            if regions[bait] != "0-0" :
                start = regions[bait].split("-")[0]
                end = regions[bait].split("-")[1]
                lenght_bait = int(end) - int(start)
                name_bait = f"{bait}_{start}_{end}"
            else :
                lenght_bait = lenght_prot[bait]
                name_bait = bait
            for prey in possible_prey :
                lenght_int = lenght_bait + lenght_prot[prey]
                if lenght_int > max_lenght_int :#remove interactions with a too big length
                    result_dict[prey] = result_dict[prey] + "Too big interactions : RF2-PPI OOM"
                else :
                    total_int += 1
                    if os.path.exists(f"{Path_Pickle_Feature}/{name_bait}_and_{prey}.a3m") == False :
                        total_int = fusioned_MSA(name_bait+".a3m",prey+".a3m", Path_Pickle_Feature, "PPI", total_int)
                    save_lenght_line[f"{name_bait}_and_{prey}.a3m"] = [lenght_bait + lenght_prot[prey], f"{Path_Pickle_Feature}/{name_bait}_and_{prey}.a3m {lenght_bait}\n"]
                    new_possible_prey.append(prey)
    sorted_list = list(sorted(save_lenght_line.items(), key=lambda item: item[1][0], reverse = True)) #sorted in function of interaction lenght
    all_lines = str()
    for interactions in sorted_list :
        all_lines += interactions[1][1]
    with open(f"{Path_Pickle_Feature}/{Interaction}.txt",'w') as file_int :
        file_int.write(all_lines)
        print(all_lines)
    print("Total interactions : " + str(total_int))
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)
    return total_int

def Prepare_RF2_Lite(file, Path_Pickle_Feature, possible_baits, GPU, regions) :
    """
    Prepare script for RF2-Lite.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_Pickle_Feature : string
    Interaction : string
    GPU : list
    regions : dict

    Returns:
    ----------
    """
    total_int = 0
    save_lenght_line = dict()
    possible_prey = file.get_possible_prey()
    new_possible_prey = list()
    lenght_prot = file.get_lenght_prot()
    result_dict = file.get_result_dict()
    index_file = 0
    pynvml.nvmlInit()
    for i in range(len(GPU)):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        name = pynvml.nvmlDeviceGetName(handle).decode("utf-8")
        mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        vram = (mem_info.total / 1024**2) * 0.001 # en MiB
    pynvml.nvmlShutdown()
    max_lenght_int = int(vram * 47) #need to bench RF2-Lite on data
    if os.path.exists(f"{Path_Pickle_Feature}/RF2_PPI_int.txt") :
        os.remove(f"{Path_Pickle_Feature}/RF2_PPI_int.txt")
    for bait in possible_baits :
        if regions[bait] != "0-0" :
            lenght_bait = int(regions[bait].split("-")[1]) - int(regions[bait].split("-")[0])
            name_bait = f"{bait}_{str(regions[bait].split('-')[0])}_{str(regions[bait].split('-')[1])}"
        else :
            lenght_bait = lenght_prot[bait]
            name_bait = bait
        for prey in possible_prey :
            total_lenght = lenght_bait + lenght_prot[prey]
            if total_lenght > max_lenght_int :#remove interactions with a too big length
                result_dict[prey] = result_dict[prey] + "Too big interactions : RF2-Lite OOM"
            else :
                total_int += 1
                if os.path.exists(f"{Path_Pickle_Feature}/{name_bait}_and_{prey}.a3m") == False : #if paired MSA doesn't exist
                    new_total_int = fusioned_MSA(name_bait+".a3m",prey+".a3m", Path_Pickle_Feature, "Lite", total_int)
                    if new_total_int != total_int-1 :#if fusionned MSA have sequence to merge
                        save_lenght_line[f"{name_bait}_and_{prey}.a3m"] = [total_lenght, f"{Path_Pickle_Feature}/{name_bait}_and_{prey}.a3m {lenght_bait} {lenght_prot[prey]} {Path_Pickle_Feature}/result_RF2_Lite_PPI/{name_bait}_and_{prey}\n"]
                        new_possible_prey.append(prey) #add all preys to the new possible prey list
                    else :
                        total_int = new_total_int #set new number of prey
                        result_dict[prey] = result_dict[prey] + "No MSA to merge with bait" 
                else : #if paired MSA exist
                    save_lenght_line[f"{name_bait}_and_{prey}.a3m"] = [total_lenght, f"{Path_Pickle_Feature}/{name_bait}_and_{prey}.a3m {lenght_bait} {lenght_prot[prey]} {Path_Pickle_Feature}/result_RF2_Lite_PPI/{name_bait}_and_{prey}\n"]
                    new_possible_prey.append(prey) #add all preys to the new possible prey list
    print(str(datetime.now()) + " End of MSA fusion")
    dict_split_GPU = Split_to_GPU(file, save_lenght_line, GPU)
    for GPU_i in dict_split_GPU.keys() :
        with open(f"{Path_Pickle_Feature}/RF2_PPI_int_{GPU_i}.txt",'w') as file_int :
            file_int.write(dict_split_GPU[GPU_i])
    print("Total interactions : " + str(total_int))
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)



def fusioned_MSA(a3m1, a3m2, Path_Pickle_Feature, RF2, total_int):
    """
    Fusioned MSA files.
    Adapted from https://github.com/uw-ipd/RoseTTAFold2/blob/main/input_prep/make_paired_MSA_simple.py

    Parameters:
    ----------
    a3m1 : string
    a3m2 : string
    Path_Pickle_Feature : string
    total_int : integer

    Returns:
    total_int : integer
    ----------
    """
    msa1, lab1 = read_a3m(f"{Path_Pickle_Feature}/{a3m1}")
    msa2, lab2 = read_a3m(f"{Path_Pickle_Feature}/{a3m2}")
    if len(lab1[1:]) == 0 or len(lab2[1:]) == 0: #check if we have sequences in MSA
        print(f"No sequences to merge: {a3m2}, {a3m1}")
        return total_int-1
    valid_lab1 = [id_ for id_ in lab1[1:] if len(id_) in (6, 10)] #filtred valid IDs
    valid_lab2 = [id_ for id_ in lab2[1:] if len(id_) in (6, 10)]
    msa1_valid = [msa1[i+1] for i, id_ in enumerate(lab1[1:]) if len(id_) in (6, 10)]
    msa2_valid = [msa2[i+1] for i, id_ in enumerate(lab2[1:]) if len(id_) in (6, 10)]
    hash1 = Uni_to_idx(valid_lab1) #hash in fucntion of UniprotID
    hash2 = Uni_to_idx(valid_lab2)
    idx1, idx2 = np.where(np.abs(hash1[:, None] - hash2[None, :]) < 10) #found close pairs (hash deviation of 10)
    if RF2=="Lite" and len(idx1) >= 4095:
        idx1 = idx1[:4095]
        idx2 = idx2[:4095]
    name = str(a3m1).split(".")[0] + "_and_" + str(a3m2).split(".")[0]
    with open(f'{Path_Pickle_Feature}/{name}.a3m', 'wt') as f: #write fusion file
        f.write('>query\n%s%s\n' % (msa1[0], msa2[0])) #add first query sequence
        for i, j in zip(idx1, idx2): #add found pairs
            if not all(cara == '-' for cara in msa1_valid[i]) :
                f.write(">%s_%s\n%s%s\n" % (valid_lab1[i], valid_lab2[j],msa1_valid[i], msa2_valid[j]))
            else :
                pass
    return total_int

def read_a3m(a3m):
    """
    Parse an a3m files as a dictionary {label->sequence}.
    Adapted from https://github.com/uw-ipd/RoseTTAFold2/blob/main/input_prep/make_paired_MSA_simple.py
    Parameters:
    ----------
    a3m : string

    Returns:
    ----------
    seq : list
    lab : list
    """
    seq = []
    lab = []
    is_first = True
    is_incl = False
    for line in open(a3m, "r"):
        if line[0] == '>':
            label = line.strip()[1:]
            is_incl = True
            if is_first: #include first sequence (query)
                is_first = False
                lab.append(label)
                continue
            if "UniRef" in label:
                code = label.split()[0].split('_')[-1]
                if code.startswith("UPI"): #UniParc identifier -- exclude
                    is_incl = False
                    continue
            elif label.startswith("tr|"):
                code = label.split('|')[1]
            else:
                is_incl = False
                continue
            lab.append(code)
        else:
            if is_incl:
                seq.append(line.rstrip())
            else:
                continue
    return seq, lab


def Uni_to_idx(ids):
    """
    Convert Uniprot ID in integer.
    Adapted from https://github.com/uw-ipd/RoseTTAFold2/blob/main/input_prep/make_paired_MSA_simple.py

    Parameters:
    ----------
    ids : string

    Returns:
    ----------
    table_idx : numpy table
    """
    filtered_ids = []
    for i in ids:
        if len(i) == 6:
            filtered_ids.append(i + 'AAA0')
        elif len(i) == 10:
            filtered_ids.append(i)
        else:
            print(f"[WARN] ID ignoré : {i}")
    if not filtered_ids:
        raise ValueError("No valid ID")
    ids2 = [s.ljust(10, 'A') for s in filtered_ids]
    assert all(len(i) == 10 for i in ids2)
    arr = np.array([list(s) for s in ids2], dtype='|S1').view(np.uint8)
    for i in [1, 5, 9]:
        arr[:, i] -= ord('0')
    arr[arr >= ord('A')] -= ord('A')
    arr[arr >= ord('0')] -= ord('0') - 26
    arr[:, 0][arr[:, 0] > ord('Q') - ord('A')] -= 3

    arr = arr.astype(np.int64)

    coef = np.array([23, 10, 26, 36, 36, 10, 26, 36, 36, 1], dtype=np.int64)
    coef = np.tile(coef[None, :], [len(ids2), 1])

    c1 = [i for i, id_ in enumerate(filtered_ids) if id_[0] in 'OPQ' and len(id_) == 10]
    c2 = [i for i, id_ in enumerate(filtered_ids) if id_[0] not in 'OPQ' and len(id_) == 10]

    coef[c1] = np.array([3, 10, 36, 36, 36, 1, 1, 1, 1, 1])
    coef[c2] = np.array([23, 10, 26, 36, 36, 1, 1, 1, 1, 1])

    for i in range(1, 10):
        coef[:, -i - 1] *= coef[:, -i]
    table_idx = np.sum(arr * coef, axis=-1)
    return table_idx

def Launch_RF2_PPI(Path_Pickle_Feature, Interaction, GPU_list) :
    """
    Launch RoseTTAFold2-PPI.

    Parameters:
    ----------
    Path_Pickle_Feature : string
    Interaction : string
    GPU_list : list

    Returns:
    ----------
    """
    start_time = datetime.now()
    processes = []
    GPU_index = 0
    with open(f"{Path_Pickle_Feature}/{Interaction}.txt", "r") as all_cmd :
        interactions = [line.strip() for line in all_cmd if line.strip()]
    gpu_queues = {GPU: Queue() for GPU in GPU_list}
    for idx, interaction in enumerate(interactions):
        GPU = GPU_list[idx % len(GPU_list)]
        gpu_queues[GPU].put(interaction)
    start_time = datetime.now()
    processes = []
    for GPU in GPU_list:
        p = multiprocessing.Process(target=run_parallel_jobs_RF2_PPI, args=(Path_Pickle_Feature, Interaction, GPU, gpu_queues[GPU]))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
    end_time = datetime.now()
    duration = end_time - start_time
    print("Duration RoseTTAFold2-PPI :", duration)

def run_parallel_jobs_RF2_PPI(Path_Pickle_Feature, Interaction, gpu_index, queue) :
    """
    Run parallel jobs one by one for RoseTTAFold2-PPI to avoid OOM errors. Rewrite all scores in one file. 

    Parameters:
    ----------
    Path_Pickle_Feature : string
    Interaction : string
    gpu_index : integer
    queue : multiprocessing.Queue

    Returns:
    ----------
    """
    all_result = str()
    while not queue.empty():
        interaction = queue.get()
        tmp_file = f"{Path_Pickle_Feature}/{Interaction}_GPU_{gpu_index}.txt"
        with open(tmp_file, "w") as f:
            f.write(interaction)
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_index)
        cmd = ["python",
        "/data/Rosetta-PPI/RoseTTAFold2-PPI/src/predict_list_PPI.py",
        "-list_fn", f"{Path_Pickle_Feature}/{Interaction}_GPU_{gpu_index}.txt",
        "-model_file", "/data/Rosetta-PPI/RoseTTAFold2-PPI/src/models/RF2-PPI.pt",
        "-number_seqs", "5000"]
        print(f"[GPU {gpu_index}] Launch interaction : {interaction}")
        subprocess.run(cmd, env=env)
        print(f"[GPU {gpu_index}] Done : {interaction}")
        with open(f"{Path_Pickle_Feature}/{Interaction}_GPU_{gpu_index}.txt.log", "r") as result_file:
            for line in result_file:
                all_result += str(line)
    with open(f"{Path_Pickle_Feature}/{Interaction}_GPU_{gpu_index}.txt.log", "w") as new_result_file:
        new_result_file.write(all_result)

def Launch_RF2_Lite(Path_Pickle_Feature, GPU) :
    """
    Launch RoseTTAFold2_Lite.

    Parameters:
    ----------
    Path_Pickle_Feature : string
    GPU : list

    Returns:
    ----------
    """
    start_time = datetime.now()
    processes = []
    for GPU_index in GPU:
        env = os.environ.copy()
        script_path = os.path.expanduser("~/RF2-Lite/networks/predict_complex_list.py")
        cmd = ["python",
                script_path,
                "-list", f"{Path_Pickle_Feature}/RF2_PPI_int_GPU_{GPU_index}.txt",
                "-p", f"cuda:{GPU_index}"]
        p = subprocess.Popen(cmd,env=env)
        processes.append(p)
    for p in processes:
        p.wait()
    end_time = datetime.now()
    duration = end_time - start_time
    print("Duration RoseTTAFold2-Lite :", duration)

def Class_output_RF2_PPI(file, Path_Pickle_Feature, Interaction, GPU) :
    """
    Classified PPIs and set new possible prey classed in function of scores.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_Pickle_Feature : string
    Interaction : string
    GPU : list

    Returns:
    ----------
    """
    result_dict = file.get_result_dict()
    score_dict = dict()
    possible_prey = list()
    for GPU_index in GPU :
        with open(f"{Path_Pickle_Feature}/{Interaction}_GPU_{GPU_index}.txt.log", "r") as log_file :
            for line in log_file :
                if line.strip("\n") != "done" :
                    name = line.split("\t")[0].split("_and_")[1].split(".")[0]
                    score = line.split("\t")[1]
                    if name not in score_dict.keys() :
                        score_dict[name] = score
                    else :
                        if score_dict[name] < score :
                            score_dict[name] = score
    if Interaction == "RF2_homo_int" :
        for key in score_dict.keys() :
            result_dict[key] = result_dict[key] + f"RF2_homo_int : {score_dict[key]} | "
            if float(score_dict[key]) >= 0.24 :
                possible_prey.append(key)
    if Interaction == "RF2_PPI_int" :
        sorted_list = list(sorted(score_dict.items(), key=lambda item: item[1], reverse = True))
        for prot_score in sorted_list :
            possible_prey.append(prot_score[0]) #set a sorted list with better interactions in first
    file.set_possible_prey(possible_prey)
    file.set_result_dict(result_dict)


def Use_RF2_PPI(file, Informations_dict, Interaction, GPU, regions) :
    """
    Main of RoseTTAFold2-PPI.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary
    Interaction : string
    GPU : list

    Returns:
    ----------
    """
    print(str(datetime.now()) + " Prepare RoseTTAFold2-PPI")
    total_int = Prepare_RF2_PPI(file, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"], Interaction, GPU, regions)
    nbr_line_result = 0
    for GPU_index in GPU :
        if os.path.exists(f"{Informations_dict['Path_Pickle_Feature']}/{Interaction}_GPU_{GPU_index}.txt.log") == True :
            with open(f"{Informations_dict['Path_Pickle_Feature']}/{Interaction}_GPU_{GPU_index}.txt.log", "r") as result_file :
                for line in result_file:
                    nbr_line_result += 1
    nbr_line_result = nbr_line_result // 2
    if total_int != nbr_line_result : #check if all interactions have been done
        print(str(datetime.now()) + " Launch RoseTTAFold2-PPI")
        Launch_RF2_PPI(Informations_dict["Path_Pickle_Feature"], Interaction, GPU)
        print(str(datetime.now()) + " Scoring RoseTTAFold2-PPI")
        Class_output_RF2_PPI(file, Informations_dict["Path_Pickle_Feature"], Interaction, GPU)
    else :
        Class_output_RF2_PPI(file, Informations_dict["Path_Pickle_Feature"], Interaction, GPU)
        print(f"RoseTTAFold2-PPI script already exist, skip this step. If error occur, delete {Interaction}.txt and relaunch.")

def Class_output_RF2_Lite(file, Path_Pickle_Feature, possible_baits, GPU, regions) :
    """
    Classified PPIs and set new possible prey classed in function of scores.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_Pickle_Feature : string
    possible_baits : list
    GPU : list
    regions : dict

    Returns:
    ----------
    """
    result_dict = file.get_result_dict()
    score_dict = dict()
    possible_prey = file.get_possible_prey()
    new_possible_prey = list()
    all_lines = ""
    for bait in possible_baits :
        if regions[bait] != "0-0" :
            name_bait = f"{bait}_{str(regions[bait].split('-')[0])}_{str(regions[bait].split('-')[1])}"
        else :
            name_bait = bait
        for prey in possible_prey :
            data = load(f"{Path_Pickle_Feature}/result_RF2_Lite_PPI/{name_bait}_and_{prey}_00.npz")
            int_score = np.max(data['dist'][:-10,10:]) # from https://github.com/SNU-CSSB/RF2-Lite
            all_lines += f"{name_bait}_and_{prey} {int_score}\n"
            score_dict[prey] = int_score
            result_dict[prey] += f" RF2_Lite_int : {int_score}" 
    sorted_list = list(sorted(score_dict.items(), key=lambda item: item[1], reverse = True))
    for prot_score in sorted_list :
        new_possible_prey.append(prot_score[0])
    file.set_possible_prey(new_possible_prey)
    file.set_result_dict(result_dict)

def Use_RF2_Lite(file, Informations_dict, GPU, regions) :
    """
    Main of RoseTTA_Lite.l

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary
    GPU : list

    Returns:
    ----------
    """
    print(str(datetime.now()) + " Prepare RoseTTAFold2-Lite")
    Prepare_RF2_Lite(file, Informations_dict["Path_Pickle_Feature"], Informations_dict["Interact_with"], GPU, regions)
    print(str(datetime.now()) + " Launch RoseTTAFold2-Lite")
    Launch_RF2_Lite(Informations_dict["Path_Pickle_Feature"], GPU)
    print(str(datetime.now()) + " Scoring RoseTTAFold2-Lite")
    Class_output_RF2_Lite(file, Informations_dict["Path_Pickle_Feature"],Informations_dict["Interact_with"], GPU, regions)


def Score_interaction_APD (file, Informations_dict, Interaction) :
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
    score_dict = dict()
    possible_prey = list()
    Path_Singularity_Image = Informations_dict["Path_Singularity_Image"]
    if os.path.isdir(f"./result_{Interaction}") == True :
        cmd = f"singularity exec --no-home --bind result_{Interaction}:/mnt {Path_Singularity_Image} run_get_good_pae.sh --output_dir=/mnt --cutoff=10"
        os.system(cmd)
        with open(f"result_{Interaction}/predictions_with_good_interpae.csv", "r") as file1 :
            reader = csv.DictReader(file1)

            if Interaction == "APD_PPI_int" :
                all_lines = "jobs,pi_score,iptm_ptm,pDockQ,iQ_score\n"
                for row in reader :
                    job = row['jobs']
                    if '_and_' in job :
                        if row['pi_score'] == 'No interface detected' :
                            iQ_score = float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30 #pi_score don't detect interface so is set on -2.63
                            line =f'{row["jobs"]},-2.63,{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                        else :
                            iQ_score = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                            line =f'{row["jobs"]},{row["pi_score"]},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                    result_dict[job.split("_and_")[1]] += f" | iQ_score : {iQ_score}"
                    all_lines = all_lines + line
                    score_dict[job.split("_and_")[1]] = iQ_score
                sorted_list = list(sorted(score_dict.items(), key=lambda item: item[1], reverse = True))
                for prot_score in sorted_list :
                    possible_prey.append(prot_score[0]) #set a sorted list with better interactions in first

            if Interaction == "APD_homo_int" :
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
                        if hiQ_score >= 50 :
                            possible_prey.append(prot_name)
                    else :
                        hiQ_score = (((float(all_homo[key][0])/all_homo[key][1])+2.63)/5.26)*60+float(row['iptm_ptm'])*40
                        line =f'{key},{str(float(all_homo[key][0])/all_homo[key][1])},{row["iptm_ptm"]},{str(hiQ_score)}\n'
                        all_lines += line
                        if hiQ_score >= 50 :
                            possible_prey.append(prot_name)
                    result_dict[key.split("_homo_")[0]] += f" | hiQ_score : {hiQ_score}"
            with open(f"result_{Interaction}/new_predictions_with_good_interpae.csv", "w") as file2 :
                file2.write(all_lines)
        file.set_result_dict(result_dict)
        file.set_possible_prey(possible_prey)
    else :
        print(f"result_{Interaction}/ don't exist")
    end_time = datetime.now()
    print("Time Scoring APD:", end_time - start_time)

def Split_to_GPU(file, save_lenght_line, GPU) :
    """
    Split interactions in function of their lenght in different GPU.

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
    index_GPU = 0
    dict_split_GPU = dict()
    for nbr_GPU in GPU :
        dict_split_GPU[f"GPU_{nbr_GPU}"] = ""
    for interactions in sorted_list :
        dict_split_GPU[f"GPU_{GPU[index_GPU]}"] += interactions[1][1]
        index_GPU += 1
        if index_GPU > len(GPU)-1 :
            index_GPU = 0
    return (dict_split_GPU)

def Resume_file(file) :
    """
    Create a resume file with all interactions, their scores and their status.

    Parameters:
    ----------
    file : object of class File_proteins

    Returns:
    ----------
    """
    all_lines = str()
    result_dict = file.get_result_dict()
    possible_prey = file.get_possible_prey()
    proteins = file.get_proteins()
    print("Possible prey : " + str(possible_prey))
    for prey in possible_prey :
        all_lines += prey + " : " + result_dict[prey] + "\n"
        del (result_dict[prey])
    proteins_rest = dict(sorted(result_dict.items(), key=lambda item: len(item[1]), reverse=True))
    for rest_prey in proteins_rest.keys() :
        all_lines += rest_prey + " : " + result_dict[rest_prey] + "\n"
    with open("Final_result.txt", "w") as result_file :
        result_file.write(all_lines)
