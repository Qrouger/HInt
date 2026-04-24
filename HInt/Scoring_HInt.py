""" Scoring file of HInt

    Author: Quentin Rouger
"""
import os
import csv
import logging
import multiprocessing
from multiprocessing import Pool
from tqdm import tqdm
from datetime import datetime
import HInt.get_good_inter_pae
import gc
import pandas as pd
import json
import pickle
import gzip
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import softmax
from Bio import PDB
import copy
import string
import subprocess
import math
from pathlib import Path

# Configure global logger
logging.basicConfig(
    filename="./log_file/HInt.log",  # Log file name
    level=logging.INFO,  # Log level
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

logger = logging.getLogger()

def Score_interaction (file, Informations_dict, CPU, Interaction, bait=None, multi_scoring=False) :
    """
    Compute interaction scores (PPI or homo-oligomer) from AlphaFold predictions, aggregate them into meaningful metrics (iQ_score / hiQ_score), and update the list of valid prey proteins accordingly.

    - detects which interactions still need to be scored
    - runs scoring in parallel (CPU multiprocessing)
    - merges all results into CSV files
    - filters out bad interactions based on inter-PAE / interface quality

    Parameters :
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary
    CPU : int
    Interaction : str
        Interaction type ("PPI_int" or "homo_int")
    bait : str
    multi_scoring : str
    """
    start_time = datetime.now()
    result_dict = file.get_result_dict()
    possible_prey = file.get_possible_prey()
    int_score = file.get_int_score() #saved scores
    homo_score = file.get_homo_score()
    seq_no_SP = file.get_proteins_sequence_no_SP()
    prot_lenght = file.get_lenght_prot()
    new_possible_prey = list()
    Path_ccp4 = Informations_dict["Path_ccp4"]
    regions = Informations_dict["Regions"]
    nbr_homo = Informations_dict["Homo-oligomer"]
    AF_version = Informations_dict["AlphaFold"]
    ppi_list = list()
    already_done = list()
    already_dict = dict()
    mean_iQ_score = dict()
    cv_iQ_score = dict()

    multi_scoring = False if multi_scoring == "False" else True
    result = subprocess.run(["conda", "env", "list", "--json"],capture_output=True, text=True, check=True)
    envs = json.loads(result.stdout)["envs"]
    exists = any("pi_score" in env for env in envs)
    if not exists :
        subprocess.run(["conda", "create", "-y", "-n", "pi_score","python=2.7", "scikit-learn=0.20.4", "biopython", "biopandas"], check=True)

    if bait is not None : #setup bait name
        bait_name = bait.replace(",","_and_")
        for prot in bait.split(",") :
            if regions[prot] != "0-0" :
                start = int(regions[prot].split("-")[0])
                end = int(regions[prot].split("-")[1])
                bait_name = bait_name.replace(prot,f"{prot}_{start}-{end}")

    #Multiprocessing to score all interactions
    if os.path.isdir(f"./result_{Interaction}") == True :
        if bait is None : #for homo-oligomer
            for protein in possible_prey :
                if f"hiQ_score_{nbr_homo}er" not in homo_score[protein].keys() :
                    ppi_list.append(f"./result_{Interaction}/{protein}_homo_{nbr_homo}er") #Found a solution for score only homo-oligomer without score
                else :
                    result_dict[protein][f"hiQ_score_{nbr_homo}er"] = homo_score[protein][f"hiQ_score_{nbr_homo}er"]
                    if homo_score[protein][f"hiQ_score_{nbr_homo}er"] > 0 :
                        new_possible_prey.append(protein)
                    if homo_score[protein][f"hiQ_score_{nbr_homo}er"] == 0 :
                        result_dict[protein]["Reason_for_filtering"] = f"Bad homo-oligomer PAE : AF"
        else : #for one vs all
            for protein in possible_prey : #check if protein is already score
                if f"iQ_score_vs_{bait_name}" not in int_score[protein].keys() :
                    dir_path = Path(f"./result_{Interaction}/{bait_name}_and_{protein}")
                    inv_dir = Path(f"./result_{Interaction}/{protein}_and_{bait_name}")
                    if dir_path.exists() and dir_path.is_dir() :
                        ppi_list.append(f"./result_{Interaction}/{bait_name}_and_{protein}")
                    elif inv_dir.exists() and inv_dir.is_dir() :
                        ppi_list.append(f"./result_{Interaction}/{protein}_and_{bait_name}")
                else :
                    result_dict[protein][f"iQ_score_vs_{bait_name}"] = int_score[protein][f"iQ_score_vs_{bait_name}"]
                    if int_score[protein][f"iQ_score_vs_{bait_name}"] > 0 :
                        new_possible_prey.append(protein)
                    if int_score[protein][f"iQ_score_vs_{bait_name}"] == 0 :
                        result_dict[protein]["Reason_for_filtering"] = f"Bad interactions with {bait_name} : inter PAE > 10 A"

        results = []
        with multiprocessing.Pool(CPU) as pool : #just run scoring for interactions without score
            tasks = [(ppi, file, AF_version, Path_ccp4, multi_scoring) for ppi in ppi_list]
            results_iter = pool.imap_unordered(run_scoring, tasks)
            for df in tqdm(results_iter, total=len(ppi_list), desc="Scoring interactions") :
                if df is not None and not df.empty :
                    results.append(df)
            pool.close()
            pool.join()

        merged_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame() #write an empty dataframe if no result
        merged_df.to_csv(os.path.join(f"./result_{Interaction}", "predictions_with_good_interpae.csv"), index=False)
        if ppi_list : #Add result in dict result if there is new ppi scored
            #Resume all score and set new possible prey
            with open(f"result_{Interaction}/predictions_with_good_interpae.csv", "r") as result_file :
                reader = csv.DictReader(result_file)

                #For one vs all
                if Interaction == "PPI_int" or "Compounds" : #make int_score
                    all_lines = "jobs,pi_score,iptm_ptm,pDockQ,iQ_score\n"
                    for row in reader :
                        job = row['jobs']
                        just_name = job.split("_ranked_")[0]

                        if '_and_' in job and just_name.split("_and_")[-1] in possible_prey and bait_name in job : #check if interaction is a PPI and if prey is in possible prey list
                            if "," in bait : #multimer bait
                                if row['pi_score'] == 'No interface detected' :
                                    if job in already_done : #if protein have multi interface interaction, mean of pi_score #multimeric bait
                                        all_lines = '\n'.join(all_lines.rstrip('\n').split('\n')[:-1]) #delete last line
                                        already_dict[job].append(float(-2.63))
                                        mean_pi_score = (sum(already_dict[job]) + float(-2.63)) / len(already_dict[job])
                                        iQ_score = ((mean_pi_score+2.63)/5.26)*60+float(row['iptm_ptm'])*40
                                        line =f'\n{just_name},{str(mean_pi_score)},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                                    else :
                                        already_dict[job] = [float(-2.63)]
                                        iQ_score = float(row['iptm_ptm'])*30#pi_score don't detect interface so it's set on -2.63
                                        line =f'{just_name},-2.63,{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                                        already_done.append(job)
                                else :
                                    if job in already_done : #if protein have multi interface interaction, mean of pi_score  #multimeric bait
                                        all_lines = '\n'.join(all_lines.rstrip('\n').split('\n')[:-1]) #delete last line
                                        already_dict[job].append(float(row['pi_score']))
                                        mean_pi_score = (sum(already_dict[job]) + float(row['pi_score'])) / len(already_dict[job])
                                        iQ_score = ((mean_pi_score+2.63)/5.26)*60+float(row['iptm_ptm'])*40
                                        line =f'\n{just_name},{str(mean_pi_score)},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                                    else :
                                        already_dict[job] = [float(row['pi_score'])]
                                        iQ_score = ((float(row['pi_score'])+2.63)/5.26)*60+float(row['iptm_ptm'])*40
                                        line =f'{just_name},{row["pi_score"]},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                                        already_done.append(job)
                            
                            else : #need to use hiQ_score for multimer bait and look at interface
                                if row['pi_score'] == 'No interface detected' :
                                    iQ_score = float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30 #pi_score don't detect interface so it's set on -2.63
                                    if "ranked_0" in job :
                                        line =f'{just_name},-2.63,{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                                        mean_iQ_score[just_name.split("_and_")[-1]] = []
                                else :
                                    iQ_score = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                                    if "ranked_0" in job :
                                        line =f'{just_name},-2.63,{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                                        mean_iQ_score[just_name.split("_and_")[-1]] =[]
                            mean_iQ_score[just_name.split("_and_")[-1]].append(iQ_score)
                            int_score[just_name.split("_and_")[-1]][f"iQ_score_vs_{bait_name}"] = iQ_score
                            result_dict[just_name.split("_and_")[-1]][f"iQ_score_vs_{bait_name}"] = iQ_score
                            if just_name.split("_and_")[-1] not in new_possible_prey :
                                new_possible_prey.append(just_name.split("_and_")[-1])
                            if "ranked_0" in job :
                                all_lines = all_lines + line
                            name_int = just_name.split("/")[-1]
                    for protein in possible_prey :
                        if multi_scoring == True and protein in mean_iQ_score.keys() :
                            scores = mean_iQ_score[protein]
                            print(protein, scores)
                            mean = sum(scores) / len(scores)
                            variance = sum((x - mean) ** 2 for x in scores) / len(scores)
                            std = math.sqrt(variance)
                            cv_iQ_score[protein] = std / mean
                            mean_iQ_score[protein] = sum(scores) / len(scores)
                        if protein not in new_possible_prey :
                            int_score[protein][f"iQ_score_vs_{bait_name}"] = 0
                            result_dict[protein][f"iQ_score_vs_{bait_name}"] = 0 #if prey don't have interaction, set iQ_score to 0
                            result_dict[protein]["Reason_for_filtering"] = f"Bad interactions with {bait_name} : inter PAE > 10 A"

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
                        #number_oligo = row["jobs"].split("_")[2].replace("er","") #AFPD 2.0.4 #problem with "_" in name sequence ?
                        prot_name = row["jobs"].split("_")[0]
                        if len(save_pi_score[key]) > int(nbr_homo) : #if model have more interface than number of homo-oligomerization
                            new_sum_pi_score = 0
                            save_pi_score[key].sort(reverse=True)
                            for index in range(0,int(nbr_homo)) :
                                new_sum_pi_score += save_pi_score[key][index]
                                hiQ_score = (((float(new_sum_pi_score)/int(nbr_homo))+2.63)/5.26)*60+float(row['iptm_ptm'])*40 #cause iptm_ptm are always same for each interface
                            line =f'{key},{str(float(new_sum_pi_score)/int(nbr_homo))},{row["iptm_ptm"]},{str(hiQ_score)}\n'
                            all_lines += line   
                        else :
                            hiQ_score = (((float(all_homo[key][0])/all_homo[key][1])+2.63)/5.26)*60+float(row['iptm_ptm'])*40
                            line =f'{key},{str(float(all_homo[key][0])/all_homo[key][1])},{row["iptm_ptm"]},{str(hiQ_score)}\n'
                            all_lines += line
                        if prot_name not in new_possible_prey :
                            new_possible_prey.append(prot_name)
                        result_dict[prot_name][f"hiQ_score_{nbr_homo}er"] = hiQ_score
                        homo_score[prot_name][f"hiQ_score_{nbr_homo}er"] = hiQ_score
                        name_int = key.split("/")[-1]
                    for protein in possible_prey :
                        if protein not in new_possible_prey :
                            homo_score[protein][f"hiQ_score_{nbr_homo}er"] = 0
                            result_dict[protein][f"hiQ_score_{nbr_homo}er"] = 0
                            result_dict[protein]["Reason_for_filtering"] = "Bad homo-oligomer PAE : AF"


                if len(all_lines.strip("\n")) > 1 : #if all_lines is not empty
                    with open(f"result_{Interaction}/new_predictions_with_good_interpae.csv", "w") as file2 :
                        file2.write(all_lines)
            end_time = datetime.now()
            logger.info("Time scoring interactions : %s\n", end_time - start_time)
            file.set_homo_score(homo_score)
            file.set_int_score(int_score)
        print(cv_iQ_score)
        print(mean_iQ_score)
        file.set_result_dict(result_dict)
        file.set_possible_prey(new_possible_prey)
    else :
        logger.info(f"result_{Interaction}/ don't exist")

def run_scoring (args) :
    """
    Wrapper function executed by multiprocessing workers to score a single AlphaFold interaction using inter-chain PAE metrics.

    - it unpacks arguments passed by the multiprocessing Pool
    - calls the get_good_inter_pae scoring routine
    - returns the resulting DataFrame

    Parameters :
    ----------
    args : tuple

    Returns :
    ----------
    result : pandas.DataFrame
    """
    interaction, file, AF_version, Path_ccp4, multi_scoring = args
    try :
        result = HInt.get_good_inter_pae.main(interaction, 100, 2, file, AF_version, Path_ccp4, multi_scoring) #normal PAE is 10
        return  result
    except Exception as e:
        pid = os.getpid()
        logger.error(f"ERROR in worker PID={pid}")
        logger.error(f"Interaction: {interaction}")
        raise

def Resume_file(file, Informations_dict) :
    """
    Create a resume file with all interactions, their scores and their status.

    - PPI interaction scores (iQ_score vs one or multiple baits)
    - Homo-oligomerization scores (hiQ_score)
    - Filtering reasons (OOM, size limits, etc.)
    - Protein-level annotations (DeepLoc, signal peptide)

    Two CSV files are produced:
    ------------------------
    1) All_Final_result_HInt.csv
       Detailed table containing score, SignalP and DeepLoc informations for each protein.

    2) Summary_result_HInt.csv
       Compact table listing only the protein name and its global status ("possible hit" or reason for filtering).

    Parameters :
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary

    Returns :
    ----------
    sorted_proteins : list of tuples
    """
    logger.info("Create result table for all preys")
    list_name_baits = list()
    result_dict = file.get_result_dict()
    possible_baits = Informations_dict["Multimer_bait"]
    regions = Informations_dict["Regions"]
    informations = ["DeepLoc","Signal_peptide"]
    nbr_homo = Informations_dict["Homo-oligomer"]
    big_csv_lines = "Name,DeepLoc,Signal_peptide\n"
    small_csv_lines = "Name,Reason_for_filtering\n"
    for protein in Informations_dict["Interact_with"] :
        result_dict.pop(protein, None) #remove bait from result dict
    if Informations_dict["Interact_with"] != [""] : #sorted in function of all baits
        for multimer in possible_baits :
            bait_name = multimer.replace(",","_and_")
            for prot in multimer.split(",") :
                if regions[prot] != "0-0" :
                    start = int(regions[prot].split("-")[0])
                    end = int(regions[prot].split("-")[1])
                    bait_name = bait_name.replace(prot,f"{prot}_{start}-{end}")
            informations.append(f"iQ_score_vs_{bait_name}")
            list_name_baits.append(bait_name)
            big_csv_lines = big_csv_lines.strip("\n") + f",iQ_score_vs_{bait_name}\n"
        sorted_proteins = sorted(result_dict.items(),key=lambda x: (any(f"iQ_score_vs_{bait}" in x[1] for bait in list_name_baits), sum(x[1].get(f"iQ_score_vs_{bait}", 0.0) for bait in list_name_baits)),reverse=True) #sorted in function of key of all baits iQ_score and sum of iQ_score
    if Informations_dict["Homo-oligomer"] != "1" :
        informations.append(f"hiQ_score_{nbr_homo}er")
        big_csv_lines = big_csv_lines.strip("\n") + f",hiQ_score_{nbr_homo}er\n"
    if Informations_dict["Interact_with"] == [""] and Informations_dict["Homo-oligomer"] != "1" : #if no PPI interactions filter on hiQ_score
        sorted_proteins = sorted(result_dict.items(),key=lambda x: (x[1].get(f"hiQ_score_{nbr_homo}er", 0), len(x[1]),), reverse=True)
    if Informations_dict["Homo-oligomer"] == "1" and Informations_dict["Interact_with"] == [""] : #no bait and no homo-oligomer
        sorted_proteins = result_dict

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
    with open("All_Final_result_HInt.csv", "w") as All_result_file :
        All_result_file.write(big_csv_lines)
    with open("Summary_result_HInt.csv", "w") as summary :
        summary.write(small_csv_lines)
    return sorted_proteins


#Generate figures
def Create_figures (file, Informations_dict, AF_version, sorted_proteins, CPU) :
    """
    Generate structural and sequence-level figures for validated prey proteins, parrallelized across multiple CPU cores.

    This function produces figures only for :
    - Monomeric baits (single protein only)
    - Preys that passed all filtering steps (i.e. no "Reason_for_filtering")

    - Distogram (AlphaFold v2 only)
    - Colored PDB structures highlighting interface residues
    - Interface residue tables
    - Clustering of interfaces across ranked proteins
    - Sequence-level interface visualization

    Parameters :
    ----------
    file : object of class File_proteins
    Informations_dict : dict
    AF_version : str
    sorted_proteins : list
    CPU : int
    """
    logger.info("Create figures for all validate preys")
    regions = Informations_dict["Regions"]
    possible_prey = file.get_possible_prey()
    result_dict = file.get_result_dict()
    complete_lenght_prot = file.get_lenght_prot()
    complete_seq_prot = file.get_proteins_sequence_no_SP()
    
    interface_dict = dict()
    tasks = []
    baits_seq = {}
    baits_lenght = {}
    for baits in Informations_dict["Multimer_bait"] :
        bait_file = baits
        for bait in baits.split(",") :
            if regions[bait] != "0-0" :
                start = int(regions[bait].split("-")[0])
                end = int(regions[bait].split("-")[1])
                bait_file = bait_file.replace(bait,f"{bait}_{start}-{end}")
            baits_seq [bait] = complete_seq_prot[bait]
            baits_lenght [bait] = complete_lenght_prot[bait]
        bait_file = bait_file.replace(",","_and_")
        for prey in possible_prey :
            lenght_prot = copy.deepcopy(baits_lenght)
            lenght_prot [prey] = complete_lenght_prot[prey]
            seq_prot = copy.deepcopy(baits_seq)
            seq_prot [prey] = complete_seq_prot[prey]
            if "Reason_for_filtering" not in result_dict[prey].keys() : #only for validate preys
                tasks.append((AF_version, bait_file, prey, lenght_prot, seq_prot, baits, regions))
    if tasks : #if there is interaction to process            
        with Pool(processes=CPU) as pool :
            results_res_int = pool.map(postprocess_interaction, tasks)
        for d in results_res_int :
            for key, list_int in d.items() :
                if key in interface_dict :
                    interface_dict[key] += list_int
                else :
                    interface_dict[key] = list_int.copy()

        interface_dict = cluster_interface(interface_dict, sorted_proteins)
        plot_sequence_interface(file, interface_dict)



def postprocess_interaction (args) : #maybe split first and second part of function
    """
    Post-process a single AlphaFold interaction to generate structural figures and interface residue tables.

    Parameters :
    ----------
    args : tuple

    returns :
    ----------
    interface_dict : dict
    """
    (AF_version, bait_file, prey, lenght_prot, seq_prot, baits, region) = args
    outdir = f"./result_PPI_int/{bait_file}_and_{prey}"
    interface_dict = dict()
    if AF_version == "2" :
        plot_Distogram(outdir)

    residues_at_interface, proteins, path_int, color_res = make_table_res_int(lenght_prot, seq_prot, outdir, baits, AF_version, region)

    if residues_at_interface is not None :
        color_int_residues(path_int, color_res, proteins)
        interface_dict = define_interface(residues_at_interface)

    return interface_dict

def plot_Distogram (job) :
    """
    Generate a distance map (distogram) for the best model of a given AlphaFold job.

    Only generates the distogram if it does not already exist as a PNG.
    Works with both pickled (.pkl) and gzipped (.pkl.gz) result files.

    Parameters :
    ----------
    job : str
    """
    ranking_results = json.load(open(os.path.join(f'{job}/ranking_debug.json')))
    best_model = ranking_results["order"][0]
    del ranking_results
    if os.path.exists(f'{job}/result_{best_model}.dmap.png') == False :
        if os.path.isfile(f'{job}/result_{best_model}.pkl.gz') :
            path_file = f'{job}/result_{best_model}.pkl.gz'
        if os.path.isfile(f'{job}/result_{best_model}.pkl') :
            path_file = f'{job}/result_{best_model}.pkl'
        if path_file.endswith(".gz") :
            with gzip.open(path_file, "rb") as f :
                results = pickle.load(f)
        else :
            with open(path_file, "rb") as f :
                results = pickle.load(f)
        if "distogram" in results.keys() : #avoid error from APD release 
            bin_edges = results["distogram"]["bin_edges"]
            bin_edges = np.insert(bin_edges, 0, 0)
            distogram_softmax = softmax(results["distogram"]["logits"], axis=2)
            dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2)
            np.savetxt(f"{job}/result_{best_model}.pkl.dmap", dist)
            lenght_list = []
            for seq in results["seqs"] :
                lenght_list.append(len(seq))
            logger.info(f"Generate {job.split('/')[2]} Distogram")
            initial_lenght = 0
            fig, ax = plt.subplots()
            d = ax.imshow(dist)
            plt.colorbar(d, ax=ax, fraction=0.046, pad=0.04)
            del dist
            del results
            del distogram_softmax
            del bin_edges
            del d
            gc.collect()
            ax.title.set_text("Distance map")
            for index in range(len(lenght_list)-1) :
                initial_lenght += lenght_list[index]
                ax.axhline(initial_lenght, color="black", linewidth=1.5)
                ax.axvline(initial_lenght, color="black", linewidth=1.5)
            plt.savefig(f"{job}/result_{best_model}.dmap.png", dpi=600)
            plt.close()
            os.remove(f"{job}/result_{best_model}.pkl.dmap")
        

def make_table_res_int (lenght_prot, seq_prot, path_int, baits, AF_version, regions) :
    """
    Generate a detailed table of residue-residue interactions for a protein-protein complex.

    This function analyzes the structural model of a protein complex to identify residues at the interface between the bait and prey proteins. It considers both AlphaFold2 and 
    AlphaFold3 outputs, using distance thresholds and predicted aligned error (PAE) to filter meaningful interactions. The function also prepares color-coding for interface 
    residues for downstream visualization.

    Parameters :
    ----------
    lenght_prot : dict
    seq_prot : dict
    path_int : str
    baits : str
    AF_verison : str
    regions : dict

    Returns :
    ----------
    residues_at_interface : list of lists or None
    proteins : list of str
    path_int : str
    color_res : dict

    Notes :
    ----------
    - For AlphaFold2, interactions are extracted from the predicted aligned error (PAE) matrix and the distogram. Filtered by a distance cutoff (10 Å) and PAE threshold (10 Å).
    - For AlphaFold3, interactions are extracted directly from the atomic coordinates in the PDB file, filtered by a distance cutoff (10 Å) and PAE threshold (10 Å). Only consider standard backbone and Cβ atoms for distance calculations.
    """
    parser = PDB.PDBParser(QUIET=True)
    names_int = path_int.split('/')[2]
    dict_int = dict()
    proteins = [bait for bait in baits.split(",")]
    proteins.append(names_int.split('_and_')[-1])
    color_res = dict()
    dist_k = True
    for prot in proteins :
        color_res[prot] = set()
    if AF_version == "2" :
        ranking_results = json.load(open(os.path.join(f'{path_int}/ranking_debug.json')))
        best_model = ranking_results["order"][0]
        if os.path.isfile(f'{path_int}/result_{best_model}.pkl.gz') :
            path_file = f'{path_int}/result_{best_model}.pkl.gz'
        if os.path.isfile(f'{path_int}/result_{best_model}.pkl') :
            path_file = f'{path_int}/result_{best_model}.pkl'
        with open(os.path.join(path_file), 'rb') as inf_file :
            if ".gz" in path_file :
                pickle_dict = pickle.load(gzip.open(inf_file))
            else :
                pickle_dict = pickle.load(inf_file)

        if "distogram" not in pickle_dict.keys() :
            dist_k = False
        else :
            pae_mtx = pickle_dict['predicted_aligned_error']#take PAE
            bin_edges = pickle_dict["distogram"]["bin_edges"]#take distogram for distance
            bin_edges = np.insert(bin_edges, 0, 0)
            logits = pickle_dict["distogram"]["logits"]
            dist = bin_edges[np.argmax(logits, axis=2)]
            del pickle_dict
            del logits
            del bin_edges
            gc.collect()
            #distogram_softmax = softmax(pickle_dict["distogram"]["logits"], axis=2)
            #dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2) #center of mass of the residue
            complete_lenght = 0
            max_hori_index = 0
            for bait in baits.split(",") :
                complete_lenght += lenght_prot[bait]
            for bait in baits.split(",") :
                min_hori_index = max_hori_index
                max_hori_index += lenght_prot[bait]
                bait_prey = bait +"_and_" + proteins[-1]
                dict_int[bait_prey] = [[bait," "+proteins[-1]," Distance_Ä"," PAE_score"]]
                for line in range(complete_lenght,complete_lenght+lenght_prot[proteins[-1]]) :
                    hori_index = -1
                    for distance in dist[line] :
                        hori_index += 1
                        if hori_index < max_hori_index and hori_index >= min_hori_index :
                            if distance <= 10 :  #center of mass of the residue
                                if pae_mtx[line][hori_index] <= 10 :
                                    real_hori_index = hori_index - min_hori_index
                                    res_in_tot_seq = real_hori_index
                                    if regions[bait] != "0-0" : #if region selected, need to ajust index
                                        res_in_tot_seq = hori_index - min_hori_index + int(regions[bait].split("-")[0]) - 1
                                    residue1 = seq_prot[bait][res_in_tot_seq]
                                    residue2 = seq_prot[proteins[-1]][line-complete_lenght]
                                    dict_int[bait_prey].append([residue1+":"+str(res_in_tot_seq+1)," "+residue2+":"+str(line-complete_lenght+1)," "+str(distance), " "+str(pae_mtx[line][real_hori_index])])
                                    color_res[bait].add(str(res_in_tot_seq+1))
                                    color_res[proteins[-1]].add(str(line-complete_lenght+1))
    
    if AF_version == "3" or dist_k == False : #if no distogram, use only PAE and distance from pdb
        with open(os.path.join(path_int, f'{path_int.split("/")[-1]}_confidences.json'), 'rb') as json_f :
            pae_mtx = np.array(json.load(json_f)['pae'])
        DIST_CUTOFF = 10.0      # Å (CA/CB/C)
        PAE_CUTOFF  = 10.0 #Observation: PAE value for residue at the interaciotn of AF3 model is generally lower than AF2 model
        ATOM_CONTACT = ["C","CA","CB"]

        len_chain_last = lenght_prot[proteins[-1]]
        total_len = pae_mtx.shape[0]
        int_already_know = {}
        structure = parser.get_structure('protein',os.path.join(path_int, f"{names_int}_ranked_0.pdb"))
        for model in structure :
            chains = model.get_list()
            last_chain = chains[-1]
            baits = baits.split(",")
            for i, chain1 in enumerate(chains[:-1]) :
                chain2 = last_chain
                interaction = baits[i] +"_and_"+ proteins[-1]

                if interaction not in dict_int :
                    dict_int[interaction] = [[baits[i], " " + proteins[-1], " Distance_Å", " PAE_score"]]
                for res1 in chain1:
                    if res1.id[0] != " " :
                        continue
                    for res2 in chain2:
                        if res2.id[0] != " " :
                            continue
                        for atom1 in res1:
                            if atom1.get_id() not in ATOM_CONTACT:
                                continue
                            for atom2 in res2 :
                                if atom2.get_id() not in ATOM_CONTACT:
                                    continue
                                dist = atom1 - atom2
                                if dist > DIST_CUTOFF:
                                    continue

                                r1 = res1.id[1] - 1
                                r2 = res2.id[1] - 1

                                idx1 = r1
                                idx2 = total_len - len_chain_last + r2

                                if idx1 >= total_len or idx2 >= total_len:
                                    continue

                                pae_score = float(pae_mtx[idx1, idx2])

                                if pae_score > PAE_CUTOFF :
                                    continue

                                key = (f"{chain1.get_id()}:{res1.get_resname()} {res1.id[1]}", f"{chain2.get_id()}:{res2.get_resname()} {res2.id[1]}")
                                color_res[baits[i]].add(str(res1.id[1]))
                                color_res[proteins[-1]].add(str(res2.id[1]))

                                if key in int_already_know:
                                    if dist < int_already_know[key][0] :
                                        int_already_know[key] = (dist, pae_score)
                                else:
                                    int_already_know[key] = (dist, pae_score)


                for (resA, resB), (dist, pae) in int_already_know.items() :
                   dict_int[interaction].append([f"{resA[2:5]}:{resA.split()[-1]}",f" {resB[2:5]}:{resB.split()[-1]}",f" {dist:.2f}",f" {pae:.2f}"])

    residues_at_interface = dict()
    for chains in dict_int.keys() :
        residues_at_interface[chains] = []
        fileout = chains+"_res_int.csv"
        np_table = np.array(dict_int[chains])
        with open(f"{path_int}/"+fileout, "w", newline="") as csv_table :
            mywriter = csv.writer(csv_table, delimiter=",")
            mywriter.writerows(np_table)
        del dict_int[chains][0] #delete title of each col
        for interaction in dict_int[chains] :
            if interaction not in residues_at_interface[chains] :
                residues_at_interface[chains].append(interaction)
    if residues_at_interface != dict() : #can arrive if it don't find atom with distance < 10 or PAE < 7
        return residues_at_interface,proteins,path_int,color_res
    else :
        return None,None,None,None

def color_int_residues(pdb_path, residues_to_color, names) :
    """
    Color residues involved in protein-protein interactions in a PDB file.
   
    This function modifies the B-factor (temperature factor) column of a PDB file to indicate interface residues.

    Parameters :
    ----------
    pdb_path : str
    residues_to_color : dict
    names : str
    """
    names_int = pdb_path.split('/')[2]
    save_lines = list()
    chain_to_prot = dict()
    chain_index = 0
    with open(f"{pdb_path}/{names_int}_ranked_0.pdb", "r") as file_in : 
        for line in file_in:
            if line.startswith("ATOM") :
                chain = line[21]
                if chain not in chain_to_prot :
                    chain_to_prot[chain] = names[chain_index]
                    chain_index += 1

                prot = chain_to_prot[chain]
                res_num = line[22:26].strip()

                if res_num in residues_to_color.get(prot, set()) :
                    line = line[:60] + "100.00" + line[66:]
                else:
                    line = line[:60] + "  0.00" + line[66:]

            save_lines.append(line)

    with open(f"{pdb_path}/{names_int}_ranked_0.pdb", "w") as writer:
        writer.writelines(save_lines)


def define_interface (residues_at_interface) :
    """
    Create a dictionary of interacting residues for a protein pair.

    This function parses a list of interacting residue pairs and adds them to the existing interface dictionary. Each protein in the pair gets a list of residues 
    that interact with the other protein. The second protein's name is also added to the list for reference.

    Parameters :
    ----------
    residues_at_interface : dict of lists

    Returns :
    ----------
    old_interface_dict : dict
    """
    old_interface_dict = {}
    for inter in residues_at_interface.keys() :
        all_residues_int = residues_at_interface[inter]
        if residues_at_interface[inter] == [] :
            continue
        protein1 = inter.split("_and_")[0]
        protein2 = inter.split("_and_")[1]
        list_int_protein1 = list()
        list_int_protein2 = list()
        if protein1 not in old_interface_dict.keys() :
            old_interface_dict[protein1] = []
        if protein2 not in old_interface_dict.keys() :
            old_interface_dict[protein2] = []
        for line in all_residues_int :
            line[1] = line[1].strip()
            if line[0] != inter[0] and "chain" not in line[0] :
                if line[0].split(":")[1] not in list_int_protein1 :
                    list_int_protein1.append(line[0].split(":")[1])
                if line[1].split(":")[1] not in list_int_protein2 :
                    list_int_protein2.append(line[1].split(":")[1])
        list_int_protein1.append(protein2) #last values of each list is the second proteins
        list_int_protein2.append(protein1)
        old_interface_dict[protein1].append(list_int_protein1)
        old_interface_dict[protein2].append(list_int_protein2)

    return old_interface_dict

def cluster_interface (interface_dict, sorted_proteins) :
    """
    Cluster and classify protein interfaces based on residue overlap.
    
    This function analyzes all detected interaction interfaces for each protein and assigns a letter (a-z) to represent unique interfaces. Interfaces with significant 
    overlap (Jaccard similarity ≥ 0.20) are considered the same and share the same letter. Small interfaces or those with low similarity get a new letter. The function also 
    limits the number of interfaces per protein to 27 betters.

    Parameters :
    ----------
    interface_dict : dict
    sorted_proteins : dict

    Returns :
    ----------
    interface_dict : dict
        Updated dictionary where each interface list starts with a letter (a-z) representing the interface cluster. Interfaces with similar residues share the same letter.
    """
    alphabet = string.ascii_lowercase

    for proteins in interface_dict.keys() :
        if len(interface_dict[proteins]) >= 27 : #Limit to 27 interfaces per protein
            best_preys = list()
            for prey in sorted_proteins :
                best_preys.append(prey[0])
                if len(best_preys) >= 27 :
                    break
            copy_interface = copy.deepcopy(interface_dict[proteins])
            for interface in copy_interface :
                if interface[len(interface)-1] not in best_preys :
                    interface_dict[proteins].remove(interface)
        alpha_index = 0
        already_inter = list()
        interface_dict[proteins] = sorted(interface_dict[proteins], key=lambda x : len(x)) #sorted all interface in function of number of residues
        for interface1 in range(len(interface_dict[proteins])) :
            if interface1 == 0 : #if it's the first interface, define a
                interface_dict[proteins][interface1].insert(0,alphabet[0])
                already_inter.append(alphabet[0])
            for interface2 in range(interface1+1,len(interface_dict[proteins])) :
                alpha_index += 1
                list_inter = list(set(interface_dict[proteins][interface1]).intersection(set(interface_dict[proteins][interface2])))
                simi_inter = len(list_inter)/(len(set(interface_dict[proteins][interface1]).union(set(interface_dict[proteins][interface2])))-3) #indice jaccard # -3 just to remove interface 'a' and uniprotID from .union()
                if simi_inter < 0.20: #create a new interface #Bias on small interface
                    if interface_dict[proteins][interface2][0] in already_inter : #Don't create new interface if it already has one
                        pass
                    else :
                        interface_dict[proteins][interface2].insert(0,alphabet[alpha_index])
                        already_inter.append(alphabet[alpha_index])
                else : #if interfaces got more than 0.20 of same residues, it's the same interface
                    if interface_dict[proteins][interface2][0] in alphabet :
                        interface_dict[proteins][interface2].pop(0)
                    interface_dict[proteins][interface2].insert(0,interface_dict[proteins][interface1][0]) #set same letter for same interface
                    alpha_index -= 1
    return interface_dict
    


def plot_sequence_interface (file, dict_inter) :
    """
    Generate sequence-based interface diagrams for proteins.

    This function visualizes interacting residues along a protein sequence. Residues involved in different interfaces are color-coded, and residues 
    participating in multiple interfaces can show multiple colors. The output is a PNG figure for each protein, showing the sequence with interface annotations.

    Parameters :
    ----------
    file : object of File_proteins class
    dict_inter : dict

    Notes :
    ----------
    - Residues involved in multiple interfaces are displayed with stacked color blocks.
    - Each line of the figure contains up to 150 residues (adjustable by `line_adjust`).
    """
    if not os.path.exists("./Interface_fig/") :
        os.makedirs("./Interface_fig/")
    sequence_dict = file.get_proteins_sequence_no_SP()
    all_color = ['red','green', 'blue', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'pink', 'brown','lime', 'indigo', 'violet', 'turquoise', 'teal', 'crimson', 'gold', 'salmon', 'plum', 'chartreuse']
    for uniprotID_main in dict_inter.keys() :
        sequence = sequence_dict[uniprotID_main]
        indice_color = -1
        interface_done = dict()
        index_to_color = dict()
        uniprot_id_interface = dict()
        for interaction in dict_inter[uniprotID_main] : #list of residue + interface + UniprotID in interaction
            if interaction[0] not in interface_done.keys() : #if it's a new interface
                indice_color += 1
                interface_done[interaction[0]] = all_color[indice_color]
                uniprot_id_interface[interaction[len(interaction)-1]] = all_color[indice_color]
                for aa_to_color in interaction :
                    if " " in aa_to_color :
                        if aa_to_color.split(" ")[1] not in index_to_color.keys() :
                            index_to_color[aa_to_color.split(" ")[1]] = [all_color[indice_color]]
                        if aa_to_color.split(" ")[1] in index_to_color.keys() and all_color[indice_color] not in index_to_color[aa_to_color.split(" ")[1]] : #add two colour if it's in two interface
                            index_to_color[aa_to_color.split(" ")[1]].append(all_color[indice_color])
                    else : #for seconde residue table
                        if aa_to_color not in index_to_color.keys() :
                            index_to_color[aa_to_color] = [all_color[indice_color]]
                        if aa_to_color in index_to_color.keys() and all_color[indice_color] not in index_to_color[aa_to_color] : #add two colour if it's in two interface
                            index_to_color[aa_to_color].append(all_color[indice_color])
            else :
                uniprot_id_interface[interaction[len(interaction)-1]] = interface_done[interaction[0]]
                for aa_to_color in interaction :
                    if " " in aa_to_color :
                        if aa_to_color.split(" ")[1] not in index_to_color.keys() :
                            index_to_color[aa_to_color.split(" ")[1]] = [all_color[indice_color]]
                        if aa_to_color.split(" ")[1] in index_to_color.keys() and all_color[indice_color] not in index_to_color[aa_to_color.split(" ")[1]] : #add two colour if it's in two interface
                            index_to_color[aa_to_color.split(" ")[1]].append(all_color[indice_color])
                    else : #for seconde residue table
                        if aa_to_color not in index_to_color.keys() :
                            index_to_color[aa_to_color] = [interface_done[interaction[0]]]
                        if aa_to_color in index_to_color.keys() and interface_done[interaction[0]] not in index_to_color[aa_to_color] : #add two colour if it's in two interface
                            index_to_color[aa_to_color].append(interface_done[interaction[0]])
        line_adjust = 150 #max aa per line
        dict_name = dict()
        n_lines = (len(sequence) + line_adjust - 1) // line_adjust
        fig, ax = plt.subplots(figsize=(line_adjust / 4, n_lines*1.5)) #Adjust figsize
        for line_index in range(0, len(sequence), line_adjust) :
            sub_sequence = sequence[line_index:line_index + line_adjust]
            y_pos = -line_index // line_adjust * 1.5
            for i in range(len(sub_sequence)) :
                aa = sub_sequence[i]
                total_index = line_index + i
                if str(total_index + 1) in index_to_color.keys() :
                    colors = index_to_color[str(total_index + 1)]
                    height = 0.5 / len(colors)
                    for color_index, color in enumerate(colors) :
                        ax.add_patch(plt.Rectangle((i, y_pos + color_index * height), 1, height, color=color))
                    ax.text(i + 0.5, y_pos + 0.25, aa, ha='center', va='center', color='white')
                else :
                    ax.add_patch(plt.Rectangle((i, y_pos), 1, 0.6, color="white"))
                    ax.text(i + 0.5, y_pos + 0.25, aa, ha='center', va='center', color='black')
                if (total_index+1) % 10 == 0 or i == 0:
                    ax.text(i + 0.5, y_pos + 0.5, str(total_index + 1), ha='center', va='center', color='black', fontsize=7)
        for index_neigh, neigh in enumerate(uniprot_id_interface) :
            name_neigh = f"{neigh}({dict_name[neigh]})" if neigh in dict_name else neigh
            ax.text(index_neigh * 6, - n_lines * 2, name_neigh, ha='center', va='center', color=uniprot_id_interface[neigh], fontsize=8)
        uniprotID_main_name = f"{uniprotID_main}({dict_name[uniprotID_main]})" if uniprotID_main in dict_name else uniprotID_main
        ax.text(-2, 0.25, uniprotID_main_name, ha='right', va='center', color='black', fontsize=10, fontweight='bold')
        ax.set_xlim(0, line_adjust)
        ax.set_ylim(-n_lines*2, 1)  #Adjust high
        ax.axis('off')
        plt.savefig("./Interface_fig/"+uniprotID_main+"_interface_fig.png", dpi=300, bbox_inches='tight')
        plt.close()
