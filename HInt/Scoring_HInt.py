""" Scoring file of HInt

    Author: Quentin Rouger
"""
import os
import csv
import logging
import multiprocessing
from tqdm import tqdm
from datetime import datetime
import get_good_inter_pae
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
from pathlib import Path

# Configure global logger
logging.basicConfig(
    filename="HInt.log",  # Log file name
    level=logging.INFO,  # Log level
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

logger = logging.getLogger()

def Score_interaction (file, Informations_dict, CPU, Interaction, bait=None) :
    """
    Generate scores for all interactions and set a list of possible prey based on the scores.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary
    CPU : int
    Interaction : string
    bait : string

    Returns:
    ----------
    """
    start_time = datetime.now()
    result_dict = file.get_result_dict()
    possible_prey = file.get_possible_prey()
    int_score = file.get_int_score() #save scores
    seq_no_SP = file.get_proteins_sequence_no_SP()
    new_possible_prey = list()
    Path_ccp4 = Informations_dict["Path_ccp4"]
    regions = Informations_dict["Regions"]
    AF_version = Informations_dict["AlphaFold"]
    ppi_list = list()

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
            for direc in os.listdir(f"./result_{Interaction}") :
                if "_homo_" in direc :
                    ppi_list.append(f"./result_{Interaction}/{direc}") #Found a solution for score only homo-oligomer without score
        else : #for one vs all
            for protein in possible_prey :
                if f"iQ_score_vs_{bait_name}" not in int_score[protein].keys() :
                    ppi_list.append(f"./result_{Interaction}/{bait_name}_and_{protein}")
                else :
                    result_dict[protein][f"iQ_score_vs_{bait_name}"] = int_score[protein][f"iQ_score_vs_{bait_name}"]
                    if int_score[protein][f"iQ_score_vs_{bait_name}"] > 0 :
                        new_possible_prey.append(protein)
                    if int_score[protein][f"iQ_score_vs_{bait_name}"] == 0 :
                        result_dict[protein]["Reason_for_filtering"] = f"Bad interactions with {bait_name} : inter PAE > 10 A"
        results = []
        with multiprocessing.Pool(CPU) as pool : #just run scoring for interactions without score
            tasks = [(ppi, "./", Path_ccp4, seq_no_SP, AF_version) for ppi in ppi_list]
            results_iter = pool.imap_unordered(run_scoring, tasks)
            for df in tqdm(results_iter, total=len(ppi_list), desc="Scoring interactions") :
                if df is not None and not df.empty :
                    results.append(df)
            pool.close()
            pool.join()

        if results : #Make a dataframe if there is results
            merged_df = pd.concat(results, ignore_index=True)
            merged_df.to_csv(os.path.join(f"./result_{Interaction}", "predictions_with_good_interpae.csv"), index=False)

        if ppi_list : #Add result in dict result if there is new ppi scored
            #Resume all score and set new possible prey
            with open(f"result_{Interaction}/predictions_with_good_interpae.csv", "r") as result_file :
                reader = csv.DictReader(result_file)

                #For one vs all
                if Interaction == "PPI_int" : #make int_score
                    all_lines = "jobs,pi_score,iptm_ptm,pDockQ,iQ_score\n"
                    for row in reader :
                        job = row['jobs']
                        if '_and_' in job and job.split("_and_")[-1] in possible_prey and bait_name in job : #check if interaction is a PPI and if prey is in possible prey list
                            if row['pi_score'] == 'No interface detected' :
                                iQ_score = float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30 #pi_score don't detect interface so it's set on -2.63
                                line =f'{row["jobs"]},-2.63,{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'
                            else :
                                iQ_score = ((float(row['pi_score'])+2.63)/5.26)*40+float(row['iptm_ptm'])*30+float(row['mpDockQ/pDockQ'])*30
                                line =f'{row["jobs"]},{row["pi_score"]},{row["iptm_ptm"]},{row["mpDockQ/pDockQ"]},{str(iQ_score)}\n'

                            int_score[job.split("_and_")[-1]][f"iQ_score_vs_{bait_name}"] = iQ_score
                            result_dict[job.split("_and_")[-1]][f"iQ_score_vs_{bait_name}"] = iQ_score
                            new_possible_prey.append(job.split("_and_")[-1])
                            all_lines = all_lines + line
                            name_int = job.split("/")[-1]
                            os.system(f"cp result_{Interaction}/{job}/ranked_0.pdb result_{Interaction}/{job}/{name_int}_ranked_0.pdb") #rename pdb file with explicit name
                    for protein in possible_prey :
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
                        number_oligo = row["jobs"].split("_")[2].replace("er","") #AFPD 2.0.4 #problem with "_" in name sequence ?
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
    interaction, output_dir, Path_ccp4 ,seq_no_SP ,AF_version = args
    try:
        result = get_good_inter_pae.main(interaction, output_dir, 10, 2, Path_ccp4, seq_no_SP, AF_version) #normal PAE is 10
        return result
    except Exception as e:
        pid = os.getpid()
        logger.error(f"ERROR in worker PID={pid}")
        logger.error(f"Interaction: {interaction}")
        traceback.print_exc()
        raise

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
    logger.info("Create result table for all preys")
    iQ_score_dict = dict()
    list_name_baits = list()
    result_dict = file.get_result_dict()
    possible_prey = file.get_possible_prey()
    possible_baits = Informations_dict["Multimer_bait"]
    regions = Informations_dict["Regions"]
    proteins = file.get_proteins()
    informations = ["DeepLoc","Signal_peptide"]
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
            if len(possible_baits) == 1 :
                big_csv_lines = big_csv_lines.strip("\n") + ",iQ_score\n"
            else :
                big_csv_lines = big_csv_lines.strip("\n") + f",iQ_score_vs_{bait_name}\n"
        sorted_proteins = sorted(result_dict.items(),key=lambda x: (
        any(f"iQ_score_vs_{bait}" in x[1] for bait in list_name_baits),
        sum(x[1].get(f"iQ_score_vs_{bait}", 0.0) for bait in list_name_baits)),reverse=True) #sorted in function of key of all baits iQ_score and sum of iQ_score
    if Informations_dict["Homo-oligomer"] != "1" :
        informations.append("hiQ_score")
        big_csv_lines = big_csv_lines.strip("\n") + "hiQ_score\n"
    if Informations_dict["Interact_with"] == [""] and Informations_dict["Homo-oligomer"] != "1" : #if no PPI interactions filter on hiQ_score
        sorted_proteins = sorted(result_dict.items(),key=lambda x: (x[1].get("hiQ_score", 0), len(x[1]),), reverse=True)
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
    with open("All_Final_result.csv", "w") as All_result_file :
        All_result_file.write(big_csv_lines)
    with open("Summary_result.csv", "w") as summary :
        summary.write(small_csv_lines)


#Generate figures
def Create_figures (file, Informations_dict, AF_version) :
    """
    Create figures for all validate preys.

    Parameters:
    ----------
    file : object of class File_proteins
    Informations_dict : dictionary
    AF_version : string

    Returns:
    ----------
    """
    logger.info("Create figures for all validate preys")
    regions = Informations_dict["Regions"]
    possible_prey = file.get_possible_prey()
    result_dict = file.get_result_dict()
    interface_dict = dict()
    for bait in Informations_dict["Multimer_bait"] :
        if len(bait.split(",")) > 1 or bait == "" : #No figures for multimers bait
            pass
        else :
            if regions[bait] != "0-0" :
                start = int(regions[bait].split("-")[0])
                end = int(regions[bait].split("-")[1])
                bait_file = f"{bait}_{start}-{end}"
            else :
                bait_file = bait
            for prey in possible_prey :
                if "Reason_for_filtering" not in result_dict[prey].keys() : #only for validate preys
                    if AF_version == "2" :
                        plot_Distogram (f"./result_PPI_int/{bait_file}_and_{prey}")
                    residues_at_interface,proteins,path_int,color_res = make_table_res_int(file, f"./result_PPI_int/{bait_file}_and_{prey}", bait, AF_version)
                    if residues_at_interface is not None :
                        color_int_residues(path_int, color_res, proteins, AF_version) #color residue in interaction on the pdb
                        interface_dict = define_interface(residues_at_interface, [bait,prey], interface_dict)
            interface_dict = cluster_interface(interface_dict)
            plot_sequence_interface(file,interface_dict)

def plot_Distogram (job) :
    """
    Generate distogram, only for best models.

    Parameters:
    ----------
    job : string
    
    Returns:
    ----------
    """
    ranking_results = json.load(open(os.path.join(f'{job}/ranking_debug.json')))
    best_model = ranking_results["order"][0]
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
    if "distogram" in results.keys() and os.path.exists(f'{job}/result_{best_model}.dmap.png') == False : #avoid error from APD release 
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
        ax.title.set_text("Distance map")
        for index in range(len(lenght_list)-1) :
           initial_lenght += lenght_list[index]
           ax.axhline(initial_lenght, color="black", linewidth=1.5)
           ax.axvline(initial_lenght, color="black", linewidth=1.5)
        plt.savefig(f"{job}/result_{best_model}.dmap.png", dpi=600)
        plt.close()
        os.remove(f"{job}/result_{best_model}.pkl.dmap")
        del dist
        del results
        del distogram_softmax
        del bin_edges
        del d
        gc.collect()

def make_table_res_int (file, path_int, bait, AF_version) :
    """
    Generate a table of residues in interactions.

    Parameters:
    ----------
    file : object of class File_proteins
    path_int : string
    bait : string
    AF_verison : string

    Returns:
    ----------
    residues_at_interface[names_int] : list
    proteins : string
    path_int : string
    color_res : dict
    """
    parser = PDB.PDBParser(QUIET=True)
    names_int = path_int.split('/')[2]
    dict_int = dict()
    proteins = [bait]
    proteins.append(names_int.split('_and_')[1])
    color_res = dict()
    color_res[proteins[0]] = set()
    color_res[proteins[1]] = set()
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
        lenght_prot = file.get_lenght_prot()
        seq_prot = file.get_proteins_sequence_no_SP()
        dict_int = dict()
        color_res = dict()
        color_res[proteins[0]] = set()
        color_res[proteins[1]] = set()
        pae_mtx = pickle_dict['predicted_aligned_error']#take PAE
        bin_edges = pickle_dict["distogram"]["bin_edges"]#take distogram for distance
        bin_edges = np.insert(bin_edges, 0, 0)
        distogram_softmax = softmax(pickle_dict["distogram"]["logits"], axis=2)
        dist = np.sum(np.multiply(distogram_softmax, bin_edges), axis=2) #center of mass of the residue
        dict_int[names_int] = [[proteins[0]," "+proteins[1]," Distance_Ä"," PAE_score"]]
        for line in range(lenght_prot[proteins[0]],lenght_prot[proteins[0]]+lenght_prot[proteins[1]]) :
            hori_index = -1
            for distance in dist[line] :
                hori_index += 1
                if hori_index < lenght_prot[proteins[0]] :
                    if distance <= 10 :  #center of mass of the residue
                        if pae_mtx[line][hori_index] < 7 :
                            residue1 = seq_prot[proteins[0]][hori_index]
                            residue2 = seq_prot[proteins[1]][line-lenght_prot[proteins[0]]]
                            dict_int[names_int].append([residue1+":"+str(hori_index+1)," "+residue2+":"+str(line-lenght_prot[proteins[0]]+1)," "+str(distance), " "+str(pae_mtx[line][hori_index])])
                            color_res[proteins[0]].add(str(hori_index+1))
                            color_res[proteins[1]].add(str(line-lenght_prot[proteins[0]]+1))

    if AF_version == "3" :
        with open(os.path.join(path_int, 'ranked_0_confidences.json'), 'rb') as json_f:
            pae_mtx = np.array(json.load(json_f)['pae'])
        DIST_CUTOFF = 10.0      # Å (CA/CB/C)
        PAE_CUTOFF  = 10.0 #Observation: PAE value for residue at the interaciotn of AF3 model is generally lower than AF2 model
        ATOM_CONTACT = ["C","CA","CB"]
        len_chainA = file.get_lenght_prot_no_SP()[proteins[0]]
        total_len = pae_mtx.shape[0]
        int_already_know = {}
        dict_int = {}
        structure = parser.get_structure('protein',os.path.join(path_int, f"{names_int}_ranked_0.pdb"))
        for model in structure:
            chains = model.get_list()
            for i, chain1 in enumerate(chains):
                for chain2 in chains[i+1:]:
                    interaction = chain1.get_id() + chain2.get_id()
                    if interaction not in dict_int:
                        dict_int[interaction] = [[proteins[0], " " + proteins[1], " Distance_Å", " PAE_score"]]
                    for res1 in chain1:
                        if res1.id[0] != " ":
                            continue
                        for res2 in chain2:
                            if res2.id[0] != " ":
                                continue
                            for atom1 in res1:
                                if atom1.get_id() not in ATOM_CONTACT:
                                    continue
                                for atom2 in res2:
                                    if atom2.get_id() not in ATOM_CONTACT:
                                        continue
                                    dist = atom1 - atom2
                                    if dist > DIST_CUTOFF:
                                        continue

                                    r1 = res1.id[1] - 1
                                    r2 = res2.id[1] - 1
                                    if interaction[0] == "A":  # bait = chain A
                                        idx1 = r1
                                        idx2 = len_chainA + r2
                                    else:
                                        idx1 = len_chainA + r1
                                        idx2 = r2
                                    if idx1 >= total_len or idx2 >= total_len:
                                        continue
                                    pae_score = float(pae_mtx[idx1, idx2])
                                    if pae_score > PAE_CUTOFF:
                                        continue
                                    key = (f"{chain1.get_id()}:{res1.get_resname()} {res1.id[1]}",f"{chain2.get_id()}:{res2.get_resname()} {res2.id[1]}")
                                    color_res[proteins[0]].add(str(res1.id[1]))
                                    color_res[proteins[1]].add(str(res2.id[1]))
                                    # Keep the interaction with the lowest distance
                                    if key in int_already_know:
                                        if dist < int_already_know[key][0]:
                                            int_already_know[key] = (dist, pae_score)
                                    else:
                                        int_already_know[key] = (dist, pae_score)

        for (resA, resB), (dist, pae) in int_already_know.items():
            interaction = resA[0] + resB[0]
            dict_int[interaction].append([f"{resA[2:5]}:{resA.split()[-1]}",f" {resB[2:5]}:{resB.split()[-1]}",f" {dist:.2f}",f" {pae:.2f}"])

    residues_at_interface = dict()
    residues_at_interface[names_int] = []
    for chains in dict_int.keys() :
        fileout = chains+"_res_int.csv"
        np_table = np.array(dict_int[chains])
        with open(f"{path_int}/"+fileout, "w", newline="") as csv_table :
            mywriter = csv.writer(csv_table, delimiter=",")
            mywriter.writerows(np_table)
        del dict_int[chains][0] #delete title of each col
        for interaction in dict_int[chains] :
            if interaction not in residues_at_interface[names_int] :
                residues_at_interface[names_int].append(interaction)
    if residues_at_interface[names_int] != [] : #can arrive if it don't find atom with distance < 10 or PAE < 7
        return residues_at_interface[names_int],proteins,path_int,color_res
    else :
        return None,None,None,None

def color_int_residues(pdb_path, residues_to_color, names, AF_version) :
    """
    Color residues in interaction in a PDB file.
   
    Parameters:
    ----------
    pdb_path : string
    residues_to_color : dict
    names : string
    AF_version : string
    
    Returns:
    ----------
    """
    if AF_version == "3" :
        chain1 = "A"
    if AF_version == "2" :
        chain1 = "B"
    names_int = pdb_path.split('/')[2]
    name_prot = names[0]
    save_line = str()
    with open(f'{pdb_path}/{names_int}_ranked_0.pdb', 'r') as file :
        for line in file:
            if line.startswith("ATOM") :
                chain2 = line[21]
                if chain1 != chain2 :
                   name_prot = names[1] #use new dict to color atoms
                res_num = line[22:26].strip()
                if res_num in residues_to_color[name_prot] : #change B-factor in color interaction residue
                    line = line[:60] + " 100  " + line[66:]
                else :
                    line = line[:60] + " 0    " + line[66:]
                chain1 = line[21]
            save_line += line
    with open(f'{pdb_path}/{names_int}_ranked_0.pdb', 'w') as writer:
        writer.write(save_line)

def define_interface (list_of_list_int, int, old_interface_dict) :
    """
    Create a dictionary of all interacting residues.

    Parameters:
    ----------
    list_of_list_int : list
    int : string
    old_interface_dict : dict

    Returns:
    ----------
    old_interface_dict : dict
    """
    all_residues_int = copy.deepcopy(list_of_list_int)
    protein1 = int[0]
    protein2 = int[1]
    list_int_protein1 = list()
    list_int_protein2 = list()
    if protein1 not in old_interface_dict.keys() :
        old_interface_dict[protein1] = []
    if protein2 not in old_interface_dict.keys() :
        old_interface_dict[protein2] = []
    for line in all_residues_int :
        line[1] = line[1].strip()
        if line[0] != int[0] and "chain" not in line[0] :
            if line[0].split(":")[1] not in list_int_protein1 :
                list_int_protein1.append(line[0].split(":")[1])
            if line[1].split(":")[1] not in list_int_protein2 :
                list_int_protein2.append(line[1].split(":")[1])
    list_int_protein1.append(protein2) #last values of each list is the second proteins
    list_int_protein2.append(protein1)
    old_interface_dict[protein1].append(list_int_protein1)
    old_interface_dict[protein2].append(list_int_protein2)
    return old_interface_dict

def cluster_interface (interface_dict) :
    """
    Compare interface in function of smaller interface and classify them with a letter representing the interface.
   
    Parameters:
    ----------
    interface_dict : dict

    Returns:
    ----------
    interface_dict : dict
    """
    alphabet = string.ascii_lowercase
    for proteins in interface_dict.keys() :
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
    Generated figures for interface in one sequence.

    Parameters:
    ----------
    file : object of File_proteins class
    dict_inter : dict

    Returns:
    ----------
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