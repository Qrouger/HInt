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


# Configure global logger
logging.basicConfig(
    filename="HInt.log",  # Log file name
    level=logging.INFO,  # Log level
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

logger = logging.getLogger()

def Score_interaction (file, Informations_dict, Interaction, bait=None) :
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
    seq_no_SP = file.get_proteins_sequence_no_SP()
    new_possible_prey = list()
    Path_ccp4 = Informations_dict["Path_ccp4"]
    regions = Informations_dict["Regions"]
    AF_version = Informations_dict["AlphaFold"]
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
            for protein in possible_prey :
                if f"iQ_score_vs_{bait}" not in int_score[protein].keys() :
                    ppi_list.append(f"./result_{Interaction}/{bait}_and_{protein}")
                else :
                    result_dict[protein][f"iQ_score_vs_{bait}"] = int_score[protein][f"iQ_score_vs_{bait}"]
                    if int_score[protein][f"iQ_score_vs_{bait}"] > 0 :
                        new_possible_prey.append(protein)

        results = []
        with multiprocessing.Pool(N_CPU) as pool :
            tasks = [(ppi, "./", Path_ccp4, seq_no_SP, AF_version) for ppi in ppi_list]
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
    interaction, output_dir, Path_ccp4 ,seq_no_SP ,AF_version = args
    result = get_good_inter_pae.main(interaction, output_dir, 10, 2, Path_ccp4, seq_no_SP, AF_version)
    return result

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
    small_csv_lines = "Name,Reason_for_filtering\n"
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