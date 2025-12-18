#!/usr/bin/env python3
#Adapted from get_good_inter_pae.py (https://github.com/KosinskiLab/AlphaPulldown/blob/main/alphapulldown/analysis_pipeline/alpha_analysis_jax0.4.def)

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "script_pi_score"))
from datetime import datetime
from absl import flags,app,logging
from calculate_mpdockq import *
import pickle
import json
import numpy as np
import pandas as pd
import subprocess
import gzip
from Bio.PDB import MMCIFParser, PDBIO
import shutil
import gzip

def examine_inter_pae(pae_mtx,seqs,cutoff) :
    """A function that checks inter-pae values in multimer prediction jobs"""
    lens = [len(seq) for seq in seqs]
    old_lenth=0
    for length in lens:
        new_length = old_lenth + length
        pae_mtx[old_lenth:new_length,old_lenth:new_length] = 50
        old_lenth = new_length
    check = np.where(pae_mtx<cutoff)[0].size !=0

    return check


def obtain_mpdockq(work_dir,pkl_dict=None) :
    """Returns chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path"""
    pdb_path = os.path.join(work_dir,'ranked_0.pdb')
    pdb_chains, chain_coords, chain_CA_inds, chain_CB_inds = read_pdb(pdb_path)
    if pkl_dict ==  None :
        best_plddt = extract_plddt_from_pdb(pdb_path)
    else :
        best_plddt = pkl_dict['plddt'] 
    plddt_per_chain = read_plddt(best_plddt,chain_CA_inds)
    return chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path

def obtain_mpdockq2(chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path) :
    """Returns mpDockQ if more than two chains otherwise return pDockQ"""
    complex_score,num_chains = score_complex(chain_coords,chain_CB_inds,plddt_per_chain)
    if complex_score is not None and num_chains>2:
        mpDockq_or_pdockq = calculate_mpDockQ(complex_score)
    elif complex_score is not None and num_chains==2:
        chain_coords,plddt_per_chain = read_pdb_pdockq(pdb_path)
        mpDockq_or_pdockq = calc_pdockq(chain_coords,plddt_per_chain,t=8)
    else:
        mpDockq_or_pdockq = "None"
    return mpDockq_or_pdockq


def extract_plddt_from_pdb(pdb_file):
    """
    Extract plddt from b-factor in pdb
    """
    plddt_values = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    b_factor = float(line[60:66].strip())
                    plddt_values.append(b_factor)
                except ValueError:
                    pass
    return np.array(plddt_values, dtype=float)

def run_and_summarise_pi_score(work_dir, jobs, surface_thres, ccp4_setup) :
    """
    A function to calculate all predicted models' pi_scores and make a pandas df of the results.
    Instrumented to log timing per major step.
    """
    result = subprocess.run(["conda", "env", "list", "--json"],
                            capture_output=True, text=True, check=True)
    envs = json.loads(result.stdout)["envs"]
    exists = any("pi_score" in env for env in envs)
    if not exists:
        subprocess.run([
            "conda", "create", "-y", "-n", "pi_score",
            "python=2.7", "scikit-learn=0.20.4", "biopython", "biopandas"
        ], check=True)

    direc = jobs[0].split("/")[len(jobs[0].split("/"))-1]
    if os.path.isdir("/scratch") :
        tmp_dir = f"/scratch/tmp/{direc}"
    else :
        tmp_dir = f"/tmp/{direc}"
    try:
        shutil.rmtree(f"{tmp_dir}/pi_score_outputs")
    except:
        pass
    print (f"Creating temporary directory {tmp_dir} for pi_score outputs")
    subprocess.run(f"mkdir -p {tmp_dir}/pi_score_outputs",
                   shell=True, executable='/bin/bash')
    pi_score_outputs = os.path.join(tmp_dir, "pi_score_outputs")

    for job in jobs:
        if not os.path.isfile(os.path.join(job, "ranked_0.pdb")):
            print(f"{job} failed. Cannot find ranked_0.pdb in {job}")
            sys.exit()

        pdb_path = os.path.join(job, "ranked_0.pdb")
        output_dir = os.path.join(pi_score_outputs)
        cmd = (
            f"source {ccp4_setup}/bin/ccp4.setup-sh && "
            f"conda run -n pi_score python ./script_pi_score/run_piscore_wc.py "
            f"-p {pdb_path} -o {output_dir} -s {surface_thres} -ps 10"
        )

        proc = subprocess.Popen(cmd, shell=True, executable="/bin/bash", close_fds=True)
        proc.wait()
    output_df = pd.DataFrame()
    for job in jobs:
        name_job = job.split("/")[-1]
        subdir = os.path.join(pi_score_outputs)
        csv_files = [f for f in os.listdir(subdir) if 'filter_intf_features' in f]
        pi_score_files = [f for f in os.listdir(subdir) if 'pi_score_' in f]

        if not csv_files or not pi_score_files:
            print(f"Warning: missing CSV or pi_score files for {name_job}")
            continue

        filtered_df = pd.read_csv(os.path.join(subdir, csv_files[0]))


        if filtered_df.shape[0] == 0:
            for column in filtered_df.columns:
                filtered_df[column] = ["None"]
            filtered_df['jobs'] = str(name_job)
            filtered_df['pi_score'] = "No interface detected"
        else:
            with open(os.path.join(subdir, pi_score_files[0]), 'r') as f:
                lines = [l for l in f.readlines() if "#" not in l]
            if len(lines) > 0:
                pi_score = pd.read_csv(os.path.join(subdir, pi_score_files[0]))
                pi_score['jobs'] = str(name_job)
            else:
                pi_score = pd.DataFrame.from_dict({"pi_score": ['SC:  mds: too many atoms']})

            pi_score['interface'] = pi_score.get('chains', None)
            filtered_df['jobs'] = str(name_job)

            filtered_df = pd.merge(filtered_df, pi_score, on=['jobs', 'interface'], how='left')

            try:
                filtered_df = filtered_df.drop(columns=["#PDB", "pdb", "pvalue", "chains", "predicted_class"])
            except:
                pass

        output_df = pd.concat([output_df, filtered_df])

    subprocess.run(f"rm -rf {tmp_dir}", shell=True, executable='/bin/bash')
    return output_df
    
    

def main(job, output_dir, cutoff, surface_thres, ccp4_setup, seq_no_SP ,AF_version) :
    good_jobs = []
    iptm_ptm = list()
    iptm = list()
    mpDockq_scores = list()
    logging.info(f"Scoring {job}")
    result_subdir = os.path.join(job)
    if AF_version == "3" : #for alphafold3
        interaction = job.split("/")[-1]
        bait = interaction.split("_and_")[0]
        if "_" in bait :
            bait = bait.split("_")[0]
        prey = interaction.split("_and_")[1]
        if os.path.isfile(os.path.join(result_subdir,'ranked_0.pdb')) == False : #create ranked_0.pdb for AF3
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure('model', os.path.join(result_subdir,'ranked_0_model.cif'))
            io = PDBIO()
            io.set_structure(structure)
            io.save(os.path.join(result_subdir,'ranked_0.pdb'))
        if os.path.isfile(os.path.join(job,'ranked_0_summary_confidences.json')):
            with open(os.path.join(result_subdir,'ranked_0_summary_confidences.json'),'rb') as json_sum_f :
                json_sum = json.load(json_sum_f)
            if "iptm" in json_sum.keys() and "ptm" in json_sum.keys():
                iptm_score = json_sum['iptm']
                ptm_score = json_sum['ptm']
                iptm_ptm_score = 0.8 * iptm_score + 0.2 * ptm_score
            with open(os.path.join(result_subdir, 'ranked_0_confidences.json'),'rb') as json_f :
                json_data = json.load(json_f)
            pae_list = json_data['pae']
            pae_mtx = np.array(pae_list)
            seqs = [seq_no_SP[bait],seq_no_SP[prey]]
            chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path = obtain_mpdockq(os.path.join(job))
            check = examine_inter_pae(pae_mtx,seqs,cutoff=cutoff)
            mpDockq_score = obtain_mpdockq2(chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path)
            if check:
                good_jobs.append(str(job))
                iptm_ptm.append(iptm_ptm_score)
                iptm.append(iptm_score)
                mpDockq_scores.append(mpDockq_score)
                


    if os.path.isfile(os.path.join(job,'ranking_debug.json')): #for alphafold2
        with open(os.path.join(result_subdir,'ranking_debug.json'),'rb') as json_f :
            data = json.load(json_f)
        best_model = data['order'][0]
        if "iptm" in data.keys() or "iptm+ptm" in data.keys():
            iptm_ptm_score = data['iptm+ptm'][best_model]
            if os.path.exists(os.path.join(result_subdir, f"result_{best_model}.pkl")) :
                pkl_path = os.path.join(result_subdir, f"result_{best_model}.pkl")
                with open(pkl_path, 'rb') as pkl :
                    check_dict = pickle.load(pkl)
            elif os.path.exists(os.path.join(result_subdir, f"result_{best_model}.pkl.gz")) :
                print("result pickle for the best model not found. Now search for zipped pickle.")
                pkl_path = os.path.join(result_subdir, f"result_{best_model}.pkl.gz")
                with gzip.open(pkl_path, 'rb') as pkl :
                    check_dict = pickle.load(pkl)
            else :
                logging.info(f"Cannot find result pickle for {job}, skipping.")
            seqs = check_dict['seqs']
            iptm_score = check_dict['iptm']
            pae_mtx = check_dict['predicted_aligned_error']
            chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path = obtain_mpdockq(os.path.join(job),check_dict)
            check = examine_inter_pae(pae_mtx,seqs,cutoff=cutoff)
            mpDockq_score = obtain_mpdockq2(chain_coords,chain_CB_inds,plddt_per_chain,best_plddt,pdb_path)
            if check:
                good_jobs.append(str(job))
                iptm_ptm.append(iptm_ptm_score)
                iptm.append(iptm_score)
                mpDockq_scores.append(mpDockq_score)
    other_measurements_df=pd.DataFrame.from_dict({
        "jobs":job.split("/")[-1],
        "iptm_ptm":iptm_ptm,
        "iptm":iptm,
        "mpDockQ/pDockQ":mpDockq_scores})

        
    if good_jobs!=[] :
        pi_score_df = run_and_summarise_pi_score(output_dir,good_jobs,surface_thres,ccp4_setup)
        pi_score_df = pd.merge(pi_score_df,other_measurements_df,on="jobs")
        columns = list(pi_score_df.columns.values)
        columns.pop(columns.index('jobs'))
        pi_score_df = pi_score_df[['jobs'] + columns]
        pi_score_df = pi_score_df.sort_values(by='iptm',ascending=False)
        return pi_score_df
    else :
        return None

