import os
import pickle
import subprocess
import numpy as np

def define_informations() :
    """
    Extract all paths from HInt.txt and store them in a dictionary.
    
    Parameters:
    ----------

    Returns:
    ----------
    informations_dict : dict
    """
    informations_dict = dict()
    with open("HInt.txt", "r") as file :
        for lines in file :
            if ":" in lines :
                informations_name = lines.split(":")[0].strip().strip("\n")
                informations = lines.split(":")[1].strip().strip("\n")
                informations_dict[informations_name] = informations
    for informations_key in informations_dict.keys() :
        if len(informations_dict[informations_key]) == 0 :
            print (f'Informations for {informations_key} is empty')
            if informations_key == "Path_Uniprot_ID" :
                exit()
            elif informations_key == "Path_AlphaFold_Data" :
                print("set by default on ./alphadata")
                informations_dict[informations_key] = "./alphadata"
            elif informations_key == "Path_Pickle_Feature" :
                print("set by default on ./feature")
                informations_dict[informations_key] = "./feature"
            elif informations_key == "Signal_peptide" :
                informations_dict[informations_key] = "None"
            elif informations_key == "Homo-oligomer" :
                informations_dict[informations_key] = 1
            elif informations_key == "Organism" :
                exit()
    if len(informations_dict["Signal_peptide"]) == 0 and len(informations_dict["Homo-oligomer"]) == 0 and len(informations_dict["Interact_with"]) == 0 : #no info
        print("need information to discriminate the potential homologue")
        exit()
    return(informations_dict)

def remove_SP (file, org) :
    """
    Create a new FASTA file without the signal peptide using SignalP.

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
    cmd = "signalp -fasta " + fasta_file + " -org " + org
    os.system(cmd)
    file_signalp = fasta_file.replace(".fasta","_summary.signalp5")
    with open(file_signalp,"r") as fh :
        for line in fh :
            new_line = line.split("\t")
            if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                prot_SP[new_line[0]] = new_line[len(new_line)-1].split("-")[1].split(".")[0]
    new_fasta_dict = dict()
    with open(fasta_file, "r") as fa_file :
        for line2 in fa_file :
            new_line2 = line2
            if SP_signal == 0 and line2[0] != ">" :
                new_line2 = line2
                new_fasta_dict[save_key] = line2
            if int(SP_signal) > 0 :
                new_line2 = line2[int(SP_signal)-1:len(line2)]
                new_fasta_dict[save_key] = line2[int(SP_signal)+1:len(line2)]
                SP_signal = 0
            if line2[0] == ">" :
                save_key = line2[1:len(line2)-1]
                if str(line2[1:len(line2)-1]) in prot_SP.keys() :
                    SP_signal = prot_SP[line2[1:len(line2)-1]]
            final_file = final_file + new_line2
    file.set_proteins_sequence(new_fasta_dict)
    cmd2 = "rm " + fasta_file
    os.system(cmd2)
    with open(fasta_file, "w") as new_file2 :
        new_file2.write(final_file)

def create_feature (file, Path_AlphaFold_Data, Path_Pickle_Feature) :
    """
    Launch command to generate features.

    Parameters:
    ----------
    file : object of class File_proteins
    Path_AlphaFold_Data : string
    Path_Pickle_Feature : string
    mmseq : boolean
    
    Returns:
    ----------
    """
    fasta_file = file.get_fasta_file()
    cmd = ["create_individual_features.py",
    f"--fasta_paths=./{fasta_file}",
    f"--data_dir={Path_AlphaFold_Data}",
    "--save_msa_files=True",
    f"--output_dir={Path_Pickle_Feature}",
    "--max_template_date=2024-05-02",
    "--skip_existing=True",
    "--use_mmseqs2=True"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
    for line in process.stdout:
       print(line, end="")
    process.stdout.close()
    process.wait()

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
    all_proteins = file.get_proteins()
    new_proteins = file.get_new_pickle() 
    shallow_MSA = str()
    result_dict = file.get_result_dict()
    for prot in new_proteins :
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
    for prot in all_proteins : #just write shallow_MSA.txt
        pre_feature_dict = pickle.load(open(f'{Path_Pickle_Feature}/{prot}.pkl','rb'))
        feature_dict = pre_feature_dict.feature_dict
        msa = feature_dict['msa']
        if len(msa) <= 100 :
            shallow_MSA += prot + " : " + str(len(msa)) + " sequences\n"
            result_dict[prot] = result_dict[prot] + "Shallow MSA"
    with open("shallow_MSA.txt", "w") as MSA_file :
        MSA_file.write(shallow_MSA)

def recover_prot_sequence(file, path_pkl) :
   """
   Take sequence from pickle files.

   Parameters:
   ----------
   file : object of File_proteins class
   path_pkl : string

   Returns:
   ----------
   """
   list_proteins = file.get_proteins()
   new_dict_sequence = file.get_proteins_sequence()
   for protein in list_proteins :
      with open(os.path.join(f'{path_pkl}/{protein}.pkl'), 'rb') as pkl_file :
         pickle_dict = pickle.load(pkl_file)
         new_dict_sequence[protein] = pickle_dict.sequence
   file.set_proteins_sequence(new_dict_sequence)

def filtered_signalP(file, SignalP) :
    file_signalp = file.get_fasta_file().replace(".fasta","_summary.signalp5")
    result_dict = file.get_result_dict()
    possible_prey = list()
    with open(file_signalp, "r") as SP_file :
        for line in SP_file :
            new_line = line.split("\t")
            if SignalP == "Yes" :
                if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                    possible_prey.append(new_line[0])
                if new_line[1] == "OTHER" and new_line[0][0] != "#" :
                    result_dict[new_line[0]] = "Don't have a signal peptide"
            if SignalP == "No" :
                if new_line[1] != "OTHER" and new_line[0][0] != "#" :
                    result_dict[new_line[0]] = "Have a signal peptide"
                if new_line[1] == "OTHER" and new_line[0][0] != "#" :
                    possible_prey.append(new_line[0])
    file.set_possible_prey(possible_prey)

def generate_bait_vs_prey (file, max_aa, informations_dict) :
    """
    Write one local script to use AlphaPullDown. This script should be written based on the maximum number of amino acids.

    Parameters:
    ----------
    file : object of class File_proteins
    max_aa : integer
    bait : list
    Returns:
    ----------
    """
    bait_vs_prey_script = str()
    OOM_int = str()
    bait = [prot for prot in informations_dict["Interact_with"].split(",")] #take all baits
    prey = file.get_possible_prey()
    prey = [protein for protein in prey if protein not in bait] #remove baits from prey list
    lenght_prot = file.get_lenght_prot()
    for index_protein in range(len(bait)) :
        lenght = lenght_prot[bait[index_protein]]
        for index2_protein in range(len(prey)) :
            int_lenght = lenght + lenght_prot[prey[index2_protein]]
            if int_lenght >= max_aa :
                OOM_int = OOM_int + bait[index_protein] + ";" +  prey[index2_protein]+ "\n"
            elif os.path.exists(f"./result_all_vs_all/{bait[index_protein]}_and_{prey[index2_protein]}/ranked_0.pdb") == False and os.path.exists(f"./result_all_vs_all/{prey[index2_protein]}_and_{bait[index_protein]}/ranked_0.pdb") == False: #make interaction if doesn't exist and is not too long
                bait_vs_prey_script = bait_vs_prey_script + bait[index_protein] + ";" +  prey[index2_protein]+ "\n"
            else :
                pass
    with open("bait_vs_prey.txt", "w") as all_file:
       all_file.write(bait_vs_prey_script)
    with open("OOM_int.txt", "w") as OOM_file :
       OOM_file.write(OOM_int)