# <img src="https://github.com/user-attachments/assets/f4701588-b624-4afa-aa8f-9a3352a6572c" alt="HInt logo" width="200"/><br>

HInt is an optimized and scalable pipeline designed for high-throughput identification of homologous proteins that retain conserved functional interactions despite substantial sequence and structural divergence. By combining efficient MSA reuse, parallelized structure prediction, and automated interaction scoring, HInt significantly accelerates large-scale interaction screening while maintaining high predictive accuracy. This enables the systematic discovery of conserved interaction networks that remain undetectable through sequence or structural similarity alone.

# 1.Instalations
```bash
conda create -n HInt -c conda-forge -c bioconda python==3.11 pdbfixer==1.9 mafft kalign2 hhsuite hmmer mmseqs2 git
conda activate HInt
pip install alphapulldown==2.1.4 nvidia-ml-py torch==2.4.0 ihm scipy==1.16.0 setuptools==80.9.0
pip install --no-warn-conflicts \ "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
pip install -U "jax[cuda12]"==0.5.3 numpy==1.26.4
```

## A. Download the GPU-indexed MMseqs2 database (2 hours, 1.5T)
To accelerate MSA generation, it is strongly recommended to store the databases on NVMe or SSD drives rather than on HDD storage.<br>
```bash
wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh
chmod +x setup_databases.sh 
GPU=1 ./setup_databases.sh ./mmseq_database
```

test : colabfold_search R388.fasta /data/colab_fold_data . --gpu 1 --db-load-mode 2 <br>

## B. Download AlphaFold database (3 hours, 2.7T)
```bash
sudo apt install aria2
git clone https://github.com/KosinskiLab/alphafold.git
cd alphafold
scripts/download_all_data.sh /<Database Directory> > download.log 2> download_all.log
```

## C.1. Install deeplocpro (Prokaryote)
```bash
git clone https://github.com/Jaimomar99/deeplocpro
cd deeplocpro
pip install .
```
## C.2. Install deeploc (Eukaryote)

Download deeploc2 package here : https://services.healthtech.dtu.dk/services/DeepLoc-2.0/
```bash
cd  deeploc2_package
pip install .
```

## D. Install SignalP5

Download SignalP5 here : [https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=5.0&packageversion=5.0b&platform=Darwin)<br>
```bash
tar -xvzf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/
sudo cp bin/signalp /usr/local/bin
sudo cp -r lib/* /usr/local/lib
```
## E. Install ccp4
Download ccp4 package here : https://www.ccp4.ac.uk/download/#os=linux
```bash
tar xvzf ccp4-9-setup.tar.gz
./ccp4-9-setup
```

<br>

# 2.Input parameters
## Setup HInt.txt

<br>

### The First part of HInt.txt file contains all a priori information about the query protein.

**Signal_peptide** : Indicates whether the protein has a signal peptide (Options : Yes,No or None). <br>

**DeepLoc** : Cellular localisation(s) of the protein. Multiple localizations can be specified, separated by commas. <br>

- Eukaryotes : Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgo apparatus, Peroxisome.
- Prokaryotes : Cell wall & surface, Extracellular, Cytoplasmic, Cytoplasmic Membrane, Outer Membrane, Periplasmic.

**Max_protein_lenght** : Maximum lenght of the protein you search (integer). <br>

**Min_protein_lenght** : Minimum lenght of the protein you search (integer). <br>

**Homo-oligomer** : Known homo-oligomerization state of the protein (integer : 1 to 20). <br>

**Interact_with** : Names of proteins expected to interact with the query protein (UniprotID or protein fasta name).
<details>
<summary>Advanced bait uses and examples </summary>

One bait :
```
Interact_with : UniprotID1
```

Region of a bait :
```
Interact_with : UniprotID1(20-200)
```

Multiple baits : # First has to be the principal. For now you can put a maximum of 3 differents bait
```
Interact_with : UniprotID1, UniprotID2 
```

Multimer bait : # Create a unique bait with multiple proteins
```
Interact_with : [Uniprot1, Uniprot2]
```

And you can mixed up all of theses examples ! <br>
/!\ HInt don't support multiple regions for baits proteins
</details>

**Organism** : Organism of interest for SignalP5 and DeepLoc (arch, gram+, gram-, or euk). <br>

>[!TIP]
>If you don’t know the information or want to skip it, you can leave this field blank.

<br>

### Second part of HInt.txt are paths.

**Path_AlphaFold_Data** : Path of AlphaFold databse (string).

**Path_ccp4** : Path of CCP4 package (string). Default set on /opt/xtal/ccp4-9.

**Path_MMseqs2_Data** : Path of GPU-indexed MMseqs2 database (string).
>[!NOTE]
>This Path is not mandatory. If not set also MMseqs2-GPU will no be used.

**Path_Uniprot_ID** : Path to the protein sequence file (string).

**Path_Pickle_Feature** : Path where MSA files will be saved (string).

<br>

## Setup protein file
The protein file must contain all UniProt IDs or all sequences in FASTA format for both preys and baits. <br>
This can be protein ncbi fasta file, classic fasta file, uniprotID's or a combination of all. <br>
>Protein file exemple

>[!TIP]
>The use of UniprotIDs is recommended for pipeline speed.

<br>

## Run HInt
You need to be in the directory with HInt.txt file.

```bash
HInt --cpu <Integer> --gpu <Integer(s)> --multi_job_per_gpu <Boolean>
```
--cpu : Number of CPUs available for computation. Enables CPU parallelization. By default, set to half of the available CPUs. <br>
--gpu : Index(es) of GPU(s) you want to uses. Declare multiple GPU allows GPU parallelisation. By default set on GPU 0. <br>
--multi_job_per_gpu : Allows multiple jobs to run on a single GPU, reducing time of modelisation. By default set on True. <br>

<br>

## Folder structure

# 3.Example
