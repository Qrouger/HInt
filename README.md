# <img src="https://github.com/user-attachments/assets/f4701588-b624-4afa-aa8f-9a3352a6572c" alt="HInt logo" width="200"/><br>



HInt enables the identification of homologous proteins that may exhibit substantial sequence and structural divergence, while preserving similar functional interactions. This allows researchers to uncover conserved interaction networks that are not apparent from sequence or structural similarity alone.
# 1.Instalations
```bash
conda create -n HInt -c conda-forge -c bioconda python==3.11 pdbfixer==1.9 mafft kalign2 hhsuite hmmer mmseqs2 git
conda activate HInt
pip install alphapulldown==2.1.4 nvidia-ml-py torch==2.4.0 ihm scipy==1.16.0 setuptools==80.9.0
pip install --no-warn-conflicts \ "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
pip install -U "jax[cuda12]"==0.5.3 numpy==1.26.4
```

## A. Download MMseqs2 database GPU indexed (2 hours, 1.5T)
To speed up MSA generation, it is recommended to put the databases on nvme or ssd disks.<br>
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

# 2.Input parameters
## Setup HInt.txt
**Signal_peptide** : No, Yes or None. <br>

**Homo-oligomer** : Oligomerization (integer : 1-20). <br>

**Interact_with** : Proteins names expected to be interacting with (UniprotID or protein fasta name)
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
Interact_with : [Uniprot1,Uniprot2]
```

And you can mixed up all of theses examples ! <br>
/!\ HInt don't support multiple regions for baits proteins
</details>

**Organism** : Organism of interest: arch, gram+, gram-, or euk for SignalP5 and DeepLoc. <br>

**DeepLoc** : Cellular localisation of your research. It's possible to select multiple localisation with coma. <br>

- Eukaryote : Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgo apparatus, Peroxisome.
- Prokaryote : Cell wall & surface, Extracellular, Cytoplasmic, Cytoplasmic Membrane, Outer Membrane, Periplasmic.

## Setup protein file
Protein file need to contains all UniprotID or all sequences in fasta format of preys and baits. <br>
This can be protein ncbi fasta file, classic fasta file, uniprotID's or a combination of all. <br>
>Protein file exemple

>[!TIP]
>The use of UniprotIDs is recommended for pipeline speed.

## Run HInt
You need to be in the directory with HInt.txt file.

```bash
HInt --cpu integer --gpu integer --multi_job_per_gpu Boolean
```
--cpu : Number of CPUs available for computation, allow CPU parallelisation. By default set on half of available CPU. <br>
--gpu : Index of GPUs you want to uses. Declare multiple GPU allows GPU parallelisation. By default set on GPU 0. <br>
--multi_job_per_gpu : Allows multiple launched of jobs in on GPU, reduced time modelisation. By default set on True. <br>
## Folder structure

# 4.Results
