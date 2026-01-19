# <img src="https://github.com/user-attachments/assets/f4701588-b624-4afa-aa8f-9a3352a6572c" alt="HInt logo" width="200"/><br>



HInt enables the identification of homologous proteins that may exhibit substantial sequence and structural divergence, while preserving similar functional interactions. This allows researchers to uncover conserved interaction networks that are not apparent from sequence or structural similarity alone.
# 1.Instalations
```bash
conda create -n HInt -c conda-forge -c bioconda python==3.11 pdbfixer==1.9 mafft kalign2 hhsuite hmmer mmseqs2 git
conda activate HInt
pip install --no-warn-conflicts \ "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
pip install alphapulldown==2.1.4 nvidia-ml-py torch==2.4.0 numpy==1.26.4 ihm scipy==1.16.0
pip install -U "jax[cuda12]"==0.5.3



```
### A. Download MMseqs2 database GPU indexed (2 hours, 1.5T)
To speed up MSA generation, it is recommended to put the databases on nvme or ssd disks.<br>
```bash
wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh
chmod +x setup_databases.sh 
GPU=1 ./setup_databases.sh ./mmseq_database
```
mmseqs databases UniRef90 ./UniRef90 tmp <br>
mmseqs createdb examples/DB.fasta targetDB <br>
mmseqs makepaddedseqdb targetDB targetDB_gpu <br>
mmseqs rmdb targetDB <br>
or GPU=1 ./setup_databases.sh /path/to/db_folder <br>

test : colabfold_search R388.fasta /data/colab_fold_data . --gpu 1 --db-load-mode 2 <br>

### B. Download AlphaFold database
```bash
sudo apt install aria2
git clone https://github.com/deepmind/alphafold.git
cd alphafold
scripts/download_all_data.sh /<Directory> > download.log 2> download_all.log
```

### C.1 Install deeplocpro (Prokaryote)

git clone https://github.com/Jaimomar99/deeplocpro <br>
cd deeplocpro <br>
pip install .

### C.2 Install deeploc (Eukaryote)

https://services.healthtech.dtu.dk/services/DeepLoc-2.0/

### D. SignalP5 installation

[https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=5.0&packageversion=5.0b&platform=Darwin)<br>
```bash
tar -xvzf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/
cp bin/signalp /usr/local/bin
sudo cp -r lib/* /usr/local/lib
```

# 2.Example 
## Folder structure
## Setup HInt.txt
## Setup protein file
The bait protein need to be in this file <br>
The use of UniprotIDs is recommended for pipeline speed
## Run HInt

# 3.All parameters
**Signal_peptide** : No, Yes or None. <br>

**Homo-oligomer** : Oligomerization (integer : 1-20). <br>

**Interact_with** : Proteins names expected to be interacting with (UniprotID or protein fasta name)
<details>
<summary>Advanced uses and examples </summary>

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

**Organism** : Organism of interest: arch, gram+, gram-, or euk for SignalP. <br>

**DeepLoc** : Cellular localisation of your research. It's possible to select multiple localisation with coma. <br>

- Eucaryote: Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgo apparatus, Peroxisome.
- Procaryote: Cell wall & surface, Extracellular, Cytoplasmic, Cytoplasmic Membrane, Outer Membrane, Periplasmic.


# 4.Results
