# <img src="https://github.com/user-attachments/assets/f4701588-b624-4afa-aa8f-9a3352a6572c" alt="HInt logo" width="200"/><br>



HInt allows you to find homologous proteins with significant differences in sequence and structure.
It allows you to find homologous proteins with similar interactions.
# 1.Instalations
```bash
conda create -n HInt -c conda-forge python==3.11 pdbfixer==1.9 mafft
pip install alphapulldown nvidia-ml-py torch==2.4.0 numpy==1.26.4 ihm scipy==1.16.0
pip install -U "jax[cuda12]"==0.5.3


RF2_PPI :

pip install torch==1.12.1+cu113 -f https://download.pytorch.org/whl/torch_stable.html

```
# MSA generations
MSA generation is mandatory for using this tool.

## Generated MSA with colabfold (recommended for less than 100 sequences)
Directly add to the pipeline.

## Generated MSA with MMseqs2 GPU (recommended for more than 100 proteins)

### Download MMseqs2 and localcolabfold
```bash
wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz
tar xvzf mmseqs-linux-gpu.tar.gz
export PATH=$(pwd)/mmseqs/bin/:$PATH
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh
export PATH=$(pwd)/localcolabfold/colabfold-conda/bin:$PATH
source ~/.bashrc
activate HInt
```

### Download MMseqs2 database GPU indexed
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

colabfold_search R388.fasta /data/colab_fold_data . --gpu 1 --db-load-mode 2 <br>

### Install deeplocpro

git clone https://github.com/Jaimomar99/deeplocpro <br>
cd deeplocpro <br>
pip install .
# 2.Example 
## Folder structure
## Setup HInt.txt
## Run HInt

# 3.All parameters
**Signal_peptide** : No, Yes or None. <br>

**Homo-oligomer** : Oligomerization (integer : 1-20). <br>

**Interact_with** : Proteins names expected to be interacting with (Uniprot name or fasta name)
<details>
<summary>Advanced uses and examples </summary>

One bait :
```
Interact_with : UniprotID1
```

Part of a bait :
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
- Procaryote:

DeepLocPro : pro https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/


# 4.Results
