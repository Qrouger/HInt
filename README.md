# <img src="https://github.com/user-attachments/assets/f4701588-b624-4afa-aa8f-9a3352a6572c" alt="HInt logo" width="200"/><br>



HInt allows you to find homologous proteins with significant differences in sequence and structure.
It allows you to find homologous proteins with similar interactions.
# Instalations
```bash
conda create -n HInt -c omnia -c bioconda -c conda-forge python==3.11 openmm==8.0 pdbfixer==1.9 kalign2 hhsuite hmmer
conda activate HInt
pip install HInt torchdata==0.9.0 pandas pydantic packaging opt_einsum torch-geometric matplotlib && \
pip install -U "jax[cuda12]"==0.5.3 && \
pip install  dgl -f https://data.dgl.ai/wheels/torch-2.1/cu121/repo.html
```
# MSA generations
MSA generation is mandatory for using this tool.

## Generated MSA with colabfold (recommended for less than 100 sequences)
Directly add to the pipeline.

## Generated MSA with mmseqs2 GPU (recommended for more than 100 proteins)

### Download mmseqs2 and localcolabfold
```bash
wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz <br>
tar xvzf mmseqs-linux-gpu.tar.gz <br>
export PATH=$(pwd)/mmseqs/bin/:$PATH <br>
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh <br>
bash install_colabbatch_linux.sh <br>
nano ~/.bashrc <br>
export PATH=$(pwd)/localcolabfold/colabfold-conda/bin:$PATH <br>
source ~/.bashrc
activate HInt
```

### Download mmseqs2 database
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

## HInt for eucaryote
## HInt for procaryote
```bash
git clone https://github.com/SNU-CSSB/RF2-Lite.git
cd RF2-Lite/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../networks/
wget http://files.ipd.uw.edu/pub/pathogens/weights.tar.gz
tar xfz weights.tar.gz
```


# Informations
Signal_peptide :#No, Yes or None.
Homo-oligomer :#Oligomerization (integer : 1-20).
Interact_with :#Proteins names expected to be interacting with (Uniprot name or fasta name), separate by ",".
Organism : #Organism of interest: arch, gram+, gram-, or euk for SignalP.
DeepLoc : #Cellular localisation of your research, only for euk. Possible localisation : Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgo apparatus, Peroxisome. It's possible to select multiple localisation with coma.
