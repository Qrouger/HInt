# <img width="240" height="240" alt="HInt_logo" src="https://github.com/user-attachments/assets/0d85a047-f02c-4819-9bc2-b51ca8bf0aba" />
# HInt
HInt accelerates AlphaFold by optimizing computations and parallelizing structure prediction. It is a scalable pipeline for high-throughput identification of homologous proteins and interologues—proteins that maintain functional interactions.  HInt enables the discovery of conserved interaction networks that may remain undetected using sequence or structural similarity alone.

# 1. Installation

## 1.1. HInt
Conda or Mamba must be installed on your system and must not be activated before running the installation.
```bash
wget https://raw.githubusercontent.com/Qrouger/HInt/main/Install_HInt.sh
bash Install_HInt.sh
```

## 1.2. DeepLoc2 (Eukaryote)

Download the DeepLoc 2.0 package from: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/

```bash
conda activate HInt
cd  deeploc2_package
pip install . torch==2.6.0
pip install triton==3.1.0 nvidia-cudnn-cu12==9.26.0.51
```

## 1.3. SignalP5

Download the SignalP5 package from: [https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=5.0&packageversion=5.0b&platform=Linux)<br>

```bash
tar -xvzf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/
sudo cp bin/signalp /usr/local/bin
sudo cp -r lib/* /usr/local/lib
```

## 1.4. CCP4

Download the ccp4 package from: https://www.ccp4.ac.uk/download/#os=linux

```bash
tar xvzf ccp4-9-setup.tar.gz
./ccp4-9-setup
```
<br>

# 2. Download databases (~1.5 TB)
For optimal performance, store MMseqs2 database and MSAs on SSD storage rather than HDDs.<br>

```bash
conda activate HInt
wget https://raw.githubusercontent.com/Qrouger/HInt/main/setup_databases.sh
chmod +x setup_databases.sh 
GPU=1 ./setup_databases.sh <DB_DIR>
```

# 3. Input parameters
## 3.1. Setup HInt.txt <br>
You need to download or copy HInt.txt file example. <br>
### *A priori* information

- **Signal_peptide** : Filter proteins based on the presence of a predicted signal peptide (Options : Yes, No or None).<br>

- **DeepLoc** : Cellular localisation(s) of the protein. Multiple localizations can be specified, separated by commas. All proteins predicted to localize to one of these compartments will be used.<br>
  - Eukaryotes : Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgi apparatus, Peroxisome.
  - Prokaryotes : Cell wall & surface, Extracellular, Cytoplasmic, Cytoplasmic Membrane, Outer Membrane, Periplasmic.

- **Max_protein_length** : Maximum length of the protein you search (integer). <br>

- **Min_protein_length** : Minimum length of the protein you search (integer), default set on 20 aa. <br>

- **AlphaFold** : AlphaFold version (Options : 2 or 3). <br>

- **Homo-oligomer** : Known homo-oligomerization state of the target protein (integer : 1 to 20), default set to 1 (monomer). <br>

- **Interact_with** : Names of proteins expected to interact with the query protein (Uniprot ID or protein fasta name).

<details>
<summary>Advanced Interact_with usage and examples </summary>

One bait :
```
Interact_with : UniprotID1
```

Region of a bait :
```
Interact_with : UniprotID1(20-200)
```

Multiple baits : # The first protein must correspond to the primary bait. For now, you can specify a maximum of three different baits
```
Interact_with : UniprotID1, UniprotID2 
```

Multimer bait : # Create a unique bait with multiple proteins
```
Interact_with : [Uniprot1, Uniprot2]
```

And you can mixed up all of theses examples ! <br>
 <br>
⚠️ **Warning** <br>
HInt currently does not support multiple regions for bait proteins.
</details>

- **Organism** : Organism of interest for SignalP5 and DeepLoc (arch, gram+, gram-, euk or None). Enables signal peptide prediction and cleavage. <br>

<br>

>[!TIP]
>If you don’t know the information or want to skip it, you can leave this field blank.

<br>

### Paths

- **Path_Database** : Path to the HInt database (string).

- **Path_ccp4** : Path of CCP4 package (string). Default set on /opt/xtal/ccp4-9.

>[!NOTE]
>This path is not mandatory. If not set also MMseqs2-GPU will not be used.

- **Path_Uniprot_ID** : Path to the protein sequence file (string).

- **Path_Pickle_Feature** : Path where MSA and pkl files will be saved (string).

<br>

## 3.2. Setup protein file
The protein file must contain all UniProt IDs and/or all sequences in FASTA format for both preys and baits. <br>
This can be an NCBI protein FASTA file, a standard FASTA file, UniProt identifiers, or a combination of these formats. <br>

>[!TIP]
>The use of UniprotIDs is recommended for pipeline speed.

<details>
<summary>Protein file examples</summary>
<br>
  
**Example 1**  <br>
```
A0ABD7FQG0,P18004,P15069,P33790
```

**Example 2**  <br>
```
>A0ABD7FQG0
MSGDENKLKKYRFPETLTNQSRWFGLPLDELIPAAICIGWGITTSKYLFGIGAAVLVYFGIKKLKKGRGSSWLRDLIYWYMPTALLRGIFHNVPDSCFRQWIK
>P18004
MNNPLEAVTQAVNSLVTALKLPDESAKANEVLGEMSFPQFSRLLPYRDYNQESGLFMNDTTMGFMLEAIPINGANESIVEALDHMLRTKLPRGIPLCIHLMSSQLVGDRIEYGLREFSWSGEQAERFNAITRAYYMKAAATQFPLPEGMNLPLTLRHYRVFISYCSPSKKKSRADILEMENLVKIIRASLQGASITTQTVDAQAFIDIVGEMINHNPDSLYPKRRQLDPYSDLNYQCVEDSFDLKVRADYLTLGLRENGRNSTARILNFHLARNPEIAFLWNMADNYSNLLNPELSISCPFILTLTLVVEDQVKTHSEANLKYMDLEKKSKTSYAKWFPSVEKEAKEWGELRQRLGSGQSSVVSYFLNITAFCKDNNETALEVEQDILNSFRKNGFELISPRFNHMRNFLTCLPFMAGKGLFKQLKEAGVVQRAESFNVANLMPLVADNPLTPAGLLAPTYRNQLAFIDIFFRGMNNTNYNMAVCGTSGAGKTGLIQPLIRSVLDSGGFAVVFDMGDGYKSLCENMGGVYLDGETLRFNPFANITDIDQSAERVRDQLSVMASPNGNLDEVHEGLLLQAVRASWLAKENRARIDDVVDFLKNASDSEQYAESPTIRSRLDEMIVLLDQYTANGTYGQYFNSDEPSLRDDAKMVVLELGGLEDRPSLLVAVMFSLIIYIENRMYRTPRNLKKLNVIDEGWRLLDFKNHKVGEFIEKGYRTARRHTGAYITITQNIVDFDSDKASSAARAAWGNSSYKIILKQSAKEFAKYNQLYPDQFLPLQRDMIGKFGAAKDQWFSSFLLQVENHSSWHRLFVDPLSRAMYSSDGPDFEFVQQKRKEGLSIHEAVWQLAWKKSGPEMASLEAWLEEHEKYRSVA
>P15069
MMPRIKPLLVLCAALLTVTPAASADVNSDMNQFFNKLGFASNTTQPGVWQGQAAGYAYGGSLYARTQVKNVQLISMTLPDINAGCGGIDAYLGSFSFINGEQLQRFVKQIMSNAAGYFFDLALQTTVPEIKTAKDFLQKMASDINSMNLSSCQAAQGIIGGLFPRTQVSQQKVCQDIAGESNIFADWAASRQGCTVGGKSDSVRDKASDKDKERVTKNINIMWNALSKNRMFDGNKELKEFVMTLTGSLVFGPNGEITPLSARTTDRSIIRAMMEGGTAKISHCNDSDKCLKVVADTPVTISRDNALKSQITKLLASIQNKAVSDTPLDDKEKGFISSTTIPVFKYLVDPQMLGVSNSMIYQLTDYIGYDILLQYIQELIQQARAMVATGNYDEAVIGHINDNMNDATRQIAAFQSQVQVQQDALLVVDRQMSYMRQQLSARMLSRYQNNYHFGGSTL
>P33790
MNEVYVIAGGEWLRNNLNAIAAFMGTWTWDSIEKIALTLSVLAVAVMWVQRHNVMDLLGWVAVFVLISLLVNVRTSVQIIDNSDLVKVHRVDNVPVGLAMPLSLTTRIGHAMVASYEMIFTQPDSVTYSKTGMLFGANLIVKSTDFLSRNPEIINLFQDYVQNCVLGDIYLNHKYTLEDLMASADPYTLIFSRPSPLRGVYDNNNNFITCKDASVTLKDRLNLDTKTGGKTWHYYVQQIFGGRPDPDLLFRQLVSDSYSYFYGSSQSASQIMRQNVTMNALKEGITSNAARNGDTASLVSLATTSSMEKQRLAHVSIGHVTMRNLPMVQTILTGIAIGIFPLLILAAVFNKLTLSVLKGYVFALMWLQTWPLLYAILNSAMTFYAKQNGAPVVLSELSQIQLKYSNLASTAGYLSAMIPPLSWMMVKGLGAGFSSVYSHFASSSISPTASAAGSVVDGNYSYGNMQTENVNGFSWSTNSTTSFGQMMYQTGSGATATQTRDGNMVMDASGAMSRLPVGINATRQIAAAQQEMAREASNRAESALHGFSSSIASAWNTLSQFGSNRGSSDSVTGGADSTMSAQDSMMASRMRSAVESYAKAHNISNEQATRELASRSTNASLGLYGDAYAKGHLGISVLGNGGGVGLQAGAKASIDGSDLDSHEASSGSRASHDARHDIDARATQDFKEASDYFTSRKVSESGSHTDNNADSRVDQLSAALNSAKQSYDQYTTNMTRSHEYAEMASRTESMSGQMSEDLSQQFAQYVMKNAPQDVEAILTNTSSPEIAERRRAMAWSFVQEQVQPGVDNTWRESRRDIGKGMESVPSGGGSQDIIADHQGHQAIIEQRTQDSNIRNDVKHQVDNMVTEYRGNIGDTQNSIRGEENIVKGQYSELQNHHKTEALTQNNKYNEEKLAQERIPGADSPKELLEKAKSYQHKE
```

**Example 3**  <br>
```
A0ABD7FQG0,P18004
>P15069
MMPRIKPLLVLCAALLTVTPAASADVNSDMNQFFNKLGFASNTTQPGVWQGQAAGYAYGGSLYARTQVKNVQLISMTLPDINAGCGGIDAYLGSFSFINGEQLQRFVKQIMSNAAGYFFDLALQTTVPEIKTAKDFLQKMASDINSMNLSSCQAAQGIIGGLFPRTQVSQQKVCQDIAGESNIFADWAASRQGCTVGGKSDSVRDKASDKDKERVTKNINIMWNALSKNRMFDGNKELKEFVMTLTGSLVFGPNGEITPLSARTTDRSIIRAMMEGGTAKISHCNDSDKCLKVVADTPVTISRDNALKSQITKLLASIQNKAVSDTPLDDKEKGFISSTTIPVFKYLVDPQMLGVSNSMIYQLTDYIGYDILLQYIQELIQQARAMVATGNYDEAVIGHINDNMNDATRQIAAFQSQVQVQQDALLVVDRQMSYMRQQLSARMLSRYQNNYHFGGSTL
>P33790
MNEVYVIAGGEWLRNNLNAIAAFMGTWTWDSIEKIALTLSVLAVAVMWVQRHNVMDLLGWVAVFVLISLLVNVRTSVQIIDNSDLVKVHRVDNVPVGLAMPLSLTTRIGHAMVASYEMIFTQPDSVTYSKTGMLFGANLIVKSTDFLSRNPEIINLFQDYVQNCVLGDIYLNHKYTLEDLMASADPYTLIFSRPSPLRGVYDNNNNFITCKDASVTLKDRLNLDTKTGGKTWHYYVQQIFGGRPDPDLLFRQLVSDSYSYFYGSSQSASQIMRQNVTMNALKEGITSNAARNGDTASLVSLATTSSMEKQRLAHVSIGHVTMRNLPMVQTILTGIAIGIFPLLILAAVFNKLTLSVLKGYVFALMWLQTWPLLYAILNSAMTFYAKQNGAPVVLSELSQIQLKYSNLASTAGYLSAMIPPLSWMMVKGLGAGFSSVYSHFASSSISPTASAAGSVVDGNYSYGNMQTENVNGFSWSTNSTTSFGQMMYQTGSGATATQTRDGNMVMDASGAMSRLPVGINATRQIAAAQQEMAREASNRAESALHGFSSSIASAWNTLSQFGSNRGSSDSVTGGADSTMSAQDSMMASRMRSAVESYAKAHNISNEQATRELASRSTNASLGLYGDAYAKGHLGISVLGNGGGVGLQAGAKASIDGSDLDSHEASSGSRASHDARHDIDARATQDFKEASDYFTSRKVSESGSHTDNNADSRVDQLSAALNSAKQSYDQYTTNMTRSHEYAEMASRTESMSGQMSEDLSQQFAQYVMKNAPQDVEAILTNTSSPEIAERRRAMAWSFVQEQVQPGVDNTWRESRRDIGKGMESVPSGGGSQDIIADHQGHQAIIEQRTQDSNIRNDVKHQVDNMVTEYRGNIGDTQNSIRGEENIVKGQYSELQNHHKTEALTQNNKYNEEKLAQERIPGADSPKELLEKAKSYQHKE
```

<br>
</details>

<br>

# 4. Run HInt
You need to be in the directory containing the HInt.txt file and activate the HInt Conda environment.

```bash
HInt --cpu <Integer> --gpu <Integer(s)> --multi_job_per_gpu <Boolean> --send_to_database <All,Negative,None>
```
<details>
<summary>Flag descriptions </summary>

```yaml
# Number of CPUs available for computation. Enables CPU parallelization. By default, set to half of the available CPUs.
  --cpu : Integer

# Index(es) of GPU(s) you want to use. Declare multiple GPU allows GPU parallelisation. By default set on GPU 0. 
  --gpu : Integer(s)
  
# Allows multiple jobs to run on a single GPU, reducing time of modeling. By default set on True. 
  --multi_job_per_gpu : Boolean

# Allows to sending interaction generated in a public database (https://qsproteome.org/). All PPI, only negative PPI or none. By default set on None.
  --send_to_database : All, Negative or None.
```
</details>

## 4.1. Run HInt on an HPC cluster (Slurm)

`HInt-slurm` splits the pipeline into Slurm jobs chained with dependencies. Run it on a login node, from the directory containing HInt.txt, with the HInt environment activated. The screen directory, the databases and `Path_Pickle_Feature` must be on a shared file system.

```bash
HInt-slurm submit --partition <partition> --account <account> --gpu_vram <GiB> --predict_tasks 20
```

| Stage | Slurm job | Work |
|---|---|---|
| fetch | 1 job, only with `--offline_compute` | UniProt sequences, on the partition with internet access |
| prepare | 1 job (GPU only if `Organism` is set, for DeepLoc) | Sequences, DeepLoc, SignalP, filters, job list of the first round |
| afdb | 1 job, only with `--offline_compute` | MSAs of the AlphaFold database, on the partition with internet access |
| msa | 1 GPU job | MSAs (AlphaFold database + MMseqs2-GPU `colabfold_search` on all proteins at once) |
| features | job array (CPU) | AlphaFold feature pkl files |
| predict | job array (GPU), one per round | 3D models, packed on the GPUs of each task by predicted VRAM |
| score | 1 CPU job per round | Scoring, job list of the next round, final results after the last round |

A round is one bait of `Interact_with` (in the given order: the next bait is only modeled with the preys validated with the previous one, as with `HInt`) or the homo-oligomer modeling. The predict arrays of a round run in parallel on as many GPUs as the cluster gives.

<details>
<summary>HInt-slurm options </summary>

```yaml
# Cluster
  --cluster : generic or jean-zay (IDRIS profile, with --jz_project and --jz_gpu a100|h100|v100)
  --offline_compute, --internet_partition : compute nodes without internet access, partition used for the steps that need it

# Slurm
  --partition / --gpu_partition / --cpu_partition : partition of all jobs, GPU jobs, CPU jobs
  --account, --qos : Slurm account and QOS
  --gpu_type : GPU type in --gres (gpu:<type>:N), e.g. a100
  --sbatch_extra : extra #SBATCH line, repeatable, e.g. --sbatch_extra="--constraint=a100"
  --setup : shell line run in each job before HInt, repeatable, e.g. --setup "module load cuda/12.4"
  --conda_env : conda environment activated in each job (default: the current one), --no_conda to disable
  --max_parallel : maximum number of array tasks running at the same time
  --dry_run : write the scripts in slurm_scripts/ without submitting them

# Work distribution
  --features_tasks : array tasks for pkl files (default: 10)
  --predict_tasks : array tasks per modeling round (default: 10)
  --gpus_per_task : GPUs per predict task (default: 1)
  --<stage>_cpus, --<stage>_mem, --<stage>_time : resources of each stage (fetch, prepare, afdb, msa, features, predict, score)
  --start_at, --start_round : resume a run from a stage (fetch, prepare, afdb, msa, features, predict) and a round

# HInt
  --gpu_vram : VRAM (GiB) of the GPUs used for modeling. Recommended: the job lists are built before the predict jobs start
  --multi_job_per_gpu, --multi_scoring, --send_to_database : same as HInt
  --max_retries : retries of a crashed AlphaFold job (default: 1)
  --no_afdb : do not download MSAs from the AlphaFold database
  --local_search_threshold : use local colabfold_search from this number of MSAs, below use the ColabFold server (default: 1)
  --jax_cache_dir : shared JAX compilation cache (default: next to Path_Database)
```
</details>

```bash
HInt-slurm status                            # progress of each stage and round
scancel $(cat log_file/slurm/job_ids.txt)    # cancel the pipeline
```

Logs of each job are written in `log_file/slurm/` (Slurm outputs, one log per stage/task, AlphaFold logs per round, `round_<N>/failed_jobs.txt`). The results are the same files as with `HInt`. Every stage skips the work already done, so a run can be submitted again after a timeout or a crash.

### Compute nodes without internet access

Use `--offline_compute --internet_partition <partition>`: the UniProt sequences (`fetch` stage), the MSAs of the AlphaFold database (`afdb` stage), the pkl files (AF2 template search uses the ColabFold server) and the scoring (`--send_to_database`) run on that partition. Create the `pi_score` conda environment on the login node (`HInt-slurm submit` warns if it is missing). Build the local species database (`uniref30_2302_species.sqlite`, done by `setup_databases.sh`) so the pkl files do not query UniProt to pair MSAs. Each features task uses its own species cache (`Path_Pickle_Feature/uniprot_species_cache_task_<N>.sqlite`), SQLite cannot be shared between nodes on a network file system.

### IDRIS Jean Zay

```bash
HInt-slurm submit --cluster jean-zay --jz_project <project> --jz_gpu a100 --predict_tasks 20
```

The `jean-zay` profile applies the IDRIS rules:
- GPU jobs (DeepLoc, MSA, modeling): `--constraint=a100|h100|v100-32g`, `--account=<project>@<gpu>`, `--gres=gpu:N`, the CPUs matching the GPU share of the node (8 per A100, 24 per H100, 10 per V100), `--hint=nomultithread`, no `--mem` (memory is proportional to the CPUs), QoS `qos_gpu_<gpu>-t3` (20h) or `-t4` above 20h, `module load arch/h100` on H100.
- Compute nodes have no internet access: sequences download, AlphaFold database MSAs, pkl files and scoring run on the `prepost` partition (internet access, hours not deducted, 20h max) with `--account=<project>@cpu` (change it with `--jz_internet_account`).
- `--gpu_vram` is set to 80 (A100/H100) or 32 (V100).

Recommendations:
- Install the conda environments in `$WORK`, not in `$HOME` (3 GB quota): `conda config --add envs_dirs $WORK/conda/envs` and `conda config --add pkgs_dirs $WORK/conda/pkgs` before running `Install_HInt.sh` on a login node (internet access through the proxy).
- Put the screen directory, the databases and `Path_Pickle_Feature` in `$WORK` (quota of 5 TB and 500k inodes per project, check with `idr_quota_project`). `$SCRATCH` is purged after 30 days without access and `$STORE` is not reachable from compute nodes.
- Keep `--predict_time` at 20h or less and increase `--predict_tasks` instead. If a predict task reaches its time limit, the scoring stops the pipeline and prints the command to resume it (models already done are skipped).
- Check the current limits on the IDRIS documentation (http://www.idris.fr/jean-zay/), they change with the machine extensions.

<details>
<summary>Initial folder structure </summary>
  
```
HInt_screen/
  HInt.txt
  sequences.txt
```

</details>

<br>

# 5. Results

<details>
<summary>Final folder structure </summary>
  
```
HInt_screen/
  HInt.txt
  sequences.txt
  All_Final_result_HInt.csv
  result_PPI_int/
    P33790_and_P15069
    ...
  msa_feature/
    P33790.pkl
    P33790.a3m
    P33790_coverage.pdf
    P33790_feature_metadata_2026-04-03.json
    P15069.pkl
    P15069.a3m
    P15069_coverage.pdf
    P15069_feature_metadata_2026-04-03.json
    ...
  Interface_fig/
  log_file/
    HInt.log
    Summary_result_HInt.csv
    HInt_report.txt
```
</details>

## Global results description

The pipeline generates three main output files:

### All_Final_result_HInt.csv
Final curated results after filtering steps.

- **Name**: Protein identifier  
- **Localization**: Predicted subcellular localization  
- **Signal_peptide**: Signal peptide prediction (presence or not)  
- **Score**: Final HInt interaction score  

---

### Summary_result_HInt.csv
Summary of filtering decisions applied to each entry.

- **Name**: Protein identifier  
- **Reason_for_filtering**: Reason why the entry was filtered or flagged  

---

### HInt_report.txt
Comprehensive report of all processed entries and timings.

Includes:
- All entry information  
- Execution timings for each stage (prediction, scoring, filtering, etc.)

## PPI specific results

`<PPI>_rest_int.csv`  
Table of interface residues identified at the protein-protein interface.  
Includes interface residues identified using PAE and inter-chain distance criteria (<10 Å).

`<PPI>_ranked_0.pdb`  
Structural model of the predicted complex.  
Interface residues can be visualized by coloring the structure using the B-factor field.

# Standalone iQ-score Calculation

Compute **iQ-score** independently of the full HInt workflow..

[![GitHub](https://img.shields.io/badge/GitHub-iQ--score-black?style=for-the-badge&logo=github)](https://github.com/Qrouger/iQ-score)

<img width="480" height="480" alt="ROC_curve" src="https://github.com/user-attachments/assets/16dd98f7-b6ed-4a4c-b89f-06b79199af16" />

# HInt modelisation acceleration

Time required to predict the structures of a set of PPIs using four RTX 4500 Ada GPUs compared with the conventional pipeline.
<img width="916" height="132" alt="Screenshot from 2026-10-06 11-48-40" src="https://github.com/user-attachments/assets/ae5c4279-b95c-455c-9b08-546fa7947f8f" />

# Citations
**HInt**: interaction-based homology discovery through accelerated genome-scale AlphaFold screening
Quentin Rouger, Pierre Paillard, Manon Thomet, Emma Touquet, Gwenaël Rabut, Emmanuel Giudice, Damien F. Meyer, Kévin Macé
bioRxiv 2026.07.31.741991; doi: https://doi.org/10.64898/2026.07.31.741991
</div>
