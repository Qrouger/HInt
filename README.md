# <img width="240" height="240" alt="HInt_logo" src="https://github.com/user-attachments/assets/0d85a047-f02c-4819-9bc2-b51ca8bf0aba" />
# HInt
HInt accelerates AlphaFold by optimizing computations and parallelizing structure predictions. It is a scalable pipeline for high-throughput identification of homologous proteins and interologues—proteins that maintain functional interactions. This enables the systematic discovery of conserved interaction networks that remain undetectable through sequence or structural similarity alone.

# 1. Installation

## 1.1. HInt
Conda/Mamba must be installed on your system before running this software.
```bash
wget https://github.com/Qrouger/Install_HInt.sh
bash Install_HInt.sh
```
<details>
<summary>AlphaFold 3 (optional) </summary>

⚠️ **Warning** <br>
You need to have AlphaFold 3 model parameters in Path_AlphaFold_Data (https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)

</details>

## 1.2. DeepLoc2 (Eukaryote)

Download deeploc2 package here : https://services.healthtech.dtu.dk/services/DeepLoc-2.0/
```bash
conda activate HInt
cd  deeploc2_package
pip install . torch==2.6.0
pip install triton==3.1.0
```

## 1.3. SignalP5

Download SignalP5 here : [https://services.healthtech.dtu.dk/services/SignalP-5.0/9-Downloads.php](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=5.0&packageversion=5.0b&platform=Linux)<br>

```bash
tar -xvzf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/
sudo cp bin/signalp /usr/local/bin
sudo cp -r lib/* /usr/local/lib
```

## 1.4. CCP4

Download ccp4 package here : https://www.ccp4.ac.uk/download/#os=linux
```bash
tar xvzf ccp4-9-setup.tar.gz
./ccp4-9-setup
```
<br>

# 2. Download databases
## 2.1. Download the GPU-indexed MMseqs2 database (2 hours, 1.9T)
To accelerate MSA generation, it is strongly recommended to store the databases on NVMe or SSD drives rather than on HDD storage.<br>
```bash
wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh
chmod +x setup_databases.sh 
GPU=1 ./setup_databases.sh ./MMseqs2_GPU_database
```


## 2.2. Download AlphaFold3 database (633G)
```bash
git clone https://github.com/google-deepmind/alphafold3.git
cd alphafold3
./fetch_databases.sh <DB_DIR>
```
<br>

# 3. Input parameters
## 3.1. Setup HInt.txt <br>
You need to download or copy HInt.txt file example. <br>
### *A priori* informations

- **Signal_peptide** : Filter proteins based on the presence of a predicted signal peptide (Options : Yes,No or None).<br>

- **DeepLoc** : Cellular localisation(s) of the protein. Multiple localizations can be specified, separated by commas. All proteins predicted to be in one of these compartments will be used.<br>
  - Eukaryotes : Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgo apparatus, Peroxisome.
  - Prokaryotes : Cell wall & surface, Extracellular, Cytoplasmic, Cytoplasmic Membrane, Outer Membrane, Periplasmic.

- **Max_protein_lenght** : Maximum lenght of the protein you search (integer). <br>

- **Min_protein_lenght** : Minimum lenght of the protein you search (integer), default set on 20aa. <br>

- **AlphaFold** : AlphaFold version (Options : 2 or 3). <br>

- **Homo-oligomer** : Known homo-oligomerization state of the protein (integer : 1 to 20), default set on 1 (monomer). <br>

- **Interact_with** : Names of proteins expected to interact with the query protein (UniprotID or protein fasta name).

<details>
<summary>Advanced Interact_with uses and examples </summary>

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
 <br>
⚠️ **Warning** <br>
HInt don't support multiple regions for baits proteins
</details>

- **Organism** : Organism of interest for SignalP5 and DeepLoc (arch, gram+, gram-, euk or None). Enables signal peptide prediction and cleavage. <br>

<br>

>[!TIP]
>If you don’t know an information or want to skip it, you can leave this field blank.

<br>

### Paths

- **Path_AlphaFold_Data** : Path of AlphaFold database (string).

- **Path_ccp4** : Path of CCP4 package (string). Default set on /opt/xtal/ccp4-9.

- **Path_MMseqs2_Data** : Path of GPU-indexed MMseqs2 database (string).
>[!NOTE]
>This path is not mandatory. If not set also MMseqs2-GPU will not be used.

- **Path_Uniprot_ID** : Path to the protein sequence file (string).

- **Path_Pickle_Feature** : Path where MSA files will be saved (string).

<br>

## 3.2. Setup protein file
The protein file must contain all UniProt IDs or all sequences in FASTA format for both preys and baits. <br>
This can be protein ncbi fasta file, classic fasta file, uniprotID's or a combination of all. <br>

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
You need to be in the directory with HInt.txt file.

```bash
HInt --cpu <Integer> --gpu <Integer(s)> --multi_job_per_gpu <Boolean>
```
<details>
<summary>Flags description </summary>

```yaml
# Number of CPUs available for computation. Enables CPU parallelization. By default, set to half of the available CPUs.
  --cpu : Integer

# Index(es) of GPU(s) you want to uses. Declare multiple GPU allows GPU parallelisation. By default set on GPU 0. 
  --gpu : Integer(s)
  
# Allows multiple jobs to run on a single GPU, reducing time of modelisation. By default set on True. 
  --multi_job_per_gpu : Boolean
```
</details>

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
Includes residues selected based on PAE and distance criteria (< 10 Å).

`<PPI>_ranked_0.pdb`  
Structural model of the predicted complex.  
Interface residues can be visualized by coloring the structure using the B-factor field.

# Standalone iQ-score Calculation

Compute **iQ-score** independently from the full HInt workflow.

[![GitHub](https://img.shields.io/badge/GitHub-iQ--score-black?style=for-the-badge&logo=github)](https://github.com/Qrouger/iQ-score)

</div>
