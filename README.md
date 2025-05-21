# <img src="https://github.com/user-attachments/assets/f4701588-b624-4afa-aa8f-9a3352a6572c" alt="HInt logo" width="200"/><br>


HInt allows you to find homologous proteins with significant differences in sequence and structure.
It allows you to find homologous proteins with similar interactions.

##MSA generations
MSA generation is mandatory for using this tool.

#MSA already generated

#Generated MSA with colabfold (recommended for less than 100 sequences)

#Generated MSA with mmseq2 GPU (recommended for more than 100 proteins)


##Instalations
Colabfoldlocal
https://github.com/YoshitakaMo/localcolabfold

##Informations
Signal_peptide :#No, Yes or None
Homo-oligomer :#Oligomerization (integer : 1-20)
Interact_with :#Proteins names expected to be interacting with (Uniprot name or fasta name), separate by ","
Organism : #Organism of interest: arch, gram+, gram-, or euk for SignalP
