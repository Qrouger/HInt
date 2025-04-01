import requests

# Définir l'URL de l'API UniProt avec le taxon 2492


def find_plasmid(taxonomy_id):
   url = 'https://rest.uniprot.org/uniprotkb/search?query=taxonomy_id:{taxonomy_id}&format=json&size=500'
   response = requests.get(url)
   if response.status_code == 200:
      prot_in_plasmid = []
      data = response.json()
      print(len(data['results']))
      for entry in data['results']:
         prot_in_plasmid.append(entry['primaryAccession'])
#        print(f"ID: {entry['primaryAccession']}, Nom: {entry['proteinDescription']['submissionNames']}")
   else:
      print(f"Erreur lors de la récupération des données : {response.status_code}")


print(prot_in_plasmid)
