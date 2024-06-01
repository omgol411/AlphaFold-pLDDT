# https://alphafold.ebi.ac.uk/files/AF-O13923-F1-model_v4.pdb

import os
import urllib
import requests

def get_alphafold_download_link(uniprot_id):
	link_pattern = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'
	return link_pattern.format(uniprot_id)
	
gene_list = open("CAR_gene_list_uniprotids.txt").read().splitlines()

save_path = '.\pdbs'

for gene in gene_list:
	gene_url = get_alphafold_download_link(gene)
	r = requests.get(gene_url)
	pdb_filename = os.path.basename(r.url)
	completeName = os.path.join(save_path, pdb_filename)       
	f = open(completeName, 'wb')
	f.write(r.content)
	f.close()
