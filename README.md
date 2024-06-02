# AlphaFold-pLDDT

#### Usage

```CAR_gene_list_uniprotids.txt``` has UniProt ids of proteins involved in fission yeast contractile ring. <br />
Run ```fetch_af_predictions.py``` to fetch the AlphaFold predictions for these proteins. <br />
Then, run ```alphafold_get_plddt.py``` to obtain per residue pLDDT in a TSV file.

The script takes in input a folder with PDB files and writes a TSV file.

    python3 alphafold_get_plddt.py -i pdbs/ -o out.tsv

```separate_plddt_for_proteins_in_sheets.py``` separates the pLDDT scores of each protein in a different sheet in ```pdata.xlsx```

#### Citation

Piovesan D, Monzon AM, Tosatto SCE.<br />
Intrinsic protein disorder and conditional folding in AlphaFoldDB. 
Protein Sci. 2022 Nov;31(11):e4466.<br />
PMID: [36210722](https://pubmed.ncbi.nlm.nih.gov/36210722/)
PMCID: [PMC9601767](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9601767/).
