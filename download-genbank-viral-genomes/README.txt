The current data is from a 2018-10-19 pull of NCBI's viral
accession list. Note that the accession list in this
directory has Influenza removed (see influenza/ for downloading
Influenza data). To create the following:
  - datasets.txt
  - all data in all/data/
  - all dataset .py files in all/
  - seqs-in-201810-pull*.txt
I ran ./run_from_lineages.sh, which takes lineages.txt as
input. This creates one dataset per species. See
run_from_lineages.sh for more detail.
