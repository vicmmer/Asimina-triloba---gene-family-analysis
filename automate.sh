#Purpose
#Some steps in this pipeline (e.g., OrthoFinder, InterProScan, CAFE) can take many hours to days when run on the full dataset.
#To allow for quick troubleshooting, this repository supports a subset troubleshooting mode using a symbolic link (symlink).
#This approach avoids editing multiple scripts and ensures the full dataset is not overwritten.

#Concept
#All pipeline scripts read protein FASTA files from the same directory: protein_sequences/
#Instead of modifying scripts to handle different datasets, this directory is implemented as a symbolic link (symlink) that can point
#to either:
  # the full proteome dataset, or
  # a small subset of protein FASTAs for quick testing

#The pipeline itself does not need to change.

#Directory Structure
    #protein_sequences_full/      Full dataset (never modified)
    #protein_sequences_subset/    Small subset for troubleshooting
    #protein_sequences/           Symlink -> points to ONE of the above

- protein_sequences_full: contains the complete proteome set used for
  final analyses
- protein_sequences_subset: subsets the first 200 sequences in each FASTA file
- protein_sequences/: symlink that controls which dataset the pipeline uses
      # NOTE: Have to figure out still exactly how to work the symlink 


#AUTOMATED PIPELINE WRAPPER 
#!/bin/bash
set -e

# Purpose:
#   Run the gene family evolution pipeline sequentially.

echo " Starting Gene Family Evolution Pipeline "
echo

# 0) Data preparation 
echo "Step 0: Download accessions "
./0.download_accessions.sh
echo

# 1) Genome completeness 
echo " Step 1a: BUSCO "
./1a.busco.sh
echo

echo " Step 1b: BUSCO summaries "
./1b.busco_summaries.sh
echo

#  2) Orthogroup inference 
echo " Step 2: OrthoFinder "
./2.orthofinder.sh
echo

echo " Step 2b: Prepare InterProScan input "
./2b.prepare_interproinput.sh
echo

#  3) Functional annotation 
echo " Step 3: InterProScan "
./3.interproscan.sh
echo

#  4) Gene family overlap + CAFE prep 
echo " Step 4: UpSet + CAFE prep "
Rscript 4.Upset_PrepCafe.R
Rscript 4.make_upset_filtered.R
echo

#  5) Gene family evolution modeling 
echo " Step 5: CAFE "
./5.cafe.sh
echo

#  6) Phylogenetic dating 
echo " Step 6: treePL "
treePL 6.config.cfg
echo

#  7) GO enrichment analysis 
echo " Step 7a: topGO prep "
./7a.topgoprep.sh
echo

echo " Step 7b: topGO "
Rscript 7b.topgo.R
echo

echo " Pipeline complete "

