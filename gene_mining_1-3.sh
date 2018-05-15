#!/bin/sh

# this stops the script when a command doesn't return 0 (when it throws an error)
set -e

##SETUP INPUT VARIABLES
#FISH=$1
WORKING_DIR_PATH="/home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly"

#PROGRESS_REPORTS="/data/fg2/derezanin/species_comp/atlantic_cod/fish_paired_reads/gene_mining_results/Progressreports_2018/Progress_$FISH"
PROGRESS_REPORTS=$WORKING_DIR_PATH/gadus_CA/gene_mining_results/Progress_gadus_3

# run makedb outside of script
# makeblastdb -in $FISH.utg.fasta -dbtype nucl
# makeblastdb -in gadus_3_CA.utg.fasta -dbtype nucl


##MAKE FOLDERS
mkdir -p $WORKING_DIR_PATH/gadus_CA/9-terminator/DATABASES
cd $WORKING_DIR_PATH/gadus_CA/9-terminator/DATABASES
mkdir -p BLASTS
cd BLASTS
mkdir -p UTG
mkdir -p SIN
echo "1) Initial folders created" \
	> $PROGRESS_REPORTS

###################################### UNITIGS SEARCH ######################################
##ALIGN GENELIST TO UNITIGS (full protein seqs of 27 immune genes + 3 control genes from 10 species, merged together in 1 file)
cd UTG
# tblastn (query=protein vs. database=nucleotide seqs which get translated,so alignm. type is protein vs. protein)
# -max_target_seqs 10 max number of saved database matches for a query (genelist), assuming there are no more than 10 gene copies?
# -db ../../../$FISH.utg.fasta \
tblastn \
	-query $WORKING_DIR_PATH/gadus_CA/9-terminator/Complete_genelist.fas \
	-db $WORKING_DIR_PATH/gadus_CA/9-terminator/gadus_3_CA.utg.fasta \
	-out utg-blastout \
	-evalue 1e-10 \
 	-outfmt 6 \
	-max_target_seqs 10 \
	-num_threads 20
echo "2) Unitig BLAST finished" \
	>> $PROGRESS_REPORTS


##PARSE BEST UNITIG HIT FOR EACH GENE TO A NEW FILE
cut utg-blastout \
	-f1,2 | \
	sort -u \
	> All_unique_UTG_hits_sorted
cat All_unique_UTG_hits_sorted | \
	awk -F '_' '{print$1}' | \
	sort -u \
	> UTG_best_hits
echo "3) Found initial hits (e-value < 1e-10) in unitigs for the following genes:" \
	>> $PROGRESS_REPORTS
cat UTG_best_hits \
	>> $PROGRESS_REPORTS


	# time ~3,33 min (22.04.2018)
