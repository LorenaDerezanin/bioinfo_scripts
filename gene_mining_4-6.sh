#!/bin/sh

# script adjusted  from: https://github.com/uio-cees/teleost_genomes_immune

set -e

##SETUP INPUT VARIABLES
#FISH=$1

WORKING_DIR_PATH="/home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly"
PROGRESS_REPORTS=$WORKING_DIR_PATH/gadus_CA/gene_mining_results/Progress_gadus_3

cd $WORKING_DIR_PATH/gadus_CA/9-terminator/DATABASES/BLASTS/UTG
mkdir -p UTG_HITS
for f in `cat UTG_best_hits`; do \
	grep $f All_unique_UTG_hits_sorted \
	>> UTG_HITS/$f; done
echo "4) Gene-containing UTGs parsed to files" \
	>> $PROGRESS_REPORTS

##EXTRACT READS IN FASTA FORMAT FOR EACH GENE
for f in UTG_HITS/*; do \
	cut -f2 $f | \
	sort -u \
	> $f"_utg_reads_only" \
	; done

# use parameter -j16 ,assigns 16 threads to parallel, otherwise uses all resources
echo "Doing fastagrep in parallel"
ls UTG_HITS/*utg_reads_only | \
	time parallel -j16 '/home/derezanin/species_comp/atlantic_cod/fish_paired_reads/\
gadus_3_assembly/gadus_CA/9-terminator/\
teleost_genomes_data_descriptor/ortholog_identification/scripts/resources/\
fastagrep \
-t \
-f {} /home/derezanin/species_comp/atlantic_cod/fish_paired_reads/\
gadus_3_assembly/gadus_CA/9-terminator/gadus_3_CA.utg.fasta \
> {}"_fasta"'

#for f in UTG_HITS/*UTG_reads_only; do \
# $WORKING_DIR_PATH/gadus_CA/9-terminator/teleost_genomes_data_descriptor/\
# ortholog_identification/scripts/resources/\
# fastagrep \
# 	-t \
# 	-f $f $WORKING_DIR_PATH/gadus_CA/9-terminator/gadus_3_CA.utg.fasta \
# 	> ${f%_only}"_fasta" \
# 	; done

cd UTG_HITS
mkdir -p FASTA
mv *_fasta FASTA
echo "5) UTG-sequences extracted" \
	>> $PROGRESS_REPORTS

##GET EACH UNITIG TO SEPARATE FILES BEFORE ORF PREDICTION
for f in FASTA/*; do \
	split -l 2 $f \
	${f%_utg_reads_fasta}"_splitted_UTG_reads_" \
	; done


	# time ~38 min (23.04.2018)
	# time with parallel fastagrep, 16 CPUs assigned ~ 4,4 min (24.04.2018)
	
