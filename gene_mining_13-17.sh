#!/bin/sh

set -e

##SETUP VARIABLES

WORKING_DIR_PATH="/home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly"
PROGRESS_REPORTS=$WORKING_DIR_PATH/gadus_CA/gene_mining_results/Progress_gadus_3

##################################### SINGLETON SEARCH #####################################

##MAKE A LIST OF GENES FOUND SO FAR
mv List_of_genes_found_in_utgs ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/
cd ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/

##SELECT ONLY THE GENES THAT WERE NOT FOUND IN SCAFFOLDS OR UNITIGS
grep \
	-f List_of_genes_found_in_utgs \
	-x \
	-v \
	~/Gene_names_only \
	> Still_no_good_hits
echo "9) Looking for the following genes in the singelton reads:" \
	>> $PROGRESS_REPORTS
cat Still_no_good_hits \
	>> $PROGRESS_REPORTS

##EXSTRACT FASTA SEQUENCE FOR THE SINGLETON BLAST SEARCH
fastagrep \
	-t \
	-f Still_no_good_hits \
	~/Complete_genelist.fas \
	> Still_no_good_hits_fas

##ALIGN REDUCED GENELIST TO SINGLETON READS
tblastn \
	-query Still_no_good_hits_fas \
	-db ../../$FISH.singleton.fasta \
	-out sin-blastout \
	-evalue 1e-1 \
	-outfmt 6 \
	-max_target_seqs 10 \
	-num_threads 24
echo "10) Singleton BLAST finished" \
	>> $PROGRESS_REPORTS

##PARSE BEST HIT FOR EACH GENE TO NEW FILE
cut sin-blastout \
	-f1,2 | \
	sort -u \
	> All_uniqe_sin_hits_sorted
cat All_uniqe_sin_hits_sorted | \
	awk -F '_' '{print$1}' | \
	sort -u \
	> sin_hits
echo "11) Found initial hits in singletons (e-value < 1e-1) for the following genes:" \
	>> $PROGRESS_REPORTS
cat sin_hits \
	>> $PROGRESS_REPORTS
mkdir SIN_HITS
for f in `cat sin_hits`; do \
	grep $f All_uniqe_sin_hits_sorted \
	>> SIN_HITS/$f; done
echo "12) Gene-containing singletons parsed to files" \
	>> $PROGRESS_REPORTS
