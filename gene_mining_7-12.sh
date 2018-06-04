#!/bin/sh

set -e

WORKING_DIR_PATH="/home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly"
PROGRESS_REPORTS=$WORKING_DIR_PATH/gadus_CA/gene_mining_results/Progress_gadus_3

##PREDICT ORF FOR ALL GENES IN EACH UNITIG
for f in FASTA/*_splitted_utg_reads_*; do \
	hmmer2go getorf \
	-i $f > \
	-o $f"_" \
    --verbose \
    -n 4 \
	; done
cd FASTA
mkdir 2_ORFs
mv *_ 2_ORFs

##EXTRACT UNITIG-2_ORFs SEQUENCE
# look for file names looking like (some_string)_xx(_clean.fasta_orfs.fasta), and remember the parts in parentheses
re='(.+)_[a-z]{2}(_clean.fasta_orfs.fasta)'
for f in 2_ORFs/*; do 
	if [[ $f =~ $re ]]; then
		# construct a file name like some_string_clean.fasta_orfs.fasta based on the current file name
		# (notice we removed the middle _xx part)
		merged=${BASH_REMATCH[1]}${BASH_REMATCH[2]}"_merged_2_ORFs"
		# append the matching input files into the same output files
		cat $f >> $merged
	fi
done

cd 2_ORFs
mkdir 2_ORFs_FASTA
mv *_2_ORFs 2_ORFs_FASTA
echo "6) Utg-2_ORFs predicted and parsed" \
	>> $PROGRESS_REPORTS

# cannot trace extract_fasta command

# uniprot_sprot.fasta.gz (last update on 25.04.18) downloaded on 15.5.18
# uniprot_trembl.fasta.gz (last update on 25.04.18) downloaded on 15.5.18

# prepared merged database file:
time makeblastdb -in merged_uniprotdbs.fasta -dbtype prot -out merged_uniprotdbs 

# Error: (1431.1) FASTA-Reader: Warning: FASTA-Reader: Ignoring FASTA modifier(s) found because the input was not expected to have any
# this error/warning message is fixed in newer BLAST+ versions, appears only in our version (2.2.29), should be ignored

##ALIGN SEQUENCE TO UNIPROT DATABASE
for f in 1_ORFs_FASTA/*; do \
	blastp \
	-query $f \
	-db merged_uniprotdbs.fasta \
	-out ${f%_clean.fasta_orfs.fasta}"_reciprocal_utg_hits" \
	-evalue 1 \
	-outfmt 6 \
	-max_target_seqs 3 \
	-num_threads 12 \
	; done
cd 1_ORFs_FASTA
mkdir 1_RECIPROCAL_HITS
mv *_hits 1_RECIPROCAL_HITS
echo "7) Predicted utg-2_ORFs aligned to UniProt database" \
	>> $PROGRESS_REPORTS

#time: 218 min = 3,6 h

##DETERMINE CONFIRMED HITS AND GET ANNOTATION
for f in 1_RECIPROCAL_HITS/*; do \
	awk '$11 < 1e-10 {print $1,"\t",$11,"\t",$2}' \
	> ${f%_reciprocal_utg_hits}"_annotation" $f && \
	echo ${f%_reciprocal_utg_hits} | \
	awk -F '/' '{print$2}' \
	>> ${f%_reciprocal_utg_hits}"_annotation"\
	; done
cd 1_RECIPROCAL_HITS
mkdir ANNOTATED_1_RECIPROCAL_HITS
mv *_annotation ANNOTATED_1_RECIPROCAL_HITS

##GET ANNOTATION TO ONE LINE
for f in ANNOTATED_1_RECIPROCAL_HITS/*; do \
	VAR=`tail -n1 $f` ;
	awk -v var="$VAR" '{print var"\t"$0}' $f \
	> $f"_oneline" \
	; done
cd ANNOTATED_1_RECIPROCAL_HITS
mkdir ONELINERS
mv *_oneline ONELINERS

##REPORT CONFIRMED HITS AND CONTINUE TO SEARCH THE SINGLETONS
for f in ONELINERS/*; do \
	cat $f | \
	awk -F "|" '{$2="";print $0}' \
	> $f"_beautified" \
	; done
for f in ONELINERS/*_beautified; do \
	head $f | \
	awk '/=/{print$1}' | \
	sort -u \
	>> List_of_genes_found_in_utgs \
	; done
echo "8) Genes present in UTGs, based on reciprocal hits (e-value < 1e-10):" \
	>> $$PROGRESS_REPORTS
	
for f in ONELINERS/*_beautified; do \
	cat $f | \
	awk '/=/{print$0}' \
	>> $$PROGRESS_REPORTS \
	; done


# list of genes empty, reslove extract_fasta command before moving to the singletons

