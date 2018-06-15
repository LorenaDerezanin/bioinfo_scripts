#!/bin/sh

set -e

WORKING_DIR_PATH="/home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly"
PROGRESS_REPORTS=$WORKING_DIR_PATH/gadus_CA/gene_mining_results/Progress_gadus_3

# need to prepare unitig ids first, before running hmmer2go, otherwise it spits error and empty output files
# get https://github.com/sestaton/sesbio/blob/master/gene_annotation/clean_multifasta.pl
# run clean_multifasta.pl script to keep only the utg ID for EMBOSS getorf search

# PREPARING UTG IDs FOR HMMER2GO
for f in FASTA/*_splitted_UTG_reads_*; do \
	./clean_multifasta.pl \
	-i $f \
	-o $f"_" \


##PREDICT ORF FOR ALL GENES IN EACH UNITIG
for f in FASTA/*_*; do \
	hmmer2go getorf \
	-i $f \
	-o ${f%_splitted_UTG_reads_}"orfs.fasta" \
	--all \
	-l 100 \
    --verbose \
    -n 4 \
	; done
cd FASTA
mkdir 2_ORFs
mv *orfs.fasta 2_ORFs

##EXTRACT UNITIG-2_ORFs SEQUENCE
# look for file names looking like (some_string)_xx(_orfs.fasta), and remember the parts in parentheses
re='(.+)_[a-z]{2}(_orfs.fasta)'
for f in 2_ORFs/*; do 
	if [[ $f =~ $re ]]; then
		# construct a file name like some_string_clean.fasta_orfs.fasta based on the current file name
		# (notice we removed the middle _xx part)
		merged=${BASH_REMATCH[1]}${BASH_REMATCH[2]}"_merged_ORFs"
		# append the matching input files into the same output files
		cat $f >> $merged
	fi
done

cd 2_ORFs
mkdir 2_ORFs_FASTA
mv *_merged_ORFs 2_ORFs_FASTA
echo "6) Utg-2_ORFs predicted and parsed" \
	>> $PROGRESS_REPORTS


# uniprot_sprot.fasta.gz (last update on 25.04.18) downloaded on 15.5.18
# uniprot_trembl.fasta.gz (last update on 25.04.18) downloaded on 15.5.18

# replace all white spaces with "$" sign, to get the full name of the protein during blastp step
sed "s/ /$/g" merged_uniprotdbs.fasta > merged_uniprotdbs_nospace.fasta 

# prepared merged database file:
time makeblastdb -in merged_uniprotdbs_nospace.fasta -dbtype prot -out merged_uniprotdbs_nospace 
# time on allegro:  (20 GB, 16 CPUs)

# Error: (1431.1) FASTA-Reader: Warning: FASTA-Reader: Ignoring FASTA modifier(s) found because the input was not expected to have any
# this error/warning message is fixed in newer BLAST+ versions, appears only in the version 2.2.29, should be ignored

##ALIGN SEQUENCE TO UNIPROT DATABASE
for f in 2_ORFs_FASTA/*; do \
	blastp \
	-query $f \
	-db merged_uniprotdbs \
	-out ${f%_splitted_UTG_reads_orfs.fasta_merged_ORFs}"_reciprocal_utg_hits" \
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

# time: 218 min = 3,6 h (for only 1 ORF per unitig)
# time on allegro:  ~ 13,3 h, 20 GB, 16 threads (multiple ORFs per unitig, merged)


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

# add print 1st column if there is 2nd column present	
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



