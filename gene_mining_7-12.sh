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
mkdir ORFs
mv *_ ORFs

##EXTRACT UNITIG-ORFs SEQUENCE
for f in ORFs/*; do \
	extract_fasta $f \
	>> ${f%_splitted_utg_reads_???}"_merged_ORFs" \
	; done
cd ORFs
mkdir ORFs_FASTA
mv *_ORFs ORFs_FASTA
echo "6) Utg-ORFs predicted and parsed" \
	>> PROGRESS_REPORTS

##ALIGN SEQUENCE TO UNIPROT DATABASE
for f in ORFs_FASTA/*; do \
	blastp \
	-query $f \
	-db uniprot_complete_nospace.fasta \
	-out ${f%_merged_ORFs}"_reciprocal_utg_hits" \
	-evalue 1 \
	-outfmt 6 \
	-max_target_seqs 3 \
	-num_threads 24 \
	; done
cd ORFs_FASTA
mkdir RECIPROCAL_HITS
mv *_hits RECIPROCAL_HITS
echo "7) Predicted utg-ORFs aligned to UniProt database" \
	>> PROGRESS_REPORTS

##DETERMINE CONFIRMED HITS AND GET ANNOTATION
for f in RECIPROCAL_HITS/*; do \
	awk '$11 < 1e-10 {print $1,"\t",$11,"\t",$2}' \
	> ${f%_reciprocal_utg_hits}"_annotation" $f && \
	echo ${f%_reciprocal_utg_hits} | \
	awk -F '/' '{print$2}' \
	>> ${f%_reciprocal_utg_hits}"_annotation"\
	; done
cd RECIPROCAL_HITS
mkdir ANNOTATED_RECIPROCAL_HITS
mv *_annotation ANNOTATED_RECIPROCAL_HITS

##GET ANNOTATION TO ONE LINE
for f in ANNOTATED_RECIPROCAL_HITS/*; do \
	VAR=`tail -n1 $f` ;
	awk -v var="$VAR" '{print var"\t"$0}' $f \
	> $f"_oneline" \
	; done
cd ANNOTATED_RECIPROCAL_HITS
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
	>> PROGRESS_REPORTS
for f in ONELINERS/*_beautified; do \
	cat $f | \
	awk '/=/{print$0}' \
	>> PROGRESS_REPORTS \
	; done

