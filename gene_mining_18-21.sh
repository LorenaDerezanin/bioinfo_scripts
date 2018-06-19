#!/bin/bash

set -e 

WORKING_DIR_PATH="/home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly"
PROGRESS_REPORTS=$WORKING_DIR_PATH/gadus_CA/gene_mining_results/Progress_gadus_3


##EXTRACT READS IN FASTA FORMAT FOR EACH GENE
for f in SINHITS/*; do \
	cut -f2 $f | \
	sort -u \
	> $f"_sin_reads_only" \
	; done
for f in SINHITS/*sin_reads_only; do \
	fastagrep \
	-t \
	-f $f ../../$FISH.singleton.fasta \
	> ${f%_only}"_fasta" \
	; done
cd SINHITS
mkdir FASTA
mv *_fasta FASTA
echo "13) Singleton-sequences extracted" \
	>> ~/Progressreports_2015/Progress_$FISH

##GET EACH SINGLETON TO SEPARATE FILES BEFORE ORF PREDICTION
for f in FASTA/*; do \
	split -l 2 $f \
	${f%_sin_reads_fasta}"_splitted_singletons_" \
	; done


# PREPARING UTG IDs FOR HMMER2GO
for f in FASTA/*_splitted_singletons_*; do \
	/home/derezanin/species_comp/scripts/perl_scripts/clean_multifasta.pl \
	-i $f \
	-o $f"_" \
	; done

##PREDICT ORF FOR ALL GENES IN EACH SINGLETON
for f in 1_CLEAN_FASTA/*_*; do \
	hmmer2go getorf \
	-i $f \
	-o ${f%_splitted_singletons_}"orfs.fasta" \
	--all \
	-l 100 \
   	--verbose \
   	-n 4 \
	; done

cd 1_CLEAN_FASTA
mkdir SIN_ORFs
mv *orfs.fasta SIN_ORFs


##EXTRACT SINGLETON ORFs SEQUENCE
# look for file names looking like (some_string)_xx(_orfs.fasta), and remember the parts in parentheses
re='(.+)_[a-z]{2}(_orfs.fasta)'
for f in SIN_ORFs/*; do 
	if [[ $f =~ $re ]]; then
		# construct a file name like some_string_clean.fasta_orfs.fasta based on the current file name
		# (notice we removed the middle _xx part)
		merged=${BASH_REMATCH[1]}${BASH_REMATCH[2]}"_merged_ORFs"
		# append the matching input files into the same output files
		cat $f >> $merged
	fi
done

cd SIN_ORFs
mkdir ORFs_FASTA
mv *_ORFs ORFs_FASTA
echo "14) Singleton-ORFs predicted and parsed" \
	>> ~/Progressreports_2015/Progress_$FISH


##ALIGN SINGLETONS TO UNIPROT DATABASE
for f in ORFs_FASTA/*; do \
        blastp \
        -query $f \
        -db $WORKING_DIR_PATH/gadus_CA/9-terminator/DATABASES/BLASTS/UTG/merged_uniprotdbs_nospace \
        -out ${f%_merged_ORFs}"_reciprocal_sin_hits" \
        -evalue 1 \
        -outfmt 6 \
        -max_target_seqs 20 \
        -num_threads 16 \
        ; done

cd ORFs_FASTA
mkdir RECIPROCAL_HITS
mv *_hits RECIPROCAL_HITS
echo "15) Predicted Singleton-ORFs aligned to UniProt database" \
	>> ~/Progressreports_2015/Progress_$FISH

##DETERMINE CONFIRMED HITS AND GET ANNOTATION
for f in RECIPROCAL_HITS/*; do \
	awk '$11 < 1e-1 {print $1,"\t",$11,"\t",$2}' \
	> ${f%_reciprocal_sin_hits}"_annotation" $f && \
	echo ${f%_reciprocal_sin_hits} | \
	awk -F '/' '{print$2}' \
	>> ${f%_reciprocal_sin_hits}"_annotation"\
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

##REPORT CONFIRMED HITS 
for f in ONELINERS/*; do \
	cat $f | \
	awk -F "|" '{$2="";print $0}' \
	> $f"_beautified" \
	; done
for f in ONELINERS/*_beautified; do \
	head $f | \
	awk '/=/{print$1}' \
	>> List_of_genes_found_in_singletons \
	; done
echo "16) Genes present in singletons, based on reciprocal hits (e-value < 1e-1):" \
	>> ~/Progressreports_2015/Progress_$FISH
for f in ONELINERS/*_beautified; do \
	head $f | \
	awk '/=/{print$0}' \
	>> ~/Progressreports_2015/Progress_$FISH \
	; done

##REPORT WHICH GENES ARE ABSENT 
mv List_of_genes_found_in_singletons ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/
cd ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/
cat List_of_genes_found_in_singletons \
	List_of_genes_found_in_utgs \
	> All_genes_present_in_$FISH
grep \
	-f All_genes_present_in_$FISH \
	-x \
	-v \
	~/Gene_names_only \
	> Genes_absent
echo "17) These genes are NOT present in the genome:" \
	>> ~/Progressreports_2015/Progress_$FISH
cat Genes_absent \
>> ~/Progressreports_2015/Progress_$FISH