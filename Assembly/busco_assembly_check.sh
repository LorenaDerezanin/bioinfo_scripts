
# example of a genome assembly quality check up
# BUSCO version: BUSCO 3.0.1

# run it directly in the command line, not as a script

time run_BUSCO -i /home/derezanin/species_comp/wolverine/CYRY01.fasta \
-o wolverine_busco -l /home/derezanin/species_comp/wolverine/mammalia_odb9 \
-m geno -c 16 \
2> busco_wolverine.log

# time: 768 min (14.8.18)

time run_BUSCO -i bff_haplo_3478-KK-0002.fasta \
-o bff_busco -l /home/derezanin/species_comp/wolverine/mammalia_odb9 -m geno \
-c 16 \
2> busco_bff.log

# worked without full path to the input file
# time(bff): 608 min (12.9.18)
# time(canada lynx): 714 min (18.9.18) incomplete assembly 

time rsync -avPe "ssh -p 22220" 3478-KK-0002_S5_L005_R2_001.fastq.gz derezanin@62.141.164.6:/home/derezanin/temp_storage/BFF/raw_reads/
# continue with other reads




time blastx -query /home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly/e-value_cutoff_check/exon_MHC.fasta \
-db /home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly/e-value_cutoff_check/MHC-U_NCBI.faa -out ref_full2.out \
-outfmt 0 \
-num_threads 16 \
2> ref_full.err

time tblastn -query /home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly/e-value_cutoff_check/MHC-U_NCBI.faa \
-db /home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly/e-value_cutoff_check/full_nucl -out 2tblastn_out \
-outfmt 0 \
-num_threads 16 




2> ref_tblastn.err



-evalue 1e-5

-max_target_seqs 10 \