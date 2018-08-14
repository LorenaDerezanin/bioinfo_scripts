
# example of a genome assembly quality check up
# BUSCO version: BUSCO 3.0.1

# run it directly in the command line, not as a script

time run_BUSCO -i /home/derezanin/species_comp/wolverine/CYRY01.fasta \
-o wolverine_busco -l /home/derezanin/species_comp/wolverine/mammalia_odb9 \
-m geno -c 16 \
2> busco_wolverine.log

# time: 768 min (14.8.18)


time tblastn -query /home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly/e-value_cutoff_check/MHC-U_NCBI.faa \
-db /home/derezanin/species_comp/atlantic_cod/fish_paired_reads/gadus_3_assembly/e-value_cutoff_check/ref_g_50.fasta -out ref_g50_e5.out \
-evalue 1e-5 -outfmt 6 \
-num_threads 16 \
2> ref_g50_e5.err


-max_target_seqs 10 \