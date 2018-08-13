#!/bin/bash

# example of a genome assembly quality check up
# BUSCO version: BUSCO 3.0.1
PATH="/home/derezanin/species_comp/wolverine"

time run_BUSCO -i $PATH/CYRY01.fasta \
-o wolverine_busco -l $PATH/mammalia_odb9 \
-m geno -c 16 \
2> busco_wolverine.log
