#!/bin/bash

for file in "/data/fg2/derezanin/species_comp/tests/mhc_rel_sub_gff_cat_vs_hs/*"
do
  sed -i -e '/^tr/d' $file
  sed -i -e'/similarity/d' $file
  sed -i '/^$/d' $file
  sed -i -e 's/^gi|753572091|ref|//g' $file
  sed -i -e 's/.//12' $file
  sed -i -e '/^##gff-ver/d' $file
  #if there's different output file, pipe it
  #sed -i '/^tr/d' $file | sed '/similarity/d' - | sed '/^$/d' - | -e 's/^gi|753572091|ref|//g' - | sed -e 's/.//12' - | sed '/^##gff-ver/d' - > output_file
        done

    #or make a file containing sed commands and pass it to sed
    #sed_commands.txt:
    # /similarity/d
    # /^$/d
    # s/^gi|753572091|ref|//g
    #....

  #sed -i -f sed_commands.txt $file
