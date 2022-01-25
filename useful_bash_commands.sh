
#find and replace strings with sed (stream editor)

# -i option is used to edit in place on the file
# -e option indicates the expression/command to run, in this case s/ (substitute)
# /g stands for global, apply command to the whole file
sed -i -e 's/|/_/g' gff_no-pipes.gff3

#remove empty lines
sed -i '/^$/d' gff_no_empty_lines.gff3

#find and remove string

awk '{gsub("gi|753572091|ref|_", "");print}' gff_no-pipes.gff3 > gff_clean.gff3

#remove line with "similarity" string in it
sed '/similarity/d' gff_clean.gff3 > gff_no_similarity.gff3

# /2 -> replace for the 2nd occurrence (of search pattern .)
# & -> back reference to the whole match
sed 's/./&_/2' gff_no_similarity.gff3 > gff_repaired.gff3

sed '/^tr/d' A0A0G2JJD3_exo.gff3 >  A0A0G2JJD3_no_tr_pipes
sed -i -e 's/|/_/g' A0A0G2JJD3_no_tr_pipes.gff3
sed -i '/similarity/d' A0A0G2JJD3_no_tr_pipes.gff3
sed -i '/^$/d' A0A0G2JJD3_no_tr_pipes.gff3
sed -e 's/^gi_753572091_ref_//g'  A0A0G2JJD3_no_tr_pipes.gff3 > A0A0G2JJD3_clean.gff3
sed -e 's/.//12' A0A0G2JJD3_clean.gff3 > A0A0G2JJD3.gff3

find directory1 -type f -size +500k -exec cp {} directory2 \;
rm `find . -type f -size +309c`

#count empty lines in file
# -c print a count of matching lines
grep -c "^$" file

#count whitespace
grep -c "^\s*$" file

# remove string after first underscore
sed -e 's/_.*//' 0_List_of_genes_found_in_CTGs

# find files in current folder and add the extension
find . -type f -exec mv '{}' '{}'.fq \;


# multi line fasta to single line fasta
seqtk seq -l0 bff_3478-KK-0002.fasta | gzip > bff_3478-KK-0002_singleLines.fa.gz

# extract sequence based on its ID, first index ref genome
 samtools faidx ref_PanTig1.0_scaffolds.fa NW_006712136.1 > ~/species_comp/felids/amur_tiger/scaffold_hits.fa 

samtools faidx ref_aciJub1_scaffolds.fa NW_015132158.1 >> ~/species_comp/felids/cheetah/scaffold_hits.fa 

cheetah scaffold_hits
NW_015131015.1
NW_015131894.1
NW_015131945.1
NW_015132158.1

seqtk seq -l0 scaffold_hits.fa > scf_hits_single_lines.fa

cat scfs
gi|753572091|ref|NC_018727.2|

samtools faidx fca_ref_Felis_catus_8.0_chrB2.fa gi|753572091|ref|NC_018727.2| > ~/species_comp/felids/domestic_cat/scaffold_hits.fa
seqtk seq -l0 scaffold_hits.fa > scf_hits_single_lines.fa