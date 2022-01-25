
#comparison of human MHC protein catalogue with genomes and chromosomes of the following species:
    #amur tiger, dog, cat, leopard, giant panda, polar bear

#test 1 - 1 mhc protein vs. cat chrB2
#exonerate -m protein2genome --showalignment yes --showquerygff yes --showtargetgff yes --dpmemory 16384 -q human/hs_mhc_proteins/A0A0D5XQ77_A0A0D5XQ77.fasta  -t cat/fca_ref_Felis_catus_8.0_chrB2.fa
#run same job with -n 5
#exonerate -m protein2genome --showalignment yes --showquerygff yes --showtargetgff yes --dpmemory 16384 -n 5 -q human/hs_mhc_proteins/A0A0D5XQ77_A0A0D5XQ77.fasta  -t cat/fca_ref_Felis_catus_8.0_chrB2.fa

#test 2 - ryo argument included
#exonerate -m protein2genome --showalignment no --showquerygff yes --showtargetgff yes --fsmmemory 2048 --bestn 5 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" -Q Protein -T DNA -q human/hs_mhc_proteins/A0A0D5XQ77_A0A0D5XQ77.fasta -t cat/fca_ref_Felis_catus_8.0_chrB2.fa
#is no work

#run multi_exo.py with whole set of mhc proteins vs. cat chrB2, dog chr12,7,18,35, genomes of amur tiger, leopard, giant panda, polar bear


import sys, os

def gene_name(protein_name):
    split_name = protein_name.split(".fasta")
    return split_name[0]

path = "/data/fg2/derezanin/species_comp/human/hs_mhc_proteins/"
# os.listdir returns list of file names as strings
protein_files = os.listdir(path)
for protein_file in protein_files:
    # f_out_path contains output file path as string
    # os.path.join concatenates strings with "/"
    f_out_path = os.path.join("/data/fg2/derezanin/species_comp/tests/cat_ryo_hs_out/", gene_name(protein_file) + "_exo")
    f_in_path = os.path.join(path, protein_file)
    #remove exonerate parameter -n 5, otherwise it only writes 5 best overall matches in the output dir
    exo = "exonerate -m protein2genome --showalignment no --showquerygff yes --showtargetgff yes --fsmmemory 2048 --bestn 5 -Q Protein -T DNA --ryo \">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n\" -q " + f_in_path + " -t cat/fca_ref_Felis_catus_8.0_chrB2.fa > " + f_out_path
    # os.system('command') runs command in bash
    os.system(exo)

#run in screen
#time scripts/comparison.py 2> tests/cat_ryo_hs_test.log
