#!/usr/bin/python

# run multiple exonerate runs with each protein seq file in dir as a query

import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", dest="input", action="store", required=True)
parser.add_argument("-o", "--output", dest="output", action="store", required=True)
parser.add_argument("-n", "--bestn", dest="bestn", action="store", type=int)
parser.add_argument("-s", "--score", dest="score", action="store", type=int, default=100)
parser.add_argument("-D", "--dpmemory", dest="dpm", action="store", type=int, default=32)



args = parser.parse_args()


def gene_name(protein_name):
    split_name = protein_name.split(".fasta")
    return split_name[0]


path = args.input

protein_files = os.listdir(path)
for protein_file in protein_files:
    f_out_path = os.path.join(args.output, gene_name(protein_file) + "_exo")
    f_in_path = os.path.join(path, protein_file)
    exo = "exonerate -m protein2genome -q %s -t %s --bestn %s --dpmemory %s -s %s -V 3 -Q protein -T dna " \
          "--showalignment yes --showtargetgff yes" % (f_in_path, f_out_path, args.bestn, args.dpm, args.score)
    os.system(exo)


