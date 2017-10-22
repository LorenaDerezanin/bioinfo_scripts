#!/usr/bin/env python

# written by Lorena Derezanin 

import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", dest="input", action="store", required=True)
parser.add_argument("-o", "--output", dest="output", action="store", required=True)
parser.add_argument("-db", "--database", dest="database", action="store", required=True)
parser.add_argument("-e", "--evalue", dest="evalue", action="store", type=float)
parser.add_argument("-word_size", dest="word_size", action="store", type=int)

args = parser.parse_args()


def prot_name(protein):
    split_name = protein.split(".fasta")
    return split_name[0]


path = args.input
protein_files = os.listdir(path)
# os.listdir returns list of file names as strings
for protein_file in protein_files:
    f_out_path = os.path.join(args.output, prot_name(protein_file) + "_blastp")
    # os.path.join concatenates strings with "/"
    f_in_path = os.path.join(path, protein_file)
    blast = "blastp -query %s -out %s -db %s -evalue %s -word_size %s" % (f_in_path, f_out_path, args.database, args.evalue, args.word_size)
    os.system(blast)
# os.system('command') runs command in bash

# "blast -query(in_file) -out(out_file) -evalue 1e-2 -word_size 6"
