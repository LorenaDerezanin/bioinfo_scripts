#!/usr/bin/python
#-*- coding: utf-8 -*-

# split multifasta file into separate fasta files

import os
import sys
from optparse import OptionParser
from Bio import SeqIO
# define usage and args: prog arg[0] arg[1]
usage = "usage: %prog fasta_file_in directory_out"
parser = OptionParser(usage)
(opts, args) = parser.parse_args()
dir_out = os.getcwd()

if len(args)<1:
    print "Error: Please enter at least one argument."
    print "See program_name.py --help"
    sys.exit()
elif len(args)==2:
    dir_out = args[1]
elif len(args)>2:
    print "error: Please enter up to 2 arguments."
    print "See program_name.py --help"
    sys.exit()

# function for creating output file name
def build_name(id):
    split_id1 = id.split("|")
    split_id2 = split_id1[2].split("_")
    short_id = split_id1[1] + "_" + split_id2[0]
    return short_id

file_name = args[0]

file_in = open(file_name)
# SeqIO parses fasta files as generator, doesn't place it in memory
records = SeqIO.parse(file_in, "fasta")
for record in records:
    if "FLA" in record.description:
        f_out_path = os.path.join(dir_out, build_name(record.id) + '.fasta')
        f_out = open(f_out_path, 'w')
        SeqIO.write([record], f_out, "fasta")

# split multifasta in equal sized chunks
# import pyfasta
# from pyfasta import Fasta
# pyfasta split -n 10 proteome.fasta
#pyfasta split -n 25 proteome_chr6_hs.fasta
