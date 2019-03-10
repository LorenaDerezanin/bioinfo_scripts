#!/bin/bash

time grep "^[ACTGN]\{100\}$" -B 1 -A 2 --no-group-separator SRR5904112_1_paired_trimmed.fastq > SRR5904112_1_filtered100.fastq