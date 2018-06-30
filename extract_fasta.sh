#!/bin/bash

ID_FILE=$1
SEQ_FILE=$2

awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f'  "$1" "$2"
