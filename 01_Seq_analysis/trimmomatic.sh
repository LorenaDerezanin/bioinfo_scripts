#!/bin/sh

 FILE_PATH=/home/derezanin/temp_storage/bears/brown_bear/swedish_b_bear/raw_reads
 TRIMMOMATIC=/usr/local/bioinf/Trimmomatic-0.35/trimmomatic-0.35.jar

 java -jar $TRIMMOMATIC SE \
         $FILE_PATH/SRR1693654_2.fastq\
         $FILE_PATH/SRR1693654_2._trimmed.fastq \
         CROP:100
       


