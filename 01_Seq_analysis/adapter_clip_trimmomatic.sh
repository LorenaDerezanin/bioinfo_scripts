 #!/bin/sh

 TRIMMOMATIC=/usr/local/bioinf/Trimmomatic-0.35/trimmomatic-0.35.jar
 DIR_PATH=$1
 cd $DIR_PATH

 SRA_IDS=`ls | grep "^SRR[0-9]\+_1.fastq.gz$" | grep -o "^SRR[0-9]\+"`

for SRA_ID in SRA_IDS; do \
   java -jar $TRIMMOMATIC PE \
         -threads 5 \
         $DIR_PATH/"$SRA_ID"_1.fastq.gz \
         $DIR_PATH/"$SRA_ID"_2.fastq.gz \
         $DIR_PATH/"$SRA_ID"_1_paired.fq.gz \
         $DIR_PATH/"$SRA_ID"_1_unpaired.fq.gz \
         $DIR_PATH/"$SRA_ID"_2_paired.fq.gz \
         $DIR_PATH/"$SRA_ID"_2_unpaired.fq.gz \
         ILLUMINACLIP:TruSeq_Nextera-PE.fa:2:30:10 \
         LEADING:3 \
         TRAILING:3 \
         SLIDINGWINDOW:4:15 \
         MINLEN:30 \
         TOPHRED33 \
         ; done


