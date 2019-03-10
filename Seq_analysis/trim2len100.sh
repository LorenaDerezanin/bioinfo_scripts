 #!/bin/sh

 FILE_PATH=/home/derezanin/temp_storage/striped_hyena/raw_reads/in_size670/
 TRIMMOMATIC=/usr/local/bioinf/Trimmomatic-0.35/trimmomatic-0.35.jar

 java -jar $TRIMMOMATIC SE \
         -threads 5 \
         $FILE_PATH/adapt_clipped/SRR5904112_2_paired.fastq \
         $FILE_PATH/len100bp/SRR5904112_2_paired_trimmed.fastq\
         CROP:100
       


