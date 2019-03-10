 #!/bin/sh

 DIR_PATH=$1
 cd $DIR_PATH

for fastq_file in $DIR_PATH; do \
	fastqc $fastq_file \
	;done 

 rm *.zip 
 mkdir fastqc_reports
 mv *.html fastqc_reports/