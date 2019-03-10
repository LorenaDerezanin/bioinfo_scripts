 #!/bin/sh

  SRA_DIR_PATH=$1
  cd $SRA_DIR_PATH

  SRA_IDS=`ls | grep "^SRR[0-9]\+$"`
  

 for SRA_ID in $SRA_IDS; do \
 	fastq-dump.2.9.2 --split-files $SRA_ID && rm $SRA_ID && gzip $SRA_ID"_*.fastq" \
 	; done
