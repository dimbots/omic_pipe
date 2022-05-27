#!/bin/bash

#############################################################################################################################################

	tput setaf 6; tput bold; echo " "
	tput setaf 1; tput bold; echo "        INITIATE PIPELINE"
	tput setaf 2; tput bold; echo " "

	tput setaf 3; tput bold; echo "    Quality Control & Trimming"
	tput setaf 2; tput bold; echo " "

# Quality Control - TRIMMING

for fq in $(ls *.gz)

	do
	fastqc $fq
	THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')
	SL=$(grep "SLIDINGWINDOW" config)
	LE=$(grep "LEADING" config)
	TR=$(grep "TRAILING" config)
	MIN=$(grep "MINLEN" config)

	trimmed_fq=$(echo $fq |cut -d "." -f1 | awk '{print $1".T.fastq"}')
	tput setaf 2; tput bold; echo "    Trimming $fq"

	TrimmomaticSE -threads $THREADS $fq $trimmed_fq $SL $LE $TR $MIN > "$fq.log" 2>&1
	fastqc $trimmed_fq
	gzip $trimmed_fq
done

mkdir fastqc
mv *.zip *.html fastqc/
mkdir trimmed
mv *.T.fastq.gz *.log trimmed/

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Results are saved in fastq and trimmed directories"
	tput setaf 2; tput bold; echo " "
