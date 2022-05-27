#!/bin/bash

# Mapping

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')
ln -s trimmed/*.gz .

for trimmed_fastq in $(ls *.T.fastq.gz)

	do
	tool=$(grep "HISAT2" config | cut -d ":" -f2)
	genome=$(grep "PATH TO GENOME" config | cut -d ":" -f2)

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Processing $trimmed_fastq"
	tput setaf 2; tput bold; echo " "
	summary=$(echo $trimmed_fastq | cut -d "." -f1 | awk '{print $1".summary.txt"}')
	sam=$(echo $trimmed_fastq | cut -d "." -f1 | awk '{print $1".sam"}')
	bam=$(echo $trimmed_fastq | cut -d "." -f1 | awk '{print $1".bam"}')

		case $tool in

		TRUE)

		tput setaf 2; tput bold; echo "    Aligning sequencing reads with Hisat2"

		hisat2 --threads $THREADS --no-spliced-alignment --summary-file $summary -x $genome -U $trimmed_fastq | \
		samtools view -@ $THREADS -b -q 30 -F 4 | \
		samtools sort -@ $THREADS -o $bam 
		samtools index $bam

		;;

		FALSE)

		tput setaf 2; tput bold; echo "    Aligning sequencing reads with Bowtie2"

		bowtie2 -p $THREADS --very-sensitive -x $genome -U $trimmed_fastq -S $sam >> "log_$sam" 2>&1
		samtools view -@ $THREADS -h -S -b -q 30 -o $bam $sam
		samtools sort -@ $THREADS -o $bam $bam
		samtools index $bam
		rm $sam

		;;

		*)

		esac
done

mkdir mapping 
mv *.bam *.bai *.txt mapping
rm *.T.fastq.gz
