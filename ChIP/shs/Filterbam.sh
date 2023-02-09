#!/bin/bash

# dublicates

cd mapping
mkdir F_tmp
for bam in $(ls *.bam)

	do
	filtered=$(echo $bam |cut -d "." -f1 | awk '{print $1 ".F.bam"}')
	echo "filter $bam"
	samtools rmdup -s $bam $filtered
	samtools index $filtered
	mv *.F.bam *.F.bam.bai F_tmp	
	done
	rm *.bam *.bai
	mv F_tmp/* .
	rm -r F_tmp

# blacklist regions

	ln -s /media/dimbo/10T/TAL_LAB/Genomes/genes_info/mm10-blacklist.v2.bed .

for bam in $(ls *.bam)


        do
	echo "remove blacklist regions from $bam"
	filtered2=$(echo $bam |cut -d "." -f1 | awk '{print $1 ".bam"}')
	bedtools intersect -wa -a $bam -b mm10-blacklist.v2.bed -v > $filtered2
	samtools index $filtered2
	done

	rm *.F.bam
	rm *.F.bam.bai


# Replace number of uniquely mapped reads in summary.txt files

for bam in $(ls *.bam)

        do

        uniq=$(samtools view $bam -@ 8 | wc -l)
	summary=$(echo $bam | cut -d "." -f1 | awk '{print $1".summary.txt"}')

	sed -i "4 c \\$uniq" $summary

	done

cd ../
