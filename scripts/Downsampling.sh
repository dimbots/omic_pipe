#!/bin/bash

# Downsampling

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')
tput setaf 3; tput bold; echo " "
tput setaf 6; tput bold; echo "    Downsampling bam files"
tput setaf 2; tput bold; echo " "

downsampling=$(grep "DOWNSAMPLING" config | cut -d ":" -f2)

case $downsampling in

TRUE)

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')

ln -s mapping/*.bam .
ln -s mapping/*.txt .

for sum in $(ls *.txt)

	do

	uniq_map=$(cat $sum  | head -n4 | tail -n1 | cut -d " " -f5)
	sample=$(echo $sum |cut -d "." -f1)


	tput setaf 3; tput bold; echo "    $sample sample has $uniq_map uniquely mapped reads aligned"

#	echo $uniq_map
#	echo $sample
done

tput setaf 2; tput bold; echo " "
tput setaf 6; tput bold; echo "Type number of reads, e.g 39707152" 
tput setaf 2; tput bold; echo " "
tput setaf 6; tput bold; echo "(This number will be used to downsample all bam files!)"
read number
tput setaf 2; tput bold; echo " "


for bam in $(ls *.bam)

	do

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold;echo "    Downsampling $bam"
	tput setaf 2; tput bold; echo " "

	input="$bam"
	header="$bam.headers.tmp"

	samtools view -H $input > $header

	shuffled="$bam.shuffled.tmp"

	samtools view -@ $THREADS $input | shuf | head -n $number > $shuffled

	unsorted="$bam.downsampled.tmp"

	cat $header $shuffled > $unsorted

	sorted="$bam.downsampled.bam"

	samtools sort -@ $THREADS $unsorted -o $sorted

	samtools index -@ $THREADS $sorted

	bw=$(echo $bam |cut -d "." -f1 | awk '{print $1".bw"}')
	bamCoverage -p $THREADS -b $sorted -o $bw
done

rm *.tmp *.bai *.bam *.txt
mkdir downsampled
mv *.bw downsampled

tput setaf 2; tput bold; echo " "
tput setaf 3; tput bold; echo "    Downsampling Complete! Results are saved in downsampling directory" 

;;

FALSE)

tput setaf 3; tput bold; echo "    Abort downsampling analysis"
tput setaf 2; tput bold; echo " "

;;

*)

esac
