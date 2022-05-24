#!/bin/bash

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')
ln -s mapping/*.bam .
ln -s mapping/*.bai .
# Normalization

for bam in $(ls *.bam)

	do
		normalize=$(grep "RPKM" config | cut -d ":" -f2)

		case $normalize in

			TRUE)

			tput setaf 2; tput bold; echo " "
			tput setaf 3; tput bold; echo "    Normalize bam files with RPKM"
			tput setaf 2; tput bold; echo " "
			tput setaf 2; tput bold; echo " "
			tput setaf 3; tput bold; echo "    Processing $bam"
			tput setaf 2; tput bold; echo " "
			bw=$(echo $bam |cut -d "." -f1 | awk '{print $1 ".bw"}')

			bamCoverage --bam $bam -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPKM -p $THREADS --extendReads 200

			;;

			FALSE)

			tput setaf 2; tput bold; echo " "
 			tput setaf 3; tput bold; echo "    Normalize bam files with BPM"
			tput setaf 2; tput bold; echo " "
			tput setaf 2; tput bold; echo " "
			tput setaf 3; tput bold; echo "    Processing $bam"
			tput setaf 2; tput bold; echo " "
			bw=$(echo $bam |cut -d "." -f1 | awk '{print $1 ".bw"}')

			bamCoverage --bam $bam -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing BPM -p $THREADS --extendReads 200

			;;

			*)

			;;
		esac
done

mkdir normalized
mv *.bw normalized/
rm *.bam *.bai
