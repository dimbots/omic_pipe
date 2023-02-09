#!/bin/bash

# Peak Calling

tput setaf 2; tput bold; echo " "
tput setaf 3; tput bold; echo "    Peak Calling with MACS2"
tput setaf 2; tput bold; echo " "


peak_c=$(grep "PEAK_CALLING" config | cut -d ":" -f2)

case $peak_c in

TRUE)

ln -s mapping/*.bam .
ln -s mapping/*.bam.bai .

while [[ $treatment != "n" ]]

	do

	tput setaf 6; tput bold; echo "Type treatment bam file. else type (n)"
	tput setaf 2; tput bold; echo " "
	read treatment

		if [[ $treatment = "n" ]]

			then
			break

			else

			tput setaf 6; tput bold; echo "Type input bam file. If input file does not exist type (n)"
			tput setaf 2; tput bold; echo " "
			read input
			tput setaf 6; tput bold; echo "Type out file. e.g peaks_sample - else type (n)"
			tput setaf 2; tput bold; echo " "
			read out

			broad=$(grep "BROAD" config | cut -d ":" -f2)
			g=$(grep "GENOME_SP" config | cut -d ":" -f2)

			case $broad in

				TRUE)

				if [[ $input == "n" ]]

					then
					#broad no input
					macs2 callpeak --treatment $treatment --nomodel --broad --format BAM --gsize $g -n $out

					else
					#broad input
					macs2 callpeak --treatment $treatment --control $input --nomodel --broad --format BAM --gsize $g -n $out

				fi

				;;

				FALSE)

				if [[ $input == "n" ]]

					then
					#default no input
					macs2 callpeak --treatment $treatment --nomodel --format BAM --gsize $g -n $out

					else
					#default input
					macs2 callpeak --treatment $treatment --control $input --nomodel --format BAM --gsize $g -n $out
				fi

				;;

				*)

				;;

			esac
		fi

done

# Merge peaks from replicates

while [[ $repA != "n" ]]

	do

	tput setaf 3; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Merge Peaks From Replicates"
	tput setaf 3; tput bold; echo " "

	tput setaf 3; tput bold; echo " "
	tput setaf 6; tput bold; echo "type replicate A peak file. else type (n)"
	tput setaf 2; tput bold; echo " "
	read repA

		if [[ $repA = "n" ]]

			then
			break

			else

			tput setaf 3; tput bold; echo " "
			tput setaf 6; tput bold; echo "type replicate B peak file. else type (n)"
			tput setaf 2; tput bold; echo " "
			read repB

			tput setaf 3; tput bold; echo " "
			tput setaf 6; tput bold; echo "type peaks out file. e.g (Peaks_Merged_samples.tsv) - else type (n)"
			tput setaf 2; tput bold; echo " "
			read out_file

#			tput setaf 3; tput bold; echo " "
#			tput setaf 6; tput bold; echo "set path to blacklist_regions. else type (n)"
#			tput setaf 2; tput bold; echo " "

			bedtools intersect -a $repA -b $repB -wa > $out_file


			#final_peaks="$out_file.bed"

#			tput setaf 3; tput bold; echo " "
#			tput setaf 6; tput bold; echo "set path to blacklist_regions. else type (n)"
#			tput setaf 2; tput bold; echo " "
#			read blacklist

#			bedtools intersect -v -a $out_file -b $blacklist > $final_peaks

		fi
done

mkdir peak_calling
mv *.xls *.gappedPeak *.broadPeak peak_calling
mkdir peak_calling/merged
mv *.tsv peak_calling/merged
rm *.bam *.bai

;;

FALSE)

tput setaf 3; tput bold; echo " "
tput setaf 6; tput bold; echo "    Peak Calling Aborted"
tput setaf 2; tput bold; echo " "

;;

*)

esac
