#!/bin/bash

# Compute matrix and plot heatmaps
THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')
PLOT=$(grep "PLOT" config | cut -d ":" -f2)

if [ $PLOT == "FALSE" ]

	then
	tput setaf 3; tput bold; echo "     Heatmap will not be created"

else

ln -s normalized/*.bw .

	for bw in $(ls *.bw)

		do

		THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')
		tss=$(grep "TSS" config | cut -d ":" -f2)
 		tes=$(grep "TES" config | cut -d ":" -f2)
 		body=$(grep "BODY" config | cut -d ":" -f2)
 		gtf=$(grep "PATH TO GTF" config | cut -d ":" -f2)
		method=$(grep "REFERENCE_POINT" config | cut -d ":" -f2)

		out_matrix=$(echo $bw |cut -d "." -f1 | awk '{print $1 ".out_matrix"}')
		out_regions=$(echo $bw |cut -d "." -f1 | awk '{print $1 ".out_regions"}')

		tput setaf 2; tput bold; echo " "
		tput setaf 3; tput bold; echo "    Compute Matrix of $bw sample"
		tput setaf 2; tput bold; echo " "


		if [ $method == "TRUE" ]

			then
			tput setaf 3; tput bold; echo "     Generate Reference Point Heatmap"

			computeMatrix reference-point \
			--referencePoint TSS \
			-b $tss \
			-a $tes \
			-R $gtf \
			-S $bw \
			--skipZeros \
			-o $out_matrix \
			-p $THREADS

			out_plot=$(echo $bw |cut -d "." -f1 | awk '{print $1 ".Heatmap.png"}')
			plotHeatmap -m $out_matrix -out $out_plot --colorMap GnBu -x " " --missingDataColor 1


		else
			tput setaf 3; tput bold; echo "     Generate Scale Regions Heatmap"
			computeMatrix scale-regions \
			-S $bw \
			-R $gtf \
			--beforeRegionStartLength $tss \
			--regionBodyLength $body \
			--afterRegionStartLength $tes \
			--smartLabels --startLabel "TSS" \
			--endLabel "TES" \
			--skipZeros -o $out_matrix \
			--outFileSortedRegions $out_regions \
			-p $THREADS

			out_plot=$(echo $bw |cut -d "." -f1 | awk '{print $1 ".Heatmap.png"}')
			plotHeatmap -m $out_matrix -out $out_plot --colorMap GnBu -x " " --missingDataColor 1

		fi
	done

fi
mkdir heatmaps
mv *.png *.out_matrix *.out_regions heatmaps
rm *.bw

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Results are saved in heatmaps directory"
	tput setaf 2; tput bold; echo " "
