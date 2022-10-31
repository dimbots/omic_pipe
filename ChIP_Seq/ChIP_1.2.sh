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

#############################################################################################################################################

# Mapping

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

find -type l -delete

#############################################################################################################################################

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

#############################################################################################################################################

# QC metrics

tput setaf 2; tput bold; echo " "
tput setaf 3; tput bold; echo "    QC Metrics"
tput setaf 2; tput bold; echo " "

	bam=$(ls *.bam)
	bw=$(ls *.bw)

	multiBigwigSummary bins -b $bw -p $THREADS --smartLabels -o results.npz 
	plotCorrelation \
	-in results.npz \
	--corMethod spearman \
	--removeOutliers \
	--skipZeros \
	--colorMap YlGnBu \
	--plotHeight 11.5 \
	--plotWidth 13 \
	--whatToPlot heatmap \
	--plotNumbers \
	-o heatmap_SpearmanCor.png
	plotCorrelation \
	-in results.npz \
	--corMethod pearson \
	--removeOutliers  \
	--skipZeros --colorMap YlGnBu \
	--plotHeight 11.5 \
	--plotWidth 13 \
	--whatToPlot heatmap \
	--plotNumbers -o heatmap_PearsonCor.png
	plotCorrelation \
	-in results.npz \
	--corMethod pearson \
	--removeOutliers \
	--skipZeros \
	--whatToPlot scatterplot \
	--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
	-o scatterplot_PearsonCorr_bigwigScores.png \
	--outFileCorMatrix PearsonCorr_bigwigScores.tab
	plotCoverage \
	--bamfiles $bam \
	--smartLabels \
	--skipZeros \
	-p $THREADS \
	--verbose \
	-o coverage_plot.png
	plotFingerprint \
	-b $bam \
	--smartLabels \
	-p $THREADS \
	-plot fingerPrint_plot.png

mkdir qc_metrics | mv *.png *.npz *.tab qc_metrics
mkdir mapping | mv *.bam *.bai *.txt mapping 
mkdir normalized | mv *.bw normalized

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Results are saved in mapping and normalized and qc_directories directories"
	tput setaf 2; tput bold; echo " "

#############################################################################################################################################

# Downsampling

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

	echo $uniq_map
	echo $sample
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

rm *.tmp *.txt
mkdir downsampled
mv *.bw downsampled
mv *.downsampled.bam downsampled
mv *.bai downsampled
rm *.bam

tput setaf 2; tput bold; echo " "
tput setaf 3; tput bold; echo "    Downsampling Complete! Results are saved in downsampling directory" 

;;

FALSE)

tput setaf 3; tput bold; echo "    Abort downsampling analysis"
tput setaf 2; tput bold; echo " "

;;

*)

esac

############################################################################################################################################

# PEAK CALLING

peak_c=$(grep "PEAK_CALLING_MACS" config | cut -d ":" -f2)

case $peak_c in

TRUE)

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Peak Calling Analysis with MACS2"
	tput setaf 2; tput bold; echo " "

ln -s downsampled/*.bam .
ln -s downsampled/*.bam.bai .

	while [[ $treatment != "EXIT" ]]

        do

        tput setaf 6; tput bold; echo "Set treatment bam file. To abord peak call analysis type (EXIT)"
        tput setaf 2; tput bold; echo " "
        read treatment

		if [[ $treatment = "EXIT" ]]

			then
                       	break

			else
			tput setaf 6; tput bold; echo "Set input bam file. If there is not input file type (no_input) "
                        tput setaf 2; tput bold; echo " "
			read input

			tput setaf 6; tput bold; echo "Set output file, e.g sample name"
                        tput setaf 2; tput bold; echo " "
                        read out

			g=$(grep "GENOME_SP" config | cut -d ":" -f2)
			broad=$(grep "BROAD" config | cut -d ":" -f2)


			                        case $broad in

		                                TRUE)
						if [[ $input != "no_input" ]]

							then
                	                		# broad (with input)
                        		        	macs2 callpeak --treatment $treatment --control $input --nomodel --broad --format BAM --gsize $g -n $out
							# Filtering (pvalue cutoff = 1.00e-06 && FDR threshold of 0.1%)
							grep -v "#" "${out}_peaks.xls" | sed -e '1,2d' | awk '$7 > 6.34 {print $0}' | awk '$9 > 6.34 {print $0}' > "${out}_peaks.filtered.xls"
							else
							macs2 callpeak --treatment $treatment --nomodel --broad --format BAM --gsize $g -n $out
						fi

		                                ;;

                		                FALSE)
						if [[ $input != "no_input" ]]

                                                        then
	      	                        		# narrow (with input)
			                                macs2 callpeak --treatment $treatment --control $input --nomodel --format BAM --gsize $g -n $out
							# Filtering (pvalue cutoff = 1.00e-06 && FDR threshold of 0.1%)
                                                        grep -v "#" "${out}_peaks.xls" | sed -e '1,2d' | awk '$7 > 6.34 {print $0}' | awk '$9 > 6.34 {print $0}' > "${out}_peaks.filtered.xls"
							else
							macs2 callpeak --treatment $treatment --nomodel --format BAM --gsize $g -n $out
						fi

                		                ;;

        	                	        *)

                        	       		 ;;

			                        esac
                fi

done

mkdir peak_calling
# narrow
mv *.narrowPeak peak_calling
mv *.xls peak_calling
# broad
mv *.broadPeak peak_calling
mv *.gappedPeak peak_calling


;;

FALSE)

tput setaf 2; tput bold; echo " "
tput setaf 3; tput bold; echo "    Peak Calling Analysis with SICER"
tput setaf 2; tput bold; echo " "

ln -s downsampled/*.bam .
ln -s downsampled/*.bam.bai .
mkdir peak_calling

	while [[ $treatment != "EXIT" ]]

        do

        tput setaf 6; tput bold; echo "Set treatment bam file. To abord peak call analysis type (EXIT)"
        tput setaf 2; tput bold; echo " "
        read treatment

		if [[ $treatment = "EXIT" ]]

			then
                       	break

			else
			tput setaf 6; tput bold; echo "Set input bam file. If there is not input file type (no_input) "
                        tput setaf 2; tput bold; echo " "
                        read input

			tput setaf 6; tput bold; echo "Set output file, e.g sample name"
                        tput setaf 2; tput bold; echo " "
                        read out

				if [[ $input != "no_input" ]]

					then
                	                # (with input)
					sicer -t $treatment -c $input -s mm10 -w 200 -g 600 --false_discovery_rate 0.001 --effective_genome_fraction 0.88 -o $out -cpu 7
					mv $out peak_calling
					else
					sicer -t $treatment -s mm10 -w 200 -g 600 --false_discovery_rate 0.001 --effective_genome_fraction 0.88 -o $out -cpu 7
					mv $out peak_calling
				fi

                fi

done

;;

*)

esac

find -type l -delete


############################################################################################################################################

# Compute matrix and plot heatmaps

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')
PLOT=$(grep "PLOT" config | cut -d ":" -f2)

if [ $PLOT == "FALSE" ]

	then
	tput setaf 3; tput bold; echo "     Heatmap will not be created"

else

ln -s downsampled/*.bw .
ln -s peak_calling/* .

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

			tput setaf 6; tput bold; echo "Set peak file for sample $bw file."
                        tput setaf 2; tput bold; echo " "
                        read peak_file

			computeMatrix reference-point \
			--referencePoint center \
			-b $tss \
			-a $tes \
			-R $peak_file \
			-S $bw \
			--skipZeros \
			-o $out_matrix \
			-p $THREADS

			out_plot=$(echo $bw |cut -d "." -f1 | awk '{print $1 ".Heatmap.png"}')
			plotHeatmap -m $out_matrix -out $out_plot --interpolationMethod nearest -z "Peaks" -x " " --boxAroundHeatmaps no  --colorList white,whitesmoke,seashell,mistyrose,salmon,r  --missingDataColor 1


		else
			tput setaf 3; tput bold; echo "     Generate Scale Regions Heatmap"

			tput setaf 6; tput bold; echo "Set peak file for sample $bw file."
                        tput setaf 2; tput bold; echo " "
                        read peak_file

			computeMatrix scale-regions \
			-S $bw \
			-R $peak_file \
			--beforeRegionStartLength $tss \
			--regionBodyLength $body \
			--afterRegionStartLength $tes \
			--smartLabels --startLabel "TSS" \
			--endLabel "TES" \
			--skipZeros -o $out_matrix \
			--outFileSortedRegions $out_regions \
			-p $THREADS

			out_plot=$(echo $bw |cut -d "." -f1 | awk '{print $1 ".Heatmap.png"}')
			plotHeatmap -m $out_matrix -out $out_plot --interpolationMethod nearest -z "Peaks" -x " " --boxAroundHeatmaps no  --colorList white,whitesmoke,seashell,mistyrose,salmon,r --missingDataColor 1

		fi
	done

fi

mkdir heatmaps
mv *.png *.out_matrix *.out_regions heatmaps

find -type l -delete

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Results are saved in heatmaps directory"
	tput setaf 2; tput bold; echo " "


#############################################################################################################################################




