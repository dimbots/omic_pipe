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
		wait

mkdir fastqc
mv *.zip *.html fastqc/
mkdir trimmed
mv *.T.fastq.gz *.log trimmed/

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Results are saved in fastq and trimmed directories"
	tput setaf 2; tput bold; echo " "
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
#!/bin/bash

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')

# QC metrics
ln -s mapping/*.bam .
ln -s mapping/*.bai .
ln -s normalized/*.bw .

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

	tput setaf 2; tput bold; echo " "
	tput setaf 3; tput bold; echo "    Results are saved in mapping and normalized and qc_directories directories"
	tput setaf 2; tput bold; echo " "
rm *.bam *.bai *.bw
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
# Annotation

THREADS=$(grep "THREADS:" config | awk '{print substr($0,length,1)}')

tput setaf 3; tput bold; echo " "
tput setaf 3; tput bold; echo "    Annotating Regions in The Genome with Homer"
tput setaf 2; tput bold; echo " "

annotation=$(grep "ANNOTATION" config | cut -d ":" -f2)

case $annotation in

TRUE)

fasta=$(grep "PATH TO FASTA" config | cut -d ":" -f2)
gtf=$(grep "PATH TO GTF" config | cut -d ":" -f2)

pwd=$(pwd)

if [ -d "$pwd/peak_calling/merged" ]

	then

		ln -s peak_calling/merged/* .

		for peaks in $(ls *.tsv)

			do

			sample=$(echo $peaks | cut -d "_" -f2)
			out=$(echo "annotation_$sample")

			annotatePeaks.pl $peaks $fasta -gtf $gtf > $out

			done

			tput setaf 2; tput bold; echo " "
			tput setaf 3; tput bold; echo "    Annotation Complete!"
			tput setaf 2; tput bold; echo " "

	else

		ln -s peak_calling/*.broadPeak .

		for peaks in $(ls *.broadPeak)

                        do

                        sample=$(echo $peaks | cut -d "_" -f2)
                        out=$(echo "annotation_$sample")

                        annotatePeaks.pl $peaks $fasta -gtf $gtf > $out

                        done

                        tput setaf 2; tput bold; echo " "
                        tput setaf 3; tput bold; echo "    Annotation Complete!"
                        tput setaf 2; tput bold; echo " "

fi

;;

FALSE)

;;

*)

esac

mkdir annotated
mv annotation* annotated/
rm *.broadPeak
