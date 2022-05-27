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
