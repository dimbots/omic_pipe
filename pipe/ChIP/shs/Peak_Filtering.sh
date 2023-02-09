# Peak calling (Default)
macs2 callpeak --treatment Set8KO_SMC3.downsampled.bam -c Set8KO_Input.downsampled.bam --format BAM --gsize mm -n Set8KO_SMC3

# Filter (pvalue cutoff = 1.00e-06 && FDR threshold of 0.1%)
grep -v "#" Set8KO_SMC3_peaks.xls | sed -e '1,2d' | awk '$7 > 6.34 {print $0}' | awk '$9 > 6.34 {print $0}' | cut -f1,2,3 > Set8KO_SMC3_peaks.Filtered.bed


sort-bed Set8KO_SMC3_peaks.Filtered.bed > regions_Set8KO.bed

bedToBigBed Set8KO_SMC3_peaks.Filtered.bed mm10_chrom_sizes regions_Set8KO.bb




##########################33


# Filtering (pvalue cutoff = 1.00e-06 && FDR threshold of 0.1%)
grep -v "#" Set8KO_SMC3_peaks.xls | sed -e '1,2d' | awk '$7 > 6.34 {print $0}' | awk '$9 > 6.34 {print $0}' > Set8KO_peaks_F.xls

# Calculate the 95th percentile of peak length
awk '$11 = $3 - $2 {print $11}' Set8KO_peaks_F.xls | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}'


# Remove the 95th percentile of peaks with length 
awk '$4 < 2000 {print $0}' Set8KO_peaks_F.xls > Set8KO_peaks_Filtered.xls


# Rename
cp Set8KO_peaks_Filtered.xls Set8KO_peaks_Filtered.bed

#############################
