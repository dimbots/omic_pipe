# Keep only top 25% based on p-value
sort -V -k5 KO_CTCF_summits.bed | tail -n 5226 > KO_top25_summits.bed

# Calculate median from peak length from original file
#sort -n peak_length | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'

grep -v "#" ../peak_calling/KO_CTCF_peaks.xls | awk '{print $4}' | awk 'NR > 2 { print }' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'

# Extract bed intervals +- from 316 nt around the peak summit of peaks
awk '{print $1 "\t" ($2 - 158) "\t" $3 "\t" $4 "\t" $5}' KO_top25_summits.bed > KO_CTCF.tmp

awk '{print $1 "\t" $2 "\t" ($3 + 158) "\t" $4 "\t" $5}' KO_CTCF.tmp > KO_CTCF_Intervals.bed

bedtools getfasta -fi mm10.fa -bed KO_CTCF_Intervals.bed -fo KO_CTCF_Intervals.fa
scrambleFasta.pl KO_CTCF_Intervals.fa > control_scrambled.fa



##########################################3
# Calculate median
grep -v "#" ../peak_calling/KO_CTCF_peaks.xls | awk '{print $4}' | awk 'NR > 2 { print }' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'

# Extract top 25% percent of peaks
sort -V -k5 KO_CTCF_summits.bed | tail -n 5226 > KO_top25_summits.bed

# Run Homer
findMotifsGenome.pl KO_top25_summits.bed mm10.fa peakAnalysis -size -158,158




# for overlapped
# extract 25%
sort -V -k7 overlapped_peaks.CTCF.bed | tail -n 846 > overlapped_top25.bed

# calculate median 
awk '{print $10}' overlapped_top25.bed | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'

findMotifsGenome.pl overlapped_top25.bed mm10.fa peakAnalysis -size -65,65
