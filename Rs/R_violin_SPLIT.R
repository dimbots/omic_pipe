# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

library(tidyverse)

setwd("/media/dimbo/10T/data/TAL_LAB/Analysis/Cut_Tag/LaminAC/1889_1891_2147_2149/peak_calling/Overlap_Peaks/PEAKSKO_WITHIN_PEAKSWT/")

readfile = read.csv("Set8KO_Peaks_Within_WT_Peaks.NOCHRXMY.bed", sep = "\t")

# sample size
sample_size = readfile %>% group_by(chr1) %>% summarize(num=n())


###########################################33
# Plot
bla = readfile %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(chr1, "\n", "n=", num)) %>%
  ggplot( aes(x=fct_reorder(myaxis, parse_number(myaxis)), y=X1, fill=chr1)) +
  geom_violin(width=1.4) +
  geom_point() +
  ylim(1,5) +
  ylab(label = "") +
  geom_boxplot(width=0.1, color="white", alpha=0.2) +
  #  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values = c("chartreuse4", "dodgerblue4", "darkslateblue",  "darkslateblue", "orangered3", "orangered3", "orangered4", "orangered4", "red4", "red4",   "chartreuse4", "darkolivegreen4", "darkolivegreen4", "darkolivegreen", "darkolivegreen", "deepskyblue4", "deepskyblue4", "dodgerblue4")) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Number of Set8KO peaks within WT peaks") +
  xlab("")
