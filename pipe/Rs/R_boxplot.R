# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

# create a dataset
data <- data.frame(
  name=c( rep("A",500), rep("B",500), rep("B",500), rep("C",20), rep('D', 100)  ),
  value=c( rnorm(500, 10, 5), rnorm(500, 13, 1), rnorm(500, 18, 1), rnorm(20, 25, 4), rnorm(100, 12, 1) )
)

data1 = read.csv("/media/dimbo/10T/data/TAL_LAB/Analysis/Cut_Tag/LaminAC/1889_1891_2147_2149/peak_calling/NEW_PEAKS/peak_length/peak_length_merged_logged.csv", sep = "\t")


# Plot
pdf("Lad_Peak_Length.pdf", width = 8, height = 8)
data1 %>%
  ggplot( aes(x=Condition, y=value, fill=Condition)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#339c4f", "#9c3333")) +

  ylim(0,40) +

  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  theme_bw() +

xlab("") +
  ylab("LAD Peak Length in (kb)")
dev.off()
