# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

data1 = read.csv("/media/dimbo/10T/data/TAL_LAB/Data/Cut_Tag/LaminAC/1889_1891_2147_2149/peak_calling/peak_length/peak_length_merged_logged.csv", sep = "\t")

# sample size
sample_size = data1 %>% group_by(Condition) %>% summarize(num=n())

# Plot
pdf("peak_length_Lamin.pdf")
p = data1 %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(Condition, "\n", "n=", num)) %>%
  
  ggplot( aes(x=myaxis, y=value, fill=Condition)) +
  
  geom_violin(width=1) +
  
  geom_boxplot(width=0.06, color="black",  outlier.shape = NA) +
  
  scale_fill_manual(values=c("#339c4f", "#9c3333")) +
  
  ylim(0,40) +
  
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  
  xlab("") +
  ylab("LAD Peak Length in (kb)")
p + theme_bw()
dev.off()


