library("VennDiagram")

pdf("Venn_CTCF.pdf", width = 10, height = 10)
grid.newpage()                                       
draw.pairwise.venn(area1 = 5838,                        
                   area2 = 20905,
                   cross.area = 3387,
                   fill = c("coral1", "lightgreen"),
                   lty = "blank",
                   cex = 1.2,
                   cat.cex = 1,
                   cat.default.pos = "outer",
                   cat.pos = c(2, 2),
                   cat.dist = c(0.009, 0.009),
                   inverted = TRUE,
                   category = c("WT CTCF", "Set8KO CTCF"))
dev.off()
