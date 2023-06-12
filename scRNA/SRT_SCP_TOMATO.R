library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)
library(Signac)
library(SeuratWrappers)
library(monocle3)
library(Nebulosa)
library(BiocFileCache)
library(Matrix)
library(extrafont)
library(remotes)
library(SCP)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
#extrafont::font_import()
set.seed(7)

#----------------------------------------------------------------------------------------------------------------------------------------------------

# CREATE OBJECT Lgr5Cre_SETDB1KO

# Load dataset A (Object A)
Lgr5Cre_Setdb1KO_A.data=Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/T008/Cellranger/run_count_22L002708/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_A = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_A.data, project = "Setdb1KO A", 
                                        min.cells = 3,
                                        min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_Setdb1KO_B.data <- Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/T008/Cellranger/run_count_22L002712/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_B = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_B.data, project = "Setdb1KO B", 
                                        min.cells = 3, 
                                        min.features = 200)
# Merge objects
Lgr5Cre_Setdb1KO=merge(Lgr5Cre_Setdb1KO_A, y=Lgr5Cre_Setdb1KO_B, add.cell.ids=c("Rep_A","Rep_B"), project="Lgr5Cre_Setdb1KO")

# CREATE OBJECT LGR5Cre_WT

# Load dataset A (Object A)
Lgr5Cre_WT_A.data=Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/T008/Cellranger/run_count_22L002700/outs/filtered_feature_bc_matrix/")
Lgr5Cre_WT_A = CreateSeuratObject(counts = Lgr5Cre_WT_A.data, project = "Lgr5Cre A", 
                                  min.cells = 3, 
                                  min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_WT_B.data <- Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/T008/Cellranger/run_count_22L002704/outs/filtered_feature_bc_matrix/")
Lgr5Cre_WT_B = CreateSeuratObject(counts = Lgr5Cre_WT_B.data, project = "Lgr5Cre B", 
                                  min.cells = 3, 
                                  min.features = 200)

# Merge objects
Lgr5Cre_WT=merge(Lgr5Cre_WT_A, y=Lgr5Cre_WT_B, add.cell.ids=c("Rep_A","Rep_B"), project="Lgr5Cre_WT")

#----------------------------------------------------------------------------------------------------------------------------------------------------

#FIND DOUBLETS IN EACH CONDITION SEPARATELY

Lgr5Cre_Setdb1KO <- NormalizeData(Lgr5Cre_Setdb1KO, normalization.method = "LogNormalize", scale.factor = 10000)

Lgr5Cre_Setdb1KO <- FindVariableFeatures(Lgr5Cre_Setdb1KO, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Lgr5Cre_Setdb1KO)
Lgr5Cre_Setdb1KO <- ScaleData(Lgr5Cre_Setdb1KO, features = all.genes)

Lgr5Cre_Setdb1KO <- RunPCA(Lgr5Cre_Setdb1KO, features = VariableFeatures(object = Lgr5Cre_Setdb1KO))

Lgr5Cre_Setdb1KO <- RunUMAP(Lgr5Cre_Setdb1KO, dims = 1:10)

nExp <- round(ncol(Lgr5Cre_Setdb1KO) * 0.069)
Lgr5Cre_Setdb1KO.filt <- doubletFinder_v3(Lgr5Cre_Setdb1KO, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(Lgr5Cre_Setdb1KO.filt@meta.data)[grepl("DF.classification", colnames(Lgr5Cre_Setdb1KO.filt@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(Lgr5Cre_Setdb1KO.filt,group.by = DF.name) + NoAxes())

VlnPlot(Lgr5Cre_Setdb1KO.filt, features = "nFeature_RNA",  group.by = DF.name, pt.size = 0)

Lgr5Cre_Setdb1KO.filt = Lgr5Cre_Setdb1KO.filt[, Lgr5Cre_Setdb1KO.filt@meta.data[, DF.name] == "Singlet"]

#----------------------------------------------------------------------------------------------------------------------------------------------------

Lgr5Cre_WT <- NormalizeData(Lgr5Cre_WT, normalization.method = "LogNormalize", scale.factor = 10000)

Lgr5Cre_WT <- FindVariableFeatures(Lgr5Cre_WT, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Lgr5Cre_WT)
Lgr5Cre_WT <- ScaleData(Lgr5Cre_WT, features = all.genes)

Lgr5Cre_WT <- RunPCA(Lgr5Cre_WT, features = VariableFeatures(object = Lgr5Cre_WT))

Lgr5Cre_WT <- RunUMAP(Lgr5Cre_WT, dims = 1:10)

nExp <- round(ncol(Lgr5Cre_WT) * 0.057)
Lgr5Cre_WT.filt <- doubletFinder_v3(Lgr5Cre_WT, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(Lgr5Cre_WT.filt@meta.data)[grepl("DF.classification", colnames(Lgr5Cre_WT.filt@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(Lgr5Cre_WT.filt, group.by = DF.name) + NoAxes())

VlnPlot(Lgr5Cre_WT.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0)

Lgr5Cre_WT.filt = Lgr5Cre_WT.filt[, Lgr5Cre_WT.filt@meta.data[, DF.name] == "Singlet"]

#----------------------------------------------------------------------------------------------------------------------------------------------------

# Merge Setdb1KO with WT to new object Lgr5Cre_MERGED
Lgr5Cre_MERGED=merge(Lgr5Cre_Setdb1KO, y=Lgr5Cre_WT, 
                     add.cell.ids=c("Setdb1KO","Lgr5Cre"), 
                     project="Lgr5Cre_MERGED")

# Change replicates names within object. e.g(Setdb1KO_repA & Setdb1KO_repB -> Setdb1KO)
Lgr5Cre_MERGED$orig.ident=plyr::mapvalues(x=Lgr5Cre_MERGED$orig.ident, from = c("Setdb1KO A", "Setdb1KO B", "Lgr5Cre A", "Lgr5Cre B"), 
                                          to = c("Setdb1KO", "Setdb1KO", "Lgr5Cre", "Lgr5Cre"))

#----------------------------------------------------------------------------------------------------------------------------------------------------

# QC and selecting cells
Lgr5Cre_MERGED[["percent.mt"]] <- PercentageFeatureSet(Lgr5Cre_MERGED, pattern = "^mt-")

# Violin Plot
VlnPlot(Lgr5Cre_MERGED, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0,
        ncol = 3) 

# Violin Plot with merged replicates
plot1 = VlnPlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                group.by = "orig.ident", ncol = 3,
                pt.size = 0) 

# Feature Scatter Plot
plot1 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt",pt.size = 0.5)
plot2 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA", pt.size = 0.5)

plot1 +   theme(text=element_text(family="FreeSerif")) + plot2  +  theme(text=element_text(family="FreeSerif"))

# Feature Scatter Plot
plot1 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt", 
                        group.by = "orig.ident", pt.size = 0.5)
plot2 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA", 
                        group.by = "orig.ident", pt.size = 0.5)

plot1 + theme(text=element_text(family="FreeSerif")) + plot2 + theme(text=element_text(family="FreeSerif"))

# Filtering cells 
Lgr5Cre_MERGED <- subset(Lgr5Cre_MERGED, 
                         subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5 & nCount_RNA < 18000)

# Normalize the data
Lgr5Cre_MERGED <- NormalizeData(Lgr5Cre_MERGED, 
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000)

# Identification of highly variable features
Lgr5Cre_MERGED <- FindVariableFeatures(Lgr5Cre_MERGED, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Lgr5Cre_MERGED), 10)
top50 <- head(VariableFeatures(Lgr5Cre_MERGED), 50)
# Plot variable features
plot1 <- VariableFeaturePlot(Lgr5Cre_MERGED)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
plot2

# Scaling the data
all.genes <- rownames(Lgr5Cre_MERGED)
Lgr5Cre_MERGED <- ScaleData(Lgr5Cre_MERGED, features = all.genes)

# Performing linear dimensional reduction
Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, 
                         features = VariableFeatures(object = Lgr5Cre_MERGED))
print(Lgr5Cre_MERGED[["pca"]], dims = 1:5, nfeatures = 5)

ElbowPlot(Lgr5Cre_MERGED)
# VizDim Plot
#VizDimLoadings(Lgr5Cre_MERGED, dims = 1:2, reduction = "pca")
# Dim Heatmap
#DimHeatmap(Lgr5Cre_MERGED, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Lgr5Cre_MERGED, dims = 1:10, cells = 500, balanced = TRUE)

# Determine the dimensionality of the dataset
Lgr5Cre_MERGED <- JackStraw(Lgr5Cre_MERGED, num.replicate = 100)
Lgr5Cre_MERGED <- ScoreJackStraw(Lgr5Cre_MERGED, dims = 1:20)

# Cluster the Cells
Lgr5Cre_MERGED <- FindNeighbors(Lgr5Cre_MERGED, dims = 1:10)
Lgr5Cre_MERGED <- FindClusters(Lgr5Cre_MERGED, resolution = 0.5)

# Run non linear dimensional reduction
reticulate::py_install(packages = 'umap-learn')
Lgr5Cre_MERGED <- RunUMAP(Lgr5Cre_MERGED, dims = 1:9)

# Umap plot merged
DimPlot(Lgr5Cre_MERGED, reduction = "umap", 
        group.by = "orig.ident", 
        pt.size = 0.01) + theme(text=element_text(family="FreeSerif"))

# Split objects 
SplitObject(Lgr5Cre_MERGED, split.by = "ident")
n_cells=(FetchData(Lgr5Cre_MERGED, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

# Umap plot merged replicates per condition
DimPlot(Lgr5Cre_MERGED,label=TRUE, 
        split.by="orig.ident", 
        pt.size = 0.01) + NoLegend() + theme(text=element_text(family="FreeSerif"))

#----------------------------------------------------------------------------------------------------------------------------------------------------

# Extract metadata
md = Lgr5Cre_MERGED@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "seurat_clusters")]

cells_per_cluster = md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")
write.table(cells_per_cluster,file = "umap1_cells_per_cluster.tsv", row.names = TRUE, sep = "\t")

# Extract and save all genes
write.table(all.genes,
            file = "all_genes.tsv", 
            sep = "\t",
            row.names = TRUE)

# Rename umap based on the dotplot annoation
Lgr5Cre_MERGED_Renamed=RenameIdents(Lgr5Cre_MERGED,  `0` = "Stem I", `1` = "Stem II", `2` = "Goblet I", `3` = "Goblet III", 
                                    `4` = "Paneth II", `5` = "Progenitor/TA", `6` = "Ent.Mature", `7` = "Paneth I", 
                                    `8` = "Goblet II", `9` = "Goblet IV", `10` = "Stem III", `11` = "Ent.Immature", 
                                    `12` = "Tuft", `13` = "Enteroendocrine", `14` = "Ent.Mature")

# Set order to clusters
levels(Lgr5Cre_MERGED_Renamed) = c("Stem I", "Stem II", "Stem III","Progenitor/TA",
                                   "Ent.Immature", "Ent.Mature", "Goblet I",
                                   "Goblet II", "Goblet III", "Goblet IV","Paneth I",
                                   "Paneth II", "Tuft", "Enteroendocrine")

#----------------------------------------------------------------------------------------------------------------------------------------------------

# DotPlots per Cluster

# STEM/TA/PROGENITOR
DotPlot(Lgr5Cre_MERGED_Renamed,features = c("Ascl2","Axin2","Bex1","Fzd2","Fzd7","Gkn3",
                                            "Hes1","Igfbp4","Lgr5","Lrp5","Lrp6","Notch1","Olfm4","Prom1","Slc12a2","Smo","Yap1"   
                                            ,         "Stmn1", "Tubb5"       ,     "Aurka","Birc5","Ccna2","Ccnb1","Ccnb2",
                                            "Cdc20","Cdc25c","Cdkn2d","Cdkn3","Cenpa","Cps1","Kif22","Kif23","Melk","Nek2","Plk1",
                                            "Rbp7","Sapcd2","Tacc3","Ube2c"  ), 
        cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  + 
  theme(text=element_text(family = "FreeSerif" )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))
# Entrocyte/GOBLET
DotPlot(Lgr5Cre_MERGED_Renamed,features = c("Gsdmc4","Krt8","Prss32","Reg3g"          , 
                                            "Aldob","Alpi","Apoa1","Apoa4","Clec2h","Dpep1","Elf3","Fabp1","Fabp6","Fam151a","Hnf4a","Hnf4aos","Hnf4g",
                                            "Lct","Mep1a","Muc3","Naaladl1","Neu1","Nudt4","Phgr1","Plb1","Pmp22","Prss30","Sis","Slc34a2",
                                            "Slc51a","Slc51b","Tmigd1","Xpnpep2"        ,   
                                            "Agr2","Atoh1","Ccl6","Ccl9","Clca1","Clca3a1","Fcgbp","Foxa1",
                                            "Klf4","Klk1","Lrrc26","Muc2","Spdef","Spink4","Tff3","Tpsg1","Tspan13","Txndc5","Zg16"), 
        cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" ))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))

# Paneth/Tuft/Enteroendo
DotPlot(Lgr5Cre_MERGED_Renamed,features = c("Ang4","Defa17","Defa21","Defa22","Defa24","Defa3",
                                            "Dll4","Gfi1","Gm14851","Lyz1","Mmp7","Mptx2","Sox9"  , 
                                            "Adh1","Aldh2","Alox5ap","Avil","Cd24a","Dclk1","Fyb","Gfi1b","Gng13","Hck",
                                            "Il13ra1","Il25","Kctd12","Lrmp","Ltc4s","Ly6g6f","Ptpn6","Ptprc","Rgs13",
                                            "Sh2d6","Skap2","Tmem176a","Trpm5"       ,          
                                            "Arx","Bex2","Cck","Chga","Chgb","Cpe","Fam183b","Fev","Foxa2",
                                            "Gcg","Gch1","Gck","Gfra3","Hmgn3","Insm1","Isl1","Marcksl1","Neurod1",
                                            "Neurod2","Neurog3","Nkx2-2","Pax6","Pcsk1n","Pyy","Sst","Tac1","Tph1","Vwa5b2" ),
        cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))

# Lgr5
DotPlot(Lgr5Cre_MERGED_Renamed,features = "Lgr5",   cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))
# GFP
DotPlot(Lgr5Cre_MERGED_Renamed,features = "GFP",   cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# SCP PLOTS/ANALYSIS
# 0,1 10 stem 
# 5 TA/Progenitor
# 6  11(imm) 14 ent.mature
# 3 8 9 goblet
# 4 7 pan
# 12 tuft
# 13 entro
# 2, 

# Rename clusters
Lgr5Cre_MERGED_Renamed$seurat_clusters <- factor(x = Lgr5Cre_MERGED_Renamed$seurat_clusters,
                                                 levels = c("0","1", "10", 
                                                            "5", "11", "14", "6", "2", 
                                                            "8", "3", "9", "7", "4", 
                                                            "12", "13"),
                                                 labels = c( "Stem I", "Stem II", "Stem III",  
                                                             "Progenitor/TA", "Ent.Immature",  "Ent.Mature", "Ent.Mature","Goblet I","Goblet II",  
                                                             "Goblet III", "Goblet IV", "Paneth I", "Paneth II", "Tuft", 
                                                             "Enteroendocrine"))

colnames(Lgr5Cre_MERGED_Renamed@meta.data)[colnames(Lgr5Cre_MERGED_Renamed@meta.data) == "seurat_clusters"] <- "CellType"

# All cells
CellDimPlot(
  srt = Lgr5Cre_MERGED_Renamed, group.by = c("CellType"), 
  reduction = "UMAP", theme_use = "theme_blank", show_stat = TRUE,
  palette = "Spectral", label = TRUE,
) + theme(text=element_text(family = "FreeSerif" ))

#----------------------------------------------------------------------------------------------------------------------------------------------------

split_obj = SplitObject(Lgr5Cre_MERGED_Renamed, split.by = "orig.ident")
WT_1st = split_obj[["Lgr5Cre"]]
KO_1st = split_obj[["Setdb1KO"]]

# MERGED UMAP
CellDimPlot(
  srt = Lgr5Cre_MERGED_Renamed, group.by = c("orig.ident"), 
  reduction = "UMAP", theme_use = "theme_blank", show_stat = TRUE, palcolor = c("green3", "red2"),
  # palette = "Earth", 
) + theme(text=element_text(family = "FreeSerif" ))



CellStatPlot(Lgr5Cre_MERGED_Renamed, stat.by = "orig.ident", group.by = "CellType", 
             stat_type = "percent", position = "dodge", 
             label = TRUE, 
             palcolor = c("red3", "green3")) + theme(text=element_text(family = "FreeSerif" ))




#####################################################################################################
#####################################################################################################
#####################################################################################################
# Grouped BARPLOTS

# Create data
data <- data.frame(values = c(3,32,	 23,5,  	2,6,	2,11,	 5,3,  	16,1,	  2,18,	   3,6,	10,8,  	5,4,  	9,3,	19,1,	1,4, 1,2
                              
                              #      32	5	6	11	3	0	18	6	8	4	3	0	4
                              
                              
                              
),  # Create example data
group = rep(c("Stem I","Stem II","Stem III", "Progenitor/TA",
              "Ent.Immature", "Ent.Mature", "Goblet I", "Goblet II", "Goblet III", "Goblet IV",
              "Paneth I", "Paneth II", "Tuft", "Enteroendocrine"),
            each = 2),
subgroup = LETTERS[1:2])
data                 

data_base <- reshape(data,                        # Modify data for Base R barplot
                     idvar = "subgroup",
                     timevar = "group",
                     direction = "wide")
row.names(data_base) <- data_base$subgroup
data_base <- data_base[ , 2:ncol(data_base)]
colnames(data_base) <- c("Stem I","Stem II","Stem III", "Progenitor/TA",
                         "Ent.Immature", "Ent.Mature", "Goblet I", "Goblet II", "Goblet III", "Goblet IV",
                         "Paneth I", "Paneth II", "Tuft", "Enteroendocrine")
data_base <- as.matrix(data_base)
data_base        

level_order <- c("Stem I","Goblet I","Progenitor/TA", "Goblet II", "Stem III",
                 "Tuft", "Enteroendocrine", "Goblet IV", "Ent.Immature",  "Paneth I", "Goblet III",
                 "Ent.Mature","Paneth II", "Stem II") 


library(dplyr)

data %>% 
  mutate(subgroup = recode(subgroup, `A` = "Lgr5Cre", `B` = "Setdb1KO"))


library("ggplot2")
ggplot(data,                                    # Grouped barplot using ggplot2
       aes(x = group,
           y = values, 
           fill = subgroup)) +
  geom_bar(stat = "identity",
           position = "dodge") + 
  scale_x_discrete(limits = level_order) + 
  #theme(panel.background = element_blank()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1)) + xlab("") + ylab("% Cells Normalized") +
  #scale_fill_discrete(labels=c('LgrCre', 'Setdb1KO')) + 
  scale_fill_manual(values = c("red3","green3"), labels =c('LgrCre', 'Setdb1KO')) +
  theme(text=element_text(family = "FreeSerif" ), axis.text = element_text(size = 10, face = "bold"), legend.text = element_text(size =10)) +
  theme(legend.position="none")









