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
Lgr5Cre_Setdb1KO_A.data=Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004866/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_A = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_A.data, project = "Setdb1KO A", 
                                        min.cells = 3,
                                        min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_Setdb1KO_B.data <- Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004870/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_B = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_B.data, project = "Setdb1KO B", 
                                        min.cells = 3, 
                                        min.features = 200)
# Merge objects
Lgr5Cre_Setdb1KO=merge(Lgr5Cre_Setdb1KO_A, y=Lgr5Cre_Setdb1KO_B, add.cell.ids=c("Rep_A","Rep_B"), project="Lgr5Cre_Setdb1KO")

# CREATE OBJECT LGR5Cre_WT

# Load dataset A (Object A)
Lgr5Cre_WT_A.data=Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004858/outs/filtered_feature_bc_matrix/")
Lgr5Cre_WT_A = CreateSeuratObject(counts = Lgr5Cre_WT_A.data, project = "Lgr5Cre A", 
                                  min.cells = 3, 
                                  min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_WT_B.data <- Read10X(data.dir = "/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004862/outs/filtered_feature_bc_matrix/")
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
Lgr5Cre_MERGED=merge(Lgr5Cre_Setdb1KO.filt, y=Lgr5Cre_WT.filt, 
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
                         subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 18000)

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
Lgr5Cre_MERGED_Renamed=RenameIdents(Lgr5Cre_MERGED,  `0` = "Stem I", `1` = "Stem II", `2` = "Progenitor I", `3` = "Stem III", 
                                    `4` = "Goblet III", `5` = "Ent.Immature", `6` = "Progenitor II", `7` = "Goblet I", 
                                    `8` = "Ent.Immature", `9` = "Enteroendocrine I", `10` = "Tuft", `11` = "Paneth", 
                                    `12` = "Ent.Mature", `13` = "Enteroendocrine II", `14` = "Goblet II", `15` = "Unclassified")

# Set order to clusters
levels(Lgr5Cre_MERGED_Renamed) = c("Stem I", "Stem II", "Stem III","Progenitor I",
                                   "Progenitor II", "Ent.Immature", "Ent.Mature",
                                   "Goblet I", "Goblet II", "Goblet III", "Paneth",
                                   "Tuft", "Enteroendocrine I", "Enteroendocrine II", "Unclassified")

#----------------------------------------------------------------------------------------------------------------------------------------------------

# DotPlots per Cluster

# STEM/TA/PROGENITOR
DotPlot(Setdb1KO,features = c("Ascl2","Axin2","Bex1","Fzd2","Fzd7","Gkn3",
                                            "Hes1","Igfbp4","Lgr5","Lrp5","Lrp6","Notch1","Olfm4","Prom1","Slc12a2","Smo","Yap1"   
                                            ,         "Stmn1", "Tubb5"       ,     "Aurka","Birc5","Ccna2","Ccnb1","Ccnb2",
                                            "Cdc20","Cdc25c","Cdkn2d","Cdkn3","Cenpa","Cps1","Kif22","Kif23","Melk","Nek2","Plk1",
                                            "Rbp7","Sapcd2","Tacc3","Ube2c"  ), 
        cols = c("dodgerblue","goldenrod1"))  + xlab('') +  ylab('') + RotatedAxis()  + 
  theme(text=element_text(family = "FreeSerif" )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))
# Entrocyte/GOBLET
DotPlot(Setdb1KO,features = c("Gsdmc4","Krt8","Prss32","Reg3g"          , 
                                            "Aldob","Alpi","Apoa1","Apoa4","Clec2h","Dpep1","Elf3","Fabp1","Fabp6","Fam151a","Hnf4a","Hnf4aos","Hnf4g",
                                            "Lct","Mep1a","Muc3","Naaladl1","Neu1","Nudt4","Phgr1","Plb1","Pmp22","Prss30","Sis","Slc34a2",
                                            "Slc51a","Slc51b","Tmigd1","Xpnpep2"        ,   
                                            "Agr2","Atoh1","Ccl6","Ccl9","Clca1","Clca3a1","Fcgbp","Foxa1",
                                            "Klf4","Klk1","Lrrc26","Muc2","Spdef","Spink4","Tff3","Tpsg1","Tspan13","Txndc5","Zg16"), 
        cols = c("dodgerblue","goldenrod1"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" ))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))

# Paneth/Tuft/Enteroendo
DotPlot(Setdb1KO,features = c("Ang4","Defa17","Defa21","Defa22","Defa24","Defa3",
                                            "Dll4","Gfi1","Gm14851","Lyz1","Mmp7","Mptx2","Sox9"  , 
                                            "Adh1","Aldh2","Alox5ap","Avil","Cd24a","Dclk1","Fyb","Gfi1b","Gng13","Hck",
                                            "Il13ra1","Il25","Kctd12","Lrmp","Ltc4s","Ly6g6f","Ptpn6","Ptprc","Rgs13",
                                            "Sh2d6","Skap2","Tmem176a","Trpm5"       ,          
                                            "Arx","Bex2","Cck","Chga","Chgb","Cpe","Fam183b","Fev","Foxa2",
                                            "Gcg","Gch1","Gck","Gfra3","Hmgn3","Insm1","Isl1","Marcksl1","Neurod1",
                                            "Neurod2","Neurog3","Nkx2-2","Pax6","Pcsk1n","Pyy","Sst","Tac1","Tph1","Vwa5b2" ),
        cols = c("dodgerblue","goldenrod1"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))

# Lgr5
DotPlot(WT,features = "Lgr5",   cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))
# GFP
DotPlot(Lgr5Cre_MERGED_Renamed,features = "GFP",   cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  +  
  theme(text=element_text(family = "FreeSerif" )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))

#----------------------------------------------------------------------------------------------------------------------------------------------------

# Split Umap1 to WT & WT

split_obj = SplitObject(Lgr5Cre_MERGED_Renamed, split.by = "orig.ident")

# See how the objects are named
View(split_obj)

Setdb1KO = split_obj[["Setdb1KO"]]

WT = split_obj[["Lgr5Cre"]]

Setdb1KO <- FindVariableFeatures(Setdb1KO, selection.method = "vst", nfeatures = 2000)
WT<- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 2000)

# WT
top50WT <- head(VariableFeatures(WT), 50)
# Plot variable features
plot1 <- VariableFeaturePlot(WT)
plot2 <- LabelPoints(plot = plot1, points = top50WT, repel = TRUE)
plot2 + ggtitle("WT Most Varibale Features")+theme(text=element_text(family = "FreeSerif" ))

#----------------------------------------------------------------------------------------------------------------------------------------------------

# SCP PLOTS/ANALYSIS

# Rename clusters
Lgr5Cre_MERGED_Renamed$seurat_clusters <- factor(x = Lgr5Cre_MERGED_Renamed$seurat_clusters,
                                                 levels = c("0","1", "3", 
                                                            "2", "6", "5", "8", "12", 
                                                            "7", "14", "4", "11", "10", 
                                                            "9", "13", "15"),
                                                 labels = c( "Stem I", "Stem II", "Stem III",  
                                                             "Progenitor I", "Progenitor II",  "Ent.Immature", "Ent.Immature","Ent.Mature",  
                                                             "Goblet I", "Goblet II", "Goblet III", "Paneth", "Tuft", 
                                                             "Enteroendocrine I", "Enteroendocrine II","Unclassified"))

colnames(Lgr5Cre_MERGED_Renamed@meta.data)[colnames(Lgr5Cre_MERGED_Renamed@meta.data) == "seurat_clusters"] <- "CellType"

# dimplot highlight lgr5cre cells
CellDimPlot(
  srt = Lgr5Cre_MERGED_Renamed, group.by = c("CellType"), 
  reduction = "UMAP", theme_use = "theme_blank", show_stat = TRUE,
  palette = "Spectral", label = TRUE,
  cells.highlight = colnames(Lgr5Cre_MERGED_Renamed)[Lgr5Cre_MERGED_Renamed$CellType == c("Stem I","Stem II", "Stem III","Progenitor I", "Progenitor II")]
) + theme(text=element_text(family = "FreeSerif" ))

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

# WT Trajectory
WT_1st <- RunSlingshot(srt = WT_1st, group.by = "CellType", reduction = "UMAP", start = "Stem I")

CellDimPlot(WT_1st, group.by = "CellType", reduction = "UMAP", lineages = paste0("Lineage", 1:7), lineages_span = 0.1,
            palette = "Spectral"
) + theme(text=element_text(family = "FreeSerif" )) # 600 / 350

FeatureDimPlot(WT_1st, features = paste0("Lineage", 1:7), reduction = "UMAP", theme_use = "theme_blank") # 650 /600

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Make heatmap on lineage 7 (Loop lineage) in WT
WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage1"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage1"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

# Make heatmap on lineage 7 (Loop lineage) in WT
WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage2"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage2"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

# Make heatmap on lineage 7 (Loop lineage) in WT
WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage3"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage3"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage4"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage4"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage5"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage5"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage6"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage6"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage7"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage7"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)


# Make heatmap on lineage 7 (Loop lineage) in WT
WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage4"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage4"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

# Make heatmap on lineage 7 (Loop lineage) in WT
WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage5"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage5"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage6"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage6"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)


WT_1st <- RunDynamicFeatures(srt = WT_1st, lineages = c("Lineage7"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = WT_1st, lineages = c("Lineage7"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)



#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# KO Trajectory
KO_1st <- RunSlingshot(srt = KO_1st, group.by = "CellType", reduction = "UMAP", start = "Stem I")

CellDimPlot(KO_1st, group.by = "CellType", reduction = "UMAP", lineages = paste0("Lineage", 1:5), lineages_span = 0.1,
            palette = "Spectral"
) + theme(text=element_text(family = "FreeSerif" )) # 600 / 350

FeatureDimPlot(KO_1st, features = paste0("Lineage", 1:5), reduction = "UMAP", theme_use = "theme_blank") # 650/600

# Make heatmap on lineage 7 (Loop lineage) in WT
KO_1st <- RunDynamicFeatures(srt = KO_1st, lineages = c("Lineage1"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage1"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

# Make heatmap on lineage 7 (Loop lineage) in WT
KO_1st <- RunDynamicFeatures(srt = KO_1st, lineages = c("Lineage2"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage2"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

# Make heatmap on lineage 7 (Loop lineage) in WT
KO_1st <- RunDynamicFeatures(srt = KO_1st, lineages = c("Lineage3"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage3"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

KO_1st <- RunDynamicFeatures(srt = KO_1st, lineages = c("Lineage4"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage4"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

KO_1st <- RunDynamicFeatures(srt = KO_1st, lineages = c("Lineage5"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage5"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

KO_1st <- RunDynamicFeatures(srt = KO_1st, lineages = c("Lineage3"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage3"),
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),pseudotime_palette = "magma", nlabel = 20,
  cell_annotation = "CellType", cell_annotation_palette = "Spectral",pseudotime_label_color = "red", 
  #height = 5, width = 2
)
print(ht$plot)

#----------------------------------------------------------------------------------------------------------------------------------------------------

# DYNAMIC FEATURES ON LOOP LINEAGES OF WT AND KO OF LGR5CRE_MERGED
KO_1st <- RunDynamicFeatures(srt = KO_1st, lineages = c("Lineage1","Lineage4"), n_candidates = 200)

# KO LINEAGE1
ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage1"), heatmap_palette = "viridis", nlabel = 35, pseudotime_palette = ("RdYlBu"),
  #  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(ht$plot) + ggtitle("Setdb1KO")

# KO LINEAGE4
ht <- DynamicHeatmap(
  srt = KO_1st, lineages = c("Lineage4"), heatmap_palette = "viridis", nlabel = 35, pseudotime_palette = ("RdYlBu"),
  #  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(ht$plot) + ggtitle("Setdb1KO")

# LINEAGE 4 KO PANETH / GOBLET GENES\

# PANETH DEFA
plot = DynamicPlot(
  srt = KO_1st, lineages = c("Lineage4"), group.by = "CellType",
  features = c("Ang4","Defa17","Defa21","Defa22","Defa24","Defa3","Dll4","Gfi1","Gm14851","Lyz1","Mmp7","Mptx2","Sox9" # PANETH MARKERS EXPRESSED IN KO not IN WT
  ), point_palette = "Spectral",
  compare_lineages = TRUE, compare_features = FALSE
)
plot # 1000/900

# GOBLET
plot = DynamicPlot(
  srt = KO_1st, lineages = c("Lineage4"), group.by = "CellType",
  features = c("Agr2","Atoh1","Ccl6","Ccl9","Clca1","Clca3a1","Fcgbp","Foxa1","Klf4","Klk1","Lrrc26","Muc2","Spdef","Spink4","Tff3","Tpsg1","Tspan13","Txndc5","Zg16" # PANETH MARKERS EXPRESSED IN KO not IN WT
  ), point_palette = "Spectral",
  compare_lineages = TRUE, compare_features = FALSE
)
plot 

# Paneth markers on lineage 5  - KO / WT
#DynamicPlot(
#  srt = KO_1st, lineages = c("Lineage5"), group.by = "CellType",
#  features = c("Defa17", "Defa21", "Defa22", "Defa24", "Defa30","Lyz1"
# ),
#  compare_lineages = TRUE, compare_features = FALSE
#) + theme(text=element_text(family = "FreeSerif" ))

# Paneth markers on lineage 5  - KO / WT
#DynamicPlot(
# srt = WT_1st, lineages = c("Lineage5"), group.by = "CellType",
#  features = c("Defa17", "Defa21", "Defa22", "Defa24", "Defa30","Lyz1"
# ),
#  compare_lineages = TRUE, compare_features = FALSE
#) + theme(text=element_text(family = "FreeSerif" ))

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Subset Lgr5 clusters
# Subset clusters and recluster
Lgr5Cre_MERGED2=subset(Lgr5Cre_MERGED, idents = c(0,1,2,3,6))

top10_Lgr5Cre_MERGED2 <- head(VariableFeatures(Lgr5Cre_MERGED2), 10)

Lgr5Cre_MERGED2 <- RunPCA(Lgr5Cre_MERGED2, features = VariableFeatures(object = Lgr5Cre_MERGED2))
print(Lgr5Cre_MERGED2[["pca"]], dims = 1:5, nfeatures = 5)

Lgr5Cre_MERGED2=FindNeighbors(Lgr5Cre_MERGED2,dims = 1:10)
Lgr5Cre_MERGED2 <- FindClusters(Lgr5Cre_MERGED2, resolution = 0.5)

Lgr5Cre_MERGED2=RunUMAP(Lgr5Cre_MERGED2, dims = 1:9)
DimPlot(Lgr5Cre_MERGED2, group.by = "orig.ident") + theme(text=element_text(family = "FreeSerif" ))

# MERGED UMAP
CellDimPlot(
  srt = Lgr5Cre_MERGED2, group.by = c("orig.ident"), 
  reduction = "UMAP", theme_use = "theme_blank", show_stat = TRUE, palcolor = c("green3", "red2"),
  # palette = "Earth", 
) + theme(text=element_text(family = "FreeSerif" ))


#----------------------------------------------------------------------------------------------------------------------------------------------------

Lgr5Cre_MERGED2=RenameIdents(Lgr5Cre_MERGED2,  `0` = "Stem I", `1` = "Progenitor II", `2` = "Progenitor I", 
                             `3` = "Progenitor S", `4` = "Stem S", `5` = "Stem II", 
                             `6` = "Ent.Immature", `7` = "Goblet/Paneth")

levels(Lgr5Cre_MERGED2) = c("Stem I", "Stem II", "Stem S","Progenitor S", "Progenitor I", "Progenitor II","Ent.Immature" ,"Goblet/Paneth")

# STEM/TA/PROGENITOR DOTPLOTS
DotPlot(Lgr5Cre_MERGED2,features = c("Ascl2","Axin2","Bex1","Fzd2","Fzd7","Gkn3",
                                     "Hes1","Igfbp4","Lgr5","Lrp5","Lrp6","Notch1","Olfm4","Prom1","Slc12a2","Smo","Yap1"   
                                     ,"Stmn1", "Tubb5","Aurka","Birc5","Ccna2","Ccnb1","Ccnb2",
                                     "Cdc20","Cdc25c","Cdkn2d","Cdkn3","Cenpa","Cps1","Kif22","Kif23","Melk","Nek2","Plk1",
                                     "Rbp7","Sapcd2","Tacc3","Ube2c",
                                     "Gsdmc4","Krt8","Prss32","Reg3g", 
                                     "Agr2","Fcgbp","Spink4","Tff3","Zg16","Defa17","Defa21","Defa22","Defa24","Gm14851","Lyz1"
), 
cols = c("cadetblue1","darkred"))  + xlab('') +  ylab('') + RotatedAxis()  + 
  theme(text=element_text(family = "FreeSerif" )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"))

# Rename clusters
#Lgr5Cre_MERGED2$seurat_clusters <- factor(x = Lgr5Cre_MERGED2$seurat_clusters,
#                                         levels = c("0","1", "2", "3", "4", "5", "6", "7"),
#                                        labels = c( "Stem I",  "Progenitor IV",  "Progenitor III", 
#                                                   "Progenitor II", "Stem III", "Stem II", 
#                                                  "Progenitor I",  "Goblet/Paneth"))
#colnames(Lgr5Cre_MERGED2@meta.data)[colnames(Lgr5Cre_MERGED2@meta.data) == "seurat_clusters"] <- "CellType"

#############################################################################################################
Lgr5Cre_MERGED2$seurat_clusters <- factor(x = Lgr5Cre_MERGED2$seurat_clusters,
                                          levels = c("0","5", "4", "3", "2", "1", "6", "7"),
                                          labels = c( "Stem I",  "Stem II",  "Stem S", 
                                                      "Progenitor S", "Progenitor I", 
                                                      "Progenitor II","Ent.Immature","Goblet/Paneth"))
colnames(Lgr5Cre_MERGED2@meta.data)[colnames(Lgr5Cre_MERGED2@meta.data) == "seurat_clusters"] <- "CellType"


split_obj = SplitObject(Lgr5Cre_MERGED2,split.by = "orig.ident")

# See how the objects are named
View(split_obj)

Setdb1KO = split_obj[["Setdb1KO"]]

WT = split_obj[["Lgr5Cre"]]

#######################################################################################################

CellDimPlot(
  srt = Lgr5Cre_MERGED2, group.by = c("CellType"),
  reduction = "UMAP", theme_use = "theme_blank", show_stat = TRUE,
  palette = "Spectral"
  #palcolor = c("darkslategray3","deepskyblue3","dodgerblue4", "seagreen1", "seagreen3", "darkolivegreen3", "darkgreen" ,"firebrick4"),
  # palcolor = c("darkslategray3","darkgreen","darkolivegreen3", "seagreen3", "dodgerblue4", "deepskyblue3", "seagreen1" ,"firebrick4")
) + theme(text=element_text(family = "FreeSerif" ))

# Cell cycle
s.genes <- cc.genes$s.genes
s.genes <- tolower(s.genes)
s.genes = str_to_title(s.genes)

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- tolower(g2m.genes)
g2m.genes = str_to_title(g2m.genes)

Cell_Cycle_MERGED2 <- CellCycleScoring(Lgr5Cre_MERGED2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

CellStatPlot(Cell_Cycle_MERGED2, stat.by = "Phase", group.by = "CellType", 
             stat_type = "count", position = "dodge", 
             label = TRUE, palcolor = c("red3", "green3", "blue3")) + theme(text=element_text(family = "FreeSerif" ))
DimPlot(Cell_Cycle_MERGED2) + theme(text=element_text(family = "FreeSerif" ))

Cell_Cycle_WT <- CellCycleScoring(WT, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

CellStatPlot(Cell_Cycle_WT, stat.by = "Phase", group.by = "CellType", 
             stat_type = "count", position = "dodge", 
             label = TRUE, palcolor = c("red3", "green3", "blue3")) + theme(text=element_text(family = "FreeSerif" )) # 480 / 350
DimPlot(Cell_Cycle_WT) + theme(text=element_text(family = "FreeSerif" ))


Cell_Cycle_KO <- CellCycleScoring(Setdb1KO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

CellStatPlot(Cell_Cycle_KO, stat.by = "Phase", group.by = "CellType", 
             stat_type = "count", position = "dodge", 
             label = TRUE, palcolor = c("red3", "green3", "blue3")) + theme(text=element_text(family = "FreeSerif" )) # 480 / 350
DimPlot(Cell_Cycle_KO) + theme(text=element_text(family = "FreeSerif" ))

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# (Stem I Stem II Stem III )

# Subset 
STEMS = subset(Lgr5Cre_MERGED2, idents = c("Stem I", "Stem II", "Stem S"))

CellDimPlot(STEMS, group.by = "CellType", palette = "Spectral") + theme(text=element_text(family = "FreeSerif" )) # 370 / 250

s.genes <- cc.genes$s.genes
s.genes <- tolower(s.genes)
s.genes = str_to_title(s.genes)

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- tolower(g2m.genes)
g2m.genes = str_to_title(g2m.genes)

# Assign Cell-Cycle scores
Cell_Cycle_STEMS <- CellCycleScoring(STEMS, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

CellStatPlot(Cell_Cycle_STEMS, stat.by = "Phase", group.by = "CellType", bg.by = "CellType", palette = "Set1", stat_type = "count", position = "dodge")

# dimplot MERGED
CellDimPlot(
  srt = Cell_Cycle_STEMS, group.by = c("CellType"),
  reduction = "UMAP",  show_stat = TRUE, palette = "Spectral"
) + theme(text=element_text(family = "FreeSerif" ))


# Split STEMS_WT & STEMS_KO
SplitObject(STEMS, split.by = "ident")
n_cells=(FetchData(STEMS, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

STEMS_WT=subset(STEMS, subset = orig.ident == "Lgr5Cre")
STEMS_KO=subset(STEMS, subset = orig.ident == "Setdb1KO")


# WITH SEURAT
# Differential expression on STEMS_WT
#STEMS_WT_Markers <- FindAllMarkers(STEMS_WT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.37)
#STEMS_WT_Markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)

#STEMS_WT_Markers %>%
#  group_by(cluster) %>%
#  top_n(n = 15, wt = avg_log2FC) -> top10
#DoHeatmap(STEMS_WT, features = c("Acot1", "Cda","Gkn3","Hmgcs2", "Jaml","Onecut2","Prap1","Ttr",
#                                
#                               "Bex1","Cd74","Ceacam10","Defa21", "Defa22",
#                              "H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1", "H2-Eb1", "H2-K1",
#                             "Psmb8","Reg3b", "Reg3g",  "Tspan1", 
#                            
#                          "Atad2","Dek", "Dut", "Gmnn","Hells","Hmgb2","Lig1","Pclaf","Pcna",
#                           "Rfc4","Rpa2", "Siva1", "Slbp", "Slfn9", "Smc2", "Stmn1","Tk1","Top2a","Tyms",

#                         "Ascl2","Axin2","Bex1","Fzd2","Fzd7","Gkn3","Hes1","Igfbp4","Lgr5","Lrp5","Lrp6","Notch1",
#                        "Olfm4","Prom1","Slc12a2","Smo","Yap1"),raster = FALSE, draw.lines = TRUE, lines.width = 10, 
#group.colors = c("dodgerblue4", "lightseagreen", "darkseagreen3"), label = FALSE,
#size = 4.5, group.bar.height = 0.014 ) + NoLegend() + scale_fill_gradientn(colors = c("deepskyblue4", "black", "orange")) + theme(text = element_text(size = 12)) +  theme(text=element_text(family="FreeSerif", face = "bold"))

#write.table(STEMS_WT_Markers, 
#           file = "WT_STEMI_VS_II_VS_III_MARKERS.tsv", 
#          sep = "\t", row.names = FALSE)

################################################################################################################################

# WITH SCP
STEMS_WT$CellType <- factor(x = STEMS_WT$CellType,
                            levels = c("Stem I","Stem II", "Stem S"),
                            labels = c( "Stem I","Stem II", "Stem S"))


STEMS_WT <- RunDEtest(srt = STEMS_WT, group_by = "CellType", fc.threshold = 1, only.pos = FALSE)

DEGs_STEMS <- STEMS_WT@tools$DEtest_CellType$AllMarkers_wilcox
DEGs_STEMS <- DEGs_STEMS[with(DEGs_STEMS, avg_log2FC > 0.37 & p_val_adj < 0.01), ]

write.table(DEGs_STEMS, 
            file = "DEGs_STEMS.tsv", 
            sep = "\t", row.names = FALSE)

STEMS_WT <- AnnotateFeatures(STEMS_WT, species = "Mus_musculus", db = c("TF", "SP"))
ht <- FeatureHeatmap(
  srt = STEMS_WT, group.by = "CellType", features = c("Acot1", "Cda","Gkn3","Hmgcs2", "Jaml","Onecut2","Prap1","Ttr",
                                                      
                                                      "Bex1","Cd74","Ceacam10","Defa21", "Defa22",
                                                      "H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1", "H2-Eb1", "H2-K1",
                                                      "Psmb8","Reg3b", "Reg3g",  "Tspan1", 
                                                      
                                                      "Atad2","Dek", "Dut", "Gmnn","Hells","Hmgb2","Lig1","Pclaf","Pcna",
                                                      "Rfc4","Rpa2", "Siva1", "Slbp", "Slfn9", "Smc2", "Stmn1","Tk1","Top2a","Tyms",
                                                      
                                                      "Ascl2","Axin2","Bex1","Gkn3","Hes1","Igfbp4","Lgr5","Lrp5","Lrp6","Notch1",#"Fzd2","Fzd7" (not expressed)
                                                      "Olfm4","Prom1","Slc12a2","Smo","Yap1"), 
  #  features = DEGs_STEMS$gene,
  feature_split = DEGs_STEMS$group3,
  feature_split_palcolor = list(c("dodgerblue4", "lightseagreen", "darkseagreen3")),
  show_row_names = FALSE,show_column_names = FALSE,
  #  species = "Mus_musculus", db = c("GO_BP"), anno_terms = TRUE, 
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),
  # heatmap_palette = "viridis", 
  group_palcolor = list(c("dodgerblue4", "lightseagreen", "darkseagreen3")),
  #  feature_annotation = c("TF", "SP"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")), 
  nlabel = 27,
  # height = 6, width = 4
)
print(ht$plot) + theme(text=element_text(family = "FreeSerif" ))
################################################################################################################################


# SUBSET STEM I AND II AND CREATE VIOLIN / FEAUTRE PLOTS FOR H2 antigen MARKERS and STEM I SPECIFIC MARKERS 
STEM_I_II = subset(Lgr5Cre_MERGED2, idents = c("Stem I", "Stem II"))

# STEM II Specific
FeatureStatPlot(
  srt = STEM_I_II, group.by = "CellType",
  stat.by = c("H2-Aa","H2-Eb1", "H2-Ab1","H2-D1","H2-K1" ,"H2-DMa"), add_box = TRUE, palette = "Spectral",
  #palcolor = c("dodgerblue4", "lightseagreen"), 
  stack = TRUE, 
  box_color = "gray20", box_width = 0.03
) + theme(text=element_text(family = "FreeSerif" )) #600/400

FeatureDimPlot(
  srt = Lgr5Cre_MERGED2, features = c("H2-Aa","H2-Eb1", "H2-Ab1","H2-D1","H2-K1" ,"H2-DMa"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5
)

# STEM I Specific
FeatureStatPlot(
  srt = STEM_I_II, group.by = "CellType",
  stat.by = c("Prap1","Jaml","Gkn3", "Acot1","Hmgcs2","Ttr"
  ), add_box = TRUE, palette = "Spectral", #palcolor = c("darkslategray3", "deepskyblue3"), 
  stack = TRUE, 
  box_color = "gray20", box_width = 0.03
  
) + theme(text=element_text(family = "FreeSerif" ))

FeatureDimPlot(
  srt = Lgr5Cre_MERGED2, features = c("Prap1","Jaml","Gkn3", "Acot1","Hmgcs2","Ttr"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5
)

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# WT VS KO in STEM I AND II

# CLUSTER STEM I 
STEMI <- subset(STEMS, idents = "Stem I")

SplitObject(STEMI, split.by = "ident")
n_cells=(FetchData(STEMI, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

md = STEMI@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "CellType")]
md[, .N, by = c("orig.ident", "CellType")] %>% dcast(., orig.ident ~ CellType, value.var = "N")

Idents(STEMI) = STEMI$orig.ident

# MARKERS EXPRESSED IN WT (STEM I)
WT_EXPRESSED <- FindMarkers(STEMI, ident.1 = "Lgr5Cre", ident.2 = "Setdb1KO", logfc.threshold = 0.25, only.pos = T)
# MARKERS EXPRESSED IN KO (STEM I)
KO_EXPRESSED <- FindMarkers(STEMI, ident.1 = "Setdb1KO", ident.2 = "Lgr5Cre", logfc.threshold = 0.25, only.pos = T)

write.table(WT_EXPRESSED, file = "STEM_I_WT_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

write.table(KO_EXPRESSED, file = "STEM_I_KO_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------------------------

# CLUSTER STEM II
STEMII <- subset(STEMS, idents = "Stem II")

SplitObject(STEMII, split.by = "ident")
n_cells=(FetchData(STEMII, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

md = STEMII@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "CellType")]
md[, .N, by = c("orig.ident", "CellType")] %>% dcast(., orig.ident ~ CellType, value.var = "N")

Idents(STEMII) = STEMII$orig.ident

# MARKERS EXPRESSED IN WT (STEM I)
WT_EXPRESSED <- FindMarkers(STEMII, ident.1 = "Lgr5Cre", ident.2 = "Setdb1KO", logfc.threshold = 0.25, only.pos = T)
# MARKERS EXPRESSED IN KO (STEM I)
KO_EXPRESSED <- FindMarkers(STEMII, ident.1 = "Setdb1KO", ident.2 = "Lgr5Cre", logfc.threshold = 0.25, only.pos = T)

write.table(WT_EXPRESSED, file = "STEM_II_WT_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

write.table(KO_EXPRESSED, file = "STEM_II_KO_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Rename orig.ident to Condition for STEM I
colnames(STEMI@meta.data)[colnames(STEMI@meta.data) == "orig.ident"] <- "Condition"

# STEM I KO EXPRESSED  (DOTPLOTS)
ht <- GroupHeatmap(
  srt = STEMI,
  features = c("Defa17","Defa21","Defa22","Defa24","Defa30","Gm14851","Lyz1","Phgr1", # PANETH MARKERS EXPRESSED IN KO not IN WT
               "Spink4", "Tff3", "Zg16" # GOBLET MARKERS EXPRESSED IN KO not in WT
  ),
  group.by = c("Condition"),
  heatmap_palcolor = c("dodgerblue","goldenrod1"), 
  group_palcolor = list(c("lavenderblush4", "gray90")), 
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, 
)
print(ht$plot) + ggtitle("STEM I")  
# 350 / 350 width height
# STEM I KO EXPRESSED  (HEATMAP)
ht <- FeatureHeatmap(
  srt = STEMI, group.by = "Condition",
  features = c("1500011B03Rik","AY761184","Akr1c13","Alad","Defa17","Defa21","Defa22","Defa24","Defa26", "Defa29","Defa30", 
               "Fabp2","Gm14851", "Gm42031", "Gm49980","H2afz","Hes6","Hist1h2bc","Itln1","Kcne3","Lipt2","Ly6e",
               "Lyz1","Mrpl15", "Mrpl17","Phgr1","Plscr1","Rbm3","Rbp2","Rgcc","Rnase4","Scand1","Smim24",
               "Snhg1","Snhg18", "Spink4", "Sycn","Tff3","Tubb2b", "Zg16","Znrd2"), 
  show_row_names = FALSE,show_column_names = FALSE,
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),
  group_palcolor = list(c("lavenderblush4", "gray90")),
  nlabel = 20,
)
print(ht$plot) + theme(text=element_text(family = "FreeSerif" )) # 450 / 350

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Rename orig.ident to Condition for STEM II
colnames(STEMII@meta.data)[colnames(STEMII@meta.data) == "orig.ident"] <- "Condition"

# STEM II KO EXPRESSED
ht <- GroupHeatmap(
  srt = STEMII,
  features = c("Defa17","Defa21","Defa22","Defa24","Defa30","Gm14851","Lyz1", # PANETH MARKERS EXPRESSED IN KO not IN WT
               "Spink4", "Tff3", "Zg16" # GOBLET MARKERS EXPRESSED IN KO not in WT
  ),
  group.by = c("Condition"),
  heatmap_palcolor = c("dodgerblue","goldenrod1"), 
  group_palcolor = list(c("lavenderblush4", "gray90")), 
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, 
)
print(ht$plot) + ggtitle("STEM II")  
# 350 / 350 width height
ht <- FeatureHeatmap(
  srt = STEMII, group.by = "Condition",
  features = c("1500011B03Rik",	"2200002D01Rik",	"AY761184",	"Akr1c13",	"Alad",	"Ascl2",	"Bex1",	"Csrp2",
               "Defa17",	"Defa21",	"Defa22",	"Defa24",	"Defa29",	"Defa30",
               "Dhrs4",	"Fam32a",	"Fcgbp",	"Fos",	"Galnt12",	"Gm14851",	"Gm19696",	"Gm42031",	"Gm49980",	"Gm5485",	
               "H1f0",	"Hdhd3",	"Hes6",	"Hist1h2bc",	"Ier2",	"Ifi27l2b",	"Itln1",	"Junb",	"Kcne3",	"Krcc1",	"Lipt2",	"Ly6e",	"Lyz1",	
               "Nrn1",	"Oit1",	"Oxa1l",	"Phgr1",	"Pibf1",	"Plscr1",	"Prxl2b",	"Rbm3",	"Rgcc",	"Rnase4",	"Scand1",	"Sdhaf1",	"Serpinb1a",
               "Smim24",	"Snhg1",	"Spink4",	"Sult1d1",	"Sycn",	"Tff3",	"Tstd1",	"Ttr",	"Tubb2b",	"Uqcr10",	"Vsig10",	"Zg16"
  ), 
  show_row_names = FALSE,show_column_names = FALSE,
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),
  group_palcolor = list(c("lavenderblush4", "gray90")),
  nlabel = 20,
)
print(ht$plot) + theme(text=element_text(family = "FreeSerif" ))

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# PROGENITOR DIFFERENTIATION

PROGENITORS = subset(Lgr5Cre_MERGED2, idents = c("Progenitor S", "Progenitor I", "Progenitor II"))

CellDimPlot(PROGENITORS,group.by = "CellType", palette = "Spectral") + theme(text=element_text(family = "FreeSerif" ))   # 330 / 180


PROGENITORS$CellType <- factor(x = PROGENITORS$CellType,
                               levels = c("Progenitor S","Progenitor I", "Progenitor II"),
                               labels = c( "Progenitor S","Progenitor I", "Progenitor II"))

#colnames(Lgr5Cre_MERGED2@meta.data)[colnames(Lgr5Cre_MERGED2@meta.data) == "seurat_clusters"] <- "CellType"

# Split STEMS_WT & STEMS_KO
SplitObject(PROGENITORS, split.by = "ident")
n_cells=(FetchData(PROGENITORS, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

PROGENITORS_WT=subset(PROGENITORS, subset = orig.ident == "Lgr5Cre")
#PROGENITORS_KO=subset(PROGENITORS, subset = orig.ident == "Setdb1KO")

# Differential expression on STEMS_WT
#PROGENITORS_WT_Markers <- FindAllMarkers(PROGENITORS_WT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.37)
#PROGENITORS_WT_Markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)

#PROGENITORS_WT_Markers %>%
# group_by(cluster) %>%
#  top_n(n = 15, wt = avg_log2FC) -> top10
#DoHeatmap(PROGENITORS_WT, features = c("Acot1", "Cda","Gkn3","Hmgcs2", "Jaml","Onecut2","Prap1","Ttr",

#                                      "Bex1","Cd74","Ceacam10","Defa21", "Defa22",
#                                     "H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1", "H2-Eb1", "H2-K1",
#                                    "Psmb8","Reg3b", "Reg3g",  "Tspan1", 

#                                   "Atad2","Dek", "Dut", "Gmnn","Hells","Hmgb2","Lig1","Pclaf","Pcna",
#                                  "Rfc4","Rpa2", "Siva1", "Slbp", "Slfn9", "Smc2", "Stmn1","Tk1","Top2a","Tyms",

#                                 "Ascl2","Axin2","Bex1","Fzd2","Fzd7","Gkn3","Hes1","Igfbp4","Lgr5","Lrp5","Lrp6","Notch1",
#                                "Olfm4","Prom1","Slc12a2","Smo","Yap1"),raster = FALSE, draw.lines = TRUE, lines.width = 10, 
#  group.colors = c("dodgerblue4", "lightseagreen", "darkseagreen3"), label = FALSE,
# size = 4.5, group.bar.height = 0.014 ) + NoLegend() + scale_fill_gradientn(colors = c("deepskyblue4", "black", "orange")) + theme(text = element_text(size = 12)) +  theme(text=element_text(family="FreeSerif", face = "bold"))

#write.table(PROGENITORS_WT_Markers, 
#           file = "WT_PROGENITOR_I_VS_II_VS_III_MARKERS.tsv", 
#          sep = "\t", row.names = FALSE)


PROGENITORS_WT <- RunDEtest(srt = PROGENITORS_WT, group_by = "CellType", fc.threshold = 1, only.pos = FALSE)

DEGs_PROGENITORS <- PROGENITORS_WT@tools$DEtest_CellType$AllMarkers_wilcox
DEGs_PROGENITORS <- DEGs_PROGENITORS[with(DEGs_PROGENITORS, avg_log2FC > 0.5 & p_val_adj < 0.01), ]

write.table(DEGs_PROGENITORS, 
            file = "DEGs_PROGENITORS.tsv", 
            sep = "\t", row.names = FALSE)

ht <- FeatureHeatmap(
  srt = PROGENITORS_WT, group.by = "CellType", 
  features = c("Pcna", "Pclaf","Lig1","Dut", "Slbp", "Tyms","Siva1", "Rrm2","Hells", "Slfn9","Olfm4",
               "Ifitm3","Slc12a2", "Cenpa",  "Clca3b",  "Cenpf",  "Ccnb2", "Sox4",  "Jaml",  "Ptms",  "Tubb2b","Lgals2", 
               "Krt19",  "Pycard",  "Krt8","Smim24",  "Rbp7",  "Dmbt1",  "Phgr1",  
               "Crip1",  "Ldha","Arg2", "Aurka","Birc5","Ccna2","Ccnb1","Ccnb2","Cdc20","Cdc25c",
               "Cdkn2d","Cdkn3","Cenpa","Cps1","Kif22","Kif23","Melk","Nek2","Plk1","Rbp7","Sapcd2","Tacc3","Ube2c"),
  #features = DEGs_PROGENITORS$gene,
  feature_split = DEGs_PROGENITORS$group3,
  feature_split_palcolor = list(c("khaki", "goldenrod1", "orange2")),
  show_row_names = FALSE,show_column_names = FALSE,
  #  species = "Mus_musculus", db = c("GO_BP"), anno_terms = TRUE, 
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),
  # heatmap_palette = "viridis", 
  group_palcolor = list(c("khaki", "goldenrod1", "orange2")),
  #  feature_annotation = c("TF", "SP"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")), 
  nlabel = 27,
  # height = 6, width = 4
)
print(ht$plot) + theme(text=element_text(family = "FreeSerif" ))


#----------------------------------------------------------------------------------------------------------------------------------------------------


# SUBSET Progenitor I AND II AND CREATE VIOLIN / FEAUTRE PLOTS FOR SPECIFIC MARKERS 
PROG_I_II = subset(PROGENITORS_WT, idents = c("Progenitor I", "Progenitor II"))

PROG_I_II$CellType <- factor(x = PROG_I_II$CellType,
                             levels = c("Progenitor I", "Progenitor II"),
                             labels = c("Progenitor I", "Progenitor II"))

# PROG I Specific
FeatureStatPlot(
  srt = PROG_I_II, group.by = "CellType",
  stat.by = c("Olfm4", "Ifitm3",  "Jaml","Slc12a2", "Sox4","Ascl2"
  ), add_box = TRUE, #palette = "Spectral",Olfm4
  palcolor = c("goldenrod1", "orange2"), 
  stack = TRUE, 
  box_color = "gray20", box_width = 0.03
) + theme(text=element_text(family = "FreeSerif" )) #600/400

FeatureDimPlot(
  srt = Lgr5Cre_MERGED2, features = c("Olfm4", "Ifitm3",  "Jaml","Slc12a2", "Sox4","Ascl2"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5
)

# PROG II Specific
FeatureStatPlot(
  srt = PROG_I_II, group.by = "CellType",
  stat.by = c("Arg2","Rbp7","Dmbt1","Smim24","Phgr1","Ldha"
  ), add_box = TRUE, #palette = "Spectral", #palcolor = c("darkslategray3", "deepskyblue3"),
  palcolor = c("goldenrod1", "orange2"), 
  stack = TRUE, 
  box_color = "gray20", box_width = 0.03
  
) + theme(text=element_text(family = "FreeSerif" ))

FeatureDimPlot(
  srt = Lgr5Cre_MERGED2, features = c("Arg2","Rbp7","Dmbt1","Smim24","Phgr1","Ldha"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5
)



#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# WT VS KO in PROGENITOR I AND II

# CLUSTER PROGENITOR I 
PROGENITORI <- subset(PROGENITORS, idents = "Progenitor I")

SplitObject(PROGENITORI, split.by = "ident")
n_cells=(FetchData(PROGENITORI, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

md = PROGENITORI@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "CellType")]
md[, .N, by = c("orig.ident", "CellType")] %>% dcast(., orig.ident ~ CellType, value.var = "N")

Idents(PROGENITORI) = PROGENITORI$orig.ident

# MARKERS EXPRESSED IN WT (STEM I)
WT_EXPRESSED <- FindMarkers(PROGENITORI, ident.1 = "Lgr5Cre", ident.2 = "Setdb1KO", logfc.threshold = 0.25, only.pos = T)
# MARKERS EXPRESSED IN KO (STEM I)
KO_EXPRESSED <- FindMarkers(PROGENITORI, ident.1 = "Setdb1KO", ident.2 = "Lgr5Cre", logfc.threshold = 0.25, only.pos = T)

write.table(WT_EXPRESSED, file = "PROGENITOR_I_WT_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

write.table(KO_EXPRESSED, file = "PROGENITOR_I_KO_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------------------------

# CLUSTER PROGENITOR II
PROGENITORII <- subset(PROGENITORS, idents = "Progenitor II")

SplitObject(PROGENITORII, split.by = "ident")
n_cells=(FetchData(PROGENITORII, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

md = PROGENITORII@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "CellType")]
md[, .N, by = c("orig.ident", "CellType")] %>% dcast(., orig.ident ~ CellType, value.var = "N")

Idents(PROGENITORII) = PROGENITORII$orig.ident

# MARKERS EXPRESSED IN WT (STEM I)
WT_EXPRESSED <- FindMarkers(PROGENITORII, ident.1 = "Lgr5Cre", ident.2 = "Setdb1KO", logfc.threshold = 0.25, only.pos = T)
# MARKERS EXPRESSED IN KO (STEM I)
KO_EXPRESSED <- FindMarkers(PROGENITORII, ident.1 = "Setdb1KO", ident.2 = "Lgr5Cre", logfc.threshold = 0.25, only.pos = T)

write.table(WT_EXPRESSED, file = "PROGENITOR_II_WT_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

write.table(KO_EXPRESSED, file = "PROGENITOR_II_KO_EXPRESSED.tsv", sep = "\t", row.names = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Rename orig.ident to Condition for STEM I
colnames(PROGENITORI@meta.data)[colnames(PROGENITORI@meta.data) == "orig.ident"] <- "Condition"

# STEM I KO EXPRESSED  (DOTPLOTS)
ht <- GroupHeatmap(
  srt = PROGENITORI,
  features = c("Defa17",	"Defa21",	"Defa22",	"Defa24",	"Defa30",	"Gm14851",	"Lyz1",	# PANETH
               "Fcgbp",	"Phgr1",	"Spink4",	"Tff3",	"Zg16" # GOBLET
  ),
  group.by = c("Condition"),
  heatmap_palcolor = c("dodgerblue","goldenrod1"), 
  group_palcolor = list(c("lavenderblush4", "gray90")), 
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, 
)
print(ht$plot) + ggtitle("PROG I")  
# 350 / 350 width height
# STEM I KO EXPRESSED  (HEATMAP)
ht <- FeatureHeatmap(
  srt = PROGENITORI, group.by = "Condition",
  features = c("1500011B03Rik",	"AY761184",	"Akr1c13",	"Alad",	"Cdkn3",	"Defa17",	"Defa21",	"Defa22",	"Defa24",	"Defa29",	"Defa30",	
               "Fabp1",	"Fabp2",	"Fcgbp",	"Gm14851",	"Gm42031",	"Gm49980",	"H1f0",	"Hes6",	"Hist1h2bc",	"Itln1",	"Kcne3",	"Lipt2",
               "Ly6e",	"Lyz1",	"Phgr1",	"Prxl2b",	"Rbm3",	"Rbp2",	"Rgcc",	"Rnase4",	"Smim24",	"Snhg1",	"Snhg18",	"Spink4",	"Sult1d1",	
               "Sycn",	"Tff3",	"Tm4sf5",	"Tstd1",	"Tubb2a",	"Tubb2b",	"Uqcr10",	"Zg16"
  ), 
  show_row_names = FALSE,show_column_names = FALSE,
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),
  group_palcolor = list(c("lavenderblush4", "gray90")),
  nlabel = 20,
)
print(ht$plot) + theme(text=element_text(family = "FreeSerif" ))

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Rename orig.ident to Condition for STEM II
colnames(PROGENITORII@meta.data)[colnames(PROGENITORII@meta.data) == "orig.ident"] <- "Condition"

# STEM II KO EXPRESSED
ht <- GroupHeatmap(
  srt = PROGENITORII,
  features = c("Defa21","Defa22","Defa24","Defa30","Gm14851","Spink4","Tff3","Zg16"
               
  ),
  group.by = c("Condition"),
  heatmap_palcolor = c("dodgerblue","goldenrod1"), 
  group_palcolor = list(c("lavenderblush4", "gray90")), 
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, 
)
print(ht$plot) + ggtitle("PROG II")  
# 350 / 350 width height
ht <- FeatureHeatmap(
  srt = PROGENITORII, group.by = "Condition",
  features = c("AY761184",	"Alad",	"Defa17",	"Defa21",	"Defa22",	"Defa24",	"Defa29",	"Defa30",
               "Gm14851",	"Gm42031",	"Gm49980",	"H1f0",	"Hdhd3",	"Hist1h2bc",	"Hmgn1",	"Itln1",
               "Kcne3",	"Khk",	"Lipt2",	"Ly6e",	"Lyz1",	"Pcbd1",	"Phgr1","Qdpr",	"Rbm3",	"Rbp2",	
               "Rbp7",	"Rnase4",	"Sis",	"Smim24",	"Snhg1",	"Spink4",	"Sycn",	"Tff3",	"Ubc",	"Vsig10",	"Zg16"
  ), 
  show_row_names = FALSE,show_column_names = FALSE,
  heatmap_palcolor = c("deepskyblue4", "black", "orange"),
  group_palcolor = list(c("lavenderblush4", "gray90")),
  nlabel = 18,
)
print(ht$plot) + theme(text=element_text(family = "FreeSerif" ))

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

# Enrichment analysis on WT of STEMS-PROGENITORS ONLY
STEMS_WT <- RunEnrichment(
  srt = STEMS_WT, group_by = "CellType", db = "GO_BP", species = "Mus_musculus",
  DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05"
)
EnrichmentPlot(
  srt = STEMS_WT, group_by = "CellType", group_use = c("Stem I", "Stem II", "Stem S", "Progenitor S","Progenitor I", "Progenitor II"),
  plot_type = "comparison"
)

PROGENITORS_WT <- RunEnrichment(
  srt = PROGENITORS_WT, group_by = "CellType", db = "GO_BP", species = "Mus_musculus",
  DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05"
)
EnrichmentPlot(
  srt = PROGENITORS_WT, group_by = "CellType", group_use = c("Stem I", "Stem II", "Stem S", "Progenitor S","Progenitor I", "Progenitor II"),
  plot_type = "comparison"
)


#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

#  CELL TO CELL HETEROGENEITY mRNA (Coefficient of Variability in gene expression)

# create two matrix file for gene count per cell (KO and WT)
# WT
Lgr5Cre_MERGED2_WT=subset(Lgr5Cre_MERGED2, subset = orig.ident == "Lgr5Cre")
# KO
Lgr5Cre_MERGED2_KO=subset(Lgr5Cre_MERGED2, subset = orig.ident == "Setdb1KO")

# Save file with gene expression per cell type per condition
# WT
write.table(Lgr5Cre_MERGED2_WT@assays[["RNA"]]@counts, file='WT_Gene_Count_per_Cell.tmp', quote=FALSE, sep='\t', col.names = TRUE)
# KO
write.table(Lgr5Cre_MERGED2_KO@assays[["RNA"]]@counts, file='KO_Gene_Count_per_Cell.tmp', quote=FALSE, sep='\t', col.names = TRUE)

# Read the filtered file
KO_matrix_Filtered = read.table("/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/SRT_SCP/KO_Gene_Count_per_Cell.tmp", sep = "\t", header = TRUE, row.names = 1 )
WT_matrix_Filtered = read.table("/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/SRT_SCP/WT_Gene_Count_per_Cell.tmp", sep = "\t", header = TRUE, row.names = 1 ) 

# calculate coefficient variability and add results in the file
KO_matrix_Filtered$cv  <- apply(KO_matrix_Filtered, 1,  function(x) sd(x) / mean(x) * 100)
WT_matrix_Filtered$cv  <- apply(WT_matrix_Filtered, 1,  function(x) sd(x) / mean(x) * 100)

# send back to bash for file manipulation and filtering
write.table(KO_matrix_Filtered, file = "KO.tmp", row.names = FALSE, sep = "\t")
write.table(WT_matrix_Filtered, file = "WT.tmp", row.names = FALSE, sep = "\t")

coefficient_variation = read.table("/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/SRT_SCP/coefficient_variation.tsv", sep = "\t", header = TRUE)
colnames(coefficient_variation) = c("Setdb1KO", "Lgr5Cre")

colors <- c("green3","red3")
boxplot(coefficient_variation, ylab = "coefficient variation", 
        col = colors, main = "Lgr5-Clusters", outline = FALSE)

##################################################################################################################################################################################################################

# Mann-Whitney-Wilcoxon Test

# Read the values of coefficient of variation for KO and WT
#ko_values = read.csv(file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/Talianidis_GFP_2_fixed/analysis/Seurat/coefficient_variability/TEST/", header = FALSE)
#wt_values = read.csv(file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/Talianidis_GFP_2_fixed/analysis/Seurat/coefficient_variability/TEST/wt", header = FALSE)
# convert list to double
#ko_values = as.numeric(unlist(ko_values))
#wt_values = as.numeric(unlist(wt_values))
# Run Mann-Whitney-Wilcoxon Test
wilcox.test(coefficient_variation$Setdb1KO, coefficient_variation$Lgr5Cre)












