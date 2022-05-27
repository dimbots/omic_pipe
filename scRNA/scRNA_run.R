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
library(Matrix)

########################################################################################################################################################################################################################


# CREATE OBJECT Lgr5Cre_SETDB1KO

# Load dataset A (Object A)
Lgr5Cre_Setdb1KO_A.data=Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004866/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_A = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_A.data, project = "Setdb1KO A", min.cells = 3, min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_Setdb1KO_B.data <- Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004870/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_B = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_B.data, project = "Setdb1KO B", min.cells = 3, min.features = 200)

# Merge objects
Lgr5Cre_Setdb1KO=merge(Lgr5Cre_Setdb1KO_A, y=Lgr5Cre_Setdb1KO_B, add.cell.ids=c("Rep_A","Rep_B"), project="Lgr5Cre_Setdb1KO")


##################################################################################################################################################################################################################


# CREATE OBJECT LGR5Cre_WT

# Load dataset A (Object A)
Lgr5Cre_WT_A.data=Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004858/outs/filtered_feature_bc_matrix/")
Lgr5Cre_WT_A = CreateSeuratObject(counts = Lgr5Cre_WT_A.data, project = "Lgr5Cre A", min.cells = 3, min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_WT_B.data <- Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004862/outs/filtered_feature_bc_matrix/")
Lgr5Cre_WT_B = CreateSeuratObject(counts = Lgr5Cre_WT_B.data, project = "Lgr5Cre B", min.cells = 3, min.features = 200)

# Merge objects
Lgr5Cre_WT=merge(Lgr5Cre_WT_A, y=Lgr5Cre_WT_B, add.cell.ids=c("Rep_A","Rep_B"), project="Lgr5Cre_WT")


##################################################################################################################################################################################################################


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

pdf("CowPlot_Setdb1KO.pdf", width=15, height=8)
cowplot::plot_grid(ncol = 2, DimPlot(Lgr5Cre_Setdb1KO.filt, group.by = DF.name) + NoAxes())
dev.off()

pdf("VlnPlot_Setdb1KO.pdf", width=15, height=8)
VlnPlot(Lgr5Cre_Setdb1KO.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
dev.off()

Lgr5Cre_Setdb1KO.filt = Lgr5Cre_Setdb1KO.filt[, Lgr5Cre_Setdb1KO.filt@meta.data[, DF.name] == "Singlet"]


##################################################################################################################################################################################################################


Lgr5Cre_WT <- NormalizeData(Lgr5Cre_WT, normalization.method = "LogNormalize", scale.factor = 10000)

Lgr5Cre_WT <- FindVariableFeatures(Lgr5Cre_WT, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Lgr5Cre_WT)
Lgr5Cre_WT <- ScaleData(Lgr5Cre_WT, features = all.genes)

Lgr5Cre_WT <- RunPCA(Lgr5Cre_WT, features = VariableFeatures(object = Lgr5Cre_Setdb1KO))

Lgr5Cre_WT <- RunUMAP(Lgr5Cre_WT, dims = 1:10)

nExp <- round(ncol(Lgr5Cre_WT) * 0.057)
Lgr5Cre_WT.filt <- doubletFinder_v3(Lgr5Cre_WT, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(Lgr5Cre_WT.filt@meta.data)[grepl("DF.classification", colnames(Lgr5Cre_WT.filt@meta.data))]

pdf("CowPlot_WT.pdf", width=15, height=8)
cowplot::plot_grid(ncol = 2, DimPlot(Lgr5Cre_WT.filt, group.by = DF.name) + NoAxes())
dev.off()

pdf("VlnPlot_WT.pdf", width=15, height=8)
VlnPlot(Lgr5Cre_WT.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
dev.off()


##################################################################################################################################################################################################################


# Merge Setdb1KO with WT to new object Lgr5Cre_MERGED
Lgr5Cre_MERGED=merge(Lgr5Cre_Setdb1KO.filt, y=Lgr5Cre_WT.filt, add.cell.ids=c("Setdb1KO","Lgr5Cre"), project="Lgr5Cre_MERGED")

# Change replicates names within object. e.g(Setdb1KO_repA & Setdb1KO_repB -> Setdb1KO)
Lgr5Cre_MERGED$orig.ident=plyr::mapvalues(x=Lgr5Cre_MERGED$orig.ident, from = c("Setdb1KO A", "Setdb1KO B", "Lgr5Cre A", "Lgr5Cre B"), to = c("Setdb1KO", "Setdb1KO", "Lgr5Cre", "Lgr5Cre"))


##################################################################################################################################################################################################################


# QC and selecting cells
Lgr5Cre_MERGED[["percent.mt"]] <- PercentageFeatureSet(Lgr5Cre_MERGED, pattern = "^mt-")

# Violin Plot
pdf("violin_plot_Lgr5Cre.pdf", width=18, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Violin Plot with merged replicates
pdf("violin_plot_Lgr5Cre_merged.pdf", width=18, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
dev.off()

# Feature Scatter Plot
plot1 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("feauture_scatter.pdf", width=15, height=8)
plot1 + plot2
dev.off()

# Feature Scatter Plot
plot1 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
pdf("feauture_scatter_merged.pdf", width=15, height=8)
plot1 + plot2
dev.off()

# Filtering cells 
Lgr5Cre_MERGED <- subset(Lgr5Cre_MERGED, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# Normalize the data
Lgr5Cre_MERGED <- NormalizeData(Lgr5Cre_MERGED, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
Lgr5Cre_MERGED <- FindVariableFeatures(Lgr5Cre_MERGED, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Lgr5Cre_MERGED), 10)

# Plot variable features
pdf("variable_features_Lgr5Cre.pdf", width=15, height=8)
plot1 <- VariableFeaturePlot(Lgr5Cre_MERGED)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# Scaling the data
all.genes <- rownames(Lgr5Cre_MERGED)
Lgr5Cre_MERGED <- ScaleData(Lgr5Cre_MERGED, features = all.genes)

# Performing linear dimensional reduction
Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = VariableFeatures(object = Lgr5Cre_MERGED))
print(Lgr5Cre_MERGED[["pca"]], dims = 1:5, nfeatures = 5)

# VizDim Plot
pdf("Viz_DimPlot_PCA_Lgr5Cre.pdf", width=15, height=8)
VizDimLoadings(Lgr5Cre_MERGED, dims = 1:2, reduction = "pca")
dev.off()

# DimPlot
pdf("DimPlot_PCA_Lgr5Cre.pdf", width=15, height=8)
VizDimLoadings(Lgr5Cre_MERGED, dims = 1:2, reduction = "pca")
dev.off()

# Dim Heatmap
pdf("DimHeatmap_Lgr5Cre.pdf", width=15, height=8)
DimHeatmap(Lgr5Cre_MERGED, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# Determine the dimensionality of the dataset
Lgr5Cre_MERGED <- JackStraw(Lgr5Cre_MERGED, num.replicate = 100)
Lgr5Cre_MERGED <- ScoreJackStraw(Lgr5Cre_MERGED, dims = 1:20)

pdf("JackstrawPlot_Lgr5Cre.pdf", width=15, height=8)
JackStrawPlot(Lgr5Cre_MERGED, dims = 1:15)
dev.off()

pdf("ElbowPlot_Lgr5Cre.pdf", width=15, height=8)
ElbowPlot(Lgr5Cre_MERGED)
dev.off()

# Cluster the Cells
Lgr5Cre_MERGED <- FindNeighbors(Lgr5Cre_MERGED, dims = 1:10)
Lgr5Cre_MERGED <- FindClusters(Lgr5Cre_MERGED, resolution = 0.5)

# Run non linear dimensional reduction
reticulate::py_install(packages = 'umap-learn')
Lgr5Cre_MERGED <- RunUMAP(Lgr5Cre_MERGED, dims = 1:10)

# Umap plot
pdf("umapPlot_Lgr5Cre.pdf", width=15, height=8)
DimPlot(Lgr5Cre_MERGED, reduction = "umap")
dev.off()

# Umap plot merged
pdf("umapPlot_Lgr5Cre_merged.pdf", width=15, height=8)
DimPlot(Lgr5Cre_MERGED, reduction = "umap", group.by = "orig.ident")
dev.off()

# Split objects 
SplitObject(Lgr5Cre_MERGED, split.by = "ident")
n_cells=(FetchData(Lgr5Cre_MERGED, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

# Umap plot merged replicates per condition
pdf("umapPlot_Lgr5Cre_conditions.pdf", width=15, height=8)
DimPlot(Lgr5Cre_MERGED,label=TRUE, split.by="orig.ident") + NoLegend()
dev.off()

# Feauture plot top 10 genes
pdf("feature_plot_top_10_MostVariantGenes.pdf", width=15, height=8)
FeaturePlot(Lgr5Cre_MERGED, features = top10)
dev.off()

# Feature plot on quality metrics ()
pdf("feautre plot qc_metrics.pdf", width = 13, height = 8)
FeaturePlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = c("cornsilk1" ,"blue1"))
dev.off()

# Dot plot per top 10 genes
pdf("DotPlot_Lgr5cre.pdf", width=15, height=8)
DotPlot(Lgr5Cre_MERGED, cols = c("blue", "orange"), group.by = "orig.ident", features = top10)
dev.off()

# Dot plot per selected Markers
pdf("Top10_Markers.pdf", width=15, height=8)
DotPlot(Lgr5Cre_MERGED, features = c("Reg3g", "Gsdmc4", "Prss32", "Krt8", "Elf3", "Sis", "Fabp1", "Hnf4a", "Tmigd1", "Fabp6"), split.by = "orig.ident")
dev.off()


##################################################################################################################################################################################################################


# Extract metadata
md = Lgr5Cre_MERGED@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "seurat_clusters")]
md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")

# Extract and save all genes
write.table(all.genes, file = "all_genes.tsv", sep = "\t", row.names = TRUE)

# Find all markers distinguishing cluster 0 from the rest of the clusters
cluster0.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 0, ident.2 = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster1.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 1, ident.2 = c(0,2,3,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster2.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 2, ident.2 = c(0,1,3,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster3.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 3, ident.2 = c(0,1,2,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster4.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 4, ident.2 = c(0,1,2,3,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster5.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 5, ident.2 = c(0,1,2,3,4,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster6.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 6, ident.2 = c(0,1,2,3,4,5,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster7.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 7, ident.2 = c(0,1,2,3,4,5,6,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster8.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 8, ident.2 = c(0,1,2,3,4,5,6,7,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster9.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 9, ident.2 = c(0,1,2,3,4,5,6,7,8,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster10.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 10, ident.2 = c(0,1,2,3,4,5,6,7,8,9,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster11.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 11, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster12.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 12, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,11,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster13.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 13, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,11,12,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster14.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 14, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13), min.pct = 0.25, logfc.threshold = 0.37)

# Save files
write.table(cluster0.markers, file = "cluster0.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster1.markers, file = "cluster1.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster2.markers, file = "cluster2.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster3.markers, file = "cluster3.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster4.markers, file = "cluster4.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster5.markers, file = "cluster5.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster6.markers, file = "cluster6.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster7.markers, file = "cluster7.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster8.markers, file = "cluster8.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster9.markers, file = "cluster9.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster10.markers, file = "cluster10.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster11.markers, file = "cluster11.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster12.markers, file = "cluster12.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster13.markers, file = "cluster13.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster14.markers, file = "cluster14.markers.tsv", sep = "\t", row.names = TRUE)

# Mat code
#Lgr5Cre_MERGED.markers = FindAllMarkers(Lgr5Cre_MERGED, only.pos = F, min.pct = 0.25, logfc.threshold = 0.37)


##################################################################################################################################################################################################################


# Feature Plots based on most variable features per cluster

pdf("Feature plot cluster_0 Markers 1-5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Olfm4", "Gkn3", "Ifitm3", "Jaml", "2410006H16Rik"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_0 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Slc12a2", "Cd74", "Clca3b", "Pdgfa"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_1 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Ube2c", "Hmgb2", "Cenpa", "Cks2", "Birc5"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_1 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("H2afx", "Tubb4b", "Cenpf", "Tubb5", "Ccnb2"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_2 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Zg16", "Fcgbp", "Muc2", "Tff3", "Clca1"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_2 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Agr2", "Spink4", "Guca2a", "Ccl6", "Klk1"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_3 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Pclaf", "Stmn1", "Dut", "Lig1", "Pcna"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_3 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Tyms", "Hells", "Siva1", "Olfm4", "Dek"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_4 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Dmbt1", "Hspd1", "Hells", "Dut", "Hspe1"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_4 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("C1qbp", "Snhg4", "Ranbp1", "Pcna", "Tomm5"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_5 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Ube2c", "Birc5", "Rbp7", "Dmbt1", "H2afx"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_5 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Hmgb2", "Pclaf", "Hist1h1b", "Tubb4b", "Ube2s"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()


pdf("Feature plot cluster_6 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Agr2", "Selenom", "Spink4", "Tff3", "Guca2a"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_6 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Txndc5", "Fcgbp", "Smim6", "Fxyd3", "Muc2"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_7 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Rbp2", "Arg2", "Fabp1", "Smim24", "St3gal4"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_7 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Krt19", "Dmbt1", "Gna11", "Prap1", "Gpx1"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_8 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Defa30", "Gm14851", "AY61184", "Lyz1", "Defa24"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_8 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Defa21", "Defa29", "Itln1", "Defa36", "Defa26"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_9 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Reg1", "Rbp2", "Fabp1", "Apoa1", "Gsta1"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_9 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Fabp2", "Arg2", "Gstm3", "Adh6a", "Adh6a"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_10 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Gip", "Ghrl", "Cck", "Fabp5", "Sct"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_10 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Cdkn1c", "Rbp4", "Isl1", "Myl7", "Fam183b"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_11 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Rgs13", "Alox5ap", "Lrmp", "Cd24a", "Ltc4s"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_11 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Hck", "Fyb", "Kctd12", "Dclk1", "Adh1"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_12 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Fabp6", "Reg3b", "Reg3g", "Spink1", "Clec2h"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_12 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Slc51a", "Apoa1", "Lgals3", "2200002D01Rik", "Mgatc"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_13 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Chgb", "Chga", "Tac1", "Reg4", "Vim"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_13 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Tph1", "Gstt1", "Sct", "Pcsk1", "Hmgn3"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_14 Markers 1_5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Ccl5", "Gzma", "Cd3g", "Cd7", "Nkg7"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

pdf("Feature plot cluster_14 Markers 6_10.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Cd52", "Fcer1g", "AW112010", "Trbc2", "Tyrobp"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()

## Selected Markers
pdf("Feature plot Selected Markers.pdf", width=20, height=40)
FeaturePlot(Lgr5Cre_MERGED, features = c("Sis", "Fabp1", "Fabp6", "Muc2", "Lyz1", "Neurod1", "Avil" ), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()


##################################################################################################################################################################################################################


# DotPlots per Cluster

pdf("Enterocyte_Immature Markers.pdf", width=15, height=8)
p =  DotPlot(Lgr5Cre_MERGED,features = c("Reg3g", "Gsdmc4", "Prss32", "Krt8"),cols = c("blue","orange"))
p + ggtitle("Enterocyte Immature Markers") + RotatedAxis()
dev.off()

pdf("Enterocyte_Mature Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Elf3","Sis","Fabp1","Hnf4aos","Hnf4a", "Hnf4g", "Tmigd1","Fabp6","Slc51b","Slc51a","Mep1a","Fam151a","Naaladl1","Slc34a2","Plb1","Nudt4","Dpep1","Pmp22","Xpnpep2","Muc3","Neu1","Clec2h","Phgr1","Prss30","Aldob","Alpi","Apoa1","Apoa4","Lct"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Mature Markers")
dev.off()

pdf("Enterocyte Progenitor Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Ccnb1","Cdc20","Cenpa","Cdkn3","Cdc25c","Ccnb2","Kif22","Ube2c","Sapcd2","Rbp7","Ccna2","Aurka","Cdkn2d","Kif23","Nek2","Birc5","Plk1","Tacc3","Melk","Cps1"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Progenitor Markers")
dev.off()

pdf("Goblet Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Muc2","Spdef","Foxa1","Agr2","Spink4","Fcgbp","Tff3","Zg16","Clca1","Ccl6","Klk1","Tpsg1","Ccl9","Txndc5","Tspan13","Atoh1","Lrrc26","Clca3a1","Klf4"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Goblet Markers")
dev.off()

pdf("Paneth Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lyz1","Mmp7","Dll4","Sox9","Gfi1","Gm14851","Defa21","Defa22","Defa17","Defa24","Defa3","Mptx2","Ang4"),cols = c("blue","orange")) + RotatedAxis()
# Gene Gm15284 & Defa-rs1 excluded as there are not in the dataset 
p + ggtitle("Paneth Markers")
dev.off()

pdf("Enteroendocrine Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Gfi1","Neurog3","Neurod1","Chga","Chgb","Isl1","Arx","Pax6","Foxa2","Sst","Gck","Gcg","Tph1","Pyy","Gfra3","Cpe","Tac1","Fam183b", "Hmgn3","Cck","Fev","Gch1","Pcsk1n", "Bex2","Vwa5b2","Nkx2-2","Marcksl1","Neurod2","Insm1"),cols = c("blue","orange")) + RotatedAxis()
# Gene Ngfrap1 excluded. It can not be found in the dataset
p + ggtitle("Enteroendocrine Markers")
dev.off()

pdf("Stem Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lgr5","Ascl2","Olfm4","Prom1","Axin2","Fzd2","Fzd7","Lrp5","Lrp6","Notch1","Hes1","Smo","Yap1","Igfbp4","Bex1","Gkn3","Slc12a2"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Stem Markers")
dev.off()

pdf("Tuft Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Dclk1","Ptprc","Avil","Lrmp","Alox5ap","Rgs13","Sh2d6","Ltc4s","Hck","Cd24a","Trpm5","Kctd12","Aldh2","Il13ra1","Gng13","Tmem176a","Skap2","Ptpn6","Ly6g6f","Fyb","Adh1","Gfi1b","Il25"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Tuft Markers")
dev.off()

pdf("TA Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Stmn1", "Tubb5"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("TA Markers")
dev.off()

pdf("Reserved Stem Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lrig1", "Bmi1", "Hopx"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Reserved Stem Markers")
dev.off()

pdf("Necroptosis Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Zbp1", "Ripk1", "Ripk3", "Cxcl1", "Ccl20", "Tnf", "Csf1" ),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Necroptosis Markers")
dev.off()

pdf("Lgr5 Marker.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lgr5" ),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Lgr5 Marker")
dev.off()

pdf("GFP Marker.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("GFP" ), cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("GFP Marker")
dev.off()

pdf("Nature paper Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Mbl2","Ptprs","C4bp","Icam1","Ptpn22","Syk","Shmt2","Zfp809","Trim35","Bst2","Ltf","Mx1","Pml","Trim6","Adar","Mx2","S100a14","Gbp3","Ifitm3","Rnf135","Tmem173","Isg20","Jak3","Stat2","Oas2","Aim2","Oas1b","Dhx58","Oas1a","Eif2ak2","Zc3hav1","Ereg","Parp14","Zbp1","Cxcl16","Oasl2","Ddx58","Trim25","Rnf125","Vnn1","Pvr","Gbp4","Gbp10","Oasl1","Gbp8","Trim15","Cd55","Oas3","Isg15","Irf7","Trim12c","Ifit2","Ifit3","Trim30a","Nos2"), cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Nature Paper Markers")
dev.off()


##################################################################################################################################################################################################################


# Dot plots per condition

pdf("Enterocyte Immature Markers Grouped.pdf", width=15, height=8)
p =  DotPlot(Lgr5Cre_MERGED,features = c("Reg3g", "Gsdmc4", "Prss32", "Krt8"), group.by = "orig.ident" ,cols = c("blue","orange"))
p + ggtitle("Enterocyte Immature Markers") + RotatedAxis()
dev.off()

pdf("Enterocyte Mature Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Elf3","Sis","Fabp1","Hnf4aos","Hnf4a", "Hnf4g", "Tmigd1","Fabp6","Slc51b","Slc51a","Mep1a","Fam151a","Naaladl1","Slc34a2","Plb1","Nudt4","Dpep1","Pmp22","Xpnpep2","Muc3","Neu1","Clec2h","Phgr1","Prss30","Aldob","Alpi","Apoa1","Apoa4","Lct"), group.by = "orig.ident" ,cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Mature Markers")
dev.off()

pdf("Enterocyte Progenitor Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Ccnb1","Cdc20","Cenpa","Cdkn3","Cdc25c","Ccnb2","Kif22","Ube2c","Sapcd2","Rbp7","Ccna2","Aurka","Cdkn2d","Kif23","Nek2","Birc5","Plk1","Tacc3","Melk","Cps1"), group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Progenitor Markers")
dev.off()

pdf("Goblet Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Muc2","Spdef","Foxa1","Agr2","Spink4","Fcgbp","Tff3","Zg16","Clca1","Ccl6","Klk1","Tpsg1","Ccl9","Txndc5","Tspan13","Atoh1","Lrrc26","Clca3a1","Klf4"), group.by = "orig.ident",  cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Goblet Markers")
dev.off()

pdf("Paneth Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lyz1","Mmp7","Dll4","Sox9","Gfi1","Gm14851","Defa21","Defa22","Defa17","Defa24","Defa3","Mptx2","Ang4"), group.by = "orig.ident", group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
# Gene Gm15284 & Defa-rs1 excluded as there are not in the dataset 
p + ggtitle("Paneth Markers")
dev.off()

pdf("Enteroendocrine Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Gfi1","Neurog3","Neurod1","Chga","Chgb","Isl1","Arx","Pax6","Foxa2","Sst","Gck","Gcg","Tph1","Pyy","Gfra3","Cpe","Tac1","Fam183b", "Hmgn3","Cck","Fev","Gch1","Pcsk1n", "Bex2","Vwa5b2","Nkx2-2","Marcksl1","Neurod2","Insm1"),group.by = "orig.ident",  cols = c("blue","orange")) + RotatedAxis()
# Gene Ngfrap1 excluded. It can not be found in the dataset
p + ggtitle("Enteroendocrine Markers")
dev.off()

pdf("Stem Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lgr5","Ascl2","Olfm4","Prom1","Axin2","Fzd2","Fzd7","Lrp5","Lrp6","Notch1","Hes1","Smo","Yap1","Igfbp4","Bex1","Gkn3","Slc12a2"), group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Stem Markers")
dev.off()

pdf("Tuft Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Dclk1","Ptprc","Avil","Lrmp","Alox5ap","Rgs13","Sh2d6","Ltc4s","Hck","Cd24a","Trpm5","Kctd12","Aldh2","Il13ra1","Gng13","Tmem176a","Skap2","Ptpn6","Ly6g6f","Fyb","Adh1","Gfi1b","Il25"), group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Tuft Markers")
dev.off()

pdf("TA Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Stmn1", "Tubb5"),group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("TA Markers")
dev.off()

pdf("Reserved Stem Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lrig1", "Bmi1", "Hopx"),group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Reserved Stem Markers")
dev.off()

pdf("Necroptosis Markers Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Zbp1", "Ripk1", "Ripk3", "Cxcl1", "Ccl20", "Tnf", "Csf1" ),group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Necroptosis Markers")
dev.off()

pdf("Lgr5 Marker Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Lgr5" ), group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Lgr5 Marker")
dev.off()

pdf("GFP Marker Grouped.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("GFP"), group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("GFP Marker")
dev.off()

pdf("Nature paper Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Mbl2","Ptprs","C4bp","Icam1","Ptpn22","Syk","Shmt2","Zfp809","Trim35","Bst2","Ltf","Mx1","Pml","Trim6","Adar","Mx2","S100a14","Gbp3","Ifitm3","Rnf135","Tmem173","Isg20","Jak3","Stat2","Oas2","Aim2","Oas1b","Dhx58","Oas1a","Eif2ak2","Zc3hav1","Ereg","Parp14","Zbp1","Cxcl16","Oasl2","Ddx58","Trim25","Rnf125","Vnn1","Pvr","Gbp4","Gbp10","Oasl1","Gbp8","Trim15","Cd55","Oas3","Isg15","Irf7","Trim12c","Ifit2","Ifit3","Trim30a","Nos2"), group.by = "orig.ident", cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Nature Paper Markers")
dev.off()


##################################################################################################################################################################################################################


# Umap renamed cluster
pdf("Umap_renamed_Clusterd Markers.pdf", width=15, height=8)
Lgr5Cre_MERGED_Renamed=RenameIdents(Lgr5Cre_MERGED,  `0` = "Stem", `1` = "Enterocyte Progenitor I", `2` = "Goblet I", `3` = "Stem - TA", `4` = "?", `5` = "Enterocyte Progenitor II", `6` = "Goblet II", `7` = "Enterocyte Mature I", `8` = "Paneth", `9` = "Enterocyte Mature II", `10` = "Enderoendocrine I", `11` = "Tuft", `12` = "Enterocyte immature", `13` = "Enteroendocrine II", `14` = "Necroptosis")
DimPlot(Lgr5Cre_MERGED_Renamed, label = TRUE)
dev.off()


##################################################################################################################################################################################################################


# Split violin plots per markers
pdf("Violin Plot_Birc5.pdf", width=15, height=10)
plots <- VlnPlot(Lgr5Cre_MERGED, features = c(""), split.by = "orig.ident", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()


############################ Selected Markers


pdf("Violin Plot_Sis.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Sis"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Fabp1.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Fabp1"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Fabp6.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Fabp6"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Muc2.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Muc2"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Lyz1.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Lyz1"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Neurod1.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Neurod1"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Avil.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Avil"), split.by = "orig.ident")
dev.off()


############################ For cluster 0


pdf("Violin Plot_Olfm4.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Olfm4"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Gkn3.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Gkn3"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Ifitm3.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Ifitm3"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Jaml.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Jaml"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_2410006H16Rik.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("2410006H16Rik"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Slc12a2.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Slc12a2"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Cd74.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Cd74"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Clca3b.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Clca3b"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Pdgfa.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Pdgfa"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Gas5.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Gas5"), split.by = "orig.ident")
dev.off()


############################  For Cluster 3


pdf("Violin Plot_Pclaf.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Pclaf"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Stmn1.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Stmn1"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Dut.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Dut"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Lig1.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Lig1"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Pcna.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Pcna"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Tyms.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Tyms"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Hells.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Hells"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Siva1.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Siva1"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Olfm4.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Olfm4"), split.by = "orig.ident")
dev.off()

pdf("Violin Plot_Dek.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Dek"), split.by = "orig.ident")
dev.off()

# Split Ridge plots per markers

pdf("RidgePlot_Olfm4.pdf", width=15, height=10)
RidgePlot(Lgr5Cre_MERGED_Renamed, features = c("Olfm4"), ncol = 2)
dev.off()


##################################################################################################################################################################################################################


# Cell cycle

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Seurat/cell_cycle/", header = TRUE, as.is = TRUE, row.names = 1)

s.genes <- cc.genes$s.genes
s.genes <- tolower(s.genes)
s.genes = str_to_title(s.genes)

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- tolower(g2m.genes)
g2m.genes = str_to_title(g2m.genes)

# Assign Cell-Cycle scores
Lgr5Cre_MERGED <- CellCycleScoring(Lgr5Cre_MERGED, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf("PCA on Cell Cycle Genes.pdf", width=15, height=8)
Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = c(s.genes, g2m.genes))
DimPlot(Lgr5Cre_MERGED)
dev.off()

# Regress out cell cycle scores during data scaling
Lgr5Cre_MERGED <- ScaleData(Lgr5Cre_MERGED, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Lgr5Cre_MERGED))

Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = VariableFeatures(Lgr5Cre_MERGED), nfeatures.print = 10)

Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = c(s.genes, g2m.genes))
DimPlot(Lgr5Cre_MERGED)


##################################################################################################################################################################################################################


# Create dot plots - on specific clusters and selected markers

pdf("Setdb1KO_Enterocyte Mature top_markers.pdf", width=15, height=8)
p =  DotPlot(Setdb1KO,features = c("Smim24","St3gal4","2200002D01Rik","Krt19","Dmbt1","Gna11","Prap1","Gpx1","Serpinb6a","Reg1","Rbp2","Fabp1","Apoa1","Gsta1","Fabp2","Arg2","Gstm3","Adh6a","Sis","Apoc3","S100g"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("Setdb1KO") + RotatedAxis()
dev.off()

pdf("WT_Enterocyte Mature top_markers.pdf", width=15, height=8)
p =  DotPlot(WT,features = c("Smim24","St3gal4","2200002D01Rik","Krt19","Dmbt1","Gna11","Prap1","Gpx1","Serpinb6a","Reg1","Rbp2","Fabp1","Apoa1","Gsta1","Fabp2","Arg2","Gstm3","Adh6a","Sis","Apoc3","S100g"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("WT") + RotatedAxis()
dev.off()


pdf("Setdb1KO_Enterocyte Progenitor top_markers.pdf", width=15, height=8)
p =  DotPlot(Setdb1KO,features = c("Ube2c","Hmgb2","Cenpa","Cks2","H2afx","Tubb4b","Cenpf","Tubb5","Ccnb2","Cdca3","Cdca8","Birc5","Rbp7","Dmbt1" ), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("Setdb1KO") + RotatedAxis()
dev.off()

pdf("WT_Enterocyte Progenitor top_markers.pdf", width=15, height=8)
p =  DotPlot(WT,features = c("Ube2c","Hmgb2","Cenpa","Cks2","H2afx","Tubb4b","Cenpf","Tubb5","Ccnb2","Cdca3","Cdca8","Birc5","Rbp7","Dmbt1" ), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("WT") + RotatedAxis()
dev.off()

pdf("Setdb1KO_Tuf top_markers.pdf", width=15, height=8)
p =  DotPlot(Setdb1KO,features = c("Rgs13","Alox5ap","Lrmp","Cd24a","Ltc4s","Hck","Fyb","Kctd12","Dclk1","Adh1","Krt18","Aldh2","Avil","Reep5","Sh2d6","Tuba1a","Espn","Tmem176b","Snrnp25"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("Setdb1KO") + RotatedAxis()
dev.off()

pdf("WT_Tuf top_markers.pdf", width=15, height=8)
p =  DotPlot(WT,features = c("Rgs13","Alox5ap","Lrmp","Cd24a","Ltc4s","Hck","Fyb","Kctd12","Dclk1","Adh1","Krt18","Aldh2","Avil","Reep5","Sh2d6","Tuba1a","Espn","Tmem176b","Snrnp25"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("WT") + RotatedAxis()
dev.off()

pdf("Setdb1KO_Paneth top_markers.pdf", width=15, height=8)
p =  DotPlot(Setdb1KO,features = c("Defa30","Gm14851","AY761184","Lyz1","Defa24","Defa22","Defa21","Defa29","Itln1","Defa36","Defa26","Defa23","Ang4","Mptx2","Defa17","Clps","Defa5","Gm15292","Gm15293"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("Setdb1KO") + RotatedAxis()
dev.off()

pdf("WT_Paneth top_markers.pdf", width=15, height=8)
p =  DotPlot(WT,features = c("Defa30","Gm14851","AY761184","Lyz1","Defa24","Defa22","Defa21","Defa29","Itln1","Defa36","Defa26","Defa23","Ang4","Mptx2","Defa17","Clps","Defa5","Gm15292","Gm15293"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("WT") + RotatedAxis()
dev.off()

pdf("Setdb1KO_Goblet top_markers.pdf", width=15, height=8)
p =  DotPlot(Setdb1KO,features = c("Zg16","Fcgbp","Muc2","Tff3","Clca1","Agr2","Spink4","Guca2a","Ccl6","Klk1","S100a6","Mptx1","Rep15","Selenom","Txndc5","Fxyd3","Ramp1","Tmsb10","Tpsg1"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("Setdb1KO") + RotatedAxis()
dev.off()

pdf("WT_Goblet top_markers.pdf", width=15, height=8)
p =  DotPlot(WT,features = c("Zg16","Fcgbp","Muc2","Tff3","Clca1","Agr2","Spink4","Guca2a","Ccl6","Klk1","S100a6","Mptx1","Rep15","Selenom","Txndc5","Fxyd3","Ramp1","Tmsb10","Tpsg1"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("WT") + RotatedAxis()
dev.off()


##################################################################################################################################################################################################################


# Create metadata from main object. e.g export WT and Setdb1KO from Merged object and create two new objects

split_obj = SplitObject(Lgr5Cre_MERGED, split.by = "orig.ident")

# See how the objects are named
View(split_obj)

Setdb1KO = split_obj[["Setdb1KO"]]

WT = split_obj[["Lgr5Cre"]]


##################################################################################################################################################################################################################


# Calculate averaged expression values for each identity class
avg_exp_Setdb1KO = (AverageExpression(object = Setdb1KO))
avg_exp_WT = (AverageExpression(object = WT))

write.table(avg_exp_Setdb1KO, file = "avg_exp_Setdb1KO.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(avg_exp_WT, file = "avg_exp_WT.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


##################################################################################################################################################################################################################


# Create violin plots for setdb1KO and WT metadata (merged normalized-processed and then split) for nfeatures 
# ncount per cluster.

pdf("Setdb1KO-ViolinPlot_per_Cluster.pdf", width=18, height=10)
VlnPlot(Setdb1KO, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

pdf("WT-ViolinPlot_per_Cluster.pdf", width=18, height=10)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

pdf("Merged_ViolinPlot_per_Cluster.pdf", width=18, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA"), split.by = "orig.ident", ncol = 2)
dev.off()

# Keep cells with only 2000 reads
#Lgr5Cre_MERGED <- subset(Lgr5Cre_MERGED, subset = nCount_RNA > 2000)


##################################################################################################################################################################################################################


# Subset clusters and recluster
Lgr5Cre_0.1.3.4.5=subset(Lgr5Cre_MERGED, idents = c(0,1,3,4,5))

top10_Lgr5Cre_0.1.3.4.5 <- head(VariableFeatures(Lgr5Cre_0.1.3.4.5), 10)

Lgr5Cre_0.1.3.4.5 <- RunPCA(Lgr5Cre_0.1.3.4.5, features = VariableFeatures(object = Lgr5Cre_0.1.3.4.5))
print(Lgr5Cre_0.1.3.4.5[["pca"]], dims = 1:5, nfeatures = 5)

Lgr5Cre_0.1.3.4.5=FindNeighbors(Lgr5Cre_0.1.3.4.5,,dims = 1:10)
Lgr5Cre_0.1.3.4.5 <- FindClusters(Lgr5Cre_0.1.3.4.5, resolution = 0.5)

Lgr5Cre_0.1.3.4.5=RunUMAP(Lgr5Cre_0.1.3.4.5, dims = 1:10)


# Umap plot
pdf("umapPlot_Lgr5Cre_0.1.3.4.5.pdf", width=15, height=8)
DimPlot(Lgr5Cre_0.1.3.4.5, reduction = "umap")
dev.off()


# DOTPLOTS

pdf("Enterocyte_Immature Markers.pdf", width=15, height=8)
p =  DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Reg3g", "Gsdmc4", "Prss32", "Krt8"),cols = c("blue","orange"))
p + ggtitle("Enterocyte Immature Markers") + RotatedAxis()
dev.off()

pdf("Enterocyte_Mature Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Elf3","Sis","Fabp1","Hnf4aos","Hnf4a", "Hnf4g", "Tmigd1","Fabp6","Slc51b","Slc51a","Mep1a","Fam151a","Naaladl1","Slc34a2","Plb1","Nudt4","Dpep1","Pmp22","Xpnpep2","Muc3","Neu1","Clec2h","Phgr1","Prss30","Aldob","Alpi","Apoa1","Apoa4","Lct"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Mature Markers")
dev.off()

pdf("Enterocyte Progenitor Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Ccnb1","Cdc20","Cenpa","Cdkn3","Cdc25c","Ccnb2","Kif22","Ube2c","Sapcd2","Rbp7","Ccna2","Aurka","Cdkn2d","Kif23","Nek2","Birc5","Plk1","Tacc3","Melk","Cps1"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Progenitor Markers")
dev.off()

pdf("Goblet Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Muc2","Spdef","Foxa1","Agr2","Spink4","Fcgbp","Tff3","Zg16","Clca1","Ccl6","Klk1","Tpsg1","Ccl9","Txndc5","Tspan13","Atoh1","Lrrc26","Clca3a1","Klf4"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Goblet Markers")
dev.off()

pdf("Paneth Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Lyz1","Mmp7","Dll4","Sox9","Gfi1","Gm14851","Defa21","Defa22","Defa17","Defa24","Defa3","Mptx2","Ang4"),cols = c("blue","orange")) + RotatedAxis()
# Gene Gm15284 & Defa-rs1 excluded as there are not in the dataset 
p + ggtitle("Paneth Markers")
dev.off()

pdf("Enteroendocrine Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Gfi1","Neurog3","Neurod1","Chga","Chgb","Isl1","Arx","Pax6","Foxa2","Sst","Gck","Gcg","Tph1","Pyy","Gfra3","Cpe","Tac1","Fam183b", "Hmgn3","Cck","Fev","Gch1","Pcsk1n", "Bex2","Vwa5b2","Nkx2-2","Marcksl1","Neurod2","Insm1"),cols = c("blue","orange")) + RotatedAxis()
# Gene Ngfrap1 excluded. It can not be found in the dataset
p + ggtitle("Enteroendocrine Markers")
dev.off()

pdf("Stem Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Lgr5","Ascl2","Olfm4","Prom1","Axin2","Fzd2","Fzd7","Lrp5","Lrp6","Notch1","Hes1","Smo","Yap1","Igfbp4","Bex1","Gkn3","Slc12a2"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Stem Markers")
dev.off()

pdf("Tuft Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Dclk1","Ptprc","Avil","Lrmp","Alox5ap","Rgs13","Sh2d6","Ltc4s","Hck","Cd24a","Trpm5","Kctd12","Aldh2","Il13ra1","Gng13","Tmem176a","Skap2","Ptpn6","Ly6g6f","Fyb","Adh1","Gfi1b","Il25"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Tuft Markers")
dev.off()

pdf("TA Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Stmn1", "Tubb5"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("TA Markers")
dev.off()

pdf("Cdk2 Marker Clusters.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_0.1.3.4.5,features = c("Cdk2"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Cdk2 Marker")
dev.off()


# Rename UMAP


# Umap renamed cluster
pdf("Umap_renamed_Clusterd Markers.pdf", width=10, height=6)
Lgr5Cre_0.1.3.4.5_Renamed=RenameIdents(Lgr5Cre_0.1.3.4.5,  `0` = "Stem I", `1` = "Stem II", `2` = "TA", `3` = "Progenitor Ι", `4` = "S phase cells", `5` = "Progenitor II", `6` = "Progenitor III", `7` = "Goblet-Paneth")
DimPlot(Lgr5Cre_0.1.3.4.5_Renamed, label = TRUE)
dev.off()


# Umap plot merged
pdf("umapPlot_Lgr5Cre_0.1.3.4.5.pdf", width=15, height=8)
DimPlot(Lgr5Cre_0.1.3.4.5, reduction = "umap", group.by = "orig.ident")
dev.off()

# Split objects 
SplitObject(Lgr5Cre_0.1.3.4.5, split.by = "ident")
n_cells=(FetchData(Lgr5Cre_0.1.3.4.5, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

# Umap plot merged replicates per condition
pdf("umapPlot_Lgr5Cre_0.1.3.4.5.conditions.pdf", width=15, height=8)
DimPlot(Lgr5Cre_0.1.3.4.5,label=TRUE, split.by="orig.ident") + NoLegend()
dev.off()

md = Lgr5Cre_0.1.3.4.5@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "seurat_clusters")]
md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")


##################################################################################################################################################################################################################


# Cell cycle
s.genes <- cc.genes$s.genes
s.genes <- tolower(s.genes)
s.genes = str_to_title(s.genes)

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- tolower(g2m.genes)
g2m.genes = str_to_title(g2m.genes)


# Assign Cell-Cycle scores
Lgr5Cre_0.1.3.4.5 <- CellCycleScoring(Lgr5Cre_0.1.3.4.5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf("PCA on Cell Cycle Genes.pdf", width=10, height=6)
Lgr5Cre_0.1.3.4.5 <- RunPCA(Lgr5Cre_0.1.3.4.5, features = c(s.genes, g2m.genes))
DimPlot(Lgr5Cre_0.1.3.4.5)
dev.off()


##################################################################################################################################################################################################################

# Calculate Coefficient of Variability in gene expression

# subset cluster 0
Lgr5Cre_012345=subset(Lgr5Cre_MERGED, idents = ("0", "1", "2", "3", "4", "5"))

# create two matrix file for gene count per cell (KO and WT)
# WT
Lgr5Cre_012345_WT=subset(Lgr5Cre_012345, subset = orig.ident == "Lgr5Cre")
# KO
Lgr5Cre_012345_KO=subset(Lgr5Cre_012345, subset = orig.ident == "Setdb1KO")


# Save file with gene expression per cell type per condition
# WT
write.table(Lgr5Cre_5_WT@assays[["RNA"]]@counts, file='WT_Gene_Count_per_Cell.tmp', quote=FALSE, sep='\t', col.names = TRUE)
# KO
write.table(Lgr5Cre_5_KO@assays[["RNA"]]@counts, file='KO_Gene_Count_per_Cell.tmp', quote=FALSE, sep='\t', col.names = TRUE)

# read the filtered file
KO_matrix_Filtered = read.table("/media/dimbo/10T/data/talianidis_data/scRNA_seq/Talianidis_GFP_2_fixed/analysis/Seurat/coefficient_variability/cluster_5/KO_Gene_Count_per_Cell.tmp", sep = "\t", header = TRUE, row.names = 1 )
WT_matrix_Filtered = read.table("/media/dimbo/10T/data/talianidis_data/scRNA_seq/Talianidis_GFP_2_fixed/analysis/Seurat/coefficient_variability/cluster_5/WT_Gene_Count_per_Cell.tmp", sep = "\t", header = TRUE, row.names = 1 ) 

# calculate coefficient variability and add results in the file
KO_matrix_Filtered$cv  <- apply(KO_matrix_Filtered, 1,  function(x) sd(x) / mean(x) * 100)
WT_matrix_Filtered$cv  <- apply(WT_matrix_Filtered, 1,  function(x) sd(x) / mean(x) * 100)

# send back to bash for file manipulation and filtering
write.table(KO_matrix_Filtered, file = "KO.tmp", row.names = FALSE, sep = "\t")
write.table(WT_matrix_Filtered, file = "WT.tmp", row.names = FALSE, sep = "\t")

coefficient_variation = read.table("/media/dimbo/10T/data/talianidis_data/scRNA_seq/Talianidis_GFP_2_fixed/analysis/Seurat/coefficient_variability/cluster_5/coefficient_variation.tsv", sep = "\t", header = TRUE)
colnames(coefficient_variation) = c("Setdb1KO", "Lgr5Cre")

pdf("boxplot_cv.pdf")
colors <- c("blue","orange")
boxplot(coefficient_variation, ylab = "coefficient variation", col = colors, main = "Progenitor Cells - Cluster 5", outline = FALSE)
dev.off()


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

##################################################################################################################################################################################################################


# Trajectory analysis

monocle_Lgr5Cre_MERGED <- as.cell_data_set(Lgr5Cre_MERGED)
monocle_Lgr5Cre_MERGED <- cluster_cells(cds = monocle_Lgr5Cre_MERGED, reduction_method = "UMAP")
monocle_Lgr5Cre_MERGED <- learn_graph(monocle_Lgr5Cre_MERGED, use_partition = TRUE)
monocle_Lgr5Cre_MERGED <- order_cells(monocle_Lgr5Cre_MERGED,reduction_method = "UMAP")


pdf("trajectory_Lgr5Cre_Merged.pdf", width = 14, height = 8)
plot_cells(monocle_Lgr5Cre_MERGED,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = FALSE,
           cell_size = 0.55)
dev.off()

# SubClusters

monocle_Lgr5Cre_subset <- as.cell_data_set(Lgr5Cre_0.1.3.4.5)
monocle_Lgr5Cre_subset <- cluster_cells(cds = monocle_Lgr5Cre_subset, reduction_method = "UMAP")
monocle_Lgr5Cre_subset <- learn_graph(monocle_Lgr5Cre_subset, use_partition = TRUE)
monocle_Lgr5Cre_subset <- order_cells(monocle_Lgr5Cre_subset,reduction_method = "UMAP")


pdf("trajectory_MERGED_subset.pdf", width = 14, height = 8)
plot_cells(monocle_Lgr5Cre_subset,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = FALSE,
           cell_size = 0.55,
           )
dev.off()

###############################################

# Trajectroy analysis on split conditions
ko = subset(x = Lgr5Cre_MERGED, subset = orig.ident == "Setdb1KO")
wt = subset(x = Lgr5Cre_MERGED, subset = orig.ident == "Lgr5Cre")


monocle_ko <- as.cell_data_set(ko)
monocle_ko <- cluster_cells(cds = monocle_ko, reduction_method = "UMAP")
monocle_ko <- learn_graph(monocle_ko, use_partition = TRUE)
monocle_ko <- order_cells(monocle_ko,reduction_method = "UMAP")


pdf("trajectory_KO.pdf", width = 14, height = 8)
plot_cells(monocle_ko,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = FALSE,
           cell_size = 0.55)
dev.off()


monocle_wt <- as.cell_data_set(wt)
monocle_wt <- cluster_cells(cds = monocle_wt, reduction_method = "UMAP")
monocle_wt <- learn_graph(monocle_wt, use_partition = TRUE)
monocle_wt <- order_cells(monocle_wt,reduction_method = "UMAP")


pdf("trajectory_WT.pdf", width = 14, height = 8)
plot_cells(monocle_wt,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = FALSE,
           cell_size = 0.55)
dev.off()



###############################################################################

# SubClusters

# subset clusters for KO
KO_012345=subset(x = Lgr5Cre_0.1.3.4.5, subset = orig.ident == "Setdb1KO")
WT_012345=subset(x = Lgr5Cre_0.1.3.4.5, subset = orig.ident == "Lgr5Cre")

# Setdb1KO
monocle_ko <- as.cell_data_set(KO_012345)
monocle_ko <- cluster_cells(cds = monocle_ko, reduction_method = "UMAP")
monocle_ko <- learn_graph(monocle_ko, use_partition = TRUE)
monocle_ko <- order_cells(monocle_ko,reduction_method = "UMAP")

pdf("trajectory_Subset_KO.pdf", width = 14, height = 8)
plot_cells(monocle_ko,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = FALSE,
           cell_size = 0.55)
dev.off()

# Lgr5Cre
monocle_wt <- as.cell_data_set(WT_012345)
monocle_wt <- cluster_cells(cds = monocle_wt, reduction_method = "UMAP")
monocle_wt <- learn_graph(monocle_wt, use_partition = TRUE)
monocle_wt <- order_cells(monocle_wt,reduction_method = "UMAP")

pdf("trajectory_Subset_WT.pdf", width = 14, height = 8)
plot_cells(monocle_wt,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = FALSE,
           cell_size = 0.55)
dev.off()

############################################################################################

# Identification of highly variable features between clusters (Stem I and Stem II)

# Subset clusters and recluster
Lgr5Cre_0.1.3.4.5=subset(Lgr5Cre_MERGED, idents = c(0,1,3,4,5))

top10_Lgr5Cre_0.1.3.4.5 <- head(VariableFeatures(Lgr5Cre_0.1.3.4.5), 10)

Lgr5Cre_0.1.3.4.5 <- RunPCA(Lgr5Cre_0.1.3.4.5, features = VariableFeatures(object = Lgr5Cre_0.1.3.4.5))
print(Lgr5Cre_0.1.3.4.5[["pca"]], dims = 1:5, nfeatures = 5)

Lgr5Cre_0.1.3.4.5=FindNeighbors(Lgr5Cre_0.1.3.4.5,,dims = 1:10)
Lgr5Cre_0.1.3.4.5 <- FindClusters(Lgr5Cre_0.1.3.4.5, resolution = 0.5)

Lgr5Cre_0.1.3.4.5=RunUMAP(Lgr5Cre_0.1.3.4.5, dims = 1:10)


Lgr5Cre_0.1 = subset(Lgr5Cre_0.1.3.4.5, idents = c(0,1))

Lgr5Cre_0.1_Renamed=RenameIdents(Lgr5Cre_0.1,  `0` = "Stem II", `1` = "Stem I")
pdf("UMAP_STEM.pdf", width = 14, height = 8)
DimPlot(Lgr5Cre_0.1_Renamed)
dev.off()
# Subset Stem clusters
Lgr5Cre_0.1.3.4.5=subset(Lgr5Cre_MERGED, idents = c(0,1,3,4,5))


# Identification of highly variable features between Stem I and Stem II
Lgr5Cre_0.1 <- FindVariableFeatures(Lgr5Cre_0.1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Lgr5Cre_0.1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Lgr5Cre_0.1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("top_10_mvf_StemI_StemII.pdf", width = 14, height = 8)
plot2
dev.off()

# Plot violin in Reg3b and Reg3g
pdf("Violin_Reg3b_Reg3g.pdf", width = 14, height = 8)
VlnPlot(Lgr5Cre_0.1, features = c("Reg3b", "Reg3g"), slot = "counts", log = TRUE)
dev.off()

# Plot feature in Reg3b and Reg3g
pdf("featureplot_Reg3b_Reg3g.pdf", width = 14, height = 8)
FeaturePlot(Lgr5Cre_0.1, features = c("Reg3b", "Reg3g"), cols = c("light yellow", "blue"))
dev.off()

# Plot feature in Cd74
pdf("featureplot_Dd74.pdf", width = 14, height = 8)
FeaturePlot(Lgr5Cre_0.1, features = c("Cd74"), cols = c("light yellow", "blue"))
dev.off()


#######################################################################

# Dotplots in Stem I and Stem II from subset Lgr5Cre_0.1.3.4.5 (based on most variable features from
# stem I and stem II)


pdf("Stem I and II Markers.pdf", width=14, height=8)
p =  DotPlot(Lgr5Cre_0.1,features = c("Reg3b", "Reg3g", "Atf3", "Itln1", "S100a6", "Pclaf", "Fabp1",  "Bglap3", "Acta1") ,cols = c("blue","orange"))
p + ggtitle("Stem I vs Stem II") + RotatedAxis()
dev.off()


pdf("Vln Plot Cd74.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_0.1, features = c("Cd74"))
dev.off()




