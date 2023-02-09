library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr) 
library(data.table)
library(cowplot)
library(stringr)
set.seed(1234)

library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm9)


################################################################################################################################

# Pre-processing workflow

# WT
counts.WT <- Read10X_h5("/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/MUC23141/outs/filtered_peak_bc_matrix.h5")

metadata.WT <- read.csv(
  file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/MUC23141/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

WT_assay <- CreateChromatinAssay(
  counts = counts.WT,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/MUC23141/outs/fragments.tsv.gz',
  min.cells = 1
)

# Create Seurat object
ATAC_WT <- CreateSeuratObject(
  counts = WT_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata.WT
)

# KO
counts.KO <- Read10X_h5("/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/MUC23142/outs/filtered_peak_bc_matrix.h5")

metadata.KO <- read.csv(
  file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/MUC23142/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

KO_assay <- CreateChromatinAssay(
  counts = counts.KO,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Cellranger/02 Cellranger/cr_count/MUC23142/outs/fragments.tsv.gz',
  min.cells = 1
)

# Create Seurat object
ATAC_KO <- CreateSeuratObject(
  counts = KO_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata.KO
)


################################################################################################################################


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add information to identify dataset of origin
ATAC_WT$dataset <- 'WT'
ATAC_KO$dataset <- 'KO'

# MERGE OBJECTS
ATAC_MERGED <- merge(
  x = ATAC_WT,
  y = ATAC_KO,
  add.cell.ids = c("Lgr5Cre_ATAC", "Setdb1KO_ATAC")
)

################################################################################################################################


# add the gene information to the object
Annotation(ATAC_MERGED) <- annotations

# Computing QC metrics
ATAC_MERGED <- NucleosomeSignal(object = ATAC_MERGED)

# TSSs (Enrichment of Tn5 integration events at transcriptional start sites)
ATAC_MERGED <- TSSEnrichment(ATAC_MERGED, fast = FALSE)

pdf("TSSPlot_MERGED.pdf",
    width = 15,
    height = 8)
ATAC_MERGED$high.tss <- ifelse(ATAC_MERGED$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(ATAC_MERGED, group.by = 'high.tss') + NoLegend()
dev.off()

ATAC_MERGED$pct_reads_in_peaks <- ATAC_MERGED$peak_region_fragments / ATAC_MERGED$passed_filters * 100
ATAC_MERGED$blacklist_ratio <- ATAC_MERGED$blacklist_region_fragments / ATAC_MERGED$peak_region_fragments

pdf("VlnPlot_MERGED.pdf",
    width = 15,
    height = 8)
VlnPlot(
  object = ATAC_MERGED,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()


################################################################################################################################

# Remove cells that are outliers from QC metrics
ATAC_MERGED.filt <- subset(
  x = ATAC_MERGED,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 30 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
) # from 1714 cells dropped to 1305

# Normalization and linear dimensional reduction
ATAC_MERGED.filt <- RunTFIDF(ATAC_MERGED.filt)
ATAC_MERGED.filt <- FindTopFeatures(ATAC_MERGED.filt, min.cutoff = 'q0')
ATAC_MERGED.filt <- RunSVD(object = ATAC_MERGED.filt)

# The first LSI component often captures sequencing depth (technical variation) rather than biological variation. 
# If this is the case, the component should be removed from downstream analysis. 

pdf("DepthCor_MERGED.pdf",
    width = 15,
    height = 8)
DepthCor(ATAC_MERGED.filt)
dev.off()


################################################################################################################################


# Non-linear dimension reduction and clustering
ATAC_MERGED.filt <- RunUMAP(
  object = ATAC_MERGED.filt,
  reduction = 'lsi',
  dims = 2:25
)

ATAC_MERGED.filt <- FindNeighbors(
  object = ATAC_MERGED.filt,
  reduction = 'lsi',
  dims = 2:25
)

ATAC_MERGED.filt <- FindClusters(
  object = ATAC_MERGED.filt,
  algorithm = 3,
  resolution = 1.5,
  verbose = FALSE
)

pdf("UMAP_ATAC_MERGED.pdf",
    width = 15,
    height = 8)
DimPlot(object = ATAC_MERGED.filt, label = TRUE) + NoLegend()
dev.off()


################################################################################################################################

# Create a gene activity matrix

# compute gene activities
gene.activities <- GeneActivity(ATAC_MERGED.filt) # Compute counts per cell in gene body and promoter region

# add the gene activity matrix to the Seurat object as a new assay
ATAC_MERGED.filt[['RNA']] <- CreateAssayObject(counts = gene.activities)
ATAC_MERGED.filt <- NormalizeData(
  object = ATAC_MERGED.filt,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC_MERGED.filt$nCount_RNA)
)

pdf("FeaturePlot_Muc2.pdf",
    width = 10,
    height = 4)

DefaultAssay(ATAC_MERGED.filt) <- 'RNA'
FeaturePlot(
  object = ATAC_MERGED.filt,
  features = c("Muc2"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  col = c("cornsilk1", "blue")
)
dev.off()


################################################################################################################################

# Load the pre-processed scRNA-seq data
Lgr5Cre_MERGED <- readRDS("/media/dimbo/10T/data/talianidis_data/scRNA_seq/Talianidis_GFP_2_fixed/analysis/Seurat/Lgr5Cre_MERGED.rds")

# Classify cells based on an scRNA-seq experiment from the same biological system
Lgr5Cre_MERGED <- FindVariableFeatures(
  object = Lgr5Cre_MERGED,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = Lgr5Cre_MERGED,
  query = ATAC_MERGED.filt,
  reduction = 'cca',
  dims = 1:25
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Lgr5Cre_MERGED$seurat_clusters,
  weight.reduction = ATAC_MERGED.filt[['lsi']],
  dims = 2:25
)

# Rename in transfer labels clusters, from 0,1,2... to Stem,Progenitor etc.
predicted.labels$predicted.id[predicted.labels$predicted.id == 0] <- "Stem"
predicted.labels$predicted.id[predicted.labels$predicted.id == 1] <- "Enterocyte Progenitor I"
predicted.labels$predicted.id[predicted.labels$predicted.id == 2] <- "Goblet I"
predicted.labels$predicted.id[predicted.labels$predicted.id == 3] <- "Stem"
predicted.labels$predicted.id[predicted.labels$predicted.id == 4] <- "S phase Cells"
predicted.labels$predicted.id[predicted.labels$predicted.id == 5] <- "Enterocyte Progenitor II"
predicted.labels$predicted.id[predicted.labels$predicted.id == 6] <- "Goblet II"
predicted.labels$predicted.id[predicted.labels$predicted.id == 7] <- "Entterocyte Mature I"
predicted.labels$predicted.id[predicted.labels$predicted.id == 8] <- "Paneth"
predicted.labels$predicted.id[predicted.labels$predicted.id == 9] <- "Enterocyte Mature II"
predicted.labels$predicted.id[predicted.labels$predicted.id == 10] <- "Enderoendocrine I"
predicted.labels$predicted.id[predicted.labels$predicted.id == 11] <- "Tuft"
predicted.labels$predicted.id[predicted.labels$predicted.id == 12] <- "Enterocyte Immature"
predicted.labels$predicted.id[predicted.labels$predicted.id == 13] <- "Enteroendocrine II"
predicted.labels$predicted.id[predicted.labels$predicted.id == 14] <- "Necroptosis"



ATAC_MERGED.filt <- AddMetaData(object = ATAC_MERGED.filt, metadata = predicted.labels)

pdf("UMAP_scRNA-scATAC.pdf", 
    width = 15,
    height = 8)
plot1 <- DimPlot(Lgr5Cre_MERGED, label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(ATAC_MERGED.filt, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2
dev.off()


################################################################################################################################

#  Find differentially accessible peaks between clusters

#switch back to working with peaks instead of gene activities
DefaultAssay(ATAC_MERGED.filt) <- 'peaks'

cluster0_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("0"), 
  ident.2 = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster0_peaks, file = "cluster0_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)

cluster1_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("1"), 
  ident.2 = c("0", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster1_peaks, file = "cluster1_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)

cluster2_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("2"), 
  ident.2 = c("0", "1", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster2_peaks, file = "cluster2_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)


cluster3_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("3"), 
  ident.2 = c("0", "1", "2", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster3_peaks, file = "cluster3_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)

cluster4_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("4"), 
  ident.2 = c("0", "1", "2", "3", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster4_peaks, file = "cluster4_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)

cluster5_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("5"), 
  ident.2 = c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster5_peaks, file = "cluster5_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)



cluster6_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("6"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "7", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster6_peaks, file = "cluster6_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)


cluster7_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("7"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "8", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster7_peaks, file = "cluster7_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)


cluster8_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("8"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "7", "9", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster8_peaks, file = "cluster8_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)


cluster9_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("9"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster9_peaks, file = "cluster9_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)



cluster10_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("10"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster10_peaks, file = "cluster10_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)



cluster11_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("11"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "12", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster11_peaks, file = "cluster11_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)



cluster12_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("12"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "13", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster12_peaks, file = "cluster12_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)


cluster13_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("13"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "14"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster12_peaks, file = "cluster13_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)



cluster14_peaks <- FindMarkers(
  object = ATAC_MERGED.filt,
  ident.1 = c("14"), 
  ident.2 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.table(cluster12_peaks, file = "cluster14_peaks.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)


################################################################################################################################

# Violin Plots on accessible peaks between clusters
plot1 <- VlnPlot(
  object = ATAC_MERGED.filt,
  features = rownames(cluster0_peaks)[1],
  pt.size = 0.1,
  idents = c("0","1","3")
)
plot2 <- FeaturePlot(
  object = brain,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  max.cutoff = 'q95'
)
plot1 | plot2


################################################################################################################################

# Open - Close
open_cluster_0 = rownames(cluster0_peaks[cluster0_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_0, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_0.tsv")

open_cluster_1 = rownames(cluster1_peaks[cluster1_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_1, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_1.tsv")

open_cluster_2 = rownames(cluster2_peaks[cluster2_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_2, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_2.tsv")

open_cluster_3 = rownames(cluster3_peaks[cluster3_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_3, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_3.tsv")

open_cluster_4 = rownames(cluster4_peaks[cluster4_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_4, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_4.tsv")

open_cluster_5 = rownames(cluster5_peaks[cluster5_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_5, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_5.tsv")


open_cluster_6 = rownames(cluster6_peaks[cluster6_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_6, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_6.tsv")

open_cluster_7 = rownames(cluster7_peaks[cluster7_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_7, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_7.tsv")

open_cluster_8 = rownames(cluster8_peaks[cluster8_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_8, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_8.tsv")

open_cluster_9 = rownames(cluster9_peaks[cluster9_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_9, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_9.tsv")

open_cluster_10 = rownames(cluster10_peaks[cluster10_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_10, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_10.tsv")

open_cluster_11 = rownames(cluster11_peaks[cluster11_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_11, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_11.tsv")

open_cluster_12 = rownames(cluster12_peaks[cluster12_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_12, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_12.tsv")

open_cluster_13 = rownames(cluster13_peaks[cluster13_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_13, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_13.tsv")

open_cluster_14 = rownames(cluster14_peaks[cluster14_peaks$avg_log2FC > 0.25, ])
write.table(open_cluster_14, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/open_cluster_14.tsv")

# Closest features

closest_cluster_0 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_0)
write.table(closest_cluster_0, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_0.tsv")

closest_cluster_1 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_1)
write.table(closest_cluster_1, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_1.tsv")

closest_cluster_2 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_2)
write.table(closest_cluster_2, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_2.tsv")

closest_cluster_3 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_3)
write.table(closest_cluster_3, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_3.tsv")

closest_cluster_4 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_4)
write.table(closest_cluster_4, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_4.tsv")

closest_cluster_5 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_5)
write.table(closest_cluster_5, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_5.tsv")

closest_cluster_6 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_6)
write.table(closest_cluster_6, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_6.tsv")

closest_cluster_7 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_7)
write.table(closest_cluster_7, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_7.tsv")

closest_cluster_8 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_8)
write.table(closest_cluster_8, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_8.tsv")

closest_cluster_9 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_9)
write.table(closest_cluster_9, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_9.tsv")

closest_cluster_10 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_10)
write.table(closest_cluster_10, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_10.tsv")

closest_cluster_11 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_11)
write.table(closest_cluster_11, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_11.tsv")

closest_cluster_12 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_12)
write.table(closest_cluster_12, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_12.tsv")

closest_cluster_13 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_13)
write.table(closest_cluster_13, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_13.tsv")

closest_cluster_14 <- ClosestFeature(ATAC_MERGED.filt, open_cluster_14)
write.table(closest_cluster_14, file = "/media/dimbo/10T/data/talianidis_data/scRNA_seq/ATAC/analysis/Signac/closest_features_cluster_14.tsv")

################################################################################################################################

# Rename Clusters
ATAC_Renamed=RenameIdents(ATAC_MERGED.filt,  `0` = "Enterocyte Progenitor II", `1` = "Stem I", `2` = "S phase", `3` = "Stem II", `4` = "Enterocyte Immature", `5` = "Paneth - Stem", `6` = "Goblet I", `7` = "Goblet II", `8` = "Enterocyte Mature II", `9` = "?", `10` = "Tuft", `11` = "Enterocyte Mature I", `12` = "Enterocyte Progenitor I", `13` = "Enteroendocrine I", `14` = "Enteroendocrine II")

################################################################################################################################

# set plotting order
levels(ATAC_Renamed) <- c("Stem I", "Stem II", "Enterocyte Progenitor I", "Enterocyte Progenitor II", "Enterocyte Immature","Enterocyte Mature I", "Enterocyte Mature II", "Goblet I", "Goblet II", "S phase", "Tuft", "Enteroendocrine I", "Enteroendocrine II", "Paneth - Stem", "?")



# Test gene location
pdf("test_gene_location.pdf", width = 15, height = 8)
CoveragePlot(
  object = ATAC_Renamed,
  region = ("Olfm4"),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)
dev.off()


# Visualization of genomic regions (advanced)
covplot = CoveragePlot(
  object = ATAC_Renamed,
  region = "chr14-79998000-80022000",
  features = "Olfm4",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = TRUE,
  idents = c("Stem I", "Stem II", "Enterocyte Progenitor I", "Enterocyte Progenitor II", "Enterocyte Immature", "Enterocyte Mature I", "Enterocyte Mature II", "Enteroendocrine I", "Enteroendocrine II", "Goblet I", "Goblet II" ),
  
)

# Plotting per-cell fragment abundance
tile_plot <- TilePlot(
  object = ATAC_Renamed,
  region = "chr14-79998000-80022000",
  idents = c("Stem I", "Stem II" )
  
)



pdf("Combined_Olfm4.pdf", width = 12, height = 7)
CombineTracks(
  plotlist = list(covplot,tile_plot),
  heights = c(17, 5),
  widths = c(10, 5)
)
dev.off()


# Re cluster and integrate with scRNA subcluster
ATAC_Subet=subset(ATAC_MERGED.filt, idents = c(0,1,3,4,11,12))

ATAC_Subet <- RunTFIDF(ATAC_Subet)
ATAC_Subet <- FindTopFeatures(ATAC_Subet, min.cutoff = 'q0')
ATAC_Subet <- RunSVD(ATAC_Subet)

ATAC_Subet <- RunUMAP(object = ATAC_Subet, reduction = 'lsi', dims = 2:25)
ATAC_Subet <- FindNeighbors(object = ATAC_Subet, reduction = 'lsi', dims = 2:25)
ATAC_Subet <- FindClusters(object = ATAC_Subet, verbose = FALSE, algorithm = 3, resolution = 1.5)

pdf("UMAP_recluster.pdf", width = 15, height = 8)
DimPlot(object = ATAC_Subet, label = TRUE) + NoLegend()
dev.off()

gene.activities <- GeneActivity(ATAC_Subet)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
ATAC_Subet[['RNA']] <- CreateAssayObject(counts = gene.activities)
ATAC_Subet <- NormalizeData(
  object = ATAC_Subet,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC_Subet$nCount_RNA)
)

DefaultAssay(ATAC_Subet) <- 'RNA'


RNA_subset <- readRDS("/media/dimbo/10T/data/talianidis_data/scRNA_seq/Talianidis_GFP_2_fixed/analysis/Seurat/Lgr5Cre_0.1.3.4.5.rds")
RNA_subset_renamed=RenameIdents(RNA_subset,  `0` = "Stem I", `1` = "Stem II", `2` = "TA", `3` = "Progenitor Ι", `4` = "S phase cells", `5` = "Progenitor II", `6` = "Progenitor III", `7` = "Goblet-Paneth")

transfer.anchors <- FindTransferAnchors(
  reference = RNA_subset_renamed,
  query = ATAC_Subet,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA_subset_renamed$seurat_clusters,
  weight.reduction = ATAC_Subet[['lsi']],
  dims = 2:30
)

predicted.labels$predicted.id[predicted.labels$predicted.id == 0] <- "Stem I"
predicted.labels$predicted.id[predicted.labels$predicted.id == 1] <- "Stem II"
predicted.labels$predicted.id[predicted.labels$predicted.id == 2] <- "TA"
predicted.labels$predicted.id[predicted.labels$predicted.id == 3] <- "Progenitor I"
predicted.labels$predicted.id[predicted.labels$predicted.id == 4] <- "S phase Cells"
predicted.labels$predicted.id[predicted.labels$predicted.id == 5] <- "Progenitor II"
predicted.labels$predicted.id[predicted.labels$predicted.id == 6] <- "Progenitor III"
predicted.labels$predicted.id[predicted.labels$predicted.id == 7] <- "Goblet-Paneth"

ATAC_Subet <- AddMetaData(object = ATAC_Subet, metadata = predicted.labels)

plot1 <- DimPlot(
  object = RNA_subset_renamed,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = ATAC_Subet,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pdf("UMAP_scRNA-scATAC.pdf", 
    width = 17,
    height = 8)
plot1 + plot2
dev.off()



# Linking peaks to genes

DefaultAssay(ATAC_Renamed) <- "peaks"

  ATAC_Renamed <- RegionStats(ATAC_Renamed, genome = BSgenome.Mmusculus.UCSC.mm10)

# link peaks to genes
#  DefaultAssay(ATAC_Subet) <- "peaks"
  
ATAC_Renamed <- LinkPeaks(
  object = ATAC_Renamed,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Olfm4", "Muc2")
)

idents.plot = c("Stem I", "Stem II")

p1 <- CoveragePlot(
  object = ATAC_Renamed,
  region = "Olfm4",
  features = "Olfm4",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000
)

pdf("link_to_Gene_Olf4.pdf", width = 10, height = 8)
p1
dev.off()

#############################################


ATAC_Renamed <- LinkPeaks(
  object = ATAC_Renamed,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Lgr5", "Muc2")
)

idents.plot = c("Stem I", "Stem II")

p1 <- CoveragePlot(
  object = ATAC_Renamed,
  region = "Lgr5",
  features = "Lgr5",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000
)

pdf("link_to_Gene_Lgr5.pdf", width = 10, height = 8)
p1
dev.off()






