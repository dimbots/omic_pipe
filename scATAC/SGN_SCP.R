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
library(DoubletFinder)
library(data.table)
library(cowplot)
library(stringr)
library(SeuratWrappers)
library(monocle3)
library(Nebulosa)
library(BiocFileCache)
library(Matrix)
library(extrafont)
library(remotes)
library(cicero)
library(SCP)
library(JASPAR2020)
library(TFBSTools)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
#extrafont::font_import()
set.seed(7)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm9)

# setwd
setwd("/media/dimbo/10T/TAL_LAB/Data/Setdb1/scATAC/GFP/analysis/SGN_SCP/")

counts <- Read10X_h5("../Cellranger/02 Cellranger/cr_count/MUC23142/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../Cellranger/02 Cellranger/cr_count/MUC23142/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

WT_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '../Cellranger/02 Cellranger/cr_count/MUC23142/outs/fragments.tsv.gz',
  min.cells = 1
)

WT <- CreateSeuratObject(
  counts = WT_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(WT) <- annotations

# QC METRICS
WT <- NucleosomeSignal(object = WT)

WT$nucleosome_group <- ifelse(WT$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = WT, group.by = 'nucleosome_group', region = 'chr1-1-10000000')


WT <- TSSEnrichment(WT, fast = FALSE)

WT$high.tss <- ifelse(WT$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(WT, group.by = 'high.tss') + NoLegend()

WT$pct_reads_in_peaks <- WT$peak_region_fragments / WT$passed_filters * 100
WT$blacklist_ratio <- WT$blacklist_region_fragments / WT$peak_region_fragments

VlnPlot(
  object = WT,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

WT <- subset(
  x = WT,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
WT

# Normalization / linear dimensional reduction
WT <- RunTFIDF(WT)
WT <- FindTopFeatures(WT, min.cutoff = 'q0')
WT <- RunSVD(object = WT)
# Correlation
DepthCor(WT, n = 30)

WT <- RunUMAP(
  object = WT,
  reduction = 'lsi',
  dims = 2:30
)
WT <- FindNeighbors(
  object = WT,
  reduction = 'lsi',
  dims = 2:30
)
WT <- FindClusters(
  object = WT,
  algorithm = 3,
  resolution = 2.5,
  verbose = FALSE
)
DimPlot(object = WT, label = TRUE) + NoLegend()


# compute gene activities
gene.activities <- GeneActivity(WT)

# add the gene activity matrix to the Seurat object as a new assay
WT[['RNA']] <- CreateAssayObject(counts = gene.activities)
WT <- NormalizeData(
  object = WT,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(WT$nCount_RNA)
)

DefaultAssay(WT) <- 'RNA'
FeaturePlot(
  object = WT,
  features = c('Olfm4'),
  pt.size = 0.1,
  max.cutoff = 'q80',
  ncol = 1
)

# INTEGRATION WITH scRNA

# Load the pre-processed scRNA-seq data
rna <- readRDS("/media/dimbo/10T/TAL_LAB/Data/Setdb1/scRNA/Talianidis_GFP_2_fixed/analysis/SRT_SCP/Lgr5Cre_MERGED_Renamed.rds")

split_obj = SplitObject(rna, split.by = "orig.ident")

# See how the objects are named
View(split_obj)

rna_KO = split_obj[["Setdb1KO"]]

rna_WT = split_obj[["Lgr5Cre"]]


rna_WT <- FindVariableFeatures(
  object = rna_WT,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = rna_WT,
  query = WT,
  reduction = 'cca',
  dims = 2:30
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_WT$CellType,
  weight.reduction = WT[['lsi']],
  dims = 2:30
)

WT <- AddMetaData(object = WT, metadata = predicted.labels)


plot1 <- CellDimPlot(rna_WT, group.by = 'CellType',label = FALSE,
                     palcolor = c("darkslategray3", "deepskyblue3", "dodgerblue4", "darkolivegreen3", 
                                  "seagreen4", "navajowhite2", "navajowhite3", "sandybrown", "salmon1", 
                                  "salmon3", "firebrick2", "pink2", "plum3", "orchid3", "orangered4"), ) + NoLegend() + ggtitle('scRNA WT')
plot2 <- DimPlot(WT, group.by = 'predicted.id', label = FALSE, repel = TRUE,
                 cols = c('Stem I' = 'darkslategray3', 'Stem II' = 'deepskyblue3', 'Stem III' = 'dodgerblue4', 
                          'Progenitor I' = 'darkolivegreen3','Progenitor II' = 'seagreen4', 'Ent.Immature' = 'navajowhite2', 
                          'Ent.Mature' = 'navajowhite3', 'Goblet I' = 'sandybrown', 'Goblet II' = 'salmon1', 
                          'Goblet III' = 'salmon3','Paneth' = "firebrick2", 'Tuft' = 'pink2',
                          'Enteroendocrine I' = 'plum3', 'Enteroendocrine' = 'orchid3',
                          'Unclassified' = 'orangered4'), pt.size = 0.5
) + NoLegend() + ggtitle('scATAC WT \n nCells:753')
plot1 + plot2

# replace each label with its most likely prediction
for(i in levels(WT)) {
  cells_to_reid <- WhichCells(WT, idents = i)
  newid <- names(which.max(table(WT$predicted.id[cells_to_reid])))
  Idents(WT, cells = cells_to_reid) <- newid
}

# set plotting order
levels(WT) <- c("Stem I","Stem II","Stem III","Progenitor I","Progenitor II", "Ent.Immature","Ent.Mature","Goblet I","Goblet II","Goblet III","Paneth",
                "Tuft","Enteroendocrine I","Enteroendocrine II","Unclassified")

#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

DefaultAssay(WT) <- 'peaks'
# SUBSET STEM I
WT_STEMI =subset(WT, idents = "Stem I")

# PEAK CALLING WITH MACS
peaks <- CallPeaks(WT_STEMI, broad = TRUE,format = "BED", effective.genome.size = 1.87e+09 )
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse",species = "Mus_musculus")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(WT_STEMI),
  features = peaks,
  cells = colnames(WT_STEMI)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
WT_STEMI[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = ("/media/dimbo/10T/TAL_LAB/Data/Setdb1/scATAC/GFP/analysis/Cellranger/02 Cellranger/cr_count/MUC23142/outs/fragments.tsv.gz"),
  annotation = annotations
)

# Linking peak to genes
WT_STEMI <- RegionStats(WT_STEMI, genome = BSgenome.Mmusculus.UCSC.mm10)

# Link peaks in all genes
WT_STEMI <- LinkPeaks(
  object = WT_STEMI,
  peak.assay = "peaks",
  expression.assay = "RNA"
)

LINK_TO_GENES = Links(WT_STEMI)
write.table(LINK_TO_GENES, file='WT_LINK_TO_GENES.tsv', quote=FALSE, sep='\t', col.names = NA)

#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

# PLOTING SPECIFIC GENES

# COVPLOT
cov_plot <- CoveragePlot(
  object = WT_STEMI,
  region = "Axin2",
  annotation = FALSE,
  peaks = FALSE,links = FALSE
)
cov_plot
# PEAKPLOT
peak_plot <- PeakPlot(
  object = WT_STEMI,
  region = "chr11-108920000-108950000"
)
peak_plot
# LINKS
link_plot <- LinkPlot(
  object = WT_STEMI,
  region = "Axin2"
)
link_plot
# GENEPLOT
gene_plot <- AnnotationPlot(
  object = WT_STEMI,
  region = "Axin2"
)
gene_plot

# COMBINE
CombineTracks(
  plotlist = list(cov_plot, peak_plot,link_plot, gene_plot),
  heights = c(2,0.5, 0.8, 0.8),
  widths = c(1, 1)
)# # 480 / 320

#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

# MOTIF ANALYSIS
# Get a list of motif position frequency matrices from the JASPAR database
#pfm <- getMatrixSet(
 # x = JASPAR2020,
  #opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
#)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Mus musculus", all_versions = FALSE)
)

# add motif information (WT ALL CLUSTERS)
WT <- AddMotifs(
  object = WT,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

# Find differential accessible peaks against Ent.Immature and Goblet
da_peaks_WT <- FindMarkers(
  object = WT,
  ident.1 = 'Stem I',
  ident.2 = c("Stem II","Ent.Immature","Ent.Mature","Goblet III", "Enteroendocrine II"),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak_WT <- rownames(da_peaks_WT[da_peaks_WT$p_val < 0.005, ])

# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(WT, idents = c("Stem I"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(WT, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak_WT, ],
  n = 50000
)

# test enrichment
enriched.motifs <- FindMotifs(
  object = WT,
  features = top.da.peak_WT
)

write.table(enriched.motifs, file='enriched_motifs_WT.tsv', quote=FALSE, sep='\t', col.names = NA)


MotifPlot(
  object = WT,
  motifs = head(rownames(enriched.motifs))
) # 800/250


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

# Transcription factor footprinting

# extract position frequency matrices for the motifs
# gather the footprinting information for sets of motifs
WT <- Footprint(
  object = WT,
  motif.name = c("Rbpjl"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(WT, features = c("Rbpjl"))
p2 + patchwork::plot_layout(ncol = 1)



###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

# KO SAMPLE
counts <- Read10X_h5("../Cellranger/02 Cellranger/cr_count/MUC23141/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../Cellranger/02 Cellranger/cr_count/MUC23141/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

KO_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '../Cellranger/02 Cellranger/cr_count/MUC23141/outs/fragments.tsv.gz',
  min.cells = 1
)

KO <- CreateSeuratObject(
  counts = KO_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(KO) <- annotations

# QC METRICS
KO <- NucleosomeSignal(object = KO)

KO$nucleosome_group <- ifelse(KO$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = KO, group.by = 'nucleosome_group', region = 'chr1-1-10000000')


KO <- TSSEnrichment(KO, fast = FALSE)

KO$high.tss <- ifelse(KO$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(KO, group.by = 'high.tss') + NoLegend()

KO$pct_reads_in_peaks <- KO$peak_region_fragments / KO$passed_filters * 100
KO$blacklist_ratio <- KO$blacklist_region_fragments / KO$peak_region_fragments

VlnPlot(
  object = KO,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

KO <- subset(
  x = KO,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 30000 &
    pct_reads_in_peaks > 30 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
KO

# Normalization / linear dimensional reduction
KO <- RunTFIDF(KO)
KO <- FindTopFeatures(KO, min.cutoff = 'q0')
KO <- RunSVD(object = KO)
# Correlation
DepthCor(KO, n = 30)

KO <- RunUMAP(
  object = KO,
  reduction = 'lsi',
  dims = 2:30
)
KO <- FindNeighbors(
  object = KO,
  reduction = 'lsi',
  dims = 2:30
)
KO <- FindClusters(
  object = KO,
  algorithm = 3,
  resolution = 2.5,
  verbose = FALSE
)
DimPlot(object = KO, label = TRUE) + NoLegend()


# compute gene activities
gene.activities <- GeneActivity(KO)

# add the gene activity matrix to the Seurat object as a new assay
KO[['RNA']] <- CreateAssayObject(counts = gene.activities)
KO <- NormalizeData(
  object = KO,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(KO$nCount_RNA)
)

DefaultAssay(KO) <- 'RNA'
FeaturePlot(
  object = KO,
  features = c('Olfm4'),
  pt.size = 0.1,
  max.cutoff = 'q80',
  ncol = 1
)

# INTEGRATION WITH scRNA
rna_KO <- FindVariableFeatures(
  object = rna_KO,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = rna_KO,
  query = KO,
  reduction = 'cca',
  dims = 2:30
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_KO$CellType,
  weight.reduction = KO[['lsi']],
  dims = 2:30
)

KO <- AddMetaData(object = KO, metadata = predicted.labels)


plot1 <- CellDimPlot(rna_KO, group.by = 'CellType',label = FALSE,
                     palcolor = c("darkslategray3", "deepskyblue3", "dodgerblue4", "darkolivegreen3", 
                                  "seagreen4", "navajowhite2", "navajowhite3", "sandybrown", "salmon1", 
                                  "salmon3", "firebrick2", "pink2", "plum3", "orchid3", "orangered4"), ) + NoLegend() + ggtitle('scRNA KO')
plot2 <- DimPlot(KO, group.by = 'predicted.id', label = FALSE, repel = TRUE,
                 cols = c('Stem I' = 'darkslategray3', 'Stem II' = 'deepskyblue3', 'Stem III' = 'dodgerblue4', 
                          'Progenitor I' = 'darkolivegreen3','Progenitor II' = 'seagreen4', 'Ent.Immature' = 'navajowhite2', 
                          'Ent.Mature' = 'navajowhite3', 'Goblet I' = 'sandybrown', 'Goblet II' = 'salmon1', 
                          'Goblet III' = 'salmon3','Paneth' = "firebrick2", 'Tuft' = 'pink2',
                          'Enteroendocrine I' = 'plum3', 'Enteroendocrine' = 'orchid3',
                          'Unclassified' = 'orangered4'), pt.size = 0.5
) + NoLegend() + ggtitle('scATAC KO \n nCells:775')
plot1 + plot2

# replace each label with its most likely prediction
for(i in levels(KO)) {
  cells_to_reid <- WhichCells(KO, idents = i)
  newid <- names(which.max(table(KO$predicted.id[cells_to_reid])))
  Idents(KO, cells = cells_to_reid) <- newid
}

# set plotting order
levels(KO) <- c("Stem I","Stem II","Stem III","Progenitor I","Progenitor II", "Ent.Immature","Ent.Mature","Goblet I","Goblet II","Goblet III","Paneth",
                "Tuft","Enteroendocrine I","Enteroendocrine II","Unclassified")

#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

DefaultAssay(KO) <- 'peaks'
# SUBSET STEM I
KO_STEMI =subset(KO, idents = "Stem I")

# PEAK CALLING WITH MACS
peaks <- CallPeaks(KO_STEMI, broad = TRUE,format = "BED", effective.genome.size = 1.87e+09 )
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse",species = "Mus_musculus")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(KO_STEMI),
  features = peaks,
  cells = colnames(KO_STEMI)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
KO_STEMI[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = ("/media/dimbo/10T/TAL_LAB/Data/Setdb1/scATAC/GFP/analysis/Cellranger/02 Cellranger/cr_count/MUC23141/outs/fragments.tsv.gz"),
  annotation = annotations
)

# Linking peak to genes
KO_STEMI <- RegionStats(KO_STEMI, genome = BSgenome.Mmusculus.UCSC.mm10)

# Link peaks in all genes
KO_STEMI <- LinkPeaks(
  object = KO_STEMI,
  peak.assay = "peaks",
  expression.assay = "RNA"
)

LINK_TO_GENES = Links(KO_STEMI)
write.table(LINK_TO_GENES, file='KO_LINK_TO_GENES.tsv', quote=FALSE, sep='\t', col.names = NA)

#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

# PLOTING SPECIFIC GENES

# COVPLOT
cov_plot_tmp <- CoveragePlot(
  object = KO_STEMI,
  region = "Axin2",
  annotation = FALSE,
  peaks = FALSE,links = FALSE
)
cov_plot = cov_plot_tmp + scale_fill_brewer(type = "seq", palette = "Dark2")
cov_plot
# PEAKPLOT
peak_plot <- PeakPlot(
  object = KO_STEMI,color = "forestgreen", 
  region = "chr11-108920000-108950000"
)
peak_plot
# LINKS
#KO_STEMI <- LinkPeaks(
#object = KO_STEMI,
#peak.assay = "peaks",
#expression.assay = "RNA", genes.use = "Ascl2"
#)
link_plot <- LinkPlot(
  object = KO_STEMI,
  region = "Axin2"
)
link_plot
# GENEPLOT
gene_plot <- AnnotationPlot(
  object = KO_STEMI,
  region = "Axin2"
)
gene_plot

# COMBINE
CombineTracks(
  plotlist = list(cov_plot, peak_plot,link_plot, gene_plot),
  heights = c(2,0.5, 0.8, 0.8),
  widths = c(1, 1)
)# # 480 / 320

#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

# MOTIF ANALYSIS
# Get a list of motif position frequency matrices from the JASPAR database
#pfm <- getMatrixSet(
# x = JASPAR2020,
#opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
#)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Mus musculus", all_versions = FALSE)
)

# add motif information (KO ALL CLUSTERS)
KO <- AddMotifs(
  object = KO,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

# Find differential accessible peaks against Ent.Immature and Goblet
da_peaks_KO <- FindMarkers(
  object = KO,
  ident.1 = 'Stem I',
  ident.2 = c("Stem II","Ent.Immature","Ent.Mature","Goblet III", "Enteroendocrine I"),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak_KO <- rownames(da_peaks_KO[da_peaks_KO$p_val < 0.005, ])

# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(KO, idents = c("Stem I"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(KO, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak_KO, ],
  n = 50000
)

# test enrichment
enriched.motifs <- FindMotifs(
  object = KO,
  features = top.da.peak_KO
)

write.table(enriched.motifs, file='enriched_motifs_KO.tsv', quote=FALSE, sep='\t', col.names = NA)


MotifPlot(
  object = KO,
  motifs = head(rownames(enriched.motifs))
) # 800/250


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

# Transcription factor footprinting

# extract position frequency matrices for the motifs
# gather the footprinting information for sets of motifs
KO <- Footprint(
  object = KO,
  motif.name = c("Rbpjl"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(KO, features = c("Rbpjl"))
p2 + patchwork::plot_layout(ncol = 1)





