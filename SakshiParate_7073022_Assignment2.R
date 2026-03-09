
# Sakshi Parate - 7073022

BiocManager::install("GreenleafLab/motifmatchr")

system("/opt/homebrew/Cellar/gcc/14.2.0_1/bin/g++-14 --version")
BiocManager::install("chromVAR", force = TRUE)

BiocManager::install("ge11232002/TFBSTools")
BiocManager::install("GreenleafLab/chromVARmotifs")

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library(ArchR) #run archr
ArchR::installExtraPackages()
install.packages("Cairo")
capabilities("cairo")

remotes::install_github(
  "GreenleafLab/ArchR",
  ref = "master",
  repos = BiocManager::repositories()
)

#Week 1

#1 Preprocessing and quality control

#1.1 Set up the environment

setwd("~/Pictures/sakshi")
set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = 1) 

data_dir <- "~/Pictures/sakshi/p2_updated/project_2"
inputFiles <- list.files(data_dir, pattern = "*.tsv.gz", full.names = TRUE)

sampleNames <- gsub(pattern = ".tsv.gz", replacement = "", basename(inputFiles))
print(sampleNames)

#1.2 Read the data into an appropriate data structure and apply filtering

ArrowFiles <- createArrowFiles( # create arrowfiles from .tsv.gz
  inputFiles = inputFiles,
  sampleNames = sampleNames,  
  minTSS = 5,
  minFrags = 550,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
) 

ArrowFiles <- list.files("~/Pictures/sakshi", pattern = "arrow$", full.names = TRUE)
ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,  
  outputDirectory = "ArchRProject"   #create archr project using arrowfiles
)
proj

#1.3 Identify doublets

proj <- addDoubletScores(  #asign doublet score
  input = proj,
  k = 10,               
  knnMethod = "UMAP",   
  LSIMethod = 1         
)
proj

doubletEnrichment <- proj$DoubletEnrichment 
doubletEnrichment
doubletFilter <- doubletEnrichment < 2 # adjust thresold 
doubletFilter
filteredProj <- proj[doubletFilter, ]# filter archr project to exculde doublet
filteredProj

#1.4 Collect all samples into a joint data structure

num_cells <- nCells(proj).  # Collect project statistics
num_cells
median_TSS <- median(proj$TSSEnrichment)
median_TSS
median_fragments <- median(proj$nFrags)
median_fragments
tile_dimensions <- dim(getMatrixFromProject(proj, useMatrix = "TileMatrix", binarize = TRUE))
tile_dimensions

filtered_num_cells <- nCells(filteredProj)
filtered_num_cells
filtered_median_TSS <- median(filteredProj$TSSEnrichment)
filtered_median_TSS
filtered_median_fragments <- median(filteredProj$nFrags)
filtered_median_fragments
filtered_tile_dimensions <- dim(getMatrixFromProject(filteredProj, useMatrix = "TileMatrix", binarize = TRUE))
filtered_tile_dimensions

#1.5 Quality control
table(proj$Sample)

plotFragmentSizes(proj)

plotTSSEnrichment(ArchRProj = proj, groupBy = "Sample")

library(ggplot2)
library(dplyr)

# extract metadata
df <- as.data.frame(getCellColData(proj, select = c("Sample", "nFrags", "TSSEnrichment")))
df$FragBin <- cut( # create frag count bins
  df$nFrags,
  breaks = c(0, 1000, 3000, 6000, 10000, Inf),
  labels = c("<1k", "1–3k", "3–6k", "6–10k", ">10k"),
  include.lowest = TRUE
)

p <- ggplot(df, aes(x = FragBin, y = TSSEnrichment, fill = Sample)) + # violin plot for tss enrichment vs frag bins by sample
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.3) +
  facet_wrap(~ Sample) +
  theme_minimal() +
  labs(
    title = "Fragments vs TSS Enrichment per Sample (Violin Plot)",
    x = "Fragment Count (Binned)",
    y = "TSS Enrichment"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

df <- as.data.frame(getCellColData(proj, select = c("Sample", "nFrags", "TSSEnrichment"))) # scatter plot for number of frag vs tss enrichment
samples <- unique(df$Sample)
for (s in samples) {
  dsub <- subset(df, Sample == s)
  p <- ggplot(dsub, aes(x = nFrags, y = TSSEnrichment)) +
    geom_point(alpha = 0.35, size = 0.6) +
    scale_x_log10() +
    theme_minimal() +
    labs(
      title = paste0("Number of Fragments vs TSS Enrichment — ", s),
      x = "Number of fragments (log10)",
      y = "TSS enrichment"
    )
  print(p)                                    # display in R
  ggsave(filename = paste0(gsub("[^A-Za-z0-9_\\-]", "_", s), "_frags_vs_TSS.png"),
         plot = p, width = 6, height = 4, dpi = 300)
}

#1.5 Filter the dataset
strict_filteredProj <- filteredProj[
  filteredProj$TSSEnrichment > 8 &
    filteredProj$nFrags > 3000,
]
strict_filteredProj

nCells(strict_filteredProj) # num of cells

median(strict_filteredProj$TSSEnrichment) # median tss

median(strict_filteredProj$nFrags) # median fragments

strict_tile_dimensions <- dim( # tile set dimension
  getMatrixFromProject(strict_filteredProj, useMatrix = "TileMatrix", binarize = TRUE)
)
strict_tile_dimensions

#Week 2

#2 Peaks

#2.1 Peak calling

pathToMacs2 <- "/home/sakshi/miniconda3/envs/macs2_fix/bin/macs2"

proj_peaks <- strict_filteredProj

groupBy <- "Sample"

proj_peaks <- addGroupCoverages( #create pseudo bulk coverages
  ArchRProj = proj_peaks,
  groupBy = groupBy
)

proj_peaks <- addReproduciblePeakSet( #call peaks using mac2
  ArchRProj = proj_peaks,
  groupBy = groupBy,
  pathToMacs2 = pathToMacs2
)
proj_peaks <- addPeakMatrix(proj_peaks) #add peak matrix

head(peak_set)
table(peak_set$sample)
head(peak_set)
gene_score_matrix <- getMatrixFromProject(proj_peaks, useMatrix = "GeneScoreMatrix")
head(gene_score_matrix)  # chek data structure

proj_peaks <- addPeakMatrix(proj_peaks)
getAvailableMatrices(proj_peaks)
peak_matrix <- getMatrixFromProject(proj_peaks, useMatrix = "PeakMatrix")
head(peak_matrix)
peak_matrix_df <- as.data.frame(assay(peak_matrix, "PeakMatrix"))
head(peak_matrix_df)
print(peak_matrix_df)

length(getPeakSet(proj_peaks))

#2.1 Cluster marker peaks

proj_peaks <- addIterativeLSI( #cluster marker peaks param
  ArchRProj = proj_peaks,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30
)

proj_peaks <- addClusters(
  input = proj_peaks,
  reducedDims = "IterativeLSI",
  name = "Clusters",
  resolution = 0.8
)

proj_peaks <- addUMAP(
  ArchRProj = proj_peaks,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

p1 <- plotEmbedding(proj_peaks, colorBy = "cellColData", name = "Clusters")
p1

markerPeaks <- getMarkerFeatures( #heatmap showing accessibility
  ArchRProj = proj_peaks,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "nFrags"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(
  markerPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

heatmap <- markerHeatmap(
  seMarker = markerPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)
heatmap

library(grid) #plots around gnes

genes_to_plot <- c("CD8A", "CD14", "GATA1", "PAX5", "TBX21")

for (gene in genes_to_plot) {
  
  p <- plotBrowserTrack(
    ArchRProj = proj_peaks,
    groupBy = "Clusters",
    geneSymbol = gene,
    upstream = 50000,
    downstream = 50000,
    loops = NULL
  )
  
  grid::grid.newpage()
  grid::grid.draw(p[[1]])
}

library(grid)

plot_CD8A <- plotBrowserTrack(
  ArchRProj = proj_peaks,
  groupBy = "Clusters",
  geneSymbol = "CD8A",
  upstream = 50000,
  downstream = 50000
)

plot_CD14 <- plotBrowserTrack(
  ArchRProj = proj_peaks,
  groupBy = "Clusters",
  geneSymbol = "CD14",
  upstream = 50000,
  downstream = 50000
)

plot_GATA1 <- plotBrowserTrack(
  ArchRProj = proj_peaks,
  groupBy = "Clusters",
  geneSymbol = "GATA1",
  upstream = 50000,
  downstream = 50000
)

plot_PAX5 <- plotBrowserTrack(
  ArchRProj = proj_peaks,
  groupBy = "Clusters",
  geneSymbol = "PAX5",
  upstream = 50000,
  downstream = 50000
)

plot_TBX21 <- plotBrowserTrack(
  ArchRProj = proj_peaks,
  groupBy = "Clusters",
  geneSymbol = "TBX21",
  upstream = 50000,
  downstream = 50000
)

grid.newpage()
grid.draw(plot_CD8A[[1]])
grid.text("Gene: CD8A", x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold", col = "blue"))

grid.newpage()
grid.draw(plot_CD14[[1]])
grid.text("Gene: CD14", x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold", col = "blue"))

grid.newpage()
grid.draw(plot_GATA1[[1]])
grid.text("Gene: GATA1", x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold", col = "blue"))

grid.newpage()
grid.draw(plot_PAX5[[1]])
grid.text("Gene: PAX5", x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold", col = "blue"))

grid.newpage()
grid.draw(plot_TBX21[[1]])
grid.text("Gene: TBX21", x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold", col = "blue"))

#3 Dimensionality Reduction 

#3.1 Iterative LSI

#3.1 UMAP with sample annotation and QC metrics

proj_peaks <- addUMAP(
  ArchRProj = proj_peaks,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

plot_sample <- plotEmbedding( #color by sample
  ArchRProj = proj_peaks, 
  colorBy = "cellColData", 
  name = "Sample",
  embedding = "UMAP"
)
plot_sample

install.packages("hexbin")

plot_tss <- plotEmbedding( #color by tss enrichment
  ArchRProj = proj_peaks, 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  embedding = "UMAP"
)
plot_tss

plot_fragments <- plotEmbedding( #color by number of fragments
  ArchRProj = proj_peaks, 
  colorBy = "cellColData", 
  name = "nFrags",
  embedding = "UMAP"
)
plot_fragments

plot_sample + ggtitle("UMAP Colored by Sample")
plot_tss + ggtitle("UMAP Colored by TSS Enrichment")
plot_fragments + ggtitle("UMAP Colored by Number of Fragments")

#3.2 Dealing with batch effects

proj_peaks <- addHarmony( #deALIng with batch effects
  ArchRProj = proj_peaks,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

proj_peaks <- addUMAP(
  ArchRProj = proj_peaks,
  reducedDims = "Harmony",
  name = "UMAP_Harmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

umap_sample_harmony <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Harmony"
)
umap_sample_harmony

umap_tss_harmony <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "cellColData",
  name = "TSSEnrichment",
  embedding = "UMAP_Harmony"
)
umap_tss_harmony

umap_frag_harmony <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "cellColData",
  name = "nFrags",
  embedding = "UMAP_Harmony"
)
umap_frag_harmony

#4 Clustering 

proj_peaks <- addClusters( #clustering
  input = proj_peaks,
  reducedDims = "IterativeLSI",
  method = "Seurat",     # use seurat lovain clustering
  name = "Clusters_Louvain",
  resolution = 0.8
)

p_clusters <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "cellColData",
  name = "Clusters_Louvain",
  embedding = "UMAP"
)
p_clusters

cluster_cell_counts <- table(proj_peaks$Clusters_Louvain)
cluster_cell_counts

table(proj_peaks$Clusters_Louvain, proj_peaks$Sample)

#Week 3

#5 Gene activity

#5.1 Compute gene activity scores

if (!"GeneScoreMatrix" %in% getAvailableMatrices(proj_peaks)) { #check if GeneScoreMatrix already exists
  proj_peaks <- addGeneScoreMatrix(
    ArchRProj = proj_peaks,
    matrixName = "GeneScoreMatrix"
  )
} else {
  message("Gene activity scores already computed.")
}

GeneScoreMatrix <- getMatrixFromProject( #extract GeneScoreMatrix
  ArchRProj = proj_peaks,
  useMatrix = "GeneScoreMatrix"
)
GeneScoreMatrix

#5.2 Identify Marker Genes

markersGene <- getMarkerFeatures( #identify marker genes from GeneScoreMatrix
  ArchRProj = proj_peaks, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",  #grouping by clusters to find cluster specific markers
  bias = c("TSSEnrichment", "log10(nFrags)"),  #correcting for biases in tss enrichment n frag count
  testMethod = "wilcoxon"  #using the wilcoxon test for differential accessibility
)
markersGene #show identified marker genes

markerList <- getMarkers( #apply cutoffs to get list of significant marker genes
  seMarker = markersGene, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1"  # FDR <= 0.05 and Log2 fold change >= 1 for significant markers
)
markerList #display top markers

print(markerList$C1) #view marker genes for specific cluster
print(markerList$C2)
print(markerList$C3)
print(markerList$C4)
print(markerList$C5)
print(markerList$C6)
print(markerList$C7)

#5.3 Using MAGIC

topGenes_C1 <- markerList$C1$name[1:5] #define top 5 marker genes for each cluster 
topGenes_C2 <- markerList$C2$name[1:5]
topGenes_C3 <- markerList$C3$name[1:5]
topGenes_C4 <- markerList$C4$name[1:5]
topGenes_C5 <- markerList$C5$name[1:5]
topGenes_C6 <- markerList$C6$name[1:5]
topGenes_C7 <- markerList$C7$name[1:5]

umap_withoutMagic_C1 <- plotEmbedding( #create umap without magic
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C1,  #top genes for cluster c1
  embedding = "UMAP_Harmony"  #umap with harmony correction, no magic smoothing
)

umap_withoutMagic_C2 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C2,  #top genes for cluster c2
  embedding = "UMAP_Harmony"
)

umap_withoutMagic_C3 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C3,  #top genes for cluster c3
  embedding = "UMAP_Harmony"
)

umap_withoutMagic_C4 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C4,  #top genes for cluster c4
  embedding = "UMAP_Harmony"
)

umap_withoutMagic_C5 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C5,  #top genes for cluster c5
  embedding = "UMAP_Harmony"
)

umap_withoutMagic_C6 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C6,  #top genes for cluster c6
  embedding = "UMAP_Harmony"
)

umap_withoutMagic_C7 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C7,  #top genes for cluster c7
  embedding = "UMAP_Harmony"
)

proj_peaks <- addImputeWeights(proj_peaks) #add imputation weights for magic

umap_withMagic_C1 <- plotEmbedding( #create umap with magic
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C1,  #top genes for cluster c1
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj_peaks)
)

umap_withMagic_C2 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C2,  #top genes for cluster c2
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj_peaks)
)

umap_withMagic_C3 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C3,  #top genes for cluster c3
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj_peaks)
)

umap_withMagic_C4 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C4,  #top genes for cluster c4
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj_peaks)
)

umap_withMagic_C5 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C5,  #top genes for cluster c5
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj_peaks)
)

umap_withMagic_C6 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C6,  #top genes for cluster c6
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj_peaks)
)

umap_withMagic_C7 <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneScoreMatrix",
  name = topGenes_C7,  #top genes for cluster c7
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj_peaks)
)

umap_withoutMagic_C1  #display umap for cluster c1 without magic
umap_withMagic_C1  #display umap for cluster c1 with magic

umap_withoutMagic_C2  #display umap for cluster c2 without magic
umap_withMagic_C2  #display umap for cluster c2 with magic

umap_withoutMagic_C3  #display umap for cluster c3 without magic
umap_withMagic_C3  #display umap for cluster c3 with magic

umap_withoutMagic_C4  #display umap for cluster c4 without magic
umap_withMagic_C4  #display umap for cluster c4 with magic

umap_withoutMagic_C5  #display umap for cluster c5 without magic
umap_withMagic_C5  #display umap for cluster c5 with magic

umap_withoutMagic_C6  #display umap for cluster c6 without magic
umap_withMagic_C6  #display umap for cluster c6 with magic

umap_withoutMagic_C7  #display umap for cluster c7 without magic
umap_withMagic_C7  #display umap for cluster c7 with magic

#6 Transcription Factor motif activity

#6.1 Compute TF motif activity
proj_peaks <- addMotifAnnotations(
  ArchRProj = proj_peaks,
  motifSet = "cisbp",  
  name = "Motif"       
)

print(getPeakAnnotation(proj_peaks, "Motif"))

proj_peaks <- addBgdPeaks(proj_peaks)

proj_peaks <- addDeviationsMatrix(
  ArchRProj = proj_peaks,
  peakAnnotation = "Motif",
  force = TRUE
)

#6.2 Plot UMAP embeddings for marker TFs

var_motifs <- getVarDeviations( #get variability scores for motifs
  ArchRProj = proj_peaks, 
  name = "MotifMatrix",  #this is matrix that contains motif activity scores
  plot = FALSE  #do not plot variability of all motifs yet
)

Top_var_motifs <- head(var_motifs, 2) #identify the most variable motifs top 2
top_names <- Top_var_motifs$name

markerMotifs <- getFeatures(proj_peaks, select = paste(top_names, collapse="|"), useMatrix = "MotifMatrix") #retrieve motifs for top 2 most variable TFs

markerMotifs_filtered <- grep("^z:", markerMotifs, value = TRUE) #filter motifs to remove unwanted entries

for (motif in markerMotifs_filtered) { #ensure each plot is displayed separately
  plot <- plotEmbedding(   #plot motif activity on umap
    ArchRProj = proj_peaks,
    colorBy = "MotifMatrix",  #use motif activity scores
    name = motif,             #specify motif name to highlight
    embedding = "UMAP_Harmony" #specify embedding to use 
  )
  dev.new()  #this will open a new graphics window for each plot 
  print(plot)
}

#6.3 Motif activity

table(proj_peaks$Clusters_Louvain)

for (motif in markerMotifs_filtered) {
  plot <- plotGroups(
    ArchRProj = proj_peaks,
    groupBy = "Clusters_Louvain",
    colorBy = "MotifMatrix",
    name = motif,
    plotAs = "violin"
  )
  print(plot)
}

#7 Integration with gene expression

#7.1 Data integration

seurat_data <- readRDS( #load annotated scRNAseq data
  "~/Pictures/sakshi/blish_awilk_seu_subset.rds"
)
seurat_data

colnames(seurat_data@meta.data)

proj_peaks <- addGeneIntegrationMatrix(
  ArchRProj = proj_peaks,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seurat_data,
  groupRNA = "cell.type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)

markerGenes <- c("PAX5", "OLIG2", "EGR1")

umap_plots <- plotEmbedding(
  ArchRProj = proj_peaks,
  colorBy = "GeneIntegrationMatrix",
  name = markerGenes,
  embedding = "UMAP_Harmony"
)
umap_plots$PAX5
umap_plots$EGR1
umap_plots$OLIG2

#7.2 Correlation Coefficients

cor_results <- correlateMatrices(
  ArchRProj = proj_peaks,
  useMatrix1 = "GeneIntegrationMatrix",  # scRNA-seq expression
  useMatrix2 = "GeneScoreMatrix"         # scATAC-seq gene activity
)

head(cor_results) #inspect the result
summary(cor_results$cor)

highest_agreement <- cor_results[which.max(cor_results$cor), ] #genes with highest agreement

lowest_agreement <- cor_results[which.min(cor_results$cor), ] #genes with lowest agreement

print("Gene with highest agreement:")
highest_agreement

print("Gene with lowest agreement:")
lowest_agreement

#7.3 Cluster labels from gene expression

plotEmbedding( #visualize RNA inferred cell type labels on atac umap
  ArchRProj = proj_peaks,
  colorBy = "cellColData",
  name = "predictedGroup_Un",
  embedding = "UMAP_Harmony"
)

proj_peaks$Clusters_RNA <- proj_peaks$predictedGroup_Un #assign rna inferred cell type labels to atac cells

confusion_matrix <- table( #compute confusion matrix between atac clusters n rna labels
  ATAC_Clusters = proj_peaks$Clusters_Louvain,
  RNA_Labels = proj_peaks$predictedGroup_Un
)

print(confusion_matrix) #display confusion matrix

#8 Peak-gene linkage

proj_peaks <- addPeak2GeneLinks( #compute peak to gene links using chromatin accessibility n gene expression
  ArchRProj = proj_peaks,
  reducedDims = "Harmony"   #use batch corrected low dimensional space
)

peak2gene_links <- getPeak2GeneLinks( #retrieve significant peak gene links
  ArchRProj = proj_peaks,
  corCutOff = 0.5,        #correlation threshold
  resolution = 10000     #10 kb resolution
)
peak2gene_links[[1]] #inspect first linked peak gene pairs

peak_gene_heatmap <- plotPeak2GeneHeatmap( #plot side by side heatmaps
  ArchRProj = proj_peaks
)
peak_gene_heatmap #display heatmap

#Week 4

#9 Differential accessibility 

#9.1 Differential peak accessibility

colnames(getCellColData(proj_peaks))
condition_map <- c(
  "ATAC_555_2_fragments_fragments"        = "COVID",
  "ATAC_557_fragments_fragments"          = "COVID",
  "ATAC_EV08_fragments_fragments"         = "Healthy",
  "ATAC_HIP02_frozen_fragments_fragments" = "Healthy"
)
proj_peaks$Condition <- condition_map[proj_peaks$Sample]
table(proj_peaks$Condition)
table(proj_peaks$Condition, proj_peaks$Sample)

diff_peaks <- getMarkerFeatures( #differential accessibility betw covid and healthy
  ArchRProj = proj_peaks,
  useMatrix = "PeakMatrix",
  groupBy = "Condition",                  #covid vs healthy
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

covid_peaks_relaxed <- getMarkers(
  diff_peaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)$COVID

healthy_peaks_relaxed <- getMarkers(
  diff_peaks,
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)$Healthy

ma_plot <- markerPlot(
  seMarker = diff_peaks,
  name = "COVID",
  cutOff = "FDR <= 0.05",
  plotAs = "MA"
)
ma_plot

volcano_plot <- markerPlot(
  seMarker = diff_peaks,
  name = "COVID",
  cutOff = "FDR <= 0.05",
  plotAs = "Volcano"
)
volcano_plot

nrow(covid_peaks)
nrow(healthy_peaks)

#9.2 TF motif enrichment

getAvailableMatrices(proj_peaks) #motifMatrix must exist

enrich_healthy <- peakAnnoEnrichment( #healthy enriched peaks
  seMarker = diff_peaks,
  ArchRProj = proj_peaks,
  peakAnnotation = "Motif",
  group = "Healthy",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5" #peaks less accessible in covid, Log2FC < 0
)
enrich_healthy

enrich_covid <- peakAnnoEnrichment( #covid enriched peaks
  seMarker = diff_peaks,
  ArchRProj = proj_peaks,
  peakAnnotation = "Motif",
  group = "COVID",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5" #peaks more accessible in covid, Log2FC > 0
)
enrich_covid

top_healthy_motifs <- rownames(enrich_healthy)[1:min(5, nrow(enrich_healthy))] #extract top enriched motifs
top_covid_motifs   <- rownames(enrich_covid)[1:min(5, nrow(enrich_covid))]

top_healthy_motifs
top_covid_motifs

plot_motifs <- unique(c(top_healthy_motifs, top_covid_motifs)) #combine motifs for plotting
plot_motifs

motif_umap <- plotEmbedding( #umap visualization of motif activity
  ArchRProj = proj_peaks,
  colorBy = "MotifMatrix",
  name = plot_motifs,
  embedding = "UMAP_Harmony"
)
motif_umap

motif_heatmap <- plotGroups( #heatmap of motif activity across clusters
  ArchRProj = proj_peaks,
  groupBy = "Clusters_Louvain",
  colorBy = "MotifMatrix",
  name = plot_motifs,
  plotAs = "heatmap"
)
motif_heatmap

#10 TF footprinting

#11 Co-accessibility

getAvailableMatrices(proj_peaks)

getReducedDims(proj_peaks)

proj_peaks <- addCoAccessibility( #add coaccessibility scores
  ArchRProj = proj_peaks,
  reducedDims = "IterativeLSI"
)

coaccessibility_links <- getCoAccessibility( #extract coaccessibility links
  ArchRProj = proj_peaks,
  corCutOff = 0.5,          #correlation threshold
  resolution = 10000        #10 kb resolution
)
coaccessibility_links[[1]]

marker_genes <- c("CD14", "CD8A") #genome track plot with 2 markers

for (gene in marker_genes) {
  
  p <- plotBrowserTrack(
    ArchRProj = proj_peaks,
    groupBy = "Condition",        #covid vs healthy
    geneSymbol = gene,
    upstream = 50000,
    downstream = 50000,
    loops = coaccessibility_links #add coaccessibility arcs
  )
  print(p[[1]])
}
p

grid::grid.newpage()
grid::grid.draw(p[[1]])

#11.2 Identify potential enhancers for marker genes

marker_genes <- c(
  "CD34", "GATA1", "PAX5", "MS4A1", "MME",
  "CD14", "MPO", "IRF8", "CD3D", "CD8A",
  "CD4", "TBX21", "CD3G", "NCAM1", "FCGR3A",
  "FOXP3", "GATA3", "RORC", "PDCD1", "HLA-DRA",
  "CD28", "IL2RA", "CD69", "CD44", "CCR7"
)

proj_peaks <- addCoAccessibility(
  ArchRProj = proj_peaks,
  reducedDims = "IterativeLSI"
)

coaccessibility_links <- getCoAccessibility(
  ArchRProj = proj_peaks,
  corCutOff = 0.5,      # significant links only
  resolution = 10000   # 10 kb bins
)

browser_tracks <- plotBrowserTrack(
  ArchRProj  = proj_peaks,
  groupBy    = "Condition",     # COVID vs Healthy
  geneSymbol = marker_genes,
  upstream   = 50000,
  downstream = 50000,
  loops      = coaccessibility_links
)

grid::grid.draw(browser_tracks$CD14) 
grid::grid.draw(browser_tracks$CD8A)

plotPDF(
  plotList   = browser_tracks,
  name       = "11.2_Marker_Genes_CoAccessibility.pdf",
  ArchRProj  = proj_peaks,
  addDOC     = FALSE,
  width      = 6,
  height     = 6
)



