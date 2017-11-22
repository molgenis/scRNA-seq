###########################################################################################################################
# Authors: Harm Brugge & Dylan de Vries
# Name: seurat_clustering.R
# Function: Script to merge and cluster the scRNA-seq Data
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################

library(Seurat)
library(Matrix)
library(Matrix.utils)
library(RColorBrewer)
library(dplyr)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

# Name: get.violin.data
# Function: returns the expression data for the createn of violin plots
# Input:
#   Name      Type          Description
#   seurat    Seurat        Seurat object containing normalized expression data
#   genes     vector        Vector containing the gene names for plotting
# Output:
#   Melted dataframe containing al the info for plotting
get.violin.data <- function(seurat, genes) {
  output <- data.frame(gene = character(0), value= numeric(0), ident = character(0))
  for (gene in genes) {
    data.use = data.frame(FetchData(seurat,gene))
    data.use = t(data.use)
    data.melt=data.frame(rep(gene, length(seurat@ident)))
    colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data.use[1,1:length(seurat@ident)])
    data.melt$id=names(data.use)[1:length(seurat@ident)]
    data.melt$ident=seurat@ident
    noise <- rnorm(length(data.melt$value))/100000
    data.melt$value=as.numeric(as.character(data.melt$value))+noise
    output <- rbind(output, data.melt)
  }
  return(output)
}

# Name: add.data
# Function: Merges cellranger count matrices into one sparse matrix
# Input:
#   Name        Type            Description
#   orig.data   Sparse matrix   Sparse matrix where the data needs to be added to.
#   lane        Integer         Number of the lane which needs to be added
# Output:
#   Merged data
add.data <- function(orig.data = NULL, lane) {
  data <- readMM(paste0("../scRNA-seq/data/lane_", lane, "/matrix.mtx"))
  rownames(data) <- make.names(sapply(readLines(paste0("../scRNA-seq/data/lane_", lane, "/genes.tsv")), extract_field, 1, delim = "\\t"), unique=TRUE)
  colnames(data) <- paste0(sapply(readLines(paste0("../scRNA-seq/data/lane_", lane, "/barcodes.tsv")), extract_field, 1, delim = "-"), "_lane", lane)
  
  if (is.null(orig.data)) return(data)
  else return(merge.Matrix(orig.data, data, by.x=rownames(orig.data), by.y=rownames(data), all.x=T, all.y=T))
}

# Name: DoClusterAnalysis
# Function: Run all the Seurat function for a clustering
# Input:
#   Name        Type            Description
#   seurat      Seurat object   Seurat object to do the analysis on.
#   pc.use      Integer         Number of principal componentes to include in the clustering
#   y.cutoff    Numeric         Cut-off for the variable genes function
#   ssn.res     Numeric         Resolution used for the SNN clustering
# Output:
#   The Seurat object
DoClusterAnalysis <- function(seurat, pc.use=16, y.cutoff = 1, ssn.res = 1.2) {
  seurat <- MeanVarPlot(seurat, x.low.cutoff = 0, y.cutoff = y.cutoff, do.plot = F)
  seurat <- RegressOut(seurat, genes.regress = seurat@var.genes, latent.vars = c("percent.mito", "nUMI"))
  seurat <- PCAFast(seurat, pc.genes = seurat@var.genes, pcs.compute = pc.use, do.print = F)
  PCElbowPlot(seurat, num.pc = pc.use)
  seurat <- RunTSNE(seurat, dims.use = 1:pc.use, do.fast = T)
  seurat <- FindClusters(seurat, pc.use = 1:pc.use, resolution = ssn.res, save.SNN = T, do.sparse = T)
  seurat <- BuildClusterTree(seurat, do.reorder = T, reorder.numeric = T)
  FeaturePlot(seurat, "llkDoublet.llkSinglet")
  return(seurat)
}

###########################################################################################################################
#
# Main
#
###########################################################################################################################
pilot3.merged <- add.data(lane = 1)
pilot3.merged <- add.data(pilot3.merged, lane = 2)
pilot3.merged <- add.data(pilot3.merged, lane = 3)
pilot3.merged <- add.data(pilot3.merged, lane = 4)
pilot3.merged <- add.data(pilot3.merged, lane = 5)
pilot3.merged <- add.data(pilot3.merged, lane = 6)
pilot3.merged <- add.data(pilot3.merged, lane = 7)
pilot3.merged <- add.data(pilot3.merged, lane = 8)

dim(pilot3.merged)

##
## Combined clustering
##
combined.seurat <- new("seurat", raw.data = pilot3.merged)
combined.seurat <- Setup(combined.seurat, min.cells = 1, min.genes = 0, project = "pilot3", do.scale = F, do.center = F, names.field = 1, names.delim = "\\-")
dim(combined.seurat@data)

## Load deAnonymized data
pilot3.de.anonymize <- data.frame(read.table("clustering/deanonymize.txt"))
row.names(pilot3.de.anonymize) <- pilot3.de.anonymize[,1]
pilot3.de.anonymize[,1] <- NULL
colnames(pilot3.de.anonymize) <- c("llkDoublet.llkSinglet", "nSNPs_tested", "sample", "lane")
combined.seurat <- AddMetaData(combined.seurat, pilot3.de.anonymize)

## Mitochondrial genes
genes <- read.table("../scRNA-seq/data/lane_1/genes.tsv")
mt.ensemble <- genes[grep("^MT-", genes$V2),"V1"]
combined.percent.mito <- colSums(expm1(combined.seurat@data[mt.ensemble, ])) / colSums(expm1(combined.seurat@data))
combined.seurat <- AddMetaData(combined.seurat, combined.percent.mito, "percent.mito")

##
## Do a preliminary clustring for QC-plots
##
## Calculata variable genes
combined.seurat <- MeanVarPlot(combined.seurat, x.low.cutoff = 0, y.cutoff = 1, do.plot = F)
## Regress out % mitochondrial genes and nUMI
combined.seurat <- RegressOut(combined.seurat, genes.regress = combined.seurat@var.genes, latent.vars = c("percent.mito", "nUMI"))
#combined.seurat <- ScaleData(combined.seurat, genes.use = combined.seurat@var.genes) # scale all data if wanted
combined.seurat <- PCAFast(combined.seurat, pc.genes = combined.seurat@var.genes, pcs.compute = 40, do.print = F)
## Use the first 16 PCs for t-SNE
combined.seurat <- RunTSNE(combined.seurat, dims.use = 1:16, do.fast = T)

combined.meta.data <- FetchData(combined.seurat, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
true.doublet <- combined.meta.data$llkDoublet.llkSinglet > 0 & !combined.meta.data$lane %in% c(2,3)
is.mito.cut <- combined.meta.data$percent.mito > 0.05
gene.high <- combined.meta.data$nGene > 3500
gene.low <- combined.meta.data$nGene < 500

qc.colors <- rep("Singlet", 28855)
qc.colors[gene.high] <- ">3500 genes"
qc.colors[is.mito.cut] <- ">5% mitochondrial"
qc.colors[true.doublet] <- "Doublet"
qc.colors <- factor(qc.colors)
qc.colors <- factor(qc.colors, levels(qc.colors)[c(4, 3, 2, 1)])

colors <- c(alpha("lightgrey",0.2), "#153057","#009ddb","#e64b50")

tsne.plot <- ggplot(combined.seurat@tsne.rot, aes(x=tSNE_1,y=tSNE_2, colour=qc.colors)) +
  geom_point(size=0.8) +
  scale_color_manual(values=colors) +
  theme_minimal(base_size = 14) +
  ylab("t-SNE 1") + xlab("t-SNE 2") +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        axis.line = element_line(size = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=8, alpha=1)))
tsne.plot

combined.seurat.subset <- SubsetData(combined.seurat, subset.name = "llkDoublet.llkSinglet", accept.high = -2)
dim(combined.seurat.subset@data)

combined.meta.data <- FetchData(combined.seurat.subset, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
ggplot(combined.meta.data, aes(nUMI, percent.mito)) + geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Fraction mtDNA-encoded genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 0.06, colour="red")

combined.seurat.subset <- SubsetData(combined.seurat.subset, subset.name = "percent.mito", accept.high = 0.05)
dim(combined.seurat.subset@data)
combined.meta.data <- FetchData(combined.seurat.subset, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
# Gene / UMI hexagon plot
ggplot(combined.meta.data, aes(nUMI, nGene)) + 
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Number of genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 3500, colour="red")

combined.seurat.subset <- SubsetData(combined.seurat.subset, subset.name = "nGene", accept.high = 3500)
dim(combined.seurat.subset@data)

# Gene / UMI doublets colored
plot(combined.meta.data$nGene~combined.meta.data$nUMI, 
     ylab = "Number of genes", xlab = "Number of reads", 
     col = ifelse(combined.meta.data$llkDoublet.llkSinglet > 100 & !combined.meta.data$lane %in% c(2,3) ,'red', rgb(0,0,0,0.1)),
     cex=0.2, cex.lab=1.4, pch=20)

combined.seurat.subset <- DoClusterAnalysis(combined.seurat.subset, pc.use = 16)
combined.major <- combined.seurat.subset

save(combined.major, file = "data/subsetted_celltypes.Rda")
#load(file = "data/subsetted_celltypes.Rda")

## Rename cluster identities
levels(combined.major@ident) <- c(levels(combined.major@ident), "cMonocyte","ncMonocyte","mDC","Megakaryocyte","CD56(bright) NK", "CD56(dim) NK", "pDC","Plasma","B", "CD4+ T", "CD8+ T")

combined.major@ident[combined.major@ident %in% c(16,17,19,20,21,22)] <- "CD4+ T"
combined.major@ident[combined.major@ident %in% c(18,15,11,10,9)] <- "CD8+ T"
combined.major@ident[combined.major@ident %in% c(1,3)] <- "cMonocyte"
combined.major@ident[combined.major@ident == 4] <- "mDC"
combined.major@ident[combined.major@ident == 2] <- "ncMonocyte"
combined.major@ident[combined.major@ident == 5] <- "Megakaryocyte"
combined.major@ident[combined.major@ident == 7] <- "CD56(dim) NK"
combined.major@ident[combined.major@ident == 6] <- "CD56(dim) NK"
combined.major@ident[combined.major@ident == 8] <- "CD56(bright) NK"
combined.major@ident[combined.major@ident == 13] <- "Plasma" 
combined.major@ident[combined.major@ident == 14] <- "B" 
combined.major@ident[combined.major@ident == 12] <- "pDC"

save(combined.major, file="clustering/subsetted_celltypes_final_ensemble.Rda")
meta <- cbind(combined.major@tsne.rot, combined.major@ident, FetchData(combined.major, c("sample")))

write.table(as.matrix(combined.major@data), file = "clustering/expression_all_cells.tsv", quote = F, sep = "\t", col.names=NA)
write.table(as.matrix(meta), file = "clustering/tsne_rot_ident_all_cells.tsv", quote = F, sep = "\t", col.names=NA)





