###########################################################################################################################
# Authors: Dylan de Vries
# Name: magic_plots.R
# Function: Generates plots to determine the effect of MAGIC imputation on the gene expression in the CD4+ T cells.
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(ggplot2)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################
magic.plot <- function(gene1, gene2, gene1.hgnc, gene2.hgnc){
	imputed.exp.gene1 <- c()
	imputed.exp.gene2 <- c()
	nonImputed.exp.gene1 <- c()
	nonImputed.exp.gene2 <- c()
	for (i in 1:length(exp.matrices.nonImputed)){
		if (gene1 %in% colnames(exp.matrices[[i]])){
			imputed.exp.gene1 <- c(imputed.exp.gene1, exp.matrices[[i]][,gene1])
		} else {
			imputed.exp.gene1 <- c(imputed.exp.gene1, rep(0, nrow(exp.matrices[[i]])))
		}

		if (gene2 %in% colnames(exp.matrices[[i]])){
			imputed.exp.gene2 <- c(imputed.exp.gene2, exp.matrices[[i]][,gene2])
		} else {
			imputed.exp.gene2 <- c(imputed.exp.gene2, rep(0, nrow(exp.matrices[[i]])))
		}

		if (gene1 %in% colnames(exp.matrices.nonImputed[[i]])){
			nonImputed.exp.gene1 <- c(nonImputed.exp.gene1, exp.matrices.nonImputed[[i]][,gene1])
		} else {
			nonImputed.exp.gene1 <- c(nonImputed.exp.gene1, rep(0, nrow(exp.matrices.nonImputed[[i]])))
		}

		if (gene2 %in% colnames(exp.matrices.nonImputed[[i]])){
			nonImputed.exp.gene2 <- c(nonImputed.exp.gene2, exp.matrices.nonImputed[[i]][,gene2])
		} else {
			nonImputed.exp.gene2 <- c(nonImputed.exp.gene2, rep(0, nrow(exp.matrices.nonImputed[[i]])))
		}
	}
	plot.nonImputed <- ggplot() +
		geom_point(aes(x=nonImputed.exp.gene1, y=nonImputed.exp.gene2)) +
		theme_minimal(base_size = 64) +
		theme(plot.title = element_text(hjust = 0.5)) +
		labs(x = paste(gene1.hgnc, "expression"), y = paste(gene2.hgnc, "expression")) +
		ggtitle(paste0("Correlation measured expression\n", gene1.hgnc, " vs. ", gene2.hgnc))


	plot.imputed <- ggplot() +
		geom_point(aes(x=imputed.exp.gene1, y=imputed.exp.gene2)) +
		theme_minimal(base_size = 64) +
		theme(plot.title = element_text(hjust = 0.5)) +
		labs(x = paste(gene1.hgnc, "expression"),
			y = paste(gene2.hgnc, "expression")) +
		ggtitle(paste0("Correlation imputed expression\n", gene1.hgnc, " vs. ", gene2.hgnc))

	return(list(plot.nonImputed, plot.imputed))
}

###########################################################################################################################
#
# Main code
#
###########################################################################################################################
##
## Read the imputed data
##
dir <- list.files(path="/Users/dylandevries/Documents/work/interactionAnalysis/expression_files/", pattern="*.tsv", full.names=T, recursive=FALSE)
exp.matrices <- list()
cell.counts <- c()
sample.names <- vector()

i <- 1
for (file in dir) {
  sample.names <- c(sample.names, tools::file_path_sans_ext(basename(file)))
  sample <- read.csv(file = file, sep = "\t", row.names = 1)
  rownames(sample) <- substr(rownames(sample), start = 7, stop = 21)
  sample <- sample[apply(sample, 1, function(x){!any(x == 0)}),]
  sample <- t(sample)
  cell.counts <- c(cell.counts, nrow(sample))

  exp.matrices[[i]] <- sample
  i <- i + 1
}

##
## Read the non-imputed data
##
exp.matrices.nonImputed <- list()
sample.names <- vector()

dir.path <- "/Users/dylandevries/Documents/work/interactionAnalysis/nonImputed_expression_files/thCellsPerSample/"
dir <- list.dirs(path=dir.path, full.names=T, recursive=FALSE)

i <- 1
for (folder in dir) {
  print(i)
  print(paste0(dir[[i]], "/matrix.mtx"))
  sample.names <- c(sample.names, tools::file_path_sans_ext(basename(folder)))
  sample.raw <- readMM(paste0(dir[[i]], "/matrix.mtx"))
  rownames(sample.raw) <- read.table(paste0(dir[[i]], "/genes.tsv"), stringsAsFactors = F)$V1
  colnames(sample.raw) <- read.table(paste0(dir[[i]], "/barcodes.tsv"), stringsAsFactors = F)$V1
  exp.matrices.nonImputed[[i]] <- t(as.matrix(sample.raw))

  i <- i + 1
}

##
## Make the plots
##
plots.CD3D.CD3E <- magic.plot("ENSG00000167286", "ENSG00000198851", "CD3D", "CD3E")
png("./plots/CD3D_CD3E_nonImputed.png", width=1500, height=1500)
print(plots.CD3D.CD3E[[1]])
dev.off()
png("./plots/CD3D_CD3E_imputed.png", width=1500, height=1500)
print(plots.CD3D.CD3E[[2]])
dev.off()

plots.CD3D.CD4 <- magic.plot("ENSG00000167286", "ENSG00000010610", "CD3D", "CD4")
png("./plots/CD3D_CD4_nonImputed.png", width=1500, height=1500)
print(plots.CD3D.CD4[[1]])
dev.off()
png("./plots/CD3D_CD4_imputed.png", width=1500, height=1500)
print(plots.CD3D.CD4[[2]])
dev.off()

plots.CD3D.CDMS4A1 <- magic.plot("ENSG00000167286", "ENSG00000156738", "CD3D", "MS4A1")
png("./plots/CD3D_MS4A1_nonImputed.png", width=1500, height=1500)
print(plots.CD3D.CDMS4A1[[1]])
dev.off()
png("./plots/CD3D_MS4A1_imputed.png", width=1500, height=1500)
print(plots.CD3D.CDMS4A1[[2]])
dev.off()

plots.CD3D.CD14 <- magic.plot("ENSG00000167286", "ENSG00000170458", "CD3D", "CD14")
png("./plots/CD3D_CD14_nonImputed.png", width=1500, height=1500)
print(plots.CD3D.CD14[[1]])
dev.off()
png("./plots/CD3D_CD14_imputed.png", width=1500, height=1500)
print(plots.CD3D.CD14[[2]])
dev.off()
