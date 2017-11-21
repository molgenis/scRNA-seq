###########################################################################################################################
# Authors: Dylan de Vries
# Name: interactionOtherCelltypes.R
# Function: Calculate and plot all interactions for the non CD4+ cell types
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(Seurat)
library(broom)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################
# Name: get.snp
# Function: returns the genotype of the input SNP for every person and makes the genotype order alphabetically consistent
# Input:
#   Name      Type          Description
#   snp.id    character     The SNP ID of the eQTL SNP
#
# Output:
#   A factor with the genotype information for the input SNP per person

get.snp <- function(snp.id) {
  snp <- unlist(genotypes[snp.id,])

  #Refactor the genotypes to be consistent in all samples.
  snp[snp == "C/A"] <- "A/C"
  snp[snp == "T/A"] <- "A/T"
  snp[snp == "T/C"] <- "C/T"
  snp[snp == "G/A"] <- "A/G"
  snp[snp == "G/C"] <- "C/G"
  snp[snp == "G/T"] <- "T/G"
  return(as.factor(snp))
}

###########################################################################################################################
#
# Main code
#
###########################################################################################################################

load("/groups/umcg-wijmenga/tmp02/projects/scRNAseq_10X_pilot/pilot3/clustering/pilot3_subsetted_celltypes_final_ensemble.Rda")


genotypes <- read.table("/groups/umcg-bios/tmp04/users/umcg-mdam/LLD_genotypes/LLDVCF_GoNLimputed/trityper/table/maf_10.calls.txt", check.names=F)
cellTypes <- combined.major@ident
cd4TCells <- which(cellTypes == "CD4+ T")

samples <- FetchData(combined.major, "sample")[-cd4TCells,]
sample.names <- levels(droplevels(samples))
expression.data <- combined.major@data[,-cd4TCells]

genotypes <- genotypes[,match(sample.names, colnames(genotypes))]
snp.genotype <- as.numeric(get.snp("rs7297175"))

correlations <- c()
cell.counts <- c()


for (sample in sample.names){
	sample.indices <- which(samples == sample)
	eQTL.gene.expression <- expression.data["ENSG00000197728",sample.indices]
	interaction.gene.expression <- expression.data["ENSG00000122026",sample.indices]
	correlations <- c(correlations, cor(eQTL.gene.expression, interaction.gene.expression, method="spearman"))
	cell.counts <- c(cell.counts, length(sample.indices))
}

model <- lm(formula = correlations~snp.genotype, weights = sqrt(cell.counts))
png("~/interactionTcells_RPS26.png")
ggplot() + 
		geom_boxplot(aes(x=snp.genotype, y=correlations, colour=as.factor(snp.genotype)), outlier.shape = NA) +
		geom_point(aes(x=snp.genotype, y=correlations, colour=as.factor(snp.genotype)), position = position_jitter(width = 0.2), size=2) +
	    scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
	    ylab("Spearman correlation") +
	    xlab("Genotype") +
	    ggtitle(paste0("CD4+ T cells", gene, "\np-value: ", formatC(tidy(model)[2,"p.value"], format = "e", digits = 2), ", R:", round(tidy(model)[2,"statistic"] / sqrt(length(snp.genotype) - 2 + tidy(model)[2,"statistic"] ** 2), 2)))
dev.off()


pdf("~/cellType_interactions.pdf", onefile=T)
for (cellType in unique(cellTypes)){
	targetCells <- which(cellTypes ==cellType)

	samples <- FetchData(combined.major, "sample")[targetCells,]
	sample.names <- levels(droplevels(samples))
	expression.data <- combined.major@data[,targetCells]


	correlations <- c()
	cell.counts <- c()


	for (sample in sample.names){
		sample.indices <- which(samples == sample)
		eQTL.gene.expression <- expression.data["ENSG00000197728",sample.indices]
		interaction.gene.expression <- expression.data["ENSG00000122026",sample.indices]
		correlations <- c(correlations, cor(eQTL.gene.expression, interaction.gene.expression, method="spearman"))
		cell.counts <- c(cell.counts, length(sample.indices))
	}
	if (length(cell.counts) != 45){next}
	model <- lm(formula = correlations~snp.genotype, weights = sqrt(cell.counts))
	tidy(model)[2,"p.value"]
	#"#57a350", "#fd7600", "#383bfe"
	colours = rep("#57a350", length(correlations))
	colours[snp.genotype == 1] = "#fd7600"
	colours[snp.genotype == 2] = "#383bfe"
	plot(snp.genotype, correlations, col=colours, main=paste(cellType, tidy(model)[2,"p.value"]))
}

dev.off()

cellTypeList = list(c("CD8+ T"), c("cMonocyte", "ncMonocyte"), c("CD56(dim) NK", "CD56(bright) NK"))
names(cellTypeList) = c("CD8+ T cells", "Monocytes", "NK cells")
snp.genotype <- as.numeric(get.snp("rs7297175"))

pdf("~/cellType_interactions.pdf", onefile=T)
for( cellType in names(cellTypeList)){
	targetCells <- which(cellTypes %in% cellTypeList[[cellType]])

	samples <- FetchData(combined.major, "sample")[targetCells,]
	sample.names <- levels(droplevels(samples))
	expression.data <- combined.major@data[,targetCells]


	correlations <- c()
	cell.counts <- c()


	for (sample in sample.names){
		sample.indices <- which(samples == sample)
		eQTL.gene.expression <- expression.data["ENSG00000197728",sample.indices]
		interaction.gene.expression <- expression.data["ENSG00000122026",sample.indices]
		correlations <- c(correlations, cor(eQTL.gene.expression, interaction.gene.expression, method="spearman"))
		cell.counts <- c(cell.counts, length(sample.indices))
	}
	if (length(cell.counts) != 45){next}
	model <- lm(formula = correlations~snp.genotype, weights = sqrt(cell.counts))
	tidy(model)[2,"p.value"]
	#"#57a350", "#fd7600", "#383bfe"
	colours = rep("#57a350", length(correlations))
	colours[snp.genotype == 1] = "#fd7600"
	colours[snp.genotype == 2] = "#383bfe"
	# png(paste0("~/interaction_", cellType, ".png"), height=2000, width=2000)
	print(ggplot() + 
		geom_boxplot(aes(x=snp.genotype, y=correlations, colour=as.factor(snp.genotype)), outlier.shape = NA) +
		geom_point(aes(x=snp.genotype, y=correlations, colour=as.factor(snp.genotype)), position = position_jitter(width = 0.2), size=2) +
	    scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
	    ylab("Spearman correlation") +
	    xlab("Genotype") +
	    ggtitle(paste0(cellType, "\np-value: ", formatC(tidy(model)[2,"p.value"], format = "e", digits = 2), ", R:", round(tidy(model)[2,"statistic"] / sqrt(length(snp.genotype) - 2 + tidy(model)[2,"statistic"] ** 2), 2))))
	# dev.off()
}

dev.off()

genderInfo <- read.table("Supplementary_Table_NEW_sample_lane_gender.txt", header=T, sep="\t")

snp.genotype.m <- snp.genotype[genderInfo[genderInfo$Sex == "Male", "LL"]]
snp.genotype.f <- snp.genotype[genderInfo[genderInfo$Sex == "Female", "LL"]]

pdf("~/gender_interactions.pdf", onefile=T)

for (gene in c("ENSG00000147403", "ENSG00000198034", "ENSG00000205542")){
	correlations.m <- c()
	correlations.f <- c()
	cell.counts.m <- c()
	cell.counts.f <- c()

	for (sample in sample.names){
		sample.indices <- which(samples == sample)
		eQTL.gene.expression <- expression.data["ENSG00000197728",sample.indices]
		interaction.gene.expression <- expression.data[gene,sample.indices]
		correlation <- cor(eQTL.gene.expression, interaction.gene.expression, method="spearman")
		cell.count <- length(sample.indices)
		if (genderInfo[genderInfo$LL == substr(sample, nchar(sample)-10, nchar(sample)), "Sex"] == "Male"){
			correlations.m <- c(correlations.m, correlation)
			cell.counts.m <- c(cell.counts.m, cell.count)
		} else if (genderInfo[genderInfo$LL == substr(sample, nchar(sample)-10, nchar(sample)), "Sex"] == "Female"){
			correlations.f <- c(correlations.f, correlation)
			cell.counts.f <- c(cell.counts.f, cell.count)
		}
	}

	model.m <- lm(formula = correlations.m~snp.genotype.m, weights = sqrt(cell.counts.m))
	model.f <- lm(formula = correlations.f~snp.genotype.f, weights = sqrt(cell.counts.f))


	print(ggplot() + 
		geom_boxplot(aes(x=snp.genotype.m, y=correlations.m, colour=as.factor(snp.genotype.m)), outlier.shape = NA) +
		geom_point(aes(x=snp.genotype.m, y=correlations.m, colour=as.factor(snp.genotype.m)), position = position_jitter(width = 0.2), size=2) +
	    scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
	    ylab("Spearman correlation") +
	    xlab("Genotype") +
	    ggtitle(paste0("CD4+ T cells - Male ", gene, "\np-value: ", formatC(tidy(model.m)[2,"p.value"], format = "e", digits = 2), ", R:", round(tidy(model.m)[2,"statistic"] / sqrt(length(snp.genotype.m) - 2 + tidy(model.m)[2,"statistic"] ** 2), 2))))

		print(ggplot() + 
		geom_boxplot(aes(x=snp.genotype.f, y=correlations.f, colour=as.factor(snp.genotype.f)), outlier.shape = NA) +
		geom_point(aes(x=snp.genotype.f, y=correlations.f, colour=as.factor(snp.genotype.f)), position = position_jitter(width = 0.2), size=2) +
	    scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
	    ylab("Spearman correlation") +
	    xlab("Genotype") +
	    ggtitle(paste0("CD4+ T cells - Female\np-value: ", formatC(tidy(model.f)[2,"p.value"], format = "e", digits = 2), ", R:", round(tidy(model.f)[2,"statistic"] / sqrt(length(snp.genotype.f) - 2 + tidy(model.f)[2,"statistic"] ** 2), 2))))
}

dev.off()



