###########################################################################################################################
# Authors: Harm Brugge
# Name: concordance_eqtls
# Function: Determine concordance RNA-seq datasets
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(Matrix)
library(Matrix.utils)

###########################################################################################################################
#
# Load data
#
###########################################################################################################################
# BIOS eQTLs independent effects
eqtls.rnaseq <- read.table("~/Documents/scRNA-seq/R/data/gene_level_eQTLs_independent_effects_interactions.txt", header = T, sep = "\t")
# DeepSAGE top eQTLS
eqtls.deepsage <- read.table("data/DeepSAGE_eqtls.csv", sep = ";", header = T, dec = ",")
# scRNA-seq eQTLs PBMCs (bulk like) confined to BIOS eQTLs
eqtls.sc.rnaseq <- read.table("~/Documents/scRNA-seq/R/data/eQTLProbesFDR0.05-ProbeLevel_RNA-seq-all-cells.txt", header = T, sep = "\t")
# scRNA-seq eQTLs PBMCs (bulk like) confined to deepSAGE eQTLs
eqtls.all.centered.deepsage <- read.table("data/eQTLProbesFDR0.05-ProbeLevel_DeepSAGE-all-cells.txt", header = T, sep = "\t")
# scRNA-seq expression data (Seurat)
load("clustering/pilot3_subsetted_celltypes_final_ensemble.Rda")
## Or a tsv 
#all.cells.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_all_log_norm.tsv", header = T, sep="\t")


###########################################################################################################################
#
# Format data
#
###########################################################################################################################
# Use SNP and ensemble ids to match the BIOS and single cell.
eqtls.rnaseq$combine = paste(eqtls.rnaseq$SNP, eqtls.rnaseq$Gene)
eqtls.sc.rnaseq$combine = paste(eqtls.sc.rnaseq$SNPName, eqtls.sc.rnaseq$ProbeName)
eqtls.deepsage$combine = paste(eqtls.deepsage$SNPName, eqtls.deepsage$HGNCName)
eqtls.all.centered.deepsage$combine = paste(eqtls.all.centered.deepsage$SNPName, eqtls.all.centered.deepsage$HGNCName)
eqtls.rnaseq.matched <- eqtls.rnaseq[match(eqtls.sc.rnaseq$combine,eqtls.rnaseq$combine, nomatch = 0),]
eqtls.matched.deepsage <- eqtls.deepsage[match(eqtls.all.centered.deepsage$combine,eqtls.deepsage$combine, nomatch = 0),]
# Flip z-scores if assesed allele is different
eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.centered.deepsage$AlleleAssessed,]$OverallZScore <- eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.centered.deepsage$AlleleAssessed,]$OverallZScore * -1
eqtls.rnaseq.matched[eqtls.rnaseq.matched$Assesed.Allele != eqtls.sc.rnaseq$AlleleAssessed,]$Z.score <- eqtls.rnaseq.matched[eqtls.rnaseq.matched$Assesed.Allele != eqtls.sc.rnaseq$AlleleAssessed,]$Z.score * -1

# Calculate average expression for every gene in the single cell data.
exp.sums <- rowMeans(combined.major@data)
# Match them to the eQTL genes
exp.sums.matched <- exp.sums[match(eqtls.rnaseq.matched$Gene, names(exp.sums), nomatch = 0)]

# calculate snp genen distance
eqtls.matched.deepsage$SNP.probe.dist <- eqtls.matched.deepsage$SNPChrPos - eqtls.matched.deepsage$ProbeCenterChrPos
#eqtls.th.centered.deepsage$SNP.probe.dist <- eqtls.th.centered.deepsage$SNPChrPos - eqtls.th.centered.deepsage$ProbeCenterChrPos


###########################################################################################################################
#
# Plots
#
###########################################################################################################################
##
## BIOS (RNA-seq) concordance
##
plot(eqtls.rnaseq.matched$Z.score, eqtls.sc.rnaseq$OverallZScore, type = "n"
     ,cex.lab=1.4,
     xlab = "RNA-seq effect size (Z-scores)",
     ylab = "scRNA-seq effect size (Z-scores)", lwd=0)
limits = par()$usr
rect(0,0,limits[2],limits[4],col="#E6EFEA", border = NA)
rect(0,0,limits[1],limits[3],col="#E6EFEA", border = NA)
box(lwd=3, col = "grey")
points(eqtls.rnaseq.matched$Z.score, eqtls.sc.rnaseq$OverallZScore, pch=21, bg=rgb(0,0,0,0.4), cex=log10(exp.sums.matched) + 1.5)
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)
title(main = "RNA-seq", adj = 0, font.main = 1, cex.main = 1.75, line = 0.7)

concordance.perc <- signif(sum(eqtls.rnaseq.matched$Z.score * eqtls.sc.rnaseq$OverallZScore > 0) / length(eqtls.sc.rnaseq$OverallZScore) * 100, 3)

##
## DeepSAGE concordance
##

# Match the expression to deepsage data
exp.sums.matched <- exp.sums[match(eqtls.all.centered.deepsage$ProbeName, names(exp.sums), nomatch = 0)]
plot(eqtls.matched.deepsage$OverallZScore, eqtls.all.centered.deepsage$OverallZScore, type = "n", 
     xlab = "DeepSAGE effect size (Z-scores)", ylab = "scRNA-seq effect size (Z-scores)", lwd=0,
     cex.lab = 1.4, cex.main=1.75, xlim = c(-10,10), ylim = c(-8,8))
limits = par()$usr
rect(0,0,limits[2],limits[4],col="#E6EFEA", border = NA)
rect(0,0,limits[1],limits[3],col="#E6EFEA", border = NA)
points(eqtls.matched.deepsage$OverallZScore, eqtls.all.centered.deepsage$OverallZScore,
       pch=21, bg=rgb(0,0,0,0.4), cex=log10(exp.sums.matched)+2)
box(lwd=3, col="grey")
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)
title(main = "DeepSAGE", adj = 0, font.main = 1, cex.main = 1.75, line = 0.7)

concordance.perc <- signif(sum(eqtls.matched.deepsage$OverallZScore * eqtls.all.centered.deepsage$OverallZScore > 0) / length(eqtls.matched.deepsage$OverallZScore) * 100, 3)

###########################################################################################################################
#
# Write tables
#
###########################################################################################################################
eqtls.sc.rnaseq$Z.score.bios <- eqtls.rnaseq.matched$Z.score
eqtls.sc.rnaseq$p.value.bios <- eqtls.rnaseq.matched$P.value
eqtls.sc.rnaseq$fdr.bios <- eqtls.rnaseq.matched$FDR

eqtls.sc.rnaseq$Concordant <- eqtls.rnaseq.matched$Z.score * eqtls.sc.rnaseq$OverallZScore > 0

write.table(eqtls.sc.rnaseq, file = "./data/concordance_eqtls_bios.csv", sep = "\t", quote = F, col.names=NA)
write.table(eqtls.sc.rnaseq[eqtls.rnaseq.matched$Z.score * eqtls.sc.rnaseq$OverallZScore > 0,], file = "./data/concordant_eqtls_bios_scrna.txt", sep = "\t", quote = F, col.names=NA)
write.table(eqtls.rnaseq.matched[eqtls.rnaseq.matched$Z.score * eqtls.sc.rnaseq$OverallZScore < 0,], file = "./data/discordant_eqtls_bios_original.txt", sep = "\t", quote = F, col.names=NA)

eqtls.all.centered.deepsage$Z.score.deepsage <- eqtls.matched.deepsage$OverallZScore
eqtls.all.centered.deepsage$Concordant <- eqtls.matched.deepsage$OverallZScore * eqtls.all.centered.deepsage$OverallZScore > 0
eqtls.all.centered.deepsage$PValue.deepsage <- eqtls.matched.deepsage$PValue
eqtls.all.centered.deepsage$FDR.deepSAGE <- eqtls.matched.deepsage$FDR

write.table(eqtls.all.centered.deepsage, file = "./data/concordance_eqtls_deepsage.csv", sep = "\t", quote = F, col.names=NA)