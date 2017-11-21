###########################################################################################################################
# Authors: Harm Brugge
# Name: eqtl_table
# Function: Creation of eQTL table of all cell-types and the replication in BIOS RNA-seq data
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(data.table)

###########################################################################################################################
#
# Load data
#
###########################################################################################################################
pbmc.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/pbmc.txt", header = T, sep ="\t", stringsAsFactors = F)
th.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/th-cells.txt", header = T, sep ="\t", stringsAsFactors = F) 
tc.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/tc-cells.txt", header = T, sep ="\t", stringsAsFactors = F) 
b.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/b-cells.txt", header = T, sep ="\t", stringsAsFactors = F) 
dend.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/dend-cells.txt", header = T, sep ="\t", stringsAsFactors = F) 
nk.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/nk-cells.txt", header = T, sep ="\t", stringsAsFactors = F) 
mono.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/monocytes.txt", header = T, sep ="\t", stringsAsFactors = F) 
c.mono.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/c-mono.txt", header = T, sep ="\t", stringsAsFactors = F) 
nc.mono.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/nc-mono.txt", header = T, sep ="\t", stringsAsFactors = F) 

#bios.eqtls <- data.frame(fread("~/Documents/datasets/eqtls/BIOS_rna-seq_gene_level_eQTLs.txt", header = T, sep ="\t", fill = T))
#bios.eqtls$combine <- paste(bios.eqtls$SNPName, bios.eqtls$ProbeName)
#save(bios.eqtls, file = "data/bios-eqtls.Rda")
load("data/bios-eqtls.Rda")

###########################################################################################################################
#
# Functions
#
###########################################################################################################################
combine.snp.gene <- function(eqtl.matrix) {
  eqtl.matrix$combine = paste(eqtl.matrix$SNPName, eqtl.matrix$ProbeName)
  return(eqtl.matrix)
}

add.to.eqtls <- function(eqtls, eqtls.to.add) {
  eqtls.to.add <- eqtls.to.add[!eqtls.to.add$combine %in% eqtls$combine,]
  eqtls.to.add <- check.replication(eqtls = eqtls.to.add, replication.table = bios.eqtls)
  return(rbind(eqtls, eqtls.to.add[,c(2,5,17,10,23:26)]))
}

check.replication <- function(eqtls, replication.table) {
  replication.matched <- replication.table[match(eqtls$combine, replication.table$combine, nomatch = NA),]
  replication.matched[!is.na(replication.matched$AlleleAssessed) & replication.matched$AlleleAssessed != eqtls$AlleleAssessed,]$OverallZScore <- replication.matched[!is.na(replication.matched$AlleleAssessed) & replication.matched$AlleleAssessed != eqtls$AlleleAssessed,]$OverallZScore * -1
  eqtls$Replication <- replication.matched$OverallZScore * eqtls$OverallZScore > 0
  eqtls <- cbind(eqtls, replication.matched[,c(11,22)])
  return(eqtls)
}

###########################################################################################################################
#
# Main code
#
###########################################################################################################################
pbmc.eqtls <- combine.snp.gene(pbmc.eqtls)
th.eqtls <- combine.snp.gene(th.eqtls)
tc.eqtls <- combine.snp.gene(tc.eqtls)
b.eqtls <- combine.snp.gene(b.eqtls)
dend.eqtls <- combine.snp.gene(dend.eqtls)
nk.eqtls <- combine.snp.gene(nk.eqtls)
mono.eqtls <- combine.snp.gene(mono.eqtls)
c.mono.eqtls <- combine.snp.gene(c.mono.eqtls)
nc.mono.eqtls <- combine.snp.gene(nc.mono.eqtls)

pbmc.eqtls <- check.replication(pbmc.eqtls, bios.eqtls)
eqtls <- pbmc.eqtls[,c(2,5,17,10,23:26)]
eqtls <- add.to.eqtls(eqtls, th.eqtls)
eqtls <- add.to.eqtls(eqtls, tc.eqtls)
eqtls <- add.to.eqtls(eqtls, b.eqtls)
eqtls <- add.to.eqtls(eqtls, dend.eqtls)
eqtls <- add.to.eqtls(eqtls, nk.eqtls)
eqtls <- add.to.eqtls(eqtls, mono.eqtls)

length(unique(pbmc.eqtls$ProbeName))
length(unique(th.eqtls$ProbeName))
length(unique(tc.eqtls$ProbeName))
length(unique(b.eqtls$ProbeName))
length(unique(dend.eqtls$ProbeName))
length(unique(nk.eqtls$ProbeName))
length(unique(mono.eqtls$ProbeName))
length(unique(eqtls$ProbeName))

eqtls <- add.to.eqtls(eqtls, c.mono.eqtls)
eqtls <- add.to.eqtls(eqtls, nc.mono.eqtls)

th.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/th-cells-all.txt", header = T, sep ="\t")
th.eqtls.all <- th.eqtls.all[,c(1,2,5,10,18,22)]
tc.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/tc-cells-all.txt", header = T, sep ="\t") 
tc.eqtls.all <- tc.eqtls.all[,c(1,2,5,10,18,22)]
b.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/b-cells-all.txt", header = T, sep ="\t") 
b.eqtls.all <- b.eqtls.all[,c(1,2,5,10,18,22)]
dend.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/dend-cells-all.txt", header = T, sep ="\t") 
dend.eqtls.all <- dend.eqtls.all[,c(1,2,5,10,18,22)]
nk.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/nk-cells-all.txt", header = T, sep ="\t") 
nk.eqtls.all <- nk.eqtls.all[,c(1,2,5,10,18,22)]
pbmc.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/pbmc-all.txt", header = T, sep ="\t")
pbmc.eqtls.all <- pbmc.eqtls.all[,c(1,2,5,10,18,22)]
mono.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/monocytes-all.txt", header = T, sep ="\t") 
mono.eqtls.all <- mono.eqtls.all[,c(1,2,5,10,18,22)]
c.mono.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/c-mono-all.txt", header = T, sep ="\t") 
c.mono.eqtls.all <- c.mono.eqtls.all[,c(1,2,5,10,18,22)]
nc.mono.eqtls.all <- read.table("data/eqtl/cis_100K_MAF_10/nc-mono-all.txt", header = T, sep ="\t") 
nc.mono.eqtls.all <- nc.mono.eqtls.all[,c(1,2,5,10,18,22)]

th.eqtls.all$combine <- paste(th.eqtls.all$SNPName, th.eqtls.all$ProbeName)
pbmc.eqtls.all$combine <- paste(pbmc.eqtls.all$SNPName, pbmc.eqtls.all$ProbeName)
tc.eqtls.all$combine <- paste(tc.eqtls.all$SNPName, tc.eqtls.all$ProbeName)
b.eqtls.all$combine <- paste(b.eqtls.all$SNPName, b.eqtls.all$ProbeName)
dend.eqtls.all$combine <- paste(dend.eqtls.all$SNPName, dend.eqtls.all$ProbeName)
nk.eqtls.all$combine <- paste(nk.eqtls.all$SNPName, nk.eqtls.all$ProbeName)
mono.eqtls.all$combine <- paste(mono.eqtls.all$SNPName, mono.eqtls.all$ProbeName)
c.mono.eqtls.all$combine <- paste(c.mono.eqtls.all$SNPName, c.mono.eqtls.all$ProbeName)
nc.mono.eqtls.all$combine <- paste(nc.mono.eqtls.all$SNPName, nc.mono.eqtls.all$ProbeName)

th.eqtls.all.matched <- th.eqtls.all[match(eqtls$combine, th.eqtls.all$combine),]
pbmc.eqtls.all.matched <- pbmc.eqtls.all[match(eqtls$combine, pbmc.eqtls.all$combine),]
tc.eqtls.all.matched <- tc.eqtls.all[match(eqtls$combine, tc.eqtls.all$combine),]
b.eqtls.all.matched <- b.eqtls.all[match(eqtls$combine, b.eqtls.all$combine),]
dend.eqtls.all.matched <- dend.eqtls.all[match(eqtls$combine, dend.eqtls.all$combine),]
nk.eqtls.all.matched <- nk.eqtls.all[match(eqtls$combine, nk.eqtls.all$combine),]
mono.eqtls.all.matched <- mono.eqtls.all[match(eqtls$combine, mono.eqtls.all$combine),]
c.mono.eqtls.all.matched <- c.mono.eqtls.all[match(eqtls$combine, c.mono.eqtls.all$combine),]
nc.mono.eqtls.all.matched <- nc.mono.eqtls.all[match(eqtls$combine, nc.mono.eqtls.all$combine),]

p.values <- cbind(pbmc.eqtls.all.matched$PValue, th.eqtls.all.matched$PValue, tc.eqtls.all.matched$PValue, nk.eqtls.all.matched$PValue, b.eqtls.all.matched$PValue, dend.eqtls.all.matched$PValue, mono.eqtls.all.matched$PValue, c.mono.eqtls.all.matched$PValue, nc.mono.eqtls.all.matched$PValue)
rownames(p.values) <- pbmc.eqtls.all.matched$combine
colnames(p.values) <- c("PBMC", "T CD4+", "T CD8+", "NK", "B", "DC", "Mono", "cMono", "ncMono")

correlations <- cbind(pbmc.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, th.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, tc.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, nk.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, b.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, dend.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, mono.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, c.mono.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient, nc.mono.eqtls.all.matched$IncludedDatasetsCorrelationCoefficient)
rownames(correlations) <- pbmc.eqtls.all.matched$combine
colnames(correlations) <- c("PBMC", "T CD4+", "T CD8+", "NK", "B", "DC", "Mono", "cMono", "ncMono")

significance <- cbind(pbmc.eqtls.all.matched$FDR < 0.05, th.eqtls.all.matched$FDR < 0.05, tc.eqtls.all.matched$FDR < 0.05, nk.eqtls.all.matched$FDR < 0.05, b.eqtls.all.matched$FDR < 0.05, dend.eqtls.all.matched$FDR < 0.05, mono.eqtls.all.matched$FDR < 0.05, c.mono.eqtls.all.matched$FDR < 0.05, nc.mono.eqtls.all.matched$FDR < 0.05)
significance[significance[T]] <- "*"
significance[significance == "FALSE"] <- ""
rownames(significance) <- pbmc.eqtls.all.matched$combine
colnames(significance) <- c("PBMC", "T CD4+", "T CD8+", "NK", "B", "DC", "Mono", "cMono", "ncMono")

eqtl.table <- cbind(eqtls[,1:4], p.values, significance, correlations)
write.table(eqtl.table, file = "./data/eqtl_table.tsv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

#eqtl.table <- read.table("data/eqtl_table.csv", header = T, sep="\t")

###
### Plot monocytes
###

mono.sign <- eqtl.table[eqtl.table$ncMono.1 == "*"  | eqtl.table$cMono.1 == "*",]
mono.non.sign <- eqtl.table[!(eqtl.table$ncMono.1 == "*"  | eqtl.table$cMono.1 == "*"),]

plot(mono.sign$ncMono.2, mono.sign$cMono.2, pch=20,
     ylab = "eQTL correlation coefficient cMonocytes",
     xlab= "eQTL correlation coefficient ncMonocytes",
     ylim = c(-0.9,0.9),
     xlim = c(-0.9,0.9),
     cex.lab = 1.5)
points(mono.non.sign$ncMono.2, mono.non.sign$cMono.2, pch=20, col=alpha("grey", 0.6))
points(mono.sign[mono.sign$HGNCName == "HLA-DQA1",]$cMono.2~mono.sign[mono.sign$HGNCName == "HLA-DQA1",]$ncMono.2, pch=20, col='red')
points(mono.sign[mono.sign$HGNCName == "CTSC",]$cMono.2~mono.sign[mono.sign$HGNCName == "CTSC",]$ncMono.2, pch=20, col='red')
points(mono.sign[mono.sign$HGNCName == "HLA-DRB5",]$cMono.2~mono.sign[mono.sign$HGNCName == "HLA-DRB5",]$ncMono.2, pch=20, col='red')
points(mono.sign[mono.sign$HGNCName == "RP11-1143G9.4",]$cMono.2~mono.sign[mono.sign$HGNCName == "RP11-1143G9.4",]$ncMono.2, pch=20, col='red')
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)
box(lwd=3, col = "grey")