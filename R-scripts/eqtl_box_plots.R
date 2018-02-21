###########################################################################################################################
# Authors: Harm Brugge
# Name: eqtl_box_plots
# Function: Creation of eQTL boxplots
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(beeswarm)
library(ggplot2)
library(ggbeeswarm)
library(data.table)

###########################################################################################################################
#
# Load data
#
###########################################################################################################################
genotypes <- read.table("data/genotypes/maf_10.genotypes.txt", check.names = F)
eqtl.table <- read.table("data/eqtl_table.csv", header = T, sep="\t")

pbmc.expression <- read.table("data/pilot3_traitfile_ensembl_all_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
b.expression <- read.table("data/pilot3_traitfile_ensembl_b_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
th.expression <- read.table("data/pilot3_traitfile_ensembl_th_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
tc.expression <- read.table("data/pilot3_traitfile_ensembl_tc_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
nk.expression <- read.table("data/pilot3_traitfile_ensembl_nk_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
dend.expression <- read.table("data/pilot3_traitfile_ensembl_dend_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
mono.expression <- read.table("data/pilot3_traitfile_ensembl_mono_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
c.mono.expression <- read.table("data/pilot3_traitfile_ensembl_c_mono_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)
nc.mono.expression <- read.table("data/pilot3_traitfile_ensembl_nc_mono_scaled.tsv", header = T, sep="\t", check.names = F, row.names = 1)

genotypes <- genotypes[,colnames(genotypes) %in% colnames(pbmc.expression)]
genotypes <- genotypes[,match(colnames(pbmc.expression), colnames(genotypes))]

genes <- read.table("../scRNA-seq/data/lane_1/genes.tsv")

###########################################################################################################################
#
# Plotting
#
###########################################################################################################################
plot.all <- function(snp, gene, flip.levels = F) {
  snp.name <- snp
  
  snp <- unlist(genotypes[snp,])
  
  snp[snp == "C/A"] <- "A/C"
  snp[snp == "T/A"] <- "A/T"
  snp[snp == "T/C"] <- "C/T"
  snp[snp == "G/A"] <- "A/G"
  snp[snp == "G/C"] <- "C/G"
  snp[snp == "G/T"] <- "T/G"
  
  snp <- droplevels(snp)
  if (flip.levels) snp = factor(snp, levels(snp)[c(3,2,1)])
  
  eqtl <- eqtl.table[eqtl.table$SNPName == snp.name & eqtl.table$ProbeName == gene,]
  correlations <- paste("r =", sprintf("%.2f", round(eqtl[23:29], digits = 2)), unlist(eqtl[14:20]))
  
  data <- data.frame(snp= snp, expression= unlist(pbmc.expression[gene,]), cell.type= "PBMC", color="black")
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(th.expression[gene,]), cell.type= "CD 4+ T", color="#153057"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(tc.expression[gene,]), cell.type= "CD 8+ T", color="#009ddb"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(nk.expression[gene,]), cell.type= "NK", color="#e64b50"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(mono.expression[gene,]), cell.type= "Monocytes", color="#edba1b"))
  data <- rbind.data.frame(data, data.frame(snp= snp[colnames(pbmc.expression) %in% colnames(b.expression)], expression= unlist(b.expression[gene,]), cell.type= "B", color="#71bc4b"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(dend.expression[gene,]), cell.type= "DC", color="#965ec8"))
  
  colors = c(rep("lightgrey", 3), rep("#153057",3), rep("#009ddb",3), rep("#e64b50",3), rep("#edba1b",3), rep("#71bc4b",3), rep("#965ec8",3))
  colors.strip = c("black","#153057","#009ddb","#e64b50","#edba1b","#71bc4b","#965ec8")
  
  gene.name <- as.character(genes[genes$V1 == gene,]$V2)

  ggplot(data, aes(x=snp, y=expression, group=snp)) +
    geom_boxplot(notch=F, color = "black", outlier.shape=NA, fill= colors, lwd=0.6, alpha=1) + 
    theme_minimal(base_family = "Helvetica") +
    theme(strip.text.x = element_text(colour = "black", size = 16, family = "Helvetica"),
          title = element_text(size = 20),
          axis.title.y = element_text(size = 16, family = "Helvetica"),
          axis.text.y = element_text(size = 12, family = "Helvetica"),
          axis.text.x = element_text(size = 12, family = "Helvetica"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill=NA, size=2)) +
    facet_wrap(~cell.type, ncol = length(levels(data$cell.type)) ) +
    geom_quasirandom(cex=2, size = 2, shape = 21, color=rgb(0,0,0,1), fill=rgb(1,1,1,0.6)) +
    ggtitle(substitute(italic(x)~ " / "~y, list(x=gene.name, y=snp.name))) +
    ylab("Expression") +
    xlab("") +
    geom_label(data=data.frame(x=2, y=max(data$expression), label=correlations,
                               cell.type=c("PBMC","CD 4+ T","CD 8+ T","NK","Monocytes","B","DC")),
               aes(x,y,label=label), inherit.aes=FALSE, colour = "white") +
    geom_text(data=data.frame(x=2, y=max(data$expression), label=correlations,
                              cell.type=c("PBMC","CD 4+ T","CD 8+ T","NK","Monocytes","B","DC")),
              aes(x,y,label=label), inherit.aes=FALSE, size=6)
}

# plot.all("rs7142883","ENSG00000166165")
# ggsave()
plot.all("rs4804315", "ENSG00000133250")
ggsave(filename = "plots/eqtl_znf414.pdf", width = 14, height = 5)
plot.all("rs2272245",	"ENSG00000106537")
ggsave(filename = "plots/eqtl_tspan13.pdf", width = 14, height = 5)

plot.mono <- function(snp, gene) {
  snp.name <- snp
  snp <- unlist(genotypes[snp,])
  
  snp[snp == "C/T"] <- "T/C"
  snp[snp == "C/G"] <- "G/C"
  snp[snp == "C/A"] <- "A/C"
  snp[snp == "T/G"] <- "G/T"
  snp[snp == "T/A"] <- "A/T"
  snp[snp == "T/C"] <- "C/T"
  snp[snp == "G/A"] <- "A/G"
  snp[snp == "G/C"] <- "C/G"
  snp[snp == "G/T"] <- "T/G"
  
  snp <- droplevels(snp)
  #snp = factor(snp, levels(snp)[c(3,2,1)])
  
  eqtl <- eqtl.table[eqtl.table$SNPName == snp.name & eqtl.table$ProbeName == gene,]
  correlations <- paste("r =",   sprintf("%.2f", round(eqtl[c(23,29:31)], digits = 2)), unlist(eqtl[ c(14,20:22) ]))
  
  data <- data.frame(snp= snp, expression= unlist(pbmc.expression[gene,]), cell.type= "PBMC", color="black")
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(mono.expression[gene,]), cell.type= "Monocytes", color="#edba1b"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(c.mono.expression[gene,]), cell.type= "cMonocytes", color="#edba1b"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(nc.mono.expression[gene,]), cell.type= "ncMonocytes", color="#edba1b"))
  
  colors = c(rep("lightgrey", 3), rep("#edba1b",9))
  
  gene.name <- as.character(genes[genes$V1 == gene,]$V2)
  
  ggplot(data, aes(x=snp, y=expression, group=snp)) +
    geom_boxplot(notch=F, color = "black", outlier.shape=NA, fill= colors, lwd=0.6, alpha=1) + 
    theme_minimal(base_family = "Helvetica") +
    theme(strip.text.x = element_text(colour = "black", size = 16, family = "Helvetica"),
          title = element_text(size = 20),
          axis.title.y = element_text(size = 16, family = "Helvetica"),
          axis.text.y = element_text(size = 12, family = "Helvetica"),
          axis.text.x = element_text(size = 12, family = "Helvetica"),
          panel.background = element_rect(fill = "white", colour = "grey"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill=NA, size=2)) +
    facet_wrap(~cell.type, ncol = length(levels(data$cell.type)) ) +
    geom_quasirandom(cex=2, size = 2, shape = 21, color=rgb(0,0,0,1), fill=rgb(1,1,1,0.6)) +
    ggtitle(substitute(italic(x)~ " / "~y, list(x=gene.name, y=snp.name))) +
    ylab("Expression") +
    xlab("") +
    geom_text(data=data.frame(x=2, y=max(data$expression, na.rm = T), label=correlations, 
                              cell.type=c("PBMC","Monocytes","cMonocytes","ncMonocytes")), 
              aes(x,y,label=label), inherit.aes=FALSE, size=6)
}

plot.mono("rs10878967", "ENSG00000257764") # RP11
ggsave(filename = "plots/eqtl_RP11.pdf", width = 8, height = 5)

plot.mono("rs61514665", "ENSG00000109861") # CTSC
ggsave(filename = "plots/eqtl_CTSC.pdf", width = 8, height = 5)

plot.mono("rs111631325", "ENSG00000198502") # HLA-DRB5
ggsave(filename = "plots/eqtl_HLA-DRB5.pdf", width = 8, height = 5)

plot.mono("rs4725361", "ENSG00000106565") # THEM176B
ggsave(filename = "plots/eqtl_THEM176B.pdf", width = 8, height = 5)

plot.mono("rs116232857", "ENSG00000196735") # HLA-DQA1
ggsave(filename = "plots/eqtl_hla-dqa1.pdf", width = 8, height = 5)
