###########################################################################################################################
# Authors: Harm Brugge & Dylan de Vries
# Name: interaction_plots.R
# Function: Functions for plotting interaction effects in single-cell data.
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################

library(Matrix)
library(scales)
library(ggplot2)
library(cowplot)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

# Name: plot.interaction
# Function: 
# Input:
#   Name      Type          Description
#   
#
# Output:
# 
plot.inter.extra <- function(exp.matrices, eqtl.gene, inter.gene, snp.name, flip.levels = F, p.value = NULL, r = NULL, sign.cutoff = NULL) {
  
  snp <- get.snp(snp.name)
  
  plot.matrix <- NULL
  
  i <- 1
  for (sample in exp.matrices) {
    if (!inter.gene %in% colnames(sample) | !eqtl.gene %in% colnames(sample)) next
    if (sum(sample[,eqtl.gene]) == 0 | sum(sample[,inter.gene]) == 0) next
    
    sample.matrix <- data.frame(eqtl=sample[,eqtl.gene], interaction=sample[,inter.gene], sample.name=sample.names[[i]], snp=snp[sample.names[[i]]])
    plot.matrix <- rbind(plot.matrix, sample.matrix)
    i <- i + 1
  }
  
  plots <- list()
  
  if (flip.levels) plot.matrix$snp <- factor(plot.matrix$snp, levels = levels(plot.matrix$snp)[3:1])
  
  for (i in 1:length(levels(plot.matrix$snp))) {
    genotype.to.plot <- plot.matrix[plot.matrix$snp == levels(plot.matrix$snp)[i],]
    
    color.high <- c("lightgreen", "orange", "lightblue")
    color.low <- c("darkgreen", "red3", "blue")
    
    plots[[i]] <- ggplot(genotype.to.plot, aes(y=eqtl, x=interaction, color=as.integer(sample.name), group = sample.name)) +
      geom_point(size=0.3, alpha = 0.2) +
      geom_smooth(method="lm", se=F, size = 0.5) +
      scale_colour_gradient(name = "sample.name", guide = F,
                            low = color.low[i], high = color.high[i]) +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white", colour = "grey"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey", fill=NA, size=2)) + 
      #scale_x_continuous(breaks = round(seq(min(genotype.to.plot$interaction), max(genotype.to.plot$interaction), by = 0.1),1)) +
      ylab("") + xlab("") +
      ggtitle(levels(plot.matrix$snp)[i])
  }
  
  plot.all <- ggplot(plot.matrix, aes(y=eqtl, x=interaction, color = snp)) + 
    geom_point(size=0.7, alpha = 0.4) +
    theme_minimal(base_size = 16) +
    theme(panel.background = element_rect(fill = "white", colour = "grey"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill=NA, size=2)) + 
    scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
    ylab(genes[genes$V1 == eqtl.gene,]$V2) +
    xlab(genes[genes$V1 == inter.gene,]$V2) +
    geom_smooth(method="lm") +
    ggtitle(paste(genes[genes$V1 == eqtl.gene,]$V2, "/", genes[genes$V1 == inter.gene,]$V2, "/", snp.name, "\n", signif(p.value, 3), "/", signif(r, 3), "/ cutoff=", sign.cutoff))
  
  if (length(plots) == 3){
    right.col <- plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=1)
    print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
  } else {
    right.col <- plot_grid(plots[[1]], plots[[2]], ncol=1)
    print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
  }
  #return(plot.matrix)
}

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
  return(droplevels(snp))
}

interaction.bios <- read.table("~/Downloads/rawData.txt", header = T)
ggplot(interaction.bios, aes(y=ENSG00000197728, x=ENSG00000122026, color = rs7297175)) + 
  geom_point(size=2, alpha = 0.6) +
  theme_minimal(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = "grey"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey", fill=NA, size=2)) + 
  scale_color_manual(values=c("#383bfe", "#fd7600", "#57a350"), guide = F) +
  ylab("RPS26 expression") +
  xlab("RPL21 expression") +
  geom_smooth(method="lm") +
  ggtitle("RNA-seq replication")
ggsave(filename = "plots/inter_rps26_rpl21_bios.pdf", width = 7, height = 7)


# Name: plot.correlations.qtl
# Function: 
# Input:
#   Name      Type          Description
#   
#
# Output:
# 
## TODO: Add zero values to outputData if absent
plot.correlations.qtl <- function(gene.1, gene.2, snp.id, exp.matrices, xlim = NULL, p.value, r) {
  snp <- get.snp(snp.id = snp.id)
  
  i <- 1
  outputData = data.frame(x=numeric(0), y=numeric(0), genotype=character(0))
  for (sample in samples) {
    if (gene.2 %in% colnames(sample) & gene.1 %in% colnames(sample)) {
      outputData = rbind(outputData, data.frame(x=sample[,gene.2], y=sample[,gene.1], genotype=snp[sample.names[[i]]]))
    }
    
    i <- i + 1
  }
  
  plot <- ggplot(outputData, aes(x=x, y=y, colour=genotype)) + 
    geom_point(alpha = 0.4) + 
    geom_smooth(method = "lm", fullrange = T, se=F) + 
    ggtitle(paste(genes[genes$V1 == gene.1,]$V2, snp.id, "R-squared:", r, "P-value:", p.value)) + 
    xlab(gene.2) + ylab(gene.1)
  
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim)
  }
  plot
}

# Name: plot.inter.per.sample
# Function: 
# Input:
#   Name      Type          Description
#   
#
# Output:
# 
plot.inter.per.sample <- function(exp.matrices, eqtl.gene, inter.gene, snp.name, output.path) {
  
  snp <- get.snp(snp.name)
  
  colors <- c("firebrick1", "deepskyblue3", "limegreen")
  
  par(mfrow=c(1,1))
  filename <- paste0(output.path, "/", eqtl.gene, "_", inter.gene, "_", snp.name, ".pdf")
  ggsave(filename, plot.correlations.qtl(gene.1 = eqtl.gene, gene.2 = inter.gene, snp.id = snp.name, exp.matrices = exp.matrices, p.value = NA, r = NA))
  
  pdf(paste0(output.path, "/", eqtl.gene, "_", inter.gene, "_", snp.name, "_per_sample.pdf"), onefile = T)
  for (level in 1:3) {
    par(mfrow=c(3,4))
    
    genotype.to.plot <- levels(snp)[level]
    i <- 1
    for (sample in exp.matrices) {
      if (!class(sample) == "matrix") {
        sample <- t(as.matrix(sample))
        sample[is.na(sample)] <- 0
      }
      
      if (!inter.gene %in% colnames(sample) | !eqtl.gene %in% colnames(sample)) next
      if (sum(sample[,eqtl.gene]) == 0 | sum(sample[,inter.gene]) == 0) next
      
      if (snp[sample.names[[i]]] == genotype.to.plot) {
        
        plot(sample[,eqtl.gene]~sample[,inter.gene], main = paste(sample.names[[i]], snp[sample.names[[i]]]),
             pch=20, cex=0.5, ylab = eqtl.gene, xlab = inter.gene, col=colors[level])
        mod <- lm(sample[,eqtl.gene]~sample[,inter.gene])
        
        abline(mod, col=colors[level])
      }
      i <- i + 1
    }
  }
  dev.off()
}

# Name: plot.inter.per.sample
# Function: 
# Input:
#   Name      Type          Description
#   
#
# Output:
# 
plot.all.inter.per.sample <- function(inter.matrix, exp.matrices, eqtl.file, output.path, max.items = 10, threshold = 0.4) {
  for (i in 1:ncol(inter.matrix)) {
    eqtl.gene <- colnames(inter.matrix)[i]
    snp.name <- eqtl.file[eqtl.file$ProbeName == eqtl.gene,]$SNPName
    
    
    interactions <- inter.matrix[,i]
    names(interactions) <- rownames(inter.matrix)
    interactions <- sort(interactions[interactions >= threshold],decreasing = T)
    
    if (length(interactions) > max.items) {
      interactions <- interactions[1:max.items]
    }
    
    for(interaction in names(interactions)) {
      plot.inter.per.sample(exp.matrices = exp.matrices, eqtl.gene = eqtl.gene, inter.gene = interaction, snp.name = snp.name, output.path = output.path)
    }
  }
}


###########################################################################################################################
#
# Main
#
###########################################################################################################################
##
## Load data.
## Change paths if needed.
##
dir <- list.files(path="../scRNA-seq/magic_imputation_cd4_t/", pattern="*.tsv", full.names=T, recursive=FALSE)
genotypes <- read.table("data/genotypes/maf_10.genotypes.txt", check.names = F)
genes <- read.table("~/Documents/scRNA-seq/data/lane_1/genes.tsv")
output.path <- "data/interactions/plots"
th.eqtls <- read.table("data/eqtl/cis_100K_MAF_10/th-cells.txt", stringsAsFactors = F, header = T)
load("data/interactionsNonImputed.Rda")

##
## Read expression data.
##
samples.imputed <- list()
sample.names <- vector()
i <- 1
for (file in dir) {
  sample <- read.csv(file = file, sep = "\t", row.names = 1)
  rownames(sample) <- substr(rownames(sample), start = 7, stop = 21)
  sample <- t(sample)
  
  samples.imputed[[i]] <- sample
  
  sample.names <- c(sample.names, tools::file_path_sans_ext(basename(file)))
  
  i <- i + 1
}

# Match genotypes sample names to expression data
genotypes <- genotypes[,match(sample.names, colnames(genotypes))]


###
### Raw 
###

##
## Read in data.
##
samples <- list()
cell.counts <- c()
sample.names <- vector()

dir.path <- "../scRNA-seq/data/th_cells_per_sample/samples/"
dir.samples <- list.dirs(path=dir.path, full.names=T, recursive=FALSE)

i <- 1
for (folder in dir.samples) {
  sample.names <- c(sample.names, tools::file_path_sans_ext(basename(folder)))
  sample.raw <- readMM(paste0(dir.samples[[i]], "/matrix.mtx"))
  rownames(sample.raw) <- read.table(paste0(dir.samples[[i]], "/genes.tsv"), stringsAsFactors = F)$V1
  colnames(sample.raw) <- read.table(paste0(dir.samples[[i]], "/barcodes.tsv"), stringsAsFactors = F)$V1
  samples[[i]] <- t(as.matrix(sample.raw))
  cell.counts <- c(cell.counts, ncol(sample.raw))
  
  i <- i + 1
}
genotypes <- genotypes[,match(sample.names, colnames(genotypes))]

plot.inter.extra(eqtl.gene = "ENSG00000197728", inter.gene = "ENSG00000122026", snp.name ="rs7297175", exp.matrices = samples, p.value = 0.00000, r =0.000000, flip.levels = T)
plot.inter.extra(eqtl.gene = "ENSG00000100376", inter.gene = "ENSG00000135018", snp.name ="rs5764704", exp.matrices = samples, p.value = 0.00000, r =0.000000, flip.levels = T)
plot.inter.extra(eqtl.gene = "ENSG00000183172", inter.gene = "ENSG00000145425", snp.name ="rs4147641", exp.matrices = samples, p.value = 3.774866e-07, r = 0.000000)
plot.inter.extra(eqtl.gene = "ENSG00000123219", inter.gene = "ENSG00000172349", snp.name ="rs152061", exp.matrices = samples, p.value = 0.00000, r =0.000000, flip.levels = T)
plot.inter.extra(eqtl.gene = "ENSG00000123219", inter.gene = "ENSG00000067596", snp.name ="rs152061", exp.matrices = samples, p.value = 0.00000, r =0.000000)

##
## Plot correlation boxplots per genotype
##
cor.matrix <- create.cor.matrix(exp.matrices = samples, eqtl.gene = "ENSG00000197728")
snp <- get.snp("rs7297175")
model <- lm(formula = cor.matrix["ENSG00000122026",]~snp, weights = sqrt(cell.counts))

cor.matrix["ENSG00000122026",]
ggplot(NULL, aes(x=snp, y=cor.matrix["ENSG00000122026",], group=snp)) + 
  theme_minimal(base_family = "Helvetica") +
  theme(strip.text.x = element_text(colour = "black", size = 16, family = "Helvetica"),
  title = element_text(size = 20),
  axis.title.y = element_text(size = 16, family = "Helvetica"),
  axis.text.y = element_text(size = 12, family = "Helvetica"),
  axis.text.x = element_text(size = 12, family = "Helvetica"),
  panel.background = element_rect(fill = "white", colour = "grey"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "grey", fill=NA, size=2)) +
  geom_boxplot(fill=c("#383bfe", "#fd7600", "#57a350"), outlier.shape = NA,
               notch=F, color = "black", lwd=0.6, alpha=1) +
  geom_beeswarm(cex=2.5, size = 4, shape = 21, color=rgb(0,0,0,1), fill=rgb(1,1,1,0.6)) +
  scale_color_manual(values=c("#383bfe", "#fd7600", "#57a350"), guide = F) +
  ylab("Spearman correlation") +
  xlab("Genotype") +
  ggtitle("RPS26 / RPL21 / rs7297175")
