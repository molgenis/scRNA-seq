###########################################################################################################################
# Authors: Harm Brugge & Dylan de Vries
# Name: co-expression_QTL.R
# Function: Calculate all interactions between a list of eQTL genes and all other expressed genes, using MAGIC imputed
#           expression data
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(pbapply)
require(broom)

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
  snp[snp == "C/T"] <- "T/C"
  snp[snp == "C/G"] <- "G/C"
  snp[snp == "C/A"] <- "A/C"
  snp[snp == "T/G"] <- "G/T"
  snp[snp == "T/A"] <- "A/T"
  snp[snp == "T/C"] <- "C/T"
  snp[snp == "G/A"] <- "A/G"
  snp[snp == "G/C"] <- "C/G"
  snp[snp == "G/T"] <- "T/G"
  return(as.factor(snp))
}

# Name: create.cor.matrix
# Function: calculate the correlations between the input gene and all other genes, for every person
# Input:
#   Name            Type          Description
#   exp.matrices    list          A list with an expression matrix for every person
#   eqtl.gene       character     The gene for which you want to calculate the correlations
#   cor.method      character     The method to use in the cor function. Default = "spearman"
#
# Output:
#   A correlation matrix with the correlation value between all genes and the input gene, for every person

create.cor.matrix <- function(exp.matrices, eqtl.gene, cor.method = "spearman") {
  cor.vectors <- list()
  for (i in 1:length(exp.matrices)) {
    #Check if the gene is in the expression matrix, since MAGIC does not include genes in its output that have no expression
    if (eqtl.gene %in% colnames(exp.matrices[[i]])) {
      samp.var <- apply(exp.matrices[[i]], 2, var)
      samp.cor <- cor(exp.matrices[[i]][,eqtl.gene], exp.matrices[[i]][,samp.var != 0], method = cor.method)
      cor.vectors[[i]] <- t(samp.cor)
    } else {
      #If there is no expression available for the gene, there is also no correlation and it is set to 0
      cor.vectors[[i]] <- matrix(NA, nrow = ncol(exp.matrices[[i]]), ncol = 1, dimnames = list(colnames(exp.matrices[[i]]), NA))
    }
  }
  #Combine all correlations within a single matrix. Since order of genes is different we can't just cbind them
  genes <- unique(unlist(lapply(cor.vectors, rownames)))
  cor.matrix <- matrix(NA, nrow = length(genes), ncol = length(cor.vectors), 
                       dimnames = list(genes, sample.names))
  for (i in seq_along(cor.vectors)) {
    cor.matrix[rownames(cor.vectors[[i]]), i] <- cor.vectors[[i]]
  }
  # #Remove the target gene itself
  # cor.matrix <- cor.matrix[-which(rownames(cor.matrix) == eqtl.gene),]
  print(table(apply(cor.matrix, 1, function(x){length(which(any(is.na(x))))})))
  return(cor.matrix)
}

# Name: create.cor.matrices
# Function: write a correlation matrix for every eQTL gene to a separate file
# Input:
#   Name            Type          Description
#   eqtl.data       data.frame    The eQTL pipeline output
#   exp.matrices    list          A list with an expression matrix for every person
#   output.dir      character     Path to the output directory
#   cor.method      character     The method to use in the cor function. Default = "spearman"
#
# Output:
#   A file with the correlation matrix, for every eQTL separately

create.cor.matrices <- function(eqtl.data, exp.matrices, output.dir, cor.method = "spearman"){
  pb <- txtProgressBar(min = 0, max = nrow(eqtl.data), style = 3)
  setTxtProgressBar(pb, 0)
  
  for (i in 1:nrow(eqtl.data)) {
    cor.matrix <- create.cor.matrix(exp.matrices = exp.matrices, eqtl.gene = eqtl.data[i,"ProbeName"])
    write.table(cor.matrix, file=paste0(output.dir, "/correlation_matrix_", eqtl.data[i,"ProbeName"], ".txt"), row.names=T, col.names=T, quote=F)
    setTxtProgressBar(pb, i)
  }
}  

# Name: interaction.regression
# Function: calculate the interaction for every gene in the correlation matrix with the eQTL gene
# Input:
#   Name            Type          Description
#   cor.matrix      matrix        A matrix with a correlation value per person, for the eQTL gene
#   eqtl.gene       character     The gene for which you want to calculate the interactions
#   snp             factor        The genotype for the eQTL SNP for every person
#
# Output:
#   A matrix with the interaction statistics for every gene with the eQTL gene

interaction.regression <- function(cor.matrix, eqtl.gene, snp) {
  interaction.statistics <- do.call("rbind", apply(cor.matrix, 1, function(x) {
    model <- lm(formula = x~snp, weights = sqrt(cell.counts))
    return(tidy(model)[2,])
  }))
  
  return(interaction.statistics)
}

# Name: do.interaction.analysis
# Function: calculate the interaction for every gene in the correlation matrix with the eQTL gene
# Input:
#   Name            Type          Description
#   eqtl.data       data.frame    The eQTL pipeline output
#   exp.matrices    list          A list with an expression matrix for every person
#   output.dir      character     Path to the output directory
#   cor.dir         character     Path to the directory with the correlation matrix files
#   permuations     logical       Parameter to indicate whether to do permutations and calculate the FDR. Default = FALSE
#   n.perm          numeric       Number of permutations to do if permutations parameter is true. Default = 20
#   fdr.thresh      numeric       The FDR significance threshold
#   perm.type       character     Indicates to do permuations per gene separately or combined, for values "gene" and "all" respectively
#
# Output:
#   A list with two matrices. The first matrix has the R values of the interaction model for every eQTL gene with every other gene.
#   The second matrix has the p-values for every interaction per eQTL gene and the p-value significance threshold to get an FDR of 0.05

do.interaction.analysis <- function(eqtl.data, exp.matrices, output.dir, cor.dir, permutations = F, n.perm = 20, fdr.thresh = 0.05, perm.type = "gene") {
  if (permutations) {
    #dir.create(paste0(output.dir, "/permutations"))
    perm.sample.orders <- list()
    for (i in 1:n.perm) {
      perm.sample.orders[[i]] <- sample(1:length(exp.matrices), length(exp.matrices), replace = F)
    }
  }
  
  r.matrix <- NULL
  p.value.matrix <- NULL
  p.value.permuted <- list()
  p.value.thresholds <- NULL
  
  eqtl.genes <- NULL
  pb <- txtProgressBar(min = 0, max = nrow(eqtl.data), style = 3)
  setTxtProgressBar(pb, 0)
  
  for (i in 1:nrow(eqtl.data)) {
    eqtl <- eqtl.data[i,]
    eqtl <- unlist(eqtl)

    snp <- as.numeric(get.snp(eqtl["SNPName"]))
    eqtl.name <- eqtl["ProbeName"]
    cor.matrix <- read.table(paste0(cor.dir, "/correlation_matrix_", eqtl.name, ".txt"), row.names=1, header=T, stringsAsFactors=F)

    #Remove all rows for which a correlation cannot be calculated within 1 or more individuals
    cor.matrix <- cor.matrix[apply(cor.matrix, 1, function(x){!any(is.na(x))}),]
    
    print(dim(cor.matrix))

    if (nrow(cor.matrix) == 0){
      p.value.permuted[[i]] <- matrix(NA, ncol=n.perm)
      next
    } else {
      eqtl.genes <- c(eqtl.genes, eqtl.name)
    }

    interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl.name, snp = snp)
    #Calculate the R value from the T statistic
    r.matrix <- cbind(r.matrix, interaction.statistics$statistic / sqrt(length(snp) - 2 + interaction.statistics$statistic ** 2))
    p.value.matrix <- cbind(p.value.matrix, interaction.statistics$p.value)
    
    if (permutations) {
      for (current.perm in 1:n.perm) {
        permuted.snp <- snp[perm.sample.orders[[current.perm]]]
        perm.interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl.name, snp = permuted.snp)
        if (current.perm == 1){
          p.value.permuted[[i]] <- perm.interaction.statistics$p.value
        } else {
          p.value.permuted[[i]] <- cbind(p.value.permuted[[i]], perm.interaction.statistics$p.value)
        }
      }
      if (perm.type == "gene"){
        write.table(p.value.permuted[[i]], file=paste0(output.dir, "/permutations/", eqtl.name, "_permutations.txt"))
        p.value.thresh <- 0
        p.values <- unique(sort(interaction.statistics$p.value, decreasing=F))
        for (current.p.value.thresh in p.values){
          signif.interactions <- length(which(interaction.statistics$p.value <= current.p.value.thresh))
          permuted.signif.interactions <- c()
          for (current.perm in 1:n.perm){
            permuted.signif.interactions <- c(permuted.signif.interactions, length(which(p.value.permuted[[i]][,current.perm] <= current.p.value.thresh)))
          }
          if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
            break
          }
          p.value.thresh <- current.p.value.thresh
        }
        p.value.thresholds <- c(p.value.thresholds, p.value.thresh)
      }
    }
    setTxtProgressBar(pb, i)
    
  }

  rownames(p.value.matrix) <- rownames(cor.matrix)
  rownames(r.matrix) <- rownames(cor.matrix)

  colnames(p.value.matrix) <- eqtl.genes
  colnames(r.matrix) <- eqtl.genes

  if (permutations & perm.type == "gene"){
    p.value.matrix <- rbind(p.value.thresholds, p.value.matrix)
    rownames(p.value.matrix)[1] = "significance_threshold"
    interaction.list <- list(r.matrix, p.value.matrix)
  } else if (permutations & perm.type == "all"){
    save(p.value.permuted, file=paste0(output.dir, "/permutations/permutedPValue.Rda"))
    p.value.thresh <- 0
    for (current.p.value.thresh in unique(sort(p.value.matrix, decreasing=F))){
      signif.interactions <- length(which(p.value.matrix <= current.p.value.thresh))
      permuted.signif.interactions <- c()
      for (current.perm in 1:n.perm){
        current.perm.p.values <- unlist(lapply(p.value.permuted, function(x){return(x[,current.perm])}))
        permuted.signif.interactions <- c(permuted.signif.interactions, length(which(current.perm.p.values <= current.p.value.thresh)))
      }
      if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
        break
      }
      p.value.thresh <- current.p.value.thresh
    }
    interaction.list <- list(r.matrix, p.value.matrix, p.value.thresh)
  } else {
    interaction.list <- list(r.matrix, p.value.matrix)
  }
  
  
  close(pb)
  
  return(interaction.list)
}

###########################################################################################################################
#
# Main code
#
###########################################################################################################################
##
## Read in data (change files).
##
genes <- read.table("/Users/dylandevries/Documents/work/interactionAnalysis/calls/genes.tsv")
genotypes <- read.table("/Users/dylandevries/Documents/work/interactionAnalysis/calls/maf_10.calls.txt", check.names = F)
eqtl.data <- read.table("/Users/dylandevries/Documents/work/interactionAnalysis/calls/eQTLProbesFDR0.05-ProbeLevel.txt", stringsAsFactors = F, header = T)

##
## Read in the expression data.
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


dir.path <- "/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/data/expression/thCellsPerSample/"
dir <- list.dirs(path=dir.path, full.names=T, recursive=FALSE)

i <- 1
for (folder in dir) {
  print(i)
  print(paste0(dir[[i]], "/matrix.mtx"))
  sample.names <- c(sample.names, tools::file_path_sans_ext(basename(folder)))
  sample.raw <- readMM(paste0(dir[[i]], "/matrix.mtx"))
  rownames(sample.raw) <- read.table(paste0(dir[[i]], "/genes.tsv"), stringsAsFactors = F)$V1
  colnames(sample.raw) <- read.table(paste0(dir[[i]], "/barcodes.tsv"), stringsAsFactors = F)$V1
  exp.matrices[[i]] <- t(as.matrix(sample.raw))
  cell.counts <- c(cell.counts, ncol(sample.raw))
  
  i <- i + 1
}



genotypes <- genotypes[,match(sample.names, colnames(genotypes))]

cor.dir = "/Users/dylandevries/Documents/work/interactionAnalysis/permutedPerGene/correlationMatricesNew/"
output.dir = "/Users/dylandevries/Documents/work/interactionAnalysis/permutedPerGene/"

if (!file.exists(paste0(output.dir, "/permutations"))){
    dir.create(paste0(output.dir, "/permutations"))
}
if (!file.exists(cor.dir)){
    dir.create(cor.dir)
}

create.cor.matrices(eqtl.data = eqtl.data, exp.matrices = exp.matrices, output.dir = cor.dir)
interaction.output <- do.interaction.analysis(eqtl.data = eqtl.data, exp.matrices = exp.matrices, output.dir = output.dir, cor.dir = cor.dir, permutations = T, n.perm=10)
save(interaction.output, file=paste0(output.dir, "/interactionOutput.Rda"))





