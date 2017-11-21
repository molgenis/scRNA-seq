###########################################################################################################################
# Authors: Harm Brugge & Dylan de Vries
# Name: co-expression_QTL.R
# Function: Calculate all interactions between a list of eQTL genes and all other expressed genes, using MAGIC imputed
#           expression data
###########################################################################################################################
#
# Main code
#
###########################################################################################################################
load("interactionOutput.Rda")

for (type in c("imputed", "nonImputed")){
  setwd(paste0("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/interactionAnalysis/", type))
  load("/interactionOutput.Rda")

  files = dir("./permutations/")
  pvalMatrix = cbind(apply(interaction.output[[2]], 2, function(x){min(x[-1])}))
  permPvalues = NULL
  for (file in files){
    if (!strsplit(file, "_")[[1]][1] %in% colnames(interaction.output[[2]])){next}
    permMatrix = read.table(paste0("./permutations/", file), header=T)
    permPvalues = rbind(permPvalues, apply(permMatrix, 2, min))
    print(file)
  }

  pvalMatrix = cbind(pvalMatrix, permPvalues)


  for (i in 1:ncol(pvalMatrix)){
    if (i == 1){
      minOrder = order(pvalMatrix[,1])
      rownames(pvalMatrix) = rownames(pvalMatrix)[minOrder]
      pvalMatrix[,1] = pvalMatrix[minOrder,1]
    } else {
      pvalMatrix[,i] = sort(pvalMatrix[,i])
    }
  }

  colnames(pvalMatrix) = c("real", paste0("perm", 1:length(files)))

  write.table(pvalMatrix, file=paste0(type, "_pval_perm_matrix.txt"), row.names=T, col.names=T, quote=F, sep="\t")
}
