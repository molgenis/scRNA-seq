###########################################################################################################################
# Authors: Dylan de Vries
# Name: make_co-expression_QTL_tables.R
# Function: Make the co-expression QTL table with all top interactions and all significant interactions
###########################################################################################################################
#
# Main code
#
###########################################################################################################################
##
## Read in data.
##
load("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/interactionAnalysis/nonImputed/interactionOutput.Rda")
nonImputedOutput <- interaction.output
load("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/interactionAnalysis/imputed/interactionOutput.Rda")
imputedOutput <- interaction.output

pvalMatrix <- read.table("/Users/dylandevries/Documents/work/interactionAnalysis/highPerm/nonImputed_pval_perm_matrix.txt", header=T)

##
## Find the p-value threshold for the FDR threshold
##

fdr.thresh = 0.1
pval.thresh = 0
for (pval in sort(as.numeric(as.matrix(pvalMatrix)))){
  fdr = mean(apply(pvalMatrix[,-1], 2, function(x){length(which(x <= pval))})) / length(which(pvalMatrix[,1] <= pval))
  if (fdr <= fdr.thresh){
    pval.thresh = pval
  }
}

signif.interaction.matrix <- data.frame(SNP = character(0), 
                                        eQTL.gene = character(0), 
                                        eQTL.symbol = character(0),
                                        interaction.gene = character(0), 
                                        interaction.symbol = character(0),
                                        AlleleAssessed = character(0),
                                        P.value = numeric(0),
                                        r = numeric(0),
                                        imputed.P.value = numeric(0),
                                        imputed.r = numeric(0),
                                        RNAseq.P.value = numeric(0))

top.interaction.matrix <- data.frame(SNP = character(0),
                                      eQTL.gene = character(0),
                                      eQTL.symbol = character(0),
                                      interaction.gene = character(0),
                                      interaction.symbol = character(0),
                                      AlleleAssessed = character(0),
                                      P.value = numeric(0),
                                      FDR = numeric(0),
                                      r = numeric(0))

for (i in 1:ncol(nonImputedOutput[[2]])){
  eQTL.gene <- colnames(nonImputedOutput[[2]])[i]
  eQTL.symbol <- genes[genes[,1] == eQTL.gene,2]

  SNP <- eqtl.data[eqtl.data$ProbeName == eQTL.gene, "SNPName"]
  AlleleAssessed <- eqtl.data[eqtl.data$ProbeName == eQTL.gene, "AlleleAssessed"]

  interaction.gene <- rownames(nonImputedOutput[[1]])[which.min(nonImputedOutput[[2]][-1,i])]
  interaction.symbol <- genes[genes[,1] == interaction.gene,2]

  P.value <- signif(nonImputedOutput[[2]][interaction.gene, eQTL.gene], 10)
  FDR <- mean(apply(pvalMatrix[,-1], 2, function(x){length(which(signif(x, 10) <= P.value))})) / length(which(signif(pvalMatrix[,1], 10) <= P.value))
  r <- nonImputedOutput[[1]][interaction.gene, eQTL.gene]

  top.interaction.matrix <- rbind(top.interaction.matrix, data.frame(SNP, eQTL.gene, eQTL.symbol, interaction.gene, interaction.symbol, AlleleAssessed, P.value, FDR, r))

  for (index in which(nonImputedOutput[[2]][-1,i] <= pval.thresh)){
    interaction.gene <- rownames(nonImputedOutput[[1]])[index]
    interaction.symbol <- genes[genes[,1] == interaction.gene,2]

    P.value <- nonImputedOutput[[2]][interaction.gene, eQTL.gene]
    r <- nonImputedOutput[[1]][interaction.gene, eQTL.gene]

    imputed.P.value <- imputedOutput[[2]][interaction.gene, eQTL.gene]
    imputed.r <- imputedOutput[[1]][interaction.gene, eQTL.gene]

    RNAseq.P.value <- biosRep[biosRep$gene == eQTL.gene & biosRep$interaction == interaction.gene, "interactionPvalue"][1]
    if (is.null(RNAseq.P.value)){
      RNAseq.P.value <- NA
    }
    signif.interaction.matrix <- rbind(signif.interaction.matrix, data.frame(SNP, eQTL.gene, eQTL.symbol, interaction.gene, interaction.symbol, AlleleAssessed, P.value, r, imputed.P.value, imputed.r, RNAseq.P.value))
  }
}

write.table(signif.interaction.matrix, file="signif_interaction_matrix.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(top.interaction.matrix, file="top_interaction_matrix.txt", row.names=F, col.names=T, quote=F, sep="\t")
