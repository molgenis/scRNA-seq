###########################################################################################################################
# Authors: Dylan de Vries
# Name: eQTL_cohort_replication_make_tables.R
# Function: Make the output replication table
###########################################################################################################################
#
# Functions
#
###########################################################################################################################
# Name: getEQTL
# Function: Get the eQTL information
# Input:
#   Name      			Type          Description
#   SNP 				character 	  The SNP rs ID
#
# Output:
#   Vector with the required eQTL information
getEQTL <- function(SNP){
	data = strsplit(SNP, ",")[[1]]
	eqtl = eQTLData[eQTLData$SNPName == data[1] & eQTLData$HGNCName == data[2],]
	return(c(eqtl$FDR, eqtl$PValue, eqtl$IncludedDatasetsCorrelationCoefficient))
}


###########################################################################################################################
#
# Main code
#
###########################################################################################################################
##
## Read in data.
##

overlapEQTLsKaselaCD8 = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/KaselaCD8_overlap.txt", header=T, stringsAsFactors=F)
overlapEQTLsKaselaCD4 = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/KaselaCD4_CD8_overlap.txt", header=T, stringsAsFactors=F)
overlapEQTLsBlueprintMonocytes = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/BlueprintMonocytes_overlap.txt", header=T, stringsAsFactors=F)
overlapEQTLsBlueprintCD4 = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/BlueprintCD4_overlap.txt", header=T, stringsAsFactors=F)

kaselaCD8Table = data.frame(SNP=overlapEQTLsKaselaCD8$SNPName, Gene=overlapEQTLsKaselaCD8$HGNCName, AlleleAssessed=overlapEQTLsKaselaCD8$AlleleAssessed, PBMC.Pvalue=0, CD4.Pvalue=0, CD8.Pvalue=0, NK.Pvalue=0, Monocytes.Pvalue=0, B.Pvalue=0, DC.Pvalue=0,PBMC.FDR=0, CD4.FDR=0, CD8.FDR=0, NK.FDR=0, Monocytes.FDR=0, B.FDR=0, DC.FDR=0, PBMC.Cor=0, CD4.Cor=0, CD8.Cor=0, NK.Cor=0, Monocytes.Cor=0, B.Cor=0, DC.Cor=0, Concordance=overlapEQTLsKaselaCD8$concordance)

kaselaCD4Table = data.frame(SNP=overlapEQTLsKaselaCD4$SNPName, Gene=overlapEQTLsKaselaCD4$HGNCName, AlleleAssessed=overlapEQTLsKaselaCD4$AlleleAssessed, PBMC.Pvalue=0, CD4.Pvalue=0, CD8.Pvalue=0, NK.Pvalue=0, Monocytes.Pvalue=0, B.Pvalue=0, DC.Pvalue=0,PBMC.FDR=0, CD4.FDR=0, CD8.FDR=0, NK.FDR=0, Monocytes.FDR=0, B.FDR=0, DC.FDR=0, PBMC.Cor=0, CD4.Cor=0, CD8.Cor=0, NK.Cor=0, Monocytes.Cor=0, B.Cor=0, DC.Cor=0, Concordance=overlapEQTLsKaselaCD4$concordance)

blueprintMonocytesTable = data.frame(SNP=overlapEQTLsBlueprintMonocytes$SNPName, Gene=overlapEQTLsBlueprintMonocytes$HGNCName, AlleleAssessed=overlapEQTLsBlueprintMonocytes$AlleleAssessed, PBMC.Pvalue=0, CD4.Pvalue=0, CD8.Pvalue=0, NK.Pvalue=0, Monocytes.Pvalue=0, B.Pvalue=0, DC.Pvalue=0,PBMC.FDR=0, CD4.FDR=0, CD8.FDR=0, NK.FDR=0, Monocytes.FDR=0, B.FDR=0, DC.FDR=0, PBMC.Cor=0, CD4.Cor=0, CD8.Cor=0, NK.Cor=0, Monocytes.Cor=0, B.Cor=0, DC.Cor=0, Concordance=overlapEQTLsBlueprintMonocytes$concordance)

blueprintCD4Table = data.frame(SNP=overlapEQTLsBlueprintCD4$SNPName, Gene=overlapEQTLsBlueprintCD4$HGNCName, AlleleAssessed=overlapEQTLsBlueprintCD4$AlleleAssessed, PBMC.Pvalue=0, CD4.Pvalue=0, CD8.Pvalue=0, NK.Pvalue=0, Monocytes.Pvalue=0, B.Pvalue=0, DC.Pvalue=0,PBMC.FDR=0, CD4.FDR=0, CD8.FDR=0, NK.FDR=0, Monocytes.FDR=0, B.FDR=0, DC.FDR=0, PBMC.Cor=0, CD4.Cor=0, CD8.Cor=0, NK.Cor=0, Monocytes.Cor=0, B.Cor=0, DC.Cor=0, Concordance=overlapEQTLsBlueprintCD4$concordance)

cellTypes = c("pbmc", "b-cells", "monocytes", "nk-cells", "tc-cells", "th-cells", "dend-cells")
eQTLfiles = c(paste0("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/scrna-seq_eqtls/", cellTypes, "-all.txt"))
cellTypes = c("PBMC", "B", "Monocytes", "NK", "CD8", "CD4", "DC")


for (i in 1:length(eQTLfiles)){
	print(i)
	eQTLData = read.table(eQTLfiles[i], header=T, stringsAsFactors=F)

	kaselaCD8Info = do.call("rbind", lapply(paste(kaselaCD8Table$SNP, kaselaCD8Table$Gene, sep=","), getEQTL))
	kaselaCD8Table[,paste0(cellTypes[i], ".FDR")] = kaselaCD8Info[,1]
	kaselaCD8Table[,paste0(cellTypes[i], ".Pvalue")] = kaselaCD8Info[,2]
	kaselaCD8Table[,paste0(cellTypes[i], ".Cor")] = kaselaCD8Info[,3]

	kaselaCD4Info = do.call("rbind", lapply(paste(kaselaCD4Table$SNP, kaselaCD4Table$Gene, sep=","), getEQTL))
	kaselaCD4Table[,paste0(cellTypes[i], ".FDR")] = kaselaCD4Info[,1]
	kaselaCD4Table[,paste0(cellTypes[i], ".Pvalue")] = kaselaCD4Info[,2]
	kaselaCD4Table[,paste0(cellTypes[i], ".Cor")] = kaselaCD4Info[,3]

	blueprintMonocytesInfo = do.call("rbind", lapply(paste(blueprintMonocytesTable$SNP, blueprintMonocytesTable$Gene, sep=","), getEQTL))
	blueprintMonocytesTable[,paste0(cellTypes[i], ".FDR")] = blueprintMonocytesInfo[,1]
	blueprintMonocytesTable[,paste0(cellTypes[i], ".Pvalue")] = blueprintMonocytesInfo[,2]
	blueprintMonocytesTable[,paste0(cellTypes[i], ".Cor")] = blueprintMonocytesInfo[,3]

	blueprintCD4Info = do.call("rbind", lapply(paste(blueprintCD4Table$SNP, blueprintCD4Table$Gene, sep=","), getEQTL))
	blueprintCD4Table[,paste0(cellTypes[i], ".FDR")] = blueprintCD4Info[,1]
	blueprintCD4Table[,paste0(cellTypes[i], ".Pvalue")] = blueprintCD4Info[,2]
	blueprintCD4Table[,paste0(cellTypes[i], ".Cor")] = blueprintCD4Info[,3]
}

write.table(kaselaCD8Table, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/kaselaCD8OverlapTable.txt", col.names=T, row.names=F, quote=F)
write.table(kaselaCD4Table, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/kaselaCD4OverlapTable.txt", col.names=T, row.names=F, quote=F)
write.table(blueprintCD4Table, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/blueprintCD4OverlapTable.txt", col.names=T, row.names=F, quote=F)
write.table(blueprintMonocytesTable, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/blueprintMonocytesOverlapTable.txt", col.names=T, row.names=F, quote=F)

