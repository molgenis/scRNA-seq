###########################################################################################################################
# Authors: Dylan de Vries
# Name: eQTL_cohort_replication.R
# Function: Compare our eQTLs with the eQTLs identified by Kasela et al. and those found in the Blueprint cohort
###########################################################################################################################
#
# Functions
#
###########################################################################################################################

# Name: compareDatasets
# Function: finds the overlap between the input eQTLs and the cohort eQTLs and checks the concordance
# Input:
#   Name      			Type          Description
#   eQTLs     			data.frame    Data frame with the eQTL pipeline output
# 	cohortEQTLs 		data.frame	  Data frame with the eQTLs of the cohort to compare with. Column 1=SNPName, 2=Allele, 3=Direction and 4=Gene
#	eQTLGenes 			character 	  All eQTL genes from our eQTLs. Is required because one cohort uses HGNC IDs and the other ENSEMBL IDs
#
# Output:
#   Returns a data frame with all overlapping eQTLs between the two cohorts and an extra column indicating concordance
compareDatasets <- function(eQTLs, cohortEQTLs, eQTLGenes){
	overlapEQTLs = data.frame(eQTLs, concordance=eQTLs)
	concordance = rep("Not_found", nrow(overlapEQTLs))
	for (eQTLIndex in 1:nrow(overlapEQTLs)){
		matchingCohortEQTLs = cohortEQTLs[which(cohortEQTLs[,1] == overlapEQTLs$SNPName[eQTLIndex]),]
		matchingCohortEQTLs = matchingCohortEQTLs[which(matchingCohortEQTLs[,4] == eQTLGenes[eQTLIndex]),]
		if (nrow(matchingCohortEQTLs) == 0){next}
		for (matchingEQTLIndex in 1:max(nrow(matchingCohortEQTLs),1)){
			if (matchingCohortEQTLs[matchingEQTLIndex,2] == overlapEQTLs$AlleleAssessed[eQTLIndex]){
				if ((matchingCohortEQTLs[matchingEQTLIndex,3] * overlapEQTLs$OverallZScore[eQTLIndex]) > 0){
					concordance[eQTLIndex] = "Yes"
				} else {
					if (concordance[eQTLIndex] != "Yes"){
						concordance[eQTLIndex] = "No"
					}
				}
			} else {
				if ((matchingCohortEQTLs[matchingEQTLIndex,3] * overlapEQTLs$OverallZScore[eQTLIndex]) < 0){
					concordance[eQTLIndex] = "Yes"
				} else {
					if (concordance[eQTLIndex] != "Yes"){
						concordance[eQTLIndex] = "No"
					}
				}
			}
		}
	}
	overlapEQTLs$concordance = concordance
	return(overlapEQTLs)
}

# Name: getENSEMBL_ID
# Function: Change the format of the input gene to match the ENSEMBL ID format of our eQTLs
# Input:
#   Name      			Type          Description
#   gene 				character 	  An ENSEMBL ID with a sub ID. Example: ENSG00000000419.8
#
# Output:
#   Returns the gene in the format used in our eQTLs
getENSEMBL_ID <- function(gene){
	return(strsplit(gene, "\\.")[[1]][1])
}


###########################################################################################################################
#
# Main code
#
###########################################################################################################################
##
## Read in data.
##

cellTypes = c("pbmc", "b-cells", "monocytes", "nk-cells", "tc-cells", "th-cells", "dend-cells")
for (cellType in cellTypes){
	print(cellType)
	if (cellType == cellTypes[1]){
		eQTLs = data.frame(read.table(paste0("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/scrna-seq_top_eqtls/", cellType, ".txt"), header=T, stringsAsFactors=F), cellType=cellType)	
	} else {
		eQTLs = rbind(eQTLs, data.frame(read.table(paste0("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/scrna-seq_top_eqtls/", cellType, ".txt"), header=T, stringsAsFactors=F), cellType=cellType))
	}
}

kaselaCD4Data = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/datasets/kasela_CD4.txt", header=T, stringsAsFactors=F)
kaselaCD8Data = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/datasets/kasela_CD8.txt", header=T, stringsAsFactors=F)
blueprintCD4Data = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/datasets/blueprint_CD4.txt", stringsAsFactors=F)
blueprintMonocytesData = read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/datasets/blueprint_monocytes.txt", stringsAsFactors=F)

#Put gene IDs in required format
bpCD4Genes = unlist(lapply(blueprintCD4Data[,3], getENSEMBL_ID))
bpMonoGenes = unlist(lapply(blueprintMonocytesData[,3], getENSEMBL_ID))

#Kasela overlap check
overlapEQTLsKaselaCD4 = compareDatasets(eQTLs, data.frame(SNPName=kaselaCD4Data$SNPName, AlleleAssessed=kaselaCD4Data$AlleleAssessed, OverallZScore=kaselaCD4Data$OverallZScore, kaselaCD4Data$HGNCName), eQTLs$HGNCName)
overlapEQTLsKaselaCD8 = compareDatasets(eQTLs, data.frame(SNPName=kaselaCD8Data$SNPName, AlleleAssessed=kaselaCD8Data$AlleleAssessed, OverallZScore=kaselaCD8Data$OverallZScore, kaselaCD8Data$HGNCName), eQTLs$HGNCName)

write.table(overlapEQTLsKaselaCD4, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/KaselaCD4_CD8_overlap.txt", row.names=F, col.names=T, quote=F)
write.table(overlapEQTLsKaselaCD8, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/KaselaCD8_overlap.txt", row.names=F, col.names=T, quote=F)

#Blueprint overlap check
overlapEQTLsBlueprintCD4 = compareDatasets(eQTLs, data.frame(SNPName=blueprintCD4Data[,2], AlleleAssessed=substr(blueprintCD4Data[,1] , nchar(blueprintCD4Data[,1]), nchar(blueprintCD4Data[,1])), OverallZScore=blueprintCD4Data[,5], Genes=bpCD4Genes), eQTLs$ProbeName)
overlapEQTLsBlueprintMonocytes = compareDatasets(eQTLs, data.frame(SNPName=blueprintMonocytesData[,2], AlleleAssessed=substr(blueprintMonocytesData[,1] , nchar(blueprintMonocytesData[,1]), nchar(blueprintMonocytesData[,1])), OverallZScore=blueprintMonocytesData[,5], Genes=bpMonoGenes), eQTLs$ProbeName)

write.table(overlapEQTLsBlueprintCD4, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/BlueprintCD4_overlap.txt", row.names=F, col.names=T, quote=F)
write.table(overlapEQTLsBlueprintMonocytes, file="/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/replication/overlap/BlueprintMonocytes_overlap.txt", row.names=F, col.names=T, quote=F)


