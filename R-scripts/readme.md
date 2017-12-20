# R scripts

## Data processing

seurat_clustering.R - This script is used to process the output of Cell Ranger from all lanes of the chip into a single data object, which will be used for the rest of the analyses.
This includes several QC steps, clustering, cell type assignment and linking the DeAnonimyzer output.
It uses the output from Cell Ranger and the output from DeAnonimyzer and then writes the output to a file.

## Generate plots

eqtl_box_plots.R - This makes the eQTL box plots found in figure 1.B, 1.D and 1.E, using the eQTL, genotypes and expression per cell type files.
The plot.all function creates the plots used in 1.B and the plot.mono creates the plots used in 1.D and 1.E.
Both plotting functions require you to provide the gene and the SNP ID.

interaction_plots.R - Here we generate the plots used in figure 2.A, 2.B and 2.C in our paper.
The script uses the expression, genotype, genes and eQTL files as input.
Additionally, the script requires you to provide the two genes for which there is an interaction and the SNP ID.

magic_plots.R - Generates plots (supplementary figure 3) to determine the effect of MAGIC imputation on the gene expression in the CD4+ T cells.

### eQTL replication
concordance_eqtls.R - Script used to determine the concordance of the eQTLs found by doing an eQTL mapping confining to top eQTLs found in the RNA-seq and deepSAGE data. Corresponds to figure 1.A and supplementary table 1.

eQTL_cohort_replication.R - Compares our eQTLs with the eQTLs identified by Kasela et al., PLOS Gen. (2017) (<http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006643>) and those found by the Blueprint cohort (<ftp://ftp.ebi.ac.uk/pub/databases/blueprint/blueprint_Epivar/qtl_as/QTL_RESULTS/>).

eQTL_cohort_replication_make_tables.R - Makes the table in which our eQTLs are compared to those found by Kasela et al. and those found in the Blueprint cohort.

eqtl_table.R - Generates a table (supplementary table 2) with the top eQTLs of all cell-types and their replication in the BIOS RNA-seq data. Also used to plot figure 1.C

### co-expression QTL identification and validation
co-expression_QTL.R - Creates the correlation matrices for the eQTL genes and looks for interactions between genes. This script requires the genotype, gene ID conversion, expression and eQTL file as input. It will look for possible co-expression QTLs between every eQTL gene and all other genes for which there is variance within the expression in all samples.

interactionOtherCellTypes.R - Calculates and plots all interactions for all cell types except the CD4+ T cells. It uses the Seurat object and the genotype file.

make_co-expression_QTL_tables.R - Makes the co-expression QTL table with all top interactions and all significant interactions.
It uses the output from the co-expression_QTL.R script and the p-value matrix to make a comprehensive table with all information about each co-expression QTL.
It requires the output from the co-expression_QTL.R, the BIOS RNA-seq data and the eQTL data.

make_top_interaction_permutation_table.R - Creates the p-value matrix with the p-values of the strongest interactions for every eQTL SNP and the strongest interaction found in the permutations.
