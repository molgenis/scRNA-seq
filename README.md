# scRNA-seq

These scripts allow to do the process the scRNA-seq data, so that it can be used to identify eQTLs and co-expression QTLs.

## DeAnonymizer

This tool uses the method as described in Multiplexing droplet-based single cell RNA-sequencing using natural genetic barcodes, Hyun Min Kang e.a. 

It takes as input a file with the following (tab separated) info fields:

CHROM 	POS 	ALT 	REF	 cell_id 	A 	C	 G 	T	

Where CHROM and POS denote the chromosome and position of the SNP, respectively, and ALT and REF denote the alleles of the SNP. cell_id consists of the cell barcode as present in the bam-file, and A, C, G and T
denote how many unique UMIs are seen in the reads of the corresponding cell at the corresponding SNP having the corresponding bases. 

Contrary to the tool presented in the beforementioned paper, this implementation assumes a fixed error in the base calling, which can be set using the -e flag. 

All flags are explained when running the program without arguments. 

## Magic

These are the settings used to run Magic, as described by David van Dijk et al. (https://www.biorxiv.org/content/early/2017/02/25/111591)
1. n_pca_components=20
2. random_pca=True
3. t=4
4. k=9
5. ka=3
6. epsilon=1
7. rescale_percent=99

In addition, the Magic script requires an output location and an input directory. The input directory needs to have all samples in a separate directory, with the following 3 files (standard 10X output):
1. barcodes.tsv
2. genes.tsv
3. matrix.mtx

## R-scripts

Here we describe all the R scripts that were used to process the data, make plots, replicate the eQTLs and to identify co-expression QTLs. More detailed descriptions can be found in the R-scripts directory.

### Data processing
seurat_clustering.R

### Generate plots
eqtl_box_plots.R

interaction_plots.R

magic_plots.R

### eQTL replication
concordance_eqtls.R

eQTL_cohort_replication.R

eQTL_cohort_replication_make_tables.R

eqtl_table.R

### co-expression QTL identification and validation
co-expression_QTL.R

interactionOtherCellTypes.R

make_co-expression_QTL_tables.R

make_top_interaction_permutation_table.R

