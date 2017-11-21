# DeAnonymizer

This tool uses the method as described in Multiplexing droplet-based single cell RNA-sequencing using natural genetic barcodes, Hyun Min Kang e.a. 

It takes as input a file with the following (tab separated) info fields:

CHROM 	POS 	ALT 	REF	 cell_id 	A 	C	 G 	T	

Where CHROM and POS denote the chromosome and position of the SNP, respectively, and ALT and REF denote the alleles of the SNP. cell_id consists of the cell barcode as present in the bam-file, and A, C, G and T
denote how many unique UMIs are seen in the reads of the corresponding cell at the corresponding SNP having the corresponding bases. 

Contrary to the tool presented in the beforementioned paper, this implementation assumes a fixed error in the base calling, which can be set using the -e flag. 

All flags are explained when running the program without arguments.  
