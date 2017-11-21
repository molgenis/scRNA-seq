

import java.io.IOException;
import org.apache.commons.cli.*;

/**
 * Main class for DeAnonymizer
 *
 * @author Marion Dam
 * @version 1.1
 * 
 * The likelihood (as computed in SNPcompare class) is described in:
 * Multiplexing droplet-based single cell RNA-sequencing using natural genetic barcodes, H.M. Kang e.a., biorxiv preprint
 */



public class Main 
{
	public static void main(String[] args) throws IOException 
	{
	    Options options = new Options();
	    options.addOption("g", "genotypeVCF", true, "Genotype VCF file path.");
	    options.addOption("s", "SNPfile", true, "SNP-reads file from countUMI-script, containing the reads per SNP per cell.");
	    options.addOption("n", "cellnames", true, "Path and file containing names of all cells to be considered.");
	    options.addOption("o", "output", true, "Output path and file name.");
	    options.addOption("a", "alpha", true, "OPTIONAL. Integer from 1 to 4, where: \t1: {0.5}, \t2: {0.25, 0.5, 0.75}, \t3: {0.10, 0.20, ..., 0.90}. \t4: {0.05, 0.10, ..., 0.95}, \tDefault = 1.");
	    options.addOption("t", "tValCutoff", true, "OPTIONAL. Set cut-off for log-likelihood difference (integer). \tDefault = 1.");
	    options.addOption("e", "error", true, "OPTIONAL. Set error in base quality. \tDefault = 0.001.");
	    options.addOption("S", "samples", true, "OPTIONAL. Text file to specify which samples in the VCF to be included, one sample name per line. \tDefault = ALL samples in VCF.");

	    
	    CommandLineParser parser = new DefaultParser();
	        
	    try 
	    {
            CommandLine line = parser.parse(options, args);
            if(!line.hasOption("g") || !line.hasOption("s") || !line.hasOption("o")) 
            {
                HelpFormatter formatter = new HelpFormatter();
                String header = "DeAnonymizer. \nThe program will compute the most likely sample (or combination of two samples) to which each input cell belongs."
                		+ " Provide as input: the genotype VCF of the samples to be considered, a list of cell-ids, a file with (SNP position, cell-id, genotype), "
                		+ "as generated with count_umi-script.\n\n";
                formatter.printHelp(header, options, true);
                System.exit(1);
            }
            
            String genoVCF = line.getOptionValue("g");
            String SNPfile = line.getOptionValue("s");
            String cellNameFile = line.getOptionValue("n");
            String outputPath = line.getOptionValue("o");
            
            
            String sampleFile = null;
            if (line.hasOption("S"))
            	sampleFile = line.getOptionValue("S");
            
            int alpha = 1;
            if(line.hasOption("a"))
            	alpha = Integer.parseInt(line.getOptionValue("a"));
            
            int tVal = 1;
            if(line.hasOption("t"))
            	tVal = Integer.parseInt(line.getOptionValue("t"));
            
            double error = 0.001;
            if(line.hasOption("e"))
            	error = Double.parseDouble(line.getOptionValue("e"));

            System.out.println("Running DeAnonymizer with the following settings:\nError = "  
            	  + error + "\nLikelihood difference cut-off = " + tVal +"\n\nStarted reading files...\n");
            
            GenoTable genoTable;

	        genoTable = new GenoTable(genoVCF, SNPfile, cellNameFile, alpha, tVal, error, sampleFile);
	        genoTable.run(); 

            genoTable.writeTable(outputPath);

            System.out.println("Finished the program!");  

	    }
	    catch (ParseException exp) 
	    {
            System.out.println("Parse exception:" + exp.getMessage());
        } 
	}
}
	    


	





