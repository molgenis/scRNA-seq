import org.molgenis.genotype.Allele;

public class SNPCompare {

	private int[] d_genoCode = new int[4]; // in alphabetic order: A, C, G, T
	private int d_numbSamples;
	private double[] d_result;
	private double[] d_doublets;
	private float[][] d_snpGP;
	private Allele d_ref;
	private Allele d_alt;
	private double[] d_alphas; 
	private boolean d_goodSNP = true;
	private double d_error;
	
	public SNPCompare(int[] genoCode, float[][] snpGP, Allele ref, Allele alt, int numbSamples, double[] alphas, double error)
	{
		d_genoCode = genoCode;
		d_numbSamples = numbSamples;
		d_result = new double[d_numbSamples]; 
		d_alphas = alphas;
		d_doublets = new double[d_numbSamples * (d_numbSamples - 1) / 2 * d_alphas.length];
		d_snpGP = snpGP;
		d_ref = ref;
		d_alt = alt;
		d_error = error;
	}
	
	public double[] run()
	{		
		int refInt = alleleToInt(d_ref);
		int altInt = alleleToInt(d_alt);

		calcSinglet(refInt, altInt);
		
		return d_result;
	}
	
	
	public double[] runDoublets()
	{		
		int refInt = alleleToInt(d_ref);
		int altInt = alleleToInt(d_alt);

		calcDoublets(refInt, altInt);
		
		return d_doublets;
	}
		
	
	
	private void calcSinglet(int ref, int alt)
	{
		for (int samp = 0; samp < d_numbSamples; ++samp)
		{			
			int n_ref = d_genoCode[ref];
			int n_alt = d_genoCode[alt];
				
			if (n_ref + n_alt == 0)
				continue; 
			
			// The number represents the genotype of samp, the letter refers to the allele (r = ref, a = alt)
			// So, p2r means: P(ref allele | genotype samp = 1/1) (then the observed ref allele is supposed to be a mistake)
			double p0r = (1 - d_error) 						   * d_snpGP[samp][0];  
			double p1r = ((1 - d_error) * 0.5 + d_error / 6.0) * d_snpGP[samp][1] ;
			double p2r = d_error / 3.0 					       * d_snpGP[samp][2]; 

			double p0a = d_error / 3.0 					       * d_snpGP[samp][0];
			double p1a = ((1 - d_error) * 0.5 + d_error / 6.0) * d_snpGP[samp][1];
			double p2a = (1 - d_error) 						   * d_snpGP[samp][2];
			
			double tmp;

			tmp = Math.pow(p0a, n_alt) * Math.pow(p0r, n_ref) + Math.pow(p1a, n_alt) * Math.pow(p1r, n_ref) + Math.pow(p2a, n_alt) * Math.pow(p2r, n_ref);
	
			// Avoid taking the logarithm of too small values, as this will result in -INF. Instead, discard this SNP
			if (tmp < 1e-300)
				d_goodSNP = false; 
			
			d_result[samp] = Math.log(tmp);	
			
		}
		
	}
	

	private void calcDoublets(int ref, int alt)
	{
		int n_ref = d_genoCode[ref];
		int n_alt = d_genoCode[alt];
		
		if (n_ref + n_alt == 0)
			return;
		
		int combTested = 0;

		
		for (int samp1 = 0; samp1 < d_numbSamples; ++samp1)
		{
			for (int samp2 = samp1 + 1; samp2 < d_numbSamples; ++samp2)
			{	
				++combTested;

				for (int alphaInt = 0; alphaInt < d_alphas.length; ++alphaInt)
				{
					double alpha = d_alphas[alphaInt];

					
					// The first number represents the genotype of samp1, second number is genotype samp2, the letter refers to the observed allele (r = ref, a = alt)
					// So, p01r means: P(ref allele | genotype samp1 = 0/0 and genotype samp2 = 0/1)
					double p00r = (1 - d_error) 							                                      * d_snpGP[samp1][0] * d_snpGP[samp2][0];
					double p01r = ((1 - alpha) * (1 - d_error) + alpha * ((1 - d_error) * 0.5 + d_error / 6.0))   * d_snpGP[samp1][0] * d_snpGP[samp2][1];
					double p02r = ((1 - alpha) * (1 - d_error) + alpha * (d_error / 3.0)) 			  			  * d_snpGP[samp1][0] * d_snpGP[samp2][2];
					double p10r = ((1 - alpha) * ((1 - d_error) * 0.5 + d_error / 6.0) + alpha * (1 - d_error))   * d_snpGP[samp1][1] * d_snpGP[samp2][0];
					double p11r = (((1 - d_error) * 0.5 + d_error / 6.0)) 										  * d_snpGP[samp1][1] * d_snpGP[samp2][1];
					double p12r = ((1 - alpha) * ((1 - d_error) * 0.5 + d_error / 6.0) + alpha * (d_error / 3.0)) * d_snpGP[samp1][1] * d_snpGP[samp2][2];
					double p20r = ((1 - alpha) * (d_error / 3.0) + alpha * (1 - d_error)) 						  * d_snpGP[samp1][2] * d_snpGP[samp2][0];
					double p21r = ((1 - alpha) * (d_error / 3.0) + alpha * ((1 - d_error) * 0.5 + d_error / 6.0)) * d_snpGP[samp1][2] * d_snpGP[samp2][1];
					double p22r = (d_error / 3.0) 																  * d_snpGP[samp1][2] * d_snpGP[samp2][2];

					
					double p00a =  (d_error / 3.0) 																   * d_snpGP[samp1][0] * d_snpGP[samp2][0];
					double p01a =  ((1 - alpha) * (d_error / 3.0) + alpha * ((1 - d_error) * 0.5 + d_error / 6.0)) * d_snpGP[samp1][0] * d_snpGP[samp2][1];
					double p02a =  ((1 - alpha) * (d_error / 3.0) + alpha * (1 - d_error)) 					   	   * d_snpGP[samp1][0] * d_snpGP[samp2][2];
					double p10a =  ((1 - alpha) * ((1 - d_error) * 0.5 + d_error / 6.0) + alpha * (d_error / 3.0)) * d_snpGP[samp1][1] * d_snpGP[samp2][0];
					double p11a =  (((1 - d_error) * 0.5 + d_error / 6.0)) 									   	   * d_snpGP[samp1][1] * d_snpGP[samp2][1];
					double p12a =  ((1 - alpha) * ((1 - d_error) * 0.5 + d_error / 6.0) + alpha * (1 - d_error))   * d_snpGP[samp1][1] * d_snpGP[samp2][2];
					double p20a =  ((1 - alpha) * (1 - d_error) + alpha * (d_error / 3.0)) 					  	   * d_snpGP[samp1][2] * d_snpGP[samp2][0];
					double p21a =  ((1 - alpha) * (1 - d_error) + alpha * ((1 - d_error) * 0.5 + d_error / 6.0))   * d_snpGP[samp1][2] * d_snpGP[samp2][1];
					double p22a =  ((1 - d_error)) 																   * d_snpGP[samp1][2] * d_snpGP[samp2][2];
					
					double tmp = Math.pow(p00r, n_ref) * Math.pow(p00a, n_alt) +
							Math.pow(p01r, n_ref) * Math.pow(p01a, n_alt) +
							Math.pow(p02r, n_ref) * Math.pow(p02a, n_alt) +
							Math.pow(p10r, n_ref) * Math.pow(p10a, n_alt) +
							Math.pow(p11r, n_ref) * Math.pow(p11a, n_alt) +
							Math.pow(p12r, n_ref) * Math.pow(p12a, n_alt) +
							Math.pow(p20r, n_ref) * Math.pow(p20a, n_alt) +
							Math.pow(p21r, n_ref) * Math.pow(p21a, n_alt) +
							Math.pow(p22r, n_ref) * Math.pow(p22a, n_alt);

					
					if (tmp < 1e-300)
						d_goodSNP = false; 
					
					d_doublets[(combTested - 1) * d_alphas.length + alphaInt] = Math.log(tmp);
 				}
				
			}
		}
		
	}
	
	public boolean goodSNP()
	{
		return d_goodSNP;
	}
	
	private int alleleToInt(Allele al)
	{
		int alleleAsInt = 0;
		
		Allele a = Allele.create('A');
		if (al.equals(a))
			alleleAsInt = 0;
		Allele c = Allele.create('C');
		if (al.equals(c))
			alleleAsInt = 1;
		Allele g = Allele.create('G');
		if (al.equals(g))
			alleleAsInt = 2;
		Allele t = Allele.create('T');
		if (al.equals(t))
			alleleAsInt = 3;
		
		return alleleAsInt;
	}
	

}


