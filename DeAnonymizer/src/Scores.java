public class Scores {
	
	private double[] d_scores;
	private int d_snps;
	
	public Scores(int numbSamples)
	{
		d_scores = new double[numbSamples];
		d_snps = 0; 
	}

	public int getNSnps()
	{
		return d_snps;
	}
	
	public double[] getScores()
	{
		return d_scores;
	}
	
	public void addScores(double[] toAdd)
	{
		if (toAdd.length == d_scores.length)
		{
			for (int idx = 0; idx < d_scores.length; ++idx)
				d_scores[idx] += toAdd[idx];
		}
	}

	public void increaseNSnps()
	{
		++d_snps;
	}

	
}
