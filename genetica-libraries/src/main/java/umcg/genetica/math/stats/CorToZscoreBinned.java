/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

/**
 *
 * @author patri
 */
public class CorToZscoreBinned {
	
	private final double[] zscoreLookupTable;
	private final int numberOfBins;
	private final double halfStep;

	/**
	 * 
	 * 
	 * @param numberOfBins determines precision
	 */
	public CorToZscoreBinned(int numberOfBins) throws Exception {
		
		double power = Math.log10(numberOfBins);
		
		if(power != Math.round(power)){
			throw new Exception("Number of bins must be power of 10");
		}
		
		
		this.numberOfBins = numberOfBins;
		zscoreLookupTable = new double[numberOfBins];
				
		
		final double stepSize = 1/numberOfBins;
		halfStep = stepSize / 2;
		
		for(int i = 0 ; i < numberOfBins ; ++i){
			
			double corBinCenter = stepSize * i + halfStep;
			
		}
		
		
		
	}
	
	public double lookupZscoreForR(double r){
		
		long bin = Math.round(r * numberOfBins - halfStep);
		
		//this is needed because due to rounding a r of 1 will not fall into any bin. Also a r > 1 is accepted for imprecise r calculations
		if(bin > numberOfBins){
			bin = numberOfBins ;
		}
		
		return zscoreLookupTable[(int) bin];
		
	}
	
	
	
	
}
