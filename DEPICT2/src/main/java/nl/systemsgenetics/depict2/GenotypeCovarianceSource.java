/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

/**
 *
 * @author patri
 */
public interface GenotypeCovarianceSource {
	
	/**
	 * 
	 * 
	 * @param chr
	 * @param start
	 * @param stop
	 * @param doubleMaxR max correlation between variants in resulting covariance matrix. Highly correlated variants will be excluded.
	 * @return 
	 */
	public GenotypieCorrelationResult getCorrelationMatrixForRange(String chr, int start, int stop, double doubleMaxR);
	
}
