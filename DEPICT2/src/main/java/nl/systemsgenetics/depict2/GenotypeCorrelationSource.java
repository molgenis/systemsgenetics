/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public interface GenotypeCorrelationSource {
	
	/**
	 * 
	 * 
	 * @param chr
	 * @param start
	 * @param stop
	 * @return 
	 */
	public DoubleMatrixDataset<String, String> getCorrelationMatrixForRange(String chr, int start, int stop);
	
}
