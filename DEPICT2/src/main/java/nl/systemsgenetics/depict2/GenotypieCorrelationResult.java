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
public class GenotypieCorrelationResult {
	
	private final double[][] corMatrix;
	private final String[] includedVariants;

	public GenotypieCorrelationResult(double[][] corMatrix, String[] includedVariants) {
		this.corMatrix = corMatrix;
		this.includedVariants = includedVariants;
	}

	public double[][] getCorMatrix() {
		return corMatrix;
	}

	public String[] getIncludedVariants() {
		return includedVariants;
	}	
	
	
}
