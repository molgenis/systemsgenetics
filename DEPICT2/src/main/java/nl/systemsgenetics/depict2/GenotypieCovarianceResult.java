/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import gnu.trove.set.hash.TIntHashSet;

/**
 *
 * @author patri
 */
public class GenotypieCovarianceResult {
	
	private final double[][] cov;
	private final String[] includedVariants;

	public GenotypieCovarianceResult(double[][] cov, String[] includedVariants) {
		this.cov = cov;
		this.includedVariants = includedVariants;
	}

	public double[][] getCov() {
		return cov;
	}

	public String[] getIncludedVariants() {
		return includedVariants;
	}	
	
	
}
