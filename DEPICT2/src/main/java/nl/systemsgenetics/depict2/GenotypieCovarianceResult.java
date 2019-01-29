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
	private final TIntHashSet includedVariantsPositions;

	public GenotypieCovarianceResult(double[][] cov, TIntHashSet includedVariantsPositions) {
		this.cov = cov;
		this.includedVariantsPositions = includedVariantsPositions;
	}

	public double[][] getCov() {
		return cov;
	}

	public TIntHashSet getIncludedVariantsPositions() {
		return includedVariantsPositions;
	}
	
	
	
}
