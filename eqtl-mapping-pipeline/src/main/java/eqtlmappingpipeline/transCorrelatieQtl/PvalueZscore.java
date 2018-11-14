/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.transCorrelatieQtl;

/**
 *
 * @author patri
 */
public class PvalueZscore {
	
	private final double pvalue;
	private final double zscore;

	public PvalueZscore(double pvalue, double zscore) {
		this.pvalue = pvalue;
		this.zscore = zscore;
	}

	public double getPvalue() {
		return pvalue;
	}

	public double getZscore() {
		return zscore;
	}
	
}
