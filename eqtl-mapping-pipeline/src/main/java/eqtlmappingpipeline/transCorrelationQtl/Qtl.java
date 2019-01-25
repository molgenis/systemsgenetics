/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.transCorrelationQtl;

/**
 *
 * @author patri
 */
public class Qtl {
	
	private final String snp;
	private final String trait;
	private final double Zscore;

	public Qtl(String snp, String trait, double Zscore) {
		this.snp = snp;
		this.trait = trait;
		this.Zscore = Zscore;
	}

	public String getSnp() {
		return snp;
	}

	public String getTrait() {
		return trait;
	}

	public double getZscore() {
		return Zscore;
	}
	
	
	
}
