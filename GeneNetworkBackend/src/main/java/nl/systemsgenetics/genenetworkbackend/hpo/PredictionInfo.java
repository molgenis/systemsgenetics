/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

/**
 *
 * @author patri
 */
public class PredictionInfo {
	
	final String hpo;
	final double pValue;
	final double auc;
	final double fdr;

	public PredictionInfo(String hpo, double pValue, double auc, double fdr) {
		this.hpo = hpo;
		this.pValue = pValue;
		this.auc = auc;
		this.fdr = fdr;
	}

	public String getHpo() {
		return hpo;
	}

	public double getpValue() {
		return pValue;
	}

	public double getAuc() {
		return auc;
	}

	public double getFdr() {
		return fdr;
	}
	
}
