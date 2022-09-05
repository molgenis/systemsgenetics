/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.toolbox;

/**
 *
 * @author patri
 */
public enum Mode {
	
	OVERLAP_GWAS(true, true);
	
	private final boolean needGenotypes;
	private final boolean needGwasCatalog;

	private Mode(boolean needGenotypes, boolean needGwasCatalog) {
		this.needGenotypes = needGenotypes;
		this.needGwasCatalog = needGwasCatalog;
	}

	public boolean isNeedGenotypes() {
		return needGenotypes;
	}

	public boolean isNeedGwasCatalog() {
		return needGwasCatalog;
	}
	
	
	
}
