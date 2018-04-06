package Decon_eQTL;


public class Qtl {
	private double[] expression;
	private double[] genotypes;
	private CellCount cellCounts;
	private String qtlName;
	private double[] swappedGenotypes;
	public Qtl(){};
	public Qtl( double[] expression, double[] genotypes, CellCount cellCounts, String qtlName){
	    this.expression = expression;
	    this.genotypes = genotypes;
	    this.cellCounts = cellCounts;
	    this.qtlName = qtlName;
	}

	public void setGenotypes(double[] genotypes) {
		this.genotypes = genotypes;
		swapGenotypes();
	}
	/** 
	 * Get a list of all the celltypes given as input
	 */
	private void swapGenotypes(){
		this.swappedGenotypes = this.genotypes.clone();
		for(int i = 0; i < this.genotypes.length; i++) {
			if(this.genotypes[i] == 0){
				this.swappedGenotypes[i] = 2;
			}
			else if(this.genotypes[i] == 2){
				this.swappedGenotypes[i] = 0;
			}
		}
	}
	public void emptyGenotypes(){
		this.genotypes = null;
	}

	
	public double[] getExpressionVector(){
		return(expression);
	}
	public double[] getGenotypeVector(){
		return(genotypes);
	}
	public CellCount getCellCounts(){
		return(cellCounts);
	}
	public String getQtlName(){
		return(qtlName);
	}
	/*
	 * Get the genotypes of all the interaction models
	 */
	public double[] getGenotypes() throws IllegalAccessException {
		if(this.genotypes == null){
			throw new IllegalAccessException("genotypes not set for this model");
		}
		return this.genotypes;
	}

	/*
	 * Get the genotypes of all the interaction models
	 */
	public double[] getSwappedGenotypes() throws IllegalAccessException {
		if(this.swappedGenotypes == null){
			throw new IllegalAccessException("genotypes not set for this model");
		}
		return this.swappedGenotypes;	
	}
}
