package deconvolution;


public class Qtl {
	private double[] expressionVector;
	private double[] genotypeVector;
	private CellCount cellCounts;
	private String qtlName;
	public Qtl(){};
	public Qtl( double[] expressionVector, double[] genotypeVector, CellCount cellCounts, String qtlName){
	    this.expressionVector = expressionVector;
	    this.genotypeVector = genotypeVector;
	    this.cellCounts = cellCounts;
	    this.qtlName = qtlName;
	}
	public double[] getExpressionVector(){
		return(expressionVector);
	}
	public double[] getGenotypeVector(){
		return(genotypeVector);
	}
	public CellCount getCellCounts(){
		return(cellCounts);
	}
	public String getQtlName(){
		return(qtlName);
	}

}
