package deconvolution;

import java.util.ArrayList;
import java.util.List;

public class Qtl {
	private ArrayList<String> expressionVector;
	private ArrayList<String> genotypeVector;
	private List<List<String>> cellcountTable;
	private String qtlName;
	public Qtl(){};
	public Qtl( ArrayList<String> expressionVector, ArrayList<String> genotypeVector, List<List<String>> cellcountTable, String qtlName){
	    this.expressionVector = expressionVector;
	    this.genotypeVector = genotypeVector;
	    this.cellcountTable = cellcountTable;
	    this.qtlName = qtlName;
	}
	public ArrayList<String> getExpressionVector(){
		return(expressionVector);
	}
	public ArrayList<String> getGenotypeVector(){
		return(genotypeVector);
	}
	public List<List<String>> getCellcountTable(){
		return(cellcountTable);
	}
	public String getQtlName(){
		return(qtlName);
	}

}
