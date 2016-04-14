package deconvolution;

import java.util.ArrayList;
import java.util.List;

public class Qtl {
	public ArrayList<String> expressionVector;
	public ArrayList<String> genotypeVector;
	public List<List<String>> cellcountTable;
	public Qtl(){};
	  public Qtl( ArrayList<String> expressionVector, ArrayList<String> genotypeVector, List<List<String>> cellcountTable){
		    this.expressionVector = expressionVector;
		    this.genotypeVector = genotypeVector;
		    this.cellcountTable = cellcountTable;
		  }
}
