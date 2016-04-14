package deconvolution;

import java.util.List;

public class CellCount {
	public List<String> celltypes;
	public List<String> samplenames;
	public List<Double> cellcountPercentages;
	public CellCount(){};
	  public CellCount( List<String> celltypes, List<String> samplenames, List<Double> cellcountPercentages){
		    this.celltypes = celltypes;
		    this.samplenames = samplenames;
		    this.cellcountPercentages = cellcountPercentages;
		  }
}
