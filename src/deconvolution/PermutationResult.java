package deconvolution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class PermutationResult {
	private Map<String, List<Double>> pvaluePerCelltype = new HashMap<String, List<Double>>();
	private Set<String> celltypes = new HashSet<String>();
	public PermutationResult(){};
	public void add( String celltype, Double pvalues){
		celltypes.add(celltype);
		if (pvaluePerCelltype.get(celltype) == null) { //gets the value for an id)
			pvaluePerCelltype.put(celltype, new ArrayList<Double>()); //no ArrayList assigned, create new ArrayList
		}
		pvaluePerCelltype.get(celltype).add(pvalues); //adds value to list.
   }
	public void add( DeconvolutionResult deconvolutionResult){
		for(String celltype : deconvolutionResult.GetCelltypes()){
			celltypes.add(celltype);
			add(celltype, deconvolutionResult.pvaluePerCelltype.get(celltype));
		}
   }
	
	public Set<String> getCelltypes(){
		return(celltypes);
   }

	public List<Double> getPvalues(String celltype){
		return(pvaluePerCelltype.get(celltype));
   }
}
