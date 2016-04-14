package deconvolution;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DeconvolutionResult {
	private List<String> celltypes;
	private String qtlName;
	private List<Double> pvalues;
	private DeconvolutionResult(){};
	Map<String, Double> pvaluePerCelltype = new HashMap<String, Double>();

	public DeconvolutionResult( List<String> celltypes, String qtlName, List<Double> pvalues){
		this.celltypes = celltypes;
		this.qtlName = qtlName;
		this.pvalues = pvalues;
		for (int i = 0; i < celltypes.size(); i++){
			pvaluePerCelltype.put(celltypes.get(i), pvalues.get(i));
		}
	}
	
	public void SetQtlName(String qtlName){
		this.qtlName = qtlName;
	}

	public String GetQtlName(){
		return(this.qtlName);
	}
	
	public List<String> GetCelltypes(){
		return(this.celltypes);
	}
	
	public void SetPvalues(List<Double> pvalues){
		this.pvalues = pvalues;
	}

	public List<Double> GetPvalues(){
		return(this.pvalues);
	}

}
