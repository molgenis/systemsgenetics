package deconvolution;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DeconvolutionResult {
	private List<String> celltypes;
	private String qtlName;
	private List<Double> pvalues;
	private Map<String, Double> pvaluePerCelltype = new HashMap<String, Double>();
	private InteractionModel fullModel;
	private List<InteractionModel> ctModels;
	public DeconvolutionResult(){};
	public DeconvolutionResult( List<String> celltypes, String qtlName, List<Double> pvalues, InteractionModel fullModel, List<InteractionModel> ctModels){
		/*
		 * Set the deconvolutionResult with the InteractionModels.
		 * 
		 * @param celltypes List of celltypes in the deconvolution result
		 * 
		 * @param qtlName The name of the QTL
		 * 
		 * @param pvalues The pvalues from the deconvolution model
		 */
		this.celltypes = celltypes;
		this.qtlName = qtlName;
		this.pvalues = pvalues;
		for (int i = 0; i < celltypes.size(); i++){
			pvaluePerCelltype.put(celltypes.get(i), pvalues.get(i));
		}
		this.fullModel = fullModel;
		this.ctModels = ctModels;
	}
	
	public DeconvolutionResult( List<String> celltypes, String qtlName, List<Double> pvalues){
		/*
		 * Set the deconvolutionResult without the InteractionModels. This is for when the models are not used, e.g. when option -m is used
		 * and not all genotypes have enough samples
		 * 
		 * @param celltypes List of celltypes in the deconvolution result
		 * 
		 * @param qtlName The name of the QTL
		 * 
		 * @param pvalues The pvalues from the deconvolution model
		 */
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

	public String GetQtlName() throws IllegalAccessException{
		if(this.qtlName == null){
			throw new IllegalAccessException("QTL name not set for this model");
		}
		return(this.qtlName);
	}
	
	public List<String> GetCelltypes() throws IllegalAccessException{
		/* 
		 * Get a list of all the celltypes given as input
		 */
		if(this.celltypes == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
		return(this.celltypes);
	}
	
	public void SetPvalues(List<Double> pvalues){
		this.pvalues = pvalues;
	}

	public List<Double> GetPvalues() throws IllegalAccessException{
		if(this.pvalues == null){
			throw new IllegalAccessException("pvalues not set for this model");
		}
		return(this.pvalues);
	}
	
	public InteractionModel GetFullModel() throws IllegalAccessException{
		if(this.fullModel == null){
			throw new IllegalAccessException("fullModel not set for this model");
		}
		return(this.fullModel);
	}
	
	public List<InteractionModel> GetCtModels() throws IllegalAccessException{
		if(this.ctModels == null){
			throw new IllegalAccessException("ctModels not set for this model");
		}
		return(this.ctModels);
	}
	
}
