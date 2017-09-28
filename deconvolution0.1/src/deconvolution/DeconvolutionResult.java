package deconvolution;

import java.util.ArrayList;
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
	private List<Double> correctedPvalues;
	private Map<String,Double> correctedPvaluePerCelltype = new HashMap<String, Double>();
	private double wholeBloodQTL;
	private double wholeBloodQTLpvalue;
	
	public DeconvolutionResult(){};
	
	/**
	 * Set the deconvolutionResult with the InteractionModels.
	 * 
	 * @param celltypes List of celltypes in the deconvolution result
	 * 
	 * @param qtlName The name of the QTL
	 * 
	 * @param pvalues The pvalues from the deconvolution model
	 * 
	 * @param fullModel Interaction model containing information of all celltypes
	 * 
	 * @param ctModels List of interaction models (each without all celltypes)
	 * 
	 * @param multipleTestingMethod Method to use for multiple testing correction
	 * 
	 * @param wholeBloodQTL Spearman correlation of genotypes and expression levels
	 */
	public DeconvolutionResult( List<String> celltypes, String qtlName, List<Double> pvalues, InteractionModel fullModel, 
								List<InteractionModel> ctModels, String multipleTestingMethod, double wholeBloodQTL, 
								double wholeBloodQTLpvalue){

		this.celltypes = celltypes;
		this.qtlName = qtlName;
		this.pvalues = pvalues;
		for (int i = 0; i < celltypes.size(); i++){
			pvaluePerCelltype.put(celltypes.get(i), pvalues.get(i));
		}
		this.fullModel = fullModel;
		this.ctModels = ctModels;
		correctPvaluesForMultipeTesting(pvalues, multipleTestingMethod);
		this.wholeBloodQTL = wholeBloodQTL;
		this.wholeBloodQTLpvalue = wholeBloodQTLpvalue;
	}
	
	/**
	 * Set the deconvolutionResult with the InteractionModels.
	 * 
	 * @param celltypes List of celltypes in the deconvolution result
	 * 
	 * @param qtlName The name of the QTL
	 * 
	 * @param pvalues The pvalues from the deconvolution model
	 * 
	 * @param fullModel Interaction model containing information of all celltypes
	 */
	public DeconvolutionResult( List<String> celltypes, String qtlName, List<Double> pvalues,InteractionModel fullModel,
								double wholeBloodQTL,double wholeBloodQTLpvalue){
		this.celltypes = celltypes;
		this.qtlName = qtlName;
		this.pvalues = pvalues;
		this.correctedPvalues = pvalues;
		for (int i = 0; i < celltypes.size(); i++){
			pvaluePerCelltype.put(celltypes.get(i), pvalues.get(i));
		}
		this.fullModel = fullModel;
		for (int i = 0; i < celltypes.size(); i++){
			correctedPvaluePerCelltype.put(celltypes.get(i), correctedPvalues.get(i));
		}
		this.wholeBloodQTL = wholeBloodQTL;
		this.wholeBloodQTLpvalue = wholeBloodQTLpvalue;
	}
	
	/**
	 * Set the deconvolutionResult without the InteractionModels. This is for when the models are not used, e.g. when option -m is used
	 * and not all genotypes have enough samples
	 * 
	 * @param celltypes List of celltypes in the deconvolution result
	 * 
	 * @param qtlName The name of the QTL
	 * 
	 * @param pvalues The pvalues from the deconvolution model
	 * 
	 * @param wholeBloodQTL Spearman correlation of genotypes and expression levels
	 */
	public DeconvolutionResult( List<String> celltypes, String qtlName, List<Double> pvalues, double wholeBloodQTL,double wholeBloodQTLpvalue){
		this.celltypes = celltypes;
		this.qtlName = qtlName;
		this.pvalues = pvalues;
		correctedPvalues = pvalues;
		for (int i = 0; i < celltypes.size(); i++){
			pvaluePerCelltype.put(celltypes.get(i), pvalues.get(i));
		}
		for (int i = 0; i < celltypes.size(); i++){
			correctedPvaluePerCelltype.put(celltypes.get(i), correctedPvalues.get(i));
		}
		this.wholeBloodQTL = wholeBloodQTL;
		this.wholeBloodQTLpvalue = wholeBloodQTLpvalue;
	}
	

	private void correctPvaluesForMultipeTesting(List<Double> pvalues, String method){
		correctedPvalues = new ArrayList<Double>();
		if (method == "bonferonni"){
			for(double pvalue : pvalues){
				double correctedPvalue = pvalue*pvalues.size();
				if(pvalue*pvalues.size() > 1){
					correctedPvalue = 1;
				}
				correctedPvalues.add(correctedPvalue);
			}
			for (int i = 0; i < celltypes.size(); i++){
				correctedPvaluePerCelltype.put(celltypes.get(i), correctedPvalues.get(i));
			}
		}
		else{
			throw new RuntimeException("Only bonferonni has been implemented as multiple testing correction method");
		}
	}
	
	public void setQtlName(String qtlName){
		this.qtlName = qtlName;
	}

	public String getQtlName() throws IllegalAccessException{
		if(this.qtlName == null){
			throw new IllegalAccessException("QTL name not set for this model");
		}
		return(this.qtlName);
	}
	
	/** 
	 * Get a list of all the celltypes given as input
	 */
	public List<String> getCelltypes() throws IllegalAccessException{
		if(this.celltypes == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
		return(this.celltypes);
	}
	
	public void setPvalues(List<Double> pvalues){
		this.pvalues = pvalues;
	}

	public List<Double> getPvalues() throws IllegalAccessException{
		if(this.pvalues == null){
			throw new IllegalAccessException("pvalues not set for this model");
		}
		return(this.pvalues);
	}
	
	public InteractionModel getFullModel() throws IllegalAccessException{
		if(this.fullModel == null){
			throw new IllegalAccessException("fullModel not set for this model");
		}
		return(this.fullModel);
	}
	
	public List<InteractionModel> getCtModels() throws IllegalAccessException{
		if(this.ctModels == null){
			throw new IllegalAccessException("ctModels not set for this model");
		}
		return(this.ctModels);
	}
	
	public Map<String, Double>  getPvaluePerCelltype() throws IllegalAccessException{
		if(this.pvaluePerCelltype == null){
			throw new IllegalAccessException("pvaluePerCelltype not set for this model");
		}
		return(this.pvaluePerCelltype);
	}
	public Map<String, Double>  getCorrectedPvaluePerCelltype() throws IllegalAccessException{
		if(this.correctedPvaluePerCelltype == null){
			throw new IllegalAccessException("correctedPvaluePerCelltype not set for this model");
		}
		return(this.correctedPvaluePerCelltype);
	}
	public List<Double>  getCorrectedPvalues() throws IllegalAccessException{
		if(this.correctedPvalues == null){
			throw new IllegalAccessException("correctedPvalues not set for this model");
		}
		return(this.correctedPvalues);
	}
	public double  getWholeBloodQTL() throws IllegalAccessException{
		return(this.wholeBloodQTL);
	}
	public double  getWholeBloodQTLpvalue() throws IllegalAccessException{
		return(this.wholeBloodQTLpvalue);
	}
}
