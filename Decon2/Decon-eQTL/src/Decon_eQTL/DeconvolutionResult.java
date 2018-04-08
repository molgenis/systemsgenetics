package Decon_eQTL;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DeconvolutionResult {
	private List<String> celltypes;
	private String qtlName;
	private List<Double> pvalues = new ArrayList<Double>();
	private Map<String, Double> pvaluePerCelltype = new HashMap<String, Double>();
	private InteractionModelCollection interactionModelCollection;
	private double wholeBloodQTL;
	private double wholeBloodQTLpvalue;
	
	public DeconvolutionResult(){};
	
	/**
	 * Set the deconvolutionResult with the InteractionModels.
	 * 
	 * @param interactionModelCollection Collection of the interaction models used to get the deconvolution result
	 * 
	 * @param wholeBloodQTL Spearman correlation of genotypes and expression levels
	 * 
	 * @param wholeBloodQTLpvalue pvalue of the pearman correlation of genotypes and expression levels
	 * @throws IllegalAccessException 
	 */
	public DeconvolutionResult( InteractionModelCollection interactionModelCollection, double wholeBloodQTL, 
								double wholeBloodQTLpvalue) throws IllegalAccessException{

		celltypes = interactionModelCollection.getAllCelltypes();
		this.qtlName = interactionModelCollection.getQtlName();
		for (int i = 0; i < celltypes.size(); i++){
			String modelName = celltypes.get(i);
			Double pvalue = interactionModelCollection.getPvalue(modelName);
			this.pvalues.add(pvalue);
			pvaluePerCelltype.put(modelName, pvalue);
		}
		this.interactionModelCollection = interactionModelCollection;
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
		for (int i = 0; i < celltypes.size(); i++){
			pvaluePerCelltype.put(celltypes.get(i), pvalues.get(i));
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
		for (int i = 0; i < celltypes.size(); i++){
			pvaluePerCelltype.put(celltypes.get(i), pvalues.get(i));
		}

		this.wholeBloodQTL = wholeBloodQTL;
		this.wholeBloodQTLpvalue = wholeBloodQTLpvalue;
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
	
	public InteractionModel getModel(String modelName) throws IllegalAccessException{
		InteractionModel interationModel = this.interactionModelCollection.getInteractionModel(modelName);
		if(interationModel == null){
			throw new IllegalAccessException(String.format("model not set for %s", modelName));
		}
		return(interationModel);
	}
	
	public InteractionModelCollection getInteractionModelCollection() throws IllegalAccessException{
		if(this.interactionModelCollection == null){
			throw new IllegalAccessException("interactionModelCollection not set");
		}
		return(this.interactionModelCollection);
	}
		
		
	public Map<String, Double>  getPvaluePerCelltype() throws IllegalAccessException{
		if(this.pvaluePerCelltype == null){
			throw new IllegalAccessException("pvaluePerCelltype not set for this model");
		}
		return(this.pvaluePerCelltype);
	}

	public double  getWholeBloodQTL() throws IllegalAccessException{
		return(this.wholeBloodQTL);
	}
	public double  getWholeBloodQTLpvalue() throws IllegalAccessException{
		return(this.wholeBloodQTLpvalue);
	}
}
