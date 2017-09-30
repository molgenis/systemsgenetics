package deconvolution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/*
 *  Collection of all the interaction models and their shared data (genotypes, expression etc)
 *  There are n + 1 interaction models, where n = the number of celltypes. One full model 
 *  with all celltypes, and for each celltype one model with the interaction term for that
 *  model removed
 */
public class InteractionModelCollection {
	private double[] expressionValues;
	private double[] genotypes;
	private int sampleSize;
	private List<String> celltypes = new ArrayList<String>();
	private String qtlName;
	private HashMap<String, InteractionModel> interactionModels = new HashMap<String, InteractionModel>();
	private HashMap<String, Double> pvalues = new HashMap<String, Double>();

	/*
	 * Add interaction model to the collections
	 */
	public void addInteractionModel(InteractionModel interactionModel, String modelName){
		this.interactionModels.put(modelName, interactionModel);
	}

	/*
	 * Get interaction model at index 
	 */
	public InteractionModel getInteractionModel(String modelName) throws IllegalAccessException{
		InteractionModel interactionModel = this.interactionModels.get(modelName);
		if(interactionModel == null){
			throw new IllegalAccessException(String.format("No model with name %s found", modelName));
		}
		return(interactionModel);
	}


	/** 
	 * Set the expression values (y) for all the interaction models. 
	 */
	public void setExpressionValues(double[] expression){
		this.expressionValues = expression;
	}
	/** 
	 * Get the expression values (y) of all the interaction models. 
	 */
	public double[] getExpessionValues() throws IllegalAccessException{
		if(this.expressionValues == null){
			throw new IllegalAccessException("expressionValues not set");
		}
		return(this.expressionValues);
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

	/** 
	 * Get a list of all the celltypes given as input
	 */
	public void setCelltypes(List<String> celltypes) throws IllegalAccessException{
		this.celltypes = celltypes;
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


	public void emptyExpressionValues(){
		this.expressionValues = null;
	}

	public void setGenotypes(double[] genotypes) {
		this.genotypes = genotypes;
	}

	public void emptyGenotypes(){
		this.genotypes = null;
	}

	public int getSampleSize(){
		return this.sampleSize;
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
	
	/*
	 * Each ctModel will have a p-value from ANOVA test with fullmodel, save it per ctModel
	 */
	public void setPvalue(Double pvalue, String modelName){
	    this.pvalues.put(modelName, pvalue);
	  }
	
	public Double getPvalue(String modelName) throws IllegalAccessException{
		Double pvalue = this.pvalues.get(modelName);
		if(pvalue == null){
			throw new IllegalAccessException(String.format("Pvalue not set for model %s", modelName));
		}
		System.out.println(pvalue);
	    return(pvalue);
	  }
}
