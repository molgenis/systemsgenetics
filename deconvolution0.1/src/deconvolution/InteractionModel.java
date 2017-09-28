package deconvolution;

import java.util.ArrayList;
import java.util.List;

public class InteractionModel {
	/* independentVariables = the names of the independent variables, e.g. neut%, mono%, neut%:GT */
	private List<String> independentVariableNames = new ArrayList<String>();
	private List<int[]> celltypeVariablesIndex = new ArrayList <int[]>();
	private double[][] observedValues;
	private double[] expressionValues;
	private double[] genotypes;
	private int sampleSize;
	private Boolean noIntercept;
	private String qtlName;
	private String modelName;
	private List<String> celltypes = new ArrayList<String>();
	private double[] estimatedRegressionParameters;
	private double[] estimateRegressionParametersStandardErrors;
	public InteractionModel(){};
	// initialize with number so that we can test if it has been set
	private double pvalue = -1;
	/**
	 * Set the observed values. Per QTL for each sample the observed values are each term of the 
	 * linear model. E.g. if the model is y = mono% + neut% + mono%:GT, the observedValues are
	 * [mono%, neut%, mono% * GT]
	 */
	public void InitializeObservedValue( int sampleSize, int numberOfTerms){

	    this.observedValues = new double[sampleSize][numberOfTerms];
	    this.sampleSize = sampleSize;
	  }
	
	public void addObservedValue( double observedValue, int sampleIndex, int termIndex){
	    this.observedValues[sampleIndex][termIndex] = observedValue;
	  }
	
	/**
	 * Get the observed values. Per QTL for each sample the observed values are each term of the 
	 * linear model. E.g. if the model is y = mono% + neut% + mono%:GT, the observedValues are
	 * [mono%, neut%, mono% * GT]
	 */
	public double[][] getObservedValues() throws IllegalAccessException{
		if(this.observedValues == null){
			throw new IllegalAccessException("observedValues not set for this model");
		}
	    return(this.observedValues);
	  }
	
	/**
	 * Add the index of the celltype variables of the linear model e.g.
	 * the index of the celltype% and celltype%:GT of the model. If 
	 *    y = neut% + mono% + eos% + neut%:GT + eos%:GT
	 *    celltypeTerms = [[0,3],[1],[2,4]
	 * This can be used to sum up the Beta * variable per cell type  
	 */
	public void addCelltypeVariablesIndex(int[] values){
		celltypeVariablesIndex.add(values);
	  }
	
	/**
	 * Get the index of the celltype variables of the linear model e.g.
	 * the index of the celltype% and celltype%:GT of the model. If 
	 *    y = neut% + mono% + eos% + neut%:GT + eos%:GT
	 *    celltypeTerms = [[0,3],[1],[2,4]
	 * This can be used to sum up the Beta * variable per cell type  
	 */
	public List<int[]> getCelltypeVariablesIndex() throws IllegalAccessException{
		if(this.celltypeVariablesIndex == null){
			throw new IllegalAccessException("celltypeVariables not set for this model");
		}
	    return (this.celltypeVariablesIndex);
	  }
	
	/** 
	 * Add the name of the independent variable name at the end of the existing list
	 * of independent variables.  
	 */
	public void addIndependentVariableName(String independentVariables){
	    this.independentVariableNames.add(independentVariables);
	  }
	
	/**
	 * Add the name of the independent variable name at <index> of the existing list
	 * of independent variables.  
	 */
	public void addIndependentVariableName(int index, String independentVariables){
	    this.independentVariableNames.add(index, independentVariables);
	  }
	
	/** 
	 * Get a list of the independent variables of the interaction model e.g.
	 * 		[neut%, mono%, neut%:GT]
	 */
	public List<String> getIndependentVariableNames() throws IllegalAccessException{

		if(this.independentVariableNames == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
	    return(this.independentVariableNames);
	  }
	
	/** 
	 * Set the expression values (y) of the model. 
	 */
	public void setExpressionValues(double[] expression){
	    this.expressionValues = expression;
	  }
	
	/** 
	 * Get the expression values (y) of the model. 
	 */
	public double[] getExpessionValues() throws IllegalAccessException{
		if(this.expressionValues == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
	    return(this.expressionValues);
	  }
	public double[] getGenotypes() throws IllegalAccessException {
		if(this.genotypes == null){
			throw new IllegalAccessException("genotypes not set for this model");
		}
		return this.genotypes;
	}
	
	/** 
	 * Set No Intercept for the model. If true, the intercept will be removed
	 */
	public void setNoIntercept(Boolean noIntercept){
	    this.noIntercept = noIntercept;
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
	
	/** 
	 * Get the No Intercept value for the model. If true, the intercept will be removed
	 */
	public Boolean getNoIntercept() throws IllegalAccessException{
		if(this.noIntercept == null){
			throw new IllegalAccessException("noIntercept not set for this model");
		}
	    return(this.noIntercept);
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
	
	public void setModelName(String modelName){
	    this.modelName = modelName;
	  }
	
	public String getModelName() throws IllegalAccessException{
		if(this.modelName == null){
			throw new IllegalAccessException("modelName name not set for this model");
		}
	    return(this.modelName);
	  }
		
	public void setEstimateRegressionParametersStandardErrors(double[] estimateRegressionParametersStandardErrors){
		this.estimateRegressionParametersStandardErrors = estimateRegressionParametersStandardErrors;
	}
	public double[] getEstimateRegressionParametersStandardErrors() throws IllegalAccessException{
		if(this.estimateRegressionParametersStandardErrors == null){
			throw new IllegalAccessException("estimateRegressionParametersStandardErrors not set for this model");
		}
		return(this.estimateRegressionParametersStandardErrors);
	}
	
	public void setEstimateRegressionParameters(double[] estimatedRegressionParameters){
		this.estimatedRegressionParameters = estimatedRegressionParameters;
	}
	public double[] getEstimateRegressionParameters() throws IllegalAccessException{
		if(this.estimatedRegressionParameters == null){
			throw new IllegalAccessException("estimatedRegressionParameters not set for this model");
		}
		return(this.estimatedRegressionParameters);
	}
	
	public void setAlltIndependentVariableNames(){
		for (int celltypeIndex = 0; celltypeIndex < this.celltypes.size(); celltypeIndex++) {
			int[] index = new int[] {celltypeIndex, this.celltypes.size() + celltypeIndex};
			addCelltypeVariablesIndex(index);
			addIndependentVariableName(celltypeIndex, this.celltypes.get(celltypeIndex));
			addIndependentVariableName(this.celltypes.get(celltypeIndex)+":GT");
		}
	}
	
	public void emptyObservedValues(){
		this.observedValues = null;
	}
	public void emptyExpressionValues(){
		this.expressionValues = null;
	}
	
	public void setPvalue(double pvalue){
		this.pvalue = pvalue;
	}
	public double getPvalue() throws IllegalAccessException{
		if(this.pvalue == -1){
			throw new IllegalAccessException("pvalue not set for this model");
		}
		return this.pvalue;
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
}
