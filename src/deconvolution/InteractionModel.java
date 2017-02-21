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
	public void InitializeObservedValue( int sampleSize, int numberOfTerms){
		/* 
		 * Set the observed values. Per QTL for each sample the observed values are each term of the 
		 * linear model. E.g. if the model is y = mono% + neut% + mono%:GT, the observedValues are
		 * [mono%, neut%, mono% * GT]
		 */
	    this.observedValues = new double[sampleSize][numberOfTerms];
	    this.sampleSize = sampleSize;
	  }
	
	public void addObservedValue( double observedValue, int sampleIndex, int termIndex){
	    this.observedValues[sampleIndex][termIndex] = observedValue;
	  }
	
	public double[][] getObservedValues() throws IllegalAccessException{
		/* 
		 * Get the observed values. Per QTL for each sample the observed values are each term of the 
		 * linear model. E.g. if the model is y = mono% + neut% + mono%:GT, the observedValues are
		 * [mono%, neut%, mono% * GT]
		 */
		if(this.observedValues == null){
			throw new IllegalAccessException("observedValues not set for this model");
		}
	    return(this.observedValues);
	  }
	
	public void addCelltypeVariablesIndex(int[] values){
		/* 
		 * Add the index of the celltype variables of the linear model e.g.
		 * the index of the celltype% and celltype%:GT of the model. If 
		 *    y = neut% + mono% + eos% + neut%:GT + eos%:GT
		 *    celltypeTerms = [[0,3],[1],[2,4]
		 * This can be used to sum up the Beta * variable per cell type  
		 */
		celltypeVariablesIndex.add(values);
	  }
	
	public List<int[]> getCelltypeVariablesIndex() throws IllegalAccessException{
		/* 
		 * Get the index of the celltype variables of the linear model e.g.
		 * the index of the celltype% and celltype%:GT of the model. If 
		 *    y = neut% + mono% + eos% + neut%:GT + eos%:GT
		 *    celltypeTerms = [[0,3],[1],[2,4]
		 * This can be used to sum up the Beta * variable per cell type  
		 */
		if(this.celltypeVariablesIndex == null){
			throw new IllegalAccessException("celltypeVariables not set for this model");
		}
	    return (this.celltypeVariablesIndex);
	  }
	
	public void addIndependentVariableName(String independentVariables){
		/* 
		 * Add the name of the independent variable name at the end of the existing list
		 * of independent variables.  
		 */
	    this.independentVariableNames.add(independentVariables);
	  }
	
	public void addIndependentVariableName(int index, String independentVariables){
		/* 
		 * Add the name of the independent variable name at <index> of the existing list
		 * of independent variables.  
		 */
	    this.independentVariableNames.add(index, independentVariables);
	  }
	
	public List<String> getIndependentVariableNames() throws IllegalAccessException{
		/* 
		 * Get a list of the independent variables of the interaction model e.g.
		 * 		[neut%, mono%, neut%:GT]
		 */
		if(this.independentVariableNames == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
	    return(this.independentVariableNames);
	  }
	
	public void setExpressionValues(double[] expression){
		/* 
		 * Set the expression values (y) of the model. 
		 */
	    this.expressionValues = expression;
	  }
	
	public double[] getExpessionValues() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
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
	public void setNoIntercept(Boolean noIntercept){
		/* 
		 * Set No Intercept for the model. If true, the intercept will be removed
		 */
	    this.noIntercept = noIntercept;
	  }
	
	public void setCelltypes(List<String> celltypes) throws IllegalAccessException{
		/* 
		 * Get a list of all the celltypes given as input
		 */
		this.celltypes = celltypes;
	}
	
	public List<String> getCelltypes() throws IllegalAccessException{
		/* 
		 * Get a list of all the celltypes given as input
		 */
		if(this.celltypes == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
		return(this.celltypes);
	}
	
	public Boolean getNoIntercept() throws IllegalAccessException{
		/* 
		 * Get the No Intercept value for the model. If true, the intercept will be removed
		 */
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
