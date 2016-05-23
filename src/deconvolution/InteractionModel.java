package deconvolution;

import java.util.ArrayList;
import java.util.List;

public class InteractionModel {
	/* independentVariables = the names of the independent variables, e.g. neut%, mono%, neut%:GT */
	private List<String> independentVariableNames = new ArrayList<String>();
	private List<int[]> celltypeVariablesIndex = new ArrayList <int[]>();
	private double[][] observedValues;
	private double[] expressionValues;
	private Boolean noIntercept;
	private String qtlName;
	private String modelName;
	private List<String> celltypes = new ArrayList<String>();
	private double[] estimatedRegressionParameters;
	private double[] estimateRegressionParametersStandardErrors;
	public InteractionModel(){};

	public void SetObservedValues( double[][] observedValues){
		/* 
		 * Set the observed values. Per QTL for each sample the observed values are each term of the 
		 * linear model. E.g. if the model is y = mono% + neut% + mono%:GT, the observedValues are
		 * [mono%, neut%, mono% * GT]  
		 */
	    this.observedValues = observedValues;
	  }
	
	public double[][] GetObservedValues() throws IllegalAccessException{
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
	
	public void AddCelltypeVariablesIndex(int[] values){
		/* 
		 * Add the index of the celltype variables of the linear model e.g.
		 * the index of the celltype% and celltype%:GT of the model. If 
		 *    y = neut% + mono% + eos% + neut%:GT + eos%:GT
		 *    celltypeTerms = [[0,3],[1],[2,4]
		 * This can be used to sum up the Beta * variable per cell type  
		 */
		celltypeVariablesIndex.add(values);
	  }
	
	public List<int[]> GetCelltypeVariablesIndex() throws IllegalAccessException{
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
	
	public void AddIndependentVariableName(String independentVariables){
		/* 
		 * Add the name of the independent variable name at the end of the existing list
		 * of independent variables.  
		 */
	    this.independentVariableNames.add(independentVariables);
	  }
	
	public void AddIndependentVariableName(int index, String independentVariables){
		/* 
		 * Add the name of the independent variable name at <index> of the existing list
		 * of independent variables.  
		 */
	    this.independentVariableNames.add(index, independentVariables);
	  }
	
	public List<String> GetIndependentVariableNames() throws IllegalAccessException{
		/* 
		 * Get a list of the independent variables of the interaction model e.g.
		 * 		[neut%, mono%, neut%:GT]
		 */
		if(this.independentVariableNames == null){
			throw new IllegalAccessException("celltypes not set for this model");
		}
	    return(this.independentVariableNames);
	  }
	
	public void SetExpressionValues(double[] expressionValues){
		/* 
		 * Set the expression values (y) of the model. 
		 */
	    this.expressionValues = expressionValues;
	  }
	
	public double[] GetExpessionValues() throws IllegalAccessException{
		/* 
		 * Get the expression values (y) of the model. 
		 */
		if(this.expressionValues == null){
			throw new IllegalAccessException("expressionValues not set for this model");
		}
	    return(this.expressionValues);
	  }
	
	public void SetNoIntercept(Boolean noIntercept){
		/* 
		 * Set No Intercept for the model. If true, the intercept will be removed
		 */
	    this.noIntercept = noIntercept;
	  }
	
	public void SetCelltypes(List<String> celltypes) throws IllegalAccessException{
		/* 
		 * Get a list of all the celltypes given as input
		 */
		this.celltypes = celltypes;
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
	
	public Boolean GetNoIntercept() throws IllegalAccessException{
		/* 
		 * Get the No Intercept value for the model. If true, the intercept will be removed
		 */
		if(this.noIntercept == null){
			throw new IllegalAccessException("noIntercept not set for this model");
		}
	    return(this.noIntercept);
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
	
	public void SetModelName(String modelName){
	    this.modelName = modelName;
	  }
	
	public String GetModelName() throws IllegalAccessException{
		if(this.modelName == null){
			throw new IllegalAccessException("modelName name not set for this model");
		}
	    return(this.modelName);
	  }
		
	public void SetEstimateRegressionParametersStandardErrors(double[] estimateRegressionParametersStandardErrors){
		this.estimateRegressionParametersStandardErrors = estimateRegressionParametersStandardErrors;
	}
	public double[] GetEstimateRegressionParametersStandardErrors() throws IllegalAccessException{
		if(this.estimateRegressionParametersStandardErrors == null){
			throw new IllegalAccessException("estimateRegressionParametersStandardErrors not set for this model");
		}
		return(this.estimateRegressionParametersStandardErrors);
	}
	
	public void SetEstimateRegressionParameters(double[] estimatedRegressionParameters){
		this.estimatedRegressionParameters = estimatedRegressionParameters;
	}
	public double[] GetEstimateRegressionParameters() throws IllegalAccessException{
		if(this.estimatedRegressionParameters == null){
			throw new IllegalAccessException("estimatedRegressionParameters not set for this model");
		}
		return(this.estimatedRegressionParameters);
	}
	
	public void emptyObservedValues(){
		observedValues = null;
	}
	public void emptyExpressionValues(){
		expressionValues = null;
	}
}
