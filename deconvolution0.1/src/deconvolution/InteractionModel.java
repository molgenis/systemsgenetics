package deconvolution;

import java.util.ArrayList;
import java.util.List;

/*
 * InteractionModel contains the information (observed values, pvalue, regression parameter values etc)
 * of one single interaction model. All the InteractionModels that are used in the deconvolution are collected 
 * in InteractionModelCollection
 */
public class InteractionModel {
	/* independentVariables = the names of the independent variables, e.g. neut%, mono%, neut%:GT */
	private List<String> independentVariableNames = new ArrayList<String>();
	private List<int[]> celltypeVariablesIndex = new ArrayList <int[]>();
	private double[][] observedValues;
	private String modelName;
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

	public void setPvalue(double pvalue){
		this.pvalue = pvalue;
	}
	public double getPvalue() throws IllegalAccessException{
		if(this.pvalue == -1){
			throw new IllegalAccessException("pvalue not set for this model");
		}
		return this.pvalue;
	}
	public void setAlltIndependentVariableNames(List<String> list){
		for (int celltypeIndex = 0; celltypeIndex < list.size(); celltypeIndex++) {
			int[] index = new int[] {celltypeIndex, list.size() + celltypeIndex};
			addCelltypeVariablesIndex(index);
			addIndependentVariableName(celltypeIndex, list.get(celltypeIndex));
			addIndependentVariableName(list.get(celltypeIndex)+":GT");
		}
	}
	
}
