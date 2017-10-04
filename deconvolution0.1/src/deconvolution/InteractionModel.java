package deconvolution;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

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
	private HashMap<String, Character> genotypeOrder = new HashMap<String, Character>();
	// initialize with number so that we can test if it has been set
	private double sumOfSquares = -1;
	private double pvalue = -1;
	private int degreesOfFreedom = -1;
	private int modelLength = -1;
	private int numberOfTerms;

	/**
	 * Initialize object by setting the observed values size. Per QTL for each sample the observed values are each term of the 
	 * linear model. E.g. if the model is y = mono% + neut% + mono%:GT, the observedValues are
	 * [mono%, neut%, mono% * GT]
	 * 
	 * @ param sampleSize number of samples
	 * @ param numberOfTerms number of terms that the interaction model has
	 */
	public InteractionModel( int sampleSize, int numberOfTerms){
	    this.observedValues = new double[sampleSize][numberOfTerms];
	    this.numberOfTerms = numberOfTerms;
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
			throw new IllegalAccessException(String.format("celltypes not set for this model %s", this.modelName));
		}
	    return(this.independentVariableNames);
	  }
	
	public void setModelName(String modelName){
	    this.modelName = modelName;
	  }
	
	public String getModelName() throws IllegalAccessException{
		if(this.modelName == null){
			throw new IllegalAccessException(String.format("modelName name not set for this model %s", this.modelName));
		}
	    return(this.modelName);
	  }
		
	public void setEstimateRegressionParametersStandardErrors(double[] estimateRegressionParametersStandardErrors){
		this.estimateRegressionParametersStandardErrors = estimateRegressionParametersStandardErrors;
	}
	public double[] getEstimateRegressionParametersStandardErrors() throws IllegalAccessException{
		if(this.estimateRegressionParametersStandardErrors == null){
			throw new IllegalAccessException(String.format("estimateRegressionParametersStandardErrors not set for this model %s", this.modelName));
		}
		return(this.estimateRegressionParametersStandardErrors);
	}
	
	public void setEstimateRegressionParameters(double[] estimatedRegressionParameters){
		this.estimatedRegressionParameters = estimatedRegressionParameters;
	}
	public double[] getEstimateRegressionParameters() throws IllegalAccessException{
		if(this.estimatedRegressionParameters == null){
			throw new IllegalAccessException(String.format("estimatedRegressionParameters not set for model %s", this.modelName));
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
			throw new IllegalAccessException(String.format("pvalue not set for model %s", this.modelName));
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
	public void addGenotypeOrder(String celltype, char c){
		// 0 = dont swap genotype, 1 = swap genotype
		this.genotypeOrder.put(celltype, c);
	}
	public HashMap<String, Character> getGenotypeOrder(){
		// 0 = dont swap genotype, 1 = swap genotype
		return this.genotypeOrder;
	}

	public void setModelLength(){
		this.modelLength = this.observedValues.length;
	}
	
	public int getModelLength() throws IllegalAccessException {
		if(this.modelLength == -1){
			throw new IllegalAccessException(String.format("modelLength not set for model %s", this.modelName));
		}
		return this.modelLength;
	}

	public void setSumOfSquares(double sumOfSquares) {
		this.sumOfSquares = sumOfSquares;
	}
	
	public double getSumOfSquares() throws IllegalAccessException {
		if(this.sumOfSquares == -1){
			throw new IllegalAccessException(String.format("sumOfSquares not set for model %s", this.modelName));
		}
		return(this.sumOfSquares);
	}

	public void setDegreesOfFreedom(int degreesOfFreedom) {
		this.degreesOfFreedom  = degreesOfFreedom;
	}
	public int getDegreesOfFreedom() throws IllegalAccessException {
		if(this.degreesOfFreedom == -1){
			throw new IllegalAccessException(String.format("degreesOfFreedom not set for model %s", this.modelName));
		}
		return(this.degreesOfFreedom);
	}
	
	private int getNumberOfTerms(){
		return this.numberOfTerms;
	}
	
	/**
	 * Calculate the sum of squares, using Non-Negative Linear Regression, given a y expression vector with y ~
	 * model.
	 * Remove of the intercept (equivalent to y ~ model -1 in R) is hard-code in
	 * 
	 * This uses NNLS from Rochester Institute of technology. Documention here: https://www.cs.rit.edu/~ark/pj/doc/edu/rit/numeric/NonNegativeLeastSquares.html
	 * Can download JAR from: https://www.cs.rit.edu/~ark/pj.shtml#installed
	 * 
	 * @param model An InteractionModel object including the y vector expression values and ObservedValues (model)
	 * Such that
	 * test_trait ~ geno_A + lymph% + geno_A:geno_B it can be for one QTL
	 * [[2, 43.4, 86.8], [2, 40.3, 80.6]], for another QTL [[0, 46.7, 0],
	 * [0, 51.5, 0] [0, 48.7, 0]]
	 * @return 
	 * 
	 *  @return An nnls object
	 */
	public void calculateSumOfSquaresNNLS(double[] expressionValues) throws IOException, IllegalAccessException {
		// Use clone for a and b because the solve() method changes the matrix in place, want to keep the original values
		/*
		 * The MxN-element A matrix for the least squares
		 *
		 * problem. On input to the solve() method, a contains the
		 * matrix A. On output, a has been replaced with QA,
		 * where Q is an MxM-element orthogonal matrix
		 * generated during the solve() method's execution.
		 */
		double[][] a = this.observedValues;
		/*
		 * The M-element b vector for the least squares problem. On
		 * input to the solve() method, b contains the vector
		 * b. On output, b has been replaced with Qb, where
		 * Q is an MxM-element orthogonal matrix generated
		 * during the solve() method's execution.
		*/
		double[] b = expressionValues.clone();
		NonNegativeLeastSquares nnls = new NonNegativeLeastSquares(getModelLength(), getNumberOfTerms());
		// results contain:
		// normsqr: sqroot of the norm error vector
		// x: the parameters
		// For more, check out the Class documentation
		nnls.solve(a, b);
		setSumOfSquares(nnls.normsqr);
		setDegreesOfFreedom(expressionValues.length - (getNumberOfTerms() + 1));
		nnls = null;
	}
	
	/**
	 * Calculate the sum of squares, using Ordinary Linear Regression, given a y expression vector with y ~
	 * model. Remove of the intercept (equivalent to y ~ model -1 in R) is hard-code in
	 * 
	 * @param model An InteractionModel object including the y vector expression values and ObservedValues (model)
	 * Such that
	 * test_trait ~ geno_A + lymph% + geno_A:geno_B it can be for one QTL
	 * [[2, 43.4, 86.8], [2, 40.3, 80.6]], for another QTL [[0, 46.7, 0],
	 * [0, 51.5, 0] [0, 48.7, 0]]
	 * @param expressionValues Vector of expression values to use
	 * @param estimateRegressionParameters If true, estimate the regression parameters and save them in the InteractionModel
	 * @return A regression object
	 */
	public void calculateSumOfSquaresOLS(double[] expressionValues, Boolean estimateRegressionParameters) throws IOException, IllegalAccessException {
		// OLS = Ordinary Least Squares
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		// if GetIntercept is false, remove the intercept (Beta1) from the linear model
		regression.setNoIntercept(true);
		try{
			regression.newSampleData(expressionValues, this.getObservedValues());
		}
		catch (DimensionMismatchException e){
			DeconvolutionLogger.log.info(String.format("Length of expression and and genotype data not the same\nexpression length: %d\nobserved values length: %d\n", 
					expressionValues.length, this.getNumberOfTerms()));
			throw(e);
		}
		this.setSumOfSquares(regression.calculateResidualSumOfSquares());
		this.setDegreesOfFreedom(expressionValues.length - (this.getNumberOfTerms() + 1));
		if(estimateRegressionParameters){
			this.setEstimateRegressionParameters(regression.estimateRegressionParameters());
			//this.setEstimateRegressionParametersStandardErrors(regression.estimateRegressionParametersStandardErrors());

		}
	}
}

