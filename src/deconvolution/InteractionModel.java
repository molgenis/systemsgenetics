package deconvolution;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;

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
	public InteractionModel(){};
	private String genotypeConfiguration;
	// initialize with number so that we can test if it has been set
	private Double sumOfSquares;
	private Double pvalue;
	private Integer degreesOfFreedom;
	private Integer modelLength;
	private Integer numberOfTerms;
	private double[] estimatedRegressionParameters;
	private Double AIC;
	private double[] residuals;
	private double AICdelta;
	private double[] estimatedRegressionParametersStandardError;
	private double estimatedStandardError;
	private double[] predictedValues;
	private String celltypeName;
	
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
	    return(this.independentVariableNames);
	  }
	
	public void setModelName(String modelName){
	    this.modelName = modelName;
	  }
	
	public String getModelName() throws IllegalAccessException{
	    return(this.modelName);
	}	

	public void setPvalue(double pvalue){
		this.pvalue = pvalue;
	}
	public double getPvalue() throws IllegalAccessException{
		return this.pvalue;
	}
	public void setAlltIndependentVariableNames(List<String> list){
		for (int celltypeIndex = 0; celltypeIndex < list.size(); ++celltypeIndex) {
			int[] index = new int[] {celltypeIndex, list.size() + celltypeIndex};
			addCelltypeVariablesIndex(index);
			addIndependentVariableName(celltypeIndex, list.get(celltypeIndex));
			addIndependentVariableName(list.get(celltypeIndex)+":GT");
		}
	}
	public void setGenotypeConfiguration(String genotypeConfiguration){
		// 0 = dont swap genotype, 1 = swap genotype
		this.genotypeConfiguration = genotypeConfiguration;
	}
	public String getGenotypeConfiguration() throws IllegalAccessException{
		return this.genotypeConfiguration;
	}

	public void setModelLength(){
		this.modelLength = this.observedValues.length;
	}
	
	public int getModelLength() throws IllegalAccessException {
		return this.modelLength;
	}

	public void setSumOfSquares(double sumOfSquares) {
		this.sumOfSquares = sumOfSquares;
	}
	
	public double getSumOfSquares() throws IllegalAccessException {
		return(this.sumOfSquares);
	}
	
	public void setAIC() {
		// +1 because error term is also a parameter
		this.AIC = Statistics.AIC(this.residuals, this.getNumberOfTerms()+1);
	}
	
	public double getAIC() {
		return(this.AIC);
	}

	public void setDegreesOfFreedom(int degreesOfFreedom) {
		this.degreesOfFreedom  = degreesOfFreedom;
	}
	public int getDegreesOfFreedom() throws IllegalAccessException {
		return(this.degreesOfFreedom);
	}
	
	public int getNumberOfTerms(){
		return this.numberOfTerms;
	}
	
	
	public void setEstimatedRegressionParameters(double[] estimatedRegressionParamters){
		this.estimatedRegressionParameters = estimatedRegressionParamters;
	}
	
	public double[] getEstimateRegressionParameters() throws IllegalAccessException{
		return(this.estimatedRegressionParameters);
	}
	
	private void setResiduals(double[] residuals) {
		this.residuals = residuals;
	}
	public double[] getResiduals() {
		return(this.residuals);
	}
	
	private void setPredictedValues(double[] predictedValues) {
		this.predictedValues = predictedValues;
	}
	public double[] getPredictedValues(){
		return this.predictedValues;
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
		NonNegativeLeastSquares nnls = new NonNegativeLeastSquares();
		
		try{
			nnls.newSampleData(expressionValues, this.getObservedValues());
		}
		catch (DimensionMismatchException e){
			DeconvolutionLogger.log.info(String.format("Length of expression and genotype data not the same\nexpression length: %d\nobserved values length: %d\n", 
					expressionValues.length, this.getNumberOfTerms()));
			throw(e);
		}
		
		// results contain:
		// normsqr: sqroot of the norm error vector
		// x: the parameters
		// For more, check out the Class documentation
		double[] estimatedRegressionParameters = nnls.estimateRegressionParameters();
		setEstimatedRegressionParameters(estimatedRegressionParameters);
		
		setSumOfSquares(nnls.calculateResidualSumOfSquares());
		setDegreesOfFreedom(expressionValues.length - (getNumberOfTerms() + 1));
		
		double[] predictedValues = nnls.getPredictedExpressionValues();
		double[] residuals = nnls.estimateResiduals();
		//setEstimatedStandardError(nnls.estimateRegressionStandardError());
		setResiduals(residuals);
		setPredictedValues(predictedValues);
	}

	public void setAICdelta(double comparativeAIC) {
		// if comparateAIC < this.getAIC(): comparativeAIC = better
		this.AICdelta = comparativeAIC - this.getAIC();
	}
	public double getAICdelta() {
		return this.AICdelta;
	}

	public void setEstimatedRegressionParametersStandardErrors(double[] estimatedRegressionParametersStandardErrors) {
		this.estimatedRegressionParametersStandardError = estimatedRegressionParametersStandardErrors;
	}

	public double[] getEstimatedRegressionParametersStandardErrors() {
		return this.estimatedRegressionParametersStandardError;
	}
	
	public void setEstimatedStandardError(double estimatedStandardError){
		this.estimatedStandardError = estimatedStandardError;
	}
	public double getEstimatedStandardError(){
		return this.estimatedStandardError;
	}

	public void cleanUp(Boolean removePredictedValues) {
		this.observedValues = null;
		this.residuals = null;
		if(removePredictedValues){
			this.predictedValues = null;
		}
	}

	public void setCelltypeName(String celltypeName) {
		this.celltypeName = celltypeName;
	}
	
	public String getCelltypeName(){
		return(this.celltypeName);
	}
}

