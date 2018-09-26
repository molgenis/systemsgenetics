package main.java.decon_eQTL_simple;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/*
 *  Collection of all the interaction models and their shared data (genotypes, expression etc)
 *  There are n + 1 interaction models, where n = the number of cell types. One full model 
 *  with all cell types, and for each cell type one model with the interaction term for that
 *  model removed
 */
public class InteractionModelCollection {
	private double[] expressionValues;
	private double[] genotypes;
	private double[] swappedGenotypes;
	private String qtlName;
	private HashMap<String, InteractionModel> interactionModels = new HashMap<String, InteractionModel>();
	private HashMap<String, InteractionModel> snpModels = new HashMap<String, InteractionModel>();
	private HashMap<String, Double> pvalues = new HashMap<String, Double>();
	private ArrayList<String> modelNames = new ArrayList<String>();
	private CellCount cellCount; 
	private List<String> celltypes = new ArrayList<String>();
	private List<String> sampleNames = new ArrayList<String>();

	/*
	 * Have to initialize instance with if NNLS or OLS will be used, and for that we need cellCounts
	 */
	public InteractionModelCollection(CellCount cellCount) throws IllegalAccessException{
		setCellCount(cellCount);
	}

	public List<String> getAllCelltypes(){
		return celltypes;
	}	

	/*
	 * Get interaction model with modelName 
	 */
	public InteractionModel getInteractionModel(String modelName) throws IllegalAccessException{
		InteractionModel interactionModel = this.interactionModels.get(modelName);
		return interactionModel;
	}
	/*
	 * Get snp model with modelName 
	 */
	public InteractionModel getSnpModel(String modelName) throws IllegalAccessException{
		InteractionModel snpModel = this.snpModels.get(modelName);
		return snpModel;
	}

	/*
	 * Set CellCount
	 */
	private void setCellCount(CellCount cellCount){
		this.cellCount = cellCount;
		celltypes = cellCount.getAllCelltypes();
		sampleNames = cellCount.getSampleNames();
	}

	public CellCount getCellCount() throws IllegalAccessException{
		return this.cellCount;
	}	
	/** 
	 * Set the expression values (y) for all the interaction models. 
	 * 
	 * @param expression	Expression vector
	 */
	public void setExpressionValues(double[] expression){
		this.expressionValues = expression;
	}
	/** 
	 * Get the expression values (y) of all the interaction models. 
	 * 
	 * @return Expression vector
	 */
	public double[] getExpessionValues(){
		return this.expressionValues;
	}
	/*
	 * Get the genotypes of all the interaction models
	 */
	public double[] getGenotypes(){
		return this.genotypes;
	}

	/*
	 * Get the genotypes of all the interaction models
	 */
	public double[] getSwappedGenotypes(){
		return this.swappedGenotypes;	
	}

	public void setGenotypes(double[] genotypes) {
		this.genotypes = genotypes;
	}

	public void setQtlName(String qtlName){
		this.qtlName = qtlName;
	}

	public String getQtlName() throws IllegalAccessException{
		return this.qtlName;
	}

	/*
	 * Each ctModel will have a p-value 
	 */
	public void setPvalue(Double pvalue, String modelName){
		this.pvalues.put(modelName, pvalue);
	}

	public Double getPvalue(String modelName) throws IllegalAccessException{
		Double pvalue = this.pvalues.get(modelName);
		return pvalue;
	}

	public ArrayList<String> getModelNames() throws IllegalAccessException{
		return this.modelNames;
	}


	/*
	 * Add interaction model to the collections
	 */
	private void addInteractionModel(InteractionModel interactionModel, String modelName){
		this.interactionModels.put(modelName, interactionModel);
		modelNames.add(modelName);
	}

	/*
	 * Go through all models and calculate the regression statistics
	 */
	public void calculateOlsForAllModels() throws IllegalAccessException, IOException{
		//System.out.println(this.qtlName);
		for (String modelName : getModelNames()){
		//	System.out.println(modelName);
			InteractionModel interactionModel = getInteractionModel(modelName);
			interactionModel.calculateSumOfSquaresOLS(getExpessionValues());
			InteractionModel snpModel = snpModels.get(modelName);
			snpModel.calculateSumOfSquaresOLS(getExpessionValues());
		//	System.out.println("------");
		}
		//System.exit(0);
	}


	/**
	 * Construct the observed value matrices that are used for calculating the regression for the full model.
	 *
	 * @throws IllegalAccessException	If cell counts file can not be read
	 */
	public void createObservedValueMatrices() 
			throws IllegalAccessException{
		CellCount cellCount = getCellCount();
		int numberOfCelltypes = cellCount.getNumberOfCelltypes();
		int numberOfSamples = cellCount.getNumberOfSamples();
		// number of terms is 3 because intercept + snp + CC + snp:CC
		int numberOfTerms = 3;
		// for each cell type we have to make 1 model, e.g. y ~ snp + snp : neut, y ~ snp + snp : mono etc
		for (int modelIndex = 0; modelIndex < numberOfCelltypes; modelIndex++) {
			InteractionModel interactionModel = new InteractionModel(numberOfSamples, 
					numberOfTerms);
			// number of terms is 1 becausee only snp
			InteractionModel snpModel = new InteractionModel(numberOfSamples, 1);
			
			String modelName = String.format("model_%s",cellCount.getCelltype(modelIndex) );
			interactionModel.setModelName(modelName);
			snpModel.setModelName(modelName);
			addInteractionModel(interactionModel, modelName);
			// this does same as addInteractionModel, but fel unnecesarry to make function for it
			this.snpModels.put(modelName, snpModel);
						
			// number of terms + 1 because for full model all cell types are included
			for (int sampleIndex = 0; sampleIndex <= numberOfSamples-1; ++sampleIndex) {
				double celltypePerc = cellCount.getCellCountPercentages()[sampleIndex][modelIndex];

				interactionModel.addObservedValue(genotypes[sampleIndex], sampleIndex, 0);
				interactionModel.addObservedValue(celltypePerc, sampleIndex, 1);
				interactionModel.addObservedValue(genotypes[sampleIndex]*celltypePerc, sampleIndex, 2);
				
				snpModel.addObservedValue(genotypes[sampleIndex], sampleIndex, 0);
			}
		}
	}



	public void cleanUp() throws IllegalAccessException {
		this.expressionValues = null;
		this.genotypes = null;
		this.swappedGenotypes = null;
		for(InteractionModel interactionModel : this.interactionModels.values()){
			interactionModel.cleanUp();
		}
		this.cellCount = null;
	}

	public List<String> getSampleNames() {
		return sampleNames;
	}
}
