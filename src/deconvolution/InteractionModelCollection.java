package deconvolution;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/*
 *  Collection of all the interaction models and their shared data (genotypes, expression etc)
 *  There are n + 1 interaction models, where n = the number of celltypes. One full model 
 *  with all celltypes, and for each celltype one model with the interaction term for that
 *  model removed
 */
public class InteractionModelCollection {
	private double[] expressionValues;
	private double[] genotypes;
	private double[] swappedGenotypes;
	private String qtlName;
	private HashMap<String, InteractionModel> interactionModels = new HashMap<String, InteractionModel>();
	private HashMap<String, Double> pvalues = new HashMap<String, Double>();
	private ArrayList<String> fullModelNames = new ArrayList<String>();
	private HashMap<String, ArrayList<String>> ctModelNames = new HashMap<String, ArrayList<String>>();
	private String bestFullModel;
	private HashMap<String, String> bestCtModel = new HashMap<String, String>();
	private CellCount cellCount; 
	private boolean useNNLS;
	private ArrayList<String> genotypeConfigurationsFullModel = new ArrayList<String> ();
	private ArrayList<String> genotypeConfigurationsCtModel = new ArrayList<String> ();
	private HashMap<String, String> celltypeOfModel = new HashMap<String,String>();
	/*
	 * Have to initialize instance with if NNLS or OLS will be used, and for that we need cellCounts
	 */
	public InteractionModelCollection(CellCount cellCount, boolean useNNLS) throws IllegalAccessException{
		setCellCount(cellCount);
		this.useNNLS = useNNLS;
		makeConfigurations();
	}
	
	public boolean getUseNNLS(){
		return(this.useNNLS);
	}

	/*
	 * Add interaction model to the collections
	 */
	private void addInteractionModel(InteractionModel interactionModel, String modelName, Boolean isFullModel){
		this.interactionModels.put(modelName, interactionModel);
		if(isFullModel){
			fullModelNames.add(modelName);
		}
		else{
			String cellType = this.celltypeOfModel.get(modelName);
			ctModelNames.putIfAbsent(cellType, new ArrayList<String>());
			ctModelNames.get(cellType).add(modelName);
		}
	}

	/*
	 * Get interaction model with modelName 
	 */
	public InteractionModel getInteractionModel(String modelName) throws IllegalAccessException{
		InteractionModel interactionModel = this.interactionModels.get(modelName);
		if(interactionModel == null){
			throw new IllegalAccessException(String.format("No model with name %s found", modelName));
		}
		return(interactionModel);
	}

	/*
	 * Remove interaction model with modelName
	 */
	public void removeInteractionModel(String modelName) throws IllegalAccessException{
		this.interactionModels.remove(modelName);
	}

	public String getCelltypeOfModel(String modelName){
		return this.celltypeOfModel.get(modelName);
	}
	
	/*
	 * Set CellCount
	 */
	private void setCellCount(CellCount cellCount){
		this.cellCount = cellCount;
	}
	
	public CellCount getCellCount() throws IllegalAccessException{
		if(this.cellCount == null){
			throw new IllegalAccessException("cellCount not set");
		}
		return this.cellCount;
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

	/*
	 * Get the genotypes of all the interaction models
	 */
	public double[] getSwappedGenotypes() throws IllegalAccessException {
		if(this.swappedGenotypes == null){
			throw new IllegalAccessException("genotypes not set for this model");
		}
		return this.swappedGenotypes;	
	}

	public void emptyExpressionValues(){
		this.expressionValues = null;
	}

	public void setGenotypes(double[] genotypes) {
		this.genotypes = genotypes;
		swapGenotypes();
	}
	/** 
	 * Get a list of all the celltypes given as input
	 */
	private void swapGenotypes(){
		this.swappedGenotypes = this.genotypes.clone();
		for(int i = 0; i < this.genotypes.length; i++) {
			this.swappedGenotypes[i] = 2 - this.genotypes[i];
		}
	}
	public void emptyGenotypes(){
		this.genotypes = null;
	}

	public void setQtlName(String qtlName){
		this.qtlName = qtlName;
	}

	public String getQtlName() throws IllegalAccessException{
		if(this.qtlName == null){
			throw new IllegalAccessException("QTL name not set");
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
	    return(pvalue);
	  }
	
	public ArrayList<String> getFullModelNames() throws IllegalAccessException{
		if(this.fullModelNames.size() == 0){
			throw new IllegalAccessException("No model names added to fullModelNames");
		}
		return this.fullModelNames;
	}
	public ArrayList<String> getCtModelNames(String celltype) throws IllegalAccessException{
		if(this.ctModelNames.size() == 0){
			throw new IllegalAccessException("No model names added to ctModelNames");
		}
		return this.ctModelNames.get(celltype);
	}
	
	private void setBestFullModel(String modelName){
		this.bestFullModel = modelName;
	}
	
	private void setBestCtModel(String celltype, String modelName){
		this.bestCtModel.put(celltype, modelName);
	}
	
	public InteractionModel getBestFullModel() throws IllegalAccessException{
		if(this.bestFullModel == null){
			throw new IllegalAccessException("bestFullModel not set");
		}
		return(this.getInteractionModel(this.bestFullModel));
	}
	// per celltype there is one best Ct model
	public InteractionModel getBestCtModel(String celltype) throws IllegalAccessException{
		String modelName = this.bestCtModel.get(celltype);
		if(modelName == null){
			throw new IllegalAccessException(String.format("No bestCtModel modelname for celltype %s", celltype));
		}
		InteractionModel bestCtModel = this.getInteractionModel(modelName);
		if(bestCtModel == null){
			throw new IllegalAccessException(String.format("bestCtModel not set for celltype %s", celltype));
		}
		return(this.getInteractionModel(this.bestCtModel.get(celltype)));
	}
	
	/*
	 * Go through all full models, calculate the regression statistics and select the model with the highest R2 as the new full model
	 */
	public void findBestFullModel() throws IllegalAccessException, IOException{
		// set to -1 so that first loop can be initialized
		double sumOfSquares = -1;
		for (String modelName : getFullModelNames()){
			InteractionModel fullModel = getInteractionModel(modelName);
			if(useNNLS){
				fullModel.calculateSumOfSquaresNNLS(getExpessionValues());
			}
			else{
				fullModel.calculateSumOfSquaresOLS(getExpessionValues());
			}
			if (sumOfSquares == -1){
				sumOfSquares = fullModel.getSumOfSquares();
				setBestFullModel(fullModel.getModelName());
			}
			if (fullModel.getSumOfSquares() <= sumOfSquares){
				sumOfSquares = fullModel.getSumOfSquares();
				setBestFullModel(fullModel.getModelName());
				fullModel.emptyObservedValues();
			}
			else{
				removeInteractionModel(fullModel.getModelName());
			}
		}
	}
	
	/*
	 * Go through all ct models, calculate the regression statistics and select the model with the highest R2 as the new ct model
	 * TODO: merge with findBestFullModel()
	 */
	public void findBestCtModel() throws IllegalAccessException, IOException{
		// set to -1 so that first loop can be initialized
		for(String celltype : getCellCount().getAllCelltypes()){
			double sumOfSquares = -1;
			for (String modelName : getCtModelNames(celltype)){
				InteractionModel ctModel = getInteractionModel(modelName);
				if(useNNLS){
					ctModel.calculateSumOfSquaresNNLS(getExpessionValues());
				}
				else{
					ctModel.calculateSumOfSquaresOLS(getExpessionValues());
				}
	
				if (sumOfSquares == -1){
					sumOfSquares = ctModel.getSumOfSquares();
	
					setBestCtModel(this.celltypeOfModel.get(modelName), ctModel.getModelName());
				}
				if (ctModel.getSumOfSquares() <= sumOfSquares){
					sumOfSquares = ctModel.getSumOfSquares();
					setBestCtModel(this.celltypeOfModel.get(modelName), ctModel.getModelName());
					ctModel.emptyObservedValues();
				}
				else{
					removeInteractionModel(ctModel.getModelName());
				}
			}
		}
	}
	
	/*
	 * Make the genotype configurations that will be used for the interaction terms 
	 */
	private void makeConfigurations() throws IllegalAccessException{
		if(getUseNNLS()){
			// if we use NNLS we need to see what genotype configuration works best, so get configuration permutations for all celltypes
			this.genotypeConfigurationsFullModel = Utils.binaryPermutations("",getCellCount().getAllCelltypes().size(), new ArrayList<String>());
			this.genotypeConfigurationsCtModel = Utils.binaryPermutations("",getCellCount().getAllCelltypes().size()-1, new ArrayList<String>());
		}else{
			// if we use OLS we just use default genotype orientation (all 0's)
			this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(getCellCount().getAllCelltypes().size(), "0")));
			this.genotypeConfigurationsCtModel.add(String.join("", Collections.nCopies(getCellCount().getAllCelltypes().size()-1, "0")));
		}
	}
	private ArrayList<String> getGenotypeConfigurationsFullModel() throws IllegalAccessException{
		if(this.genotypeConfigurationsFullModel == null){
			throw new IllegalAccessException("genotypeConfigurationsFullModel not set");
		}
		return this.genotypeConfigurationsFullModel;
	}
	private ArrayList<String> getGenotypeConfigurationsCtModel() throws IllegalAccessException{
		if(this.genotypeConfigurationsCtModel == null){
			throw new IllegalAccessException("genotypeConfigurationsCtModel not set");
		}
		return this.genotypeConfigurationsCtModel;
	}
	
	/**
	 * Construct the observed value matrices that are used for calculating the regression for the full model.
	 * Add all permutations of genotypes/swappedGenotypes (swappedGenotypes -> 0=2, 2=0)
	 * 
	 * @param ctModel InteractionModel object for saving the results
	 * @param m The current model that is being evaluated (for each celltype 1 model)
	 * @param fullModel InteractionModel object that contains information on the fullModel (such as expression values)
	 * 
	 * TODO: Move this to InteractionModel class. Also, merge overlapping code with createObservedValueMatricesCtModel
	 */
	public void createObservedValueMatricesFullModel() 
			throws IllegalAccessException{
		int numberOfTerms = getCellCount().getNumberOfCelltypes() * 2;
		// Have to test which genotype combination is the best, so 2**number of celltype loops
		for (String genotypeConfiguration : getGenotypeConfigurationsFullModel()){
			// things neded for fullModel defined outside of loop because every celltype model (ctModel) has to be compared to it
			InteractionModel fullModel = new InteractionModel(getCellCount().getNumberOfSamples(), 
															  numberOfTerms);
			fullModel.setGenotypeConfiguration(genotypeConfiguration);
			String modelName = String.format("fullModel_%s",genotypeConfiguration);
			fullModel.setModelName(modelName);
			addInteractionModel(fullModel, modelName, true);

			// number of terms + 1 because for full model all cell types are included
			for (int sampleIndex = 0; sampleIndex <= getCellCount().getNumberOfSamples()-1; sampleIndex++) {
				for (int celltypeIndex = 0; celltypeIndex < getCellCount().getNumberOfCelltypes(); celltypeIndex++) {

					double celltype_perc = getCellCount().getCellcountPercentages()[sampleIndex][celltypeIndex];
					// if i (cell type index) is the same as m (model index), don't add the interaction term of celltype:GT
					fullModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
					try {
						if(sampleIndex == 0){
							/** save the index of the variables related to current celltype so that this can be used later to calculate
							 * Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
							 * but for ctModel this is easiest method
							 */
							int[] index = new int[] {celltypeIndex, getCellCount().getNumberOfCelltypes() + celltypeIndex};
							fullModel.addCelltypeVariablesIndex(index);
								// add the celltype name at position i so that it gets in front of the celltype:GT
								fullModel.addIndependentVariableName(celltypeIndex, getCellCount().getCelltype(celltypeIndex));
								fullModel.addIndependentVariableName(getCellCount().getCelltype(celltypeIndex)+":GT");
								
						}
						// Have permutation of (2**number of celltypes) as binary ( so 00, 10, 01, 11 ), when 0 do normal genotype, 1 do swapped genotype
						double[] genotypes;
						char genotypeOrderAtCelltype = genotypeConfiguration.charAt(celltypeIndex);
						// Use the binary string permutation to decide if the genotype should be swapped or not
						if(genotypeOrderAtCelltype == '0'){
							genotypes = getGenotypes();
						} else{
							genotypes = getSwappedGenotypes();
						}
						fullModel.addObservedValue(celltype_perc * genotypes[sampleIndex], 
												   sampleIndex, getCellCount().getNumberOfCelltypes() + celltypeIndex);					
					} catch (ArrayIndexOutOfBoundsException error) {
						throw new RuntimeException(
								"The counts file and expression and/or genotype file do not have equal number of samples or QTLs",
								error);
					}
				}
			}
			fullModel.setModelLength();
		}
	}

	/**
	 * Construct the observed value matrices that are used for calculating the regression
	 * @param genotypeOrder 
	 * 
	 * @param InteractionModelCollection Collection of InteractionModel objects for saving the results
	 * @param genotypeOrder The order of genotypes to use, e.g. 010 means non swapped genotypes celltype 1, swapped genotypes celltype 2, non swapped genotypes celltype 3
	 * 
	 * TODO: Move this to InteractionModel class. Also, merge overlapping code with createObservedValueMatricesFullModel
	 */
	public void createObservedValueMatricesCtModels() 
			throws IllegalAccessException{
		int genotypeCounter = getCellCount().getNumberOfCelltypes();
		// -1 because one interaction term is removed
		int numberOfTerms = (getCellCount().getNumberOfCelltypes() * 2) - 1;
		for (String genotypeConfiguration : getGenotypeConfigurationsCtModel()){
			// m = model, there are equally many models as celltypes
			for (int modelIndex = 0; modelIndex < getCellCount().getNumberOfCelltypes(); modelIndex++) {
				InteractionModel ctModel = new InteractionModel(getCellCount().getNumberOfSamples(), numberOfTerms);	
				ctModel.setGenotypeConfiguration(genotypeConfiguration);
				// calculate p-value and save it, with other information, in a ctModel object. Then, add it to a list of these models to return as decon results
				String modelName = String.format("%s_%s", getCellCount().getCelltype(modelIndex), genotypeConfiguration);
				ctModel.setModelName(modelName);
				celltypeOfModel.put(modelName, getCellCount().getCelltype(modelIndex));
				addInteractionModel(ctModel,ctModel.getModelName(), false);	
				for (int sampleIndex = 0; sampleIndex <= getCellCount().getNumberOfSamples()-1; sampleIndex++) {
					int configurationIndex = 0;
					for (int celltypeIndex = 0; celltypeIndex < getCellCount().getNumberOfCelltypes(); celltypeIndex++) {
						// There is one fullModel including all celltypes add values for celltypePerc and interaction term of
						// celltypePerc * genotypePerc so that you get [[0.3, 0.6], [0.4, 0.8], [0.2, 0.4], [0.1, 0.2]]
						// where numberOfSamples = 1 and numberOfCellTypes = 4 with celltypePerc = 0.3, 0.4, 0.2, and 0.1 and genotype = 2
						// for each cell type is 1 model, celltype% * genotype without 1 celltype.
						// j+1 because j==0 is header
						double celltype_perc = getCellCount().getCellcountPercentages()[sampleIndex][celltypeIndex];
						ctModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
						if(sampleIndex == 0){
							// add the celltype name at position i so that it gets in front of the celltype:GT, but once
							try{
								ctModel.addIndependentVariableName(celltypeIndex, getCellCount().getCelltype(celltypeIndex));
							}
							catch(NullPointerException e){
								DeconvolutionLogger.log.info(String.format("Nullpoint exception with celltype %s", celltypeIndex));
								throw e;
							}
						}
	
						// if celltypeIndex is the same as m modelIndex, don't add the interaction term of celltype:GT
						if (celltypeIndex != modelIndex) {
								// Only add IndependentVariableName once per QTL (j==0)
								if(sampleIndex == 0){
	
									// Add the interaction term of celltype:genotype
									ctModel.addIndependentVariableName(getCellCount().getCelltype(celltypeIndex)+":GT");
									// save the index of the variables related to current celltype so that this can be used later to calculate
									// Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
									// but for ctModel this is easiest method
									int[] index = new int[] {celltypeIndex, getCellCount().getNumberOfCelltypes()-1+celltypeIndex};
									ctModel.addCelltypeVariablesIndex(index);
									// add the celltype name. This could be done with less code by getting it from IndependentVariableName, but this way 
									// it is explicit. Don't know if better.
								}
						try {
								double genotype = 0;
								// because the genotype configuration is of length (number of celltypes - 1), when a model is skipped we need to 
								// adjust all celltype indices from that point forward
								char genotypeOrderAtCelltype = genotypeConfiguration.charAt(configurationIndex);
								configurationIndex++;
								if(genotypeOrderAtCelltype == '0'){
									genotype = getGenotypes()[sampleIndex];
								}
								else if(genotypeOrderAtCelltype == '1'){
									genotype = getSwappedGenotypes()[sampleIndex];
								}
								else{
									throw new RuntimeException(String.format("Genotype order should be 0 or 1, was: %s", genotypeOrderAtCelltype));
								}
								ctModel.addObservedValue(celltype_perc * genotype, sampleIndex, genotypeCounter);

							} catch (ArrayIndexOutOfBoundsException error) {
								DeconvolutionLogger.log.info("ERROR: The counts file and expression and/or genotype file do not have equal number of samples or QTLs");
								throw error;
							}
							genotypeCounter++;
						}
						// if i==m there is not celltype:GT interaction term so only one index added to CelltypeVariables
						else if (sampleIndex == 0){
							int[] index = new int[] {celltypeIndex};
							ctModel.addCelltypeVariablesIndex(index);
						}
					}
					// because 1 of numberOfCelltypes + i needs to be skipped,
					// keeping it tracked with separate value is easier
					genotypeCounter = getCellCount().getNumberOfCelltypes();
				}
			ctModel.setModelLength();	
			}
		}
	}
}
