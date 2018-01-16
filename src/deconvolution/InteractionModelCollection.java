package deconvolution;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
	private double[] swappedGenotypes;
	private String qtlName;
	private HashMap<String, InteractionModel> interactionModels = new HashMap<String, InteractionModel>();
	private HashMap<String, Double> pvalues = new HashMap<String, Double>();
	private ArrayList<String> fullModelNames = new ArrayList<String>();
	private HashMap<String, ArrayList<String>> ctModelNames = new HashMap<String, ArrayList<String>>();
	private HashMap<String, String> ctModelName = new HashMap<String, String>();
	private String bestFullModel;
	private HashMap<String, String> bestCtModel = new HashMap<String, String>();
	private CellCount cellCount; 
	private boolean useNNLS;
	private List<String> genotypeConfigurationsFullModel = new ArrayList<String> ();
	private List<String> genotypeConfigurationsCtModel = new ArrayList<String> ();
	private HashMap<String, HashMap<String,String>> ctModelByGenotypeConfiguration = new HashMap<String, HashMap<String, String>>();
	private Double fullModelAIC;
	private HashMap<String, Double> ctModelAICs = new HashMap<String, Double>();
	private HashMap<String, String> ctModelsSameGenotypeConfigurationBestFullModel = new HashMap<String, String>();
	private HashMap<String, String> modelCelltype = new HashMap<String, String>();
	private HashMap<String, ArrayList<String>> genotypeConfigMap = new HashMap<String, ArrayList<String>>();
	private List<String> celltypes = new ArrayList<String>();
	private List<String> sampleNames = new ArrayList<String>();
	/*
	 * Have to initialize instance with if NNLS or OLS will be used, and for that we need cellCounts
	 */
	public InteractionModelCollection(CellCount cellCount, boolean useNNLS, String genotypeConfigurationType) throws IllegalAccessException{
		setCellCount(cellCount);
		this.useNNLS = useNNLS;
		makeConfigurations(genotypeConfigurationType);
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
			String cellType = getCelltypeOfModel(modelName);
			ctModelNames.putIfAbsent(cellType, new ArrayList<String>());
			ctModelNames.get(cellType).add(modelName);
		}
	}

	public String getCtModelName(String celltype){
		return(ctModelName.get(celltype));
	}

	public List<String> getAllCelltypes(){
		return(celltypes);
	}	
	
	public HashMap<String, String> getCtModelsByGenotypeConfiguration(String genotypeConfiguration){
		return(this.ctModelByGenotypeConfiguration.get(genotypeConfiguration));
	}

	/*
	 * Get interaction model with modelName 
	 */
	public InteractionModel getInteractionModel(String modelName) throws IllegalAccessException{
		InteractionModel interactionModel = this.interactionModels.get(modelName);
		return(interactionModel);
	}

	/*
	 * Remove interaction model with modelName
	 */
	public void removeInteractionModel(String modelName) throws IllegalAccessException{
		this.interactionModels.remove(modelName);
	}

	public String getCelltypeOfModel(String modelName){
		return modelName.split("_")[0];
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
	 */
	public void setExpressionValues(double[] expression){
		this.expressionValues = expression;
	}
	/** 
	 * Get the expression values (y) of all the interaction models. 
	 */
	public double[] getExpessionValues() throws IllegalAccessException{
		return(this.expressionValues);
	}
	/*
	 * Get the genotypes of all the interaction models
	 */
	public double[] getGenotypes() throws IllegalAccessException {
		return this.genotypes;
	}

	/*
	 * Get the genotypes of all the interaction models
	 */
	public double[] getSwappedGenotypes() throws IllegalAccessException {
		return this.swappedGenotypes;	
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

	public void setQtlName(String qtlName){
		this.qtlName = qtlName;
	}

	public String getQtlName() throws IllegalAccessException{
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
		return(pvalue);
	}

	public ArrayList<String> getFullModelNames() throws IllegalAccessException{
		return this.fullModelNames;
	}
	public ArrayList<String> getCtModelNames(String celltype) throws IllegalAccessException{
		return this.ctModelNames.get(celltype);
	}

	private void setBestFullModel(String modelName){
		this.bestFullModel = modelName;
	}

	private void setCtModelSameGenotypeConfigurationAsBestFullModel(String celltype, String modelName){
		this.ctModelsSameGenotypeConfigurationBestFullModel.put(celltype, modelName);
	}
	public String getCtModelSameGenotypeConfigurationAsBestFullModel(String celltype) throws IllegalAccessException{
		if(!this.ctModelsSameGenotypeConfigurationBestFullModel.containsKey(celltype)){
			InteractionModel bestFullModel = this.getBestFullModel();
			setCtModelSameGenotypeConfigurationAsBestFullModel(celltype, this.getCtModelsByGenotypeConfiguration(bestFullModel.getGenotypeConfiguration()).get(celltype));
		}
		return(this.ctModelsSameGenotypeConfigurationBestFullModel.get(celltype));
	}

	private void setBestCtModel(String celltype, String modelName){
		this.bestCtModel.put(celltype, modelName);
	}

	public InteractionModel getBestFullModel() throws IllegalAccessException{
		return(this.getInteractionModel(this.bestFullModel));
	}
	// per celltype there is one best Ct model
	public InteractionModel getBestCtModel(String celltype) throws IllegalAccessException{
		return(this.getInteractionModel(this.bestCtModel.get(celltype)));
	}

	private HashMap<String, String> getCtModelsWithConfiguration(String genotypConfiguration) throws IllegalAccessException{
		return(this.ctModelByGenotypeConfiguration.get(this.getBestFullModel().getGenotypeConfiguration()));
	}

	public void setAIC() throws IllegalAccessException{
		InteractionModel bestFullModel = getBestFullModel();
		bestFullModel.setAIC();
		this.fullModelAIC = bestFullModel.getAIC();
		HashMap<String, String> cellTypeCtModel = getCtModelsWithConfiguration(bestFullModel.getGenotypeConfiguration());
		for(String celltype : cellTypeCtModel.keySet()){
			String modelName = cellTypeCtModel.get(celltype);
			InteractionModel ctModel = this.getInteractionModel(modelName);
			ctModel.setAIC();
			this.ctModelAICs.put(modelName,ctModel.getAIC());
			ctModel.setAICdelta(bestFullModel.getAIC());
		}
	}

	public double getFullModelAIC() throws IllegalAccessException{
		return(fullModelAIC);
	}

	public double getCtModelAIC(String ctModelName) throws IllegalAccessException{
		return(ctModelAICs.get(ctModelName));
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
				fullModel.calculateSumOfSquaresOLS(getExpessionValues(), true);
			}
			if (sumOfSquares == -1){
				sumOfSquares = fullModel.getSumOfSquares();
			}
			if (fullModel.getSumOfSquares() <= sumOfSquares){
				sumOfSquares = fullModel.getSumOfSquares();
				setBestFullModel(fullModel.getModelName());
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
		for(String celltype : celltypes){
			double sumOfSquares = -1;
			for (String modelName : getCtModelNames(celltype)){
				InteractionModel ctModel = getInteractionModel(modelName);
				modelCelltype.put(modelName, celltype);
				if(useNNLS){
					ctModel.calculateSumOfSquaresNNLS(getExpessionValues());
				}
				else{
					ctModel.calculateSumOfSquaresOLS(getExpessionValues(), false);
				}

				if (sumOfSquares == -1){
					sumOfSquares = ctModel.getSumOfSquares();
					setBestCtModel(getCelltypeOfModel(modelName), ctModel.getModelName());
				}
				setCtModelByGenotypeConfiguration();
				double ctSumOfSquares = ctModel.getSumOfSquares();
				if (ctSumOfSquares <= sumOfSquares){
					sumOfSquares = ctSumOfSquares;
					setBestCtModel(getCelltypeOfModel(modelName), ctModel.getModelName());
				}
				else{
					// if the interaction model name of the full model is not the same as the model name of the 
					// CT model, we want to remove the ct model data to preserve RAM. If it is the same,
					// we keep it so that AIC can be calculated
					if(!ctModel.getModelName().equals(this.getCtModelSameGenotypeConfigurationAsBestFullModel(celltype))){
						removeInteractionModel(ctModel.getModelName());
					}
				}
			}
		}
	}

	/*
	 * Make the genotype configurations that will be used for the interaction terms 
	 */
	private void makeConfigurations(String genotypeConfigurationType) throws IllegalAccessException{
		if(getUseNNLS()){
			if(genotypeConfigurationType.equals("all")){
				// this gets all possible combinations, e.g. if 3 celltypes: 000, 001, 010, 100, 011, 101, 110, 111
				this.genotypeConfigurationsFullModel = Utils.binaryPermutations("",celltypes.size(), new ArrayList<String>());
			}else if(genotypeConfigurationType.equals("two")){
				// this gets two possible combinations, e.g. if 3 celltypes: 000, 111
				this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(celltypes.size(), "0")));
				this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(celltypes.size(), "1")));
			}else if(genotypeConfigurationType.equals("one")){
				// similar to "two", but can have one different, e.g. : 000, 111, 001, 010, 100
				this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(celltypes.size(), "0")));
				this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(celltypes.size(), "1")));
				for(int i = 0; i < celltypes.size(); ++i){
					StringBuilder genotypeConfiguration = new StringBuilder(String.join("", Collections.nCopies(celltypes.size(), "0")));
					genotypeConfiguration.setCharAt(i, '1');
					this.genotypeConfigurationsFullModel.add(genotypeConfiguration.toString());
				}
				for(int i = 0; i < celltypes.size(); ++i){
					StringBuilder genotypeConfiguration = new StringBuilder(String.join("", Collections.nCopies(celltypes.size(), "1")));
					genotypeConfiguration.setCharAt(i, '0');
					this.genotypeConfigurationsFullModel.add(genotypeConfiguration.toString());
				}
			}else{
				throw new RuntimeException("configurationType should be either \"all\" or \"two\", was: "+genotypeConfigurationType);
			}
			for(String genotypeConfiguration : genotypeConfigurationsFullModel){
				genotypeConfigMap.putIfAbsent(genotypeConfiguration, new ArrayList<String>());
				for(int i = 0; i < genotypeConfiguration.length()-1; i++){
					String s = genotypeConfiguration.substring(0, i);
					String s2 = genotypeConfiguration.substring(i+1, genotypeConfiguration.length());
					String newS = s.concat(s2);
					genotypeConfigMap.get(genotypeConfiguration).add(this.getCellCount().getCelltype(i)+"_"+newS);
				}
				String newS = genotypeConfiguration.substring(0, genotypeConfiguration.length()-1);
				genotypeConfigMap.get(genotypeConfiguration).add(this.getCellCount().getCelltype(genotypeConfiguration.length()-1)+"_"+newS);
			}

			this.genotypeConfigurationsCtModel = Utils.binaryPermutations("",celltypes.size()-1, new ArrayList<String>());
		}else{
			// if we use OLS we just use default genotype orientation (all 0's)
			String fullModelGenotypeConfiguration = String.join("", Collections.nCopies(celltypes.size(), "0"));
			this.genotypeConfigurationsFullModel.add(fullModelGenotypeConfiguration);
			String ctModelGenotypeConfiguration = String.join("", Collections.nCopies(celltypes.size()-1, "0"));
			this.genotypeConfigurationsCtModel.add(ctModelGenotypeConfiguration);
			genotypeConfigMap.putIfAbsent(fullModelGenotypeConfiguration, new ArrayList<String>());
			for(int i = 0; i < this.getCellCount().getNumberOfCelltypes(); i++){
				genotypeConfigMap.get(fullModelGenotypeConfiguration).add(this.getCellCount().getCelltype(i)+"_"+ctModelGenotypeConfiguration);
			}
		}
	}
	private List<String> getGenotypeConfigurationsFullModel() throws IllegalAccessException{
		return this.genotypeConfigurationsFullModel;
	}
	private List<String> getGenotypeConfigurationsCtModel() throws IllegalAccessException{
		return this.genotypeConfigurationsCtModel;
	}

	/**
	 * Construct the observed value matrices that are used for calculating the regression for the full model.
	 * Add all permutations of genotypes/swappedGenotypes (swappedGenotypes -> 0=2, 2=0)
	 * 
	 * TODO: Move this to InteractionModel class. Also, merge overlapping code with createObservedValueMatricesCtModel
	 */
	public void createObservedValueMatricesFullModel() 
			throws IllegalAccessException{
		CellCount cellCount = getCellCount();
		int numberOfCelltypes = cellCount.getNumberOfCelltypes();
		int numberOfSamples = cellCount.getNumberOfSamples();
		int numberOfTerms = numberOfCelltypes * 2;
		// Have to test which genotype combination is the best, so 2**number of celltype loops
		for (String genotypeConfiguration : getGenotypeConfigurationsFullModel()){
			// things neded for fullModel defined outside of loop because every celltype model (ctModel) has to be compared to it
			InteractionModel fullModel = new InteractionModel(numberOfSamples, 
					numberOfTerms);
			fullModel.setGenotypeConfiguration(genotypeConfiguration);
			String modelName = String.format("fullModel_%s",genotypeConfiguration);
			fullModel.setModelName(modelName);
			addInteractionModel(fullModel, modelName, true);

			// number of terms + 1 because for full model all cell types are included
			for (int sampleIndex = 0; sampleIndex <= numberOfSamples-1; sampleIndex++) {
				for (int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; celltypeIndex++) {

					double celltype_perc = cellCount.getCellcountPercentages()[sampleIndex][celltypeIndex];
					// if i (cell type index) is the same as m (model index), don't add the interaction term of celltype:GT
					fullModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
					try {
						if(sampleIndex == 0){
							/** save the index of the variables related to current celltype so that this can be used later to calculate
							 * Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
							 * but for ctModel this is easiest method
							 */
							int[] index = new int[] {celltypeIndex, numberOfCelltypes + celltypeIndex};
							fullModel.addCelltypeVariablesIndex(index);
							// add the celltype name at position i so that it gets in front of the celltype:GT
							fullModel.addIndependentVariableName(celltypeIndex, cellCount.getCelltype(celltypeIndex));
							fullModel.addIndependentVariableName(cellCount.getCelltype(celltypeIndex)+":GT");

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
								sampleIndex, numberOfCelltypes + celltypeIndex);					
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

	public void setCtModelByGenotypeConfiguration() throws IllegalAccessException{		
		InteractionModel bestFullModel = this.getBestFullModel();
		String bestFullModelgenotypeConfiguration = bestFullModel.getGenotypeConfiguration();
		HashMap<String, String> cellTypeCtModelName = new HashMap<String, String>();
		for(String ctModelName : genotypeConfigMap.get(bestFullModelgenotypeConfiguration)){
			cellTypeCtModelName.put(ctModelName.split("_")[0], ctModelName);
			ctModelByGenotypeConfiguration.putIfAbsent(bestFullModelgenotypeConfiguration, cellTypeCtModelName);			
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
		CellCount cellCount = getCellCount();
		int numberOfCelltypes = cellCount.getNumberOfCelltypes();
		int numberOfSamples = cellCount.getNumberOfSamples();
		int genotypeCounter = numberOfCelltypes;
		// -1 because one interaction term is removed
		int numberOfTerms = (numberOfCelltypes * 2) - 1;
		for (String genotypeConfiguration : getGenotypeConfigurationsCtModel()){
			// m = model, there are equally many models as celltypes
			for (int modelIndex = 0; modelIndex < numberOfCelltypes; modelIndex++) {
				InteractionModel ctModel = new InteractionModel(numberOfSamples, numberOfTerms);	
				ctModel.setGenotypeConfiguration(genotypeConfiguration);
				// calculate p-value and save it, with other information, in a ctModel object. Then, add it to a list of these models to return as decon results
				String modelName = String.format("%s_%s", cellCount.getCelltype(modelIndex), genotypeConfiguration);
				ctModel.setModelName(modelName);
				addInteractionModel(ctModel,ctModel.getModelName(), false);	
				for (int sampleIndex = 0; sampleIndex <= numberOfSamples-1; sampleIndex++) {
					int configurationIndex = 0;
					for (int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; celltypeIndex++) {
						// There is one fullModel including all celltypes add values for celltypePerc and interaction term of
						// celltypePerc * genotypePerc so that you get [[0.3, 0.6], [0.4, 0.8], [0.2, 0.4], [0.1, 0.2]]
						// where numberOfSamples = 1 and numberOfCellTypes = 4 with celltypePerc = 0.3, 0.4, 0.2, and 0.1 and genotype = 2
						// for each cell type is 1 model, celltype% * genotype without 1 celltype.
						// j+1 because j==0 is header
						double celltype_perc = cellCount.getCellcountPercentages()[sampleIndex][celltypeIndex];
						ctModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
						if(sampleIndex == 0){
							// add the celltype name at position i so that it gets in front of the celltype:GT, but once
							try{
								ctModel.addIndependentVariableName(celltypeIndex, cellCount.getCelltype(celltypeIndex));
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
								ctModel.addIndependentVariableName(cellCount.getCelltype(celltypeIndex)+":GT");
								// save the index of the variables related to current celltype so that this can be used later to calculate
								// Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
								// but for ctModel this is easiest method
								int[] index = new int[] {celltypeIndex, numberOfCelltypes-1+celltypeIndex};
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
					genotypeCounter = cellCount.getNumberOfCelltypes();
				}
				ctModel.setModelLength();

			}
		}
	}

	public void cleanUp(Boolean removePredictedValues) throws IllegalAccessException {
		this.expressionValues = null;
		this.genotypes = null;
		this.swappedGenotypes = null;
		for(InteractionModel interactionModel : this.interactionModels.values()){
			interactionModel.cleanUp(removePredictedValues);
		}
		this.genotypeConfigurationsCtModel = null;
		this.genotypeConfigurationsFullModel = null;
		this.cellCount = null;
	}

	public List<String> getSampleNames() {
		return sampleNames;
	}
}
