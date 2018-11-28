package main.java.decon_eQTL;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
	private HashMap<String, Double> pvalues = new HashMap<String, Double>();
	private ArrayList<String> fullModelNames = new ArrayList<String>();
	private HashMap<String, ArrayList<String>> ctModelNames = new HashMap<String, ArrayList<String>>();
	private HashMap<String, ArrayList<String>> fullModelNamesByCelltype = new HashMap<String, ArrayList<String>>();
	private String bestFullModelName;
	private HashMap<String, String> bestCtModel = new HashMap<String, String>();
	private CellCount cellCount; 
	private List<String> genotypeConfigurationsFullModel = new ArrayList<String> ();
	private List<String> genotypeConfigurationsCtModel = new ArrayList<String> ();
	private HashMap<String, HashMap<String,String>> ctModelByGenotypeConfiguration = new HashMap<String, HashMap<String, String>>();
	private HashMap<String, String> modelCelltype = new HashMap<String, String>();
	private HashMap<String, ArrayList<String>> genotypeConfigMap = new HashMap<String, ArrayList<String>>();
	private List<String> celltypes = new ArrayList<String>();
	private List<String> sampleNames = new ArrayList<String>();

	/*
	 * Have to initialize instance with if NNLS or OLS will be used, and for that we need cellCounts
	 */
	public InteractionModelCollection(CellCount cellCount, String genotypeConfigurationType) throws IllegalAccessException{
		setCellCount(cellCount);
		makeConfigurations(genotypeConfigurationType);
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
	 * Remove interaction model with modelName
	 */
	public void removeInteractionModel(String modelName) throws IllegalAccessException{
		this.interactionModels.remove(modelName);
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
		return this.qtlName;
	}

	/*
	 * Each ctModel will have a p-value from ANOVA test with fullmodel, save it per ctModel
	 */
	public void setPvalue(Double pvalue, String modelName){
		this.pvalues.put(modelName, pvalue);
	}

	public Double getPvalue(String modelName) throws IllegalAccessException{
		Double pvalue = this.pvalues.get(modelName);
		return pvalue;
	}

	public ArrayList<String> getFullModelNames() throws IllegalAccessException{
		return this.fullModelNames;
	}
	public ArrayList<String> getCtModelNames(String celltype) throws IllegalAccessException{
		return this.ctModelNames.get(celltype);
	}

	private void setBestFullModelName(String modelName){
		this.bestFullModelName = modelName;
	}

	public InteractionModel getBestFullModel() throws IllegalAccessException{
		return this.getInteractionModel(this.bestFullModelName);
	}

	private void setBestCtModel(String celltype, String modelName){
		this.bestCtModel.put(celltype, modelName);
	}

	// per cell type there is one best Ct model
	public InteractionModel getBestCtModel(String celltype) throws IllegalAccessException{
		return this.getInteractionModel(this.bestCtModel.get(celltype));
	}

	/*
	 * Add interaction model to the collections
	 */
	private void addInteractionModel(InteractionModel interactionModel, String modelName, Boolean isFullModel){
		this.interactionModels.put(modelName, interactionModel);
		String cellType = interactionModel.getCelltypeName();
		if(isFullModel){
			fullModelNames.add(modelName);
			fullModelNamesByCelltype.putIfAbsent(cellType, new ArrayList<String>());
			fullModelNamesByCelltype.get(cellType).add(modelName);
		}
		else{
			ctModelNames.putIfAbsent(cellType, new ArrayList<String>());
			ctModelNames.get(cellType).add(modelName);
		}
	}

	/*
	 * Go through all full models, calculate the regression statistics and 
	 * select the model with the highest R2 as the new full model
	 */
	public void findBestFullModel() throws IllegalAccessException, IOException{
		// set to -1 so that first loop can be initialised
		double sumOfSquares = -1;
		for (String modelName : getFullModelNames()){
			InteractionModel fullModel = getInteractionModel(modelName);
			fullModel.calculateSumOfSquaresNNLS(getExpessionValues());
			if (sumOfSquares == -1){
				sumOfSquares = fullModel.getSumOfSquares();
			}
			if (fullModel.getSumOfSquares() <= sumOfSquares){
				sumOfSquares = fullModel.getSumOfSquares();
				setBestFullModelName(fullModel.getModelName());
			}
			else{
				removeInteractionModel(fullModel.getModelName());
			}
		}

	}

	/*
	 * Go through all cell type models, calculate the regression statistics and select the model with the highest R2 as the new cell type model
	 * TODO: merge with findBestFullModel()
	 */
	public void findBestCtModel() throws IllegalAccessException, IOException{
		// set to -1 so that first loop can be initialised
		for(String celltype : celltypes){
			double sumOfSquares = -1;
			for (String modelName : getCtModelNames(celltype)){
				InteractionModel ctModel = getInteractionModel(modelName);
				modelCelltype.put(modelName, celltype);

				ctModel.calculateSumOfSquaresNNLS(getExpessionValues());

				if (sumOfSquares == -1){
					sumOfSquares = ctModel.getSumOfSquares();
					setBestCtModel(ctModel.getCelltypeName(), ctModel.getModelName());
				}
				setCtModelByGenotypeConfiguration();


				double ctSumOfSquares = ctModel.getSumOfSquares();
				if (ctSumOfSquares <= sumOfSquares){
					sumOfSquares = ctSumOfSquares;
					setBestCtModel(ctModel.getCelltypeName(), ctModel.getModelName());
				}
			}
		}
	}

	/*
	 * Make the genotype configurations that will be used for the interaction terms 
	 */
	private void makeConfigurations(String genotypeConfigurationType) throws IllegalAccessException{

		if(genotypeConfigurationType.equals("all")){
			// this gets all possible combinations, e.g. if 3 celltypes: 000, 001, 010, 100, 011, 101, 110, 111
			this.genotypeConfigurationsFullModel = Utils.binaryPermutations("",celltypes.size(), new ArrayList<String>());
			this.genotypeConfigurationsCtModel = Utils.binaryPermutations("",celltypes.size()-1, new ArrayList<String>());
		}else if(genotypeConfigurationType.equals("two")){
			// this gets two possible combinations, e.g. if 3 celltypes: 000, 111
			this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(celltypes.size(), "0")));
			this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(celltypes.size(), "1")));
			
			this.genotypeConfigurationsCtModel.add(String.join("", Collections.nCopies(celltypes.size()-1, "0")));
			this.genotypeConfigurationsCtModel.add(String.join("", Collections.nCopies(celltypes.size()-1, "1")));
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
			
			
			this.genotypeConfigurationsCtModel.add(String.join("", Collections.nCopies(celltypes.size()-1, "0")));
			this.genotypeConfigurationsCtModel.add(String.join("", Collections.nCopies(celltypes.size()-1, "1")));
			for(int i = 0; i < celltypes.size()-1; ++i){
				StringBuilder genotypeConfiguration = new StringBuilder(String.join("", Collections.nCopies(celltypes.size()-1, "0")));
				genotypeConfiguration.setCharAt(i, '1');
				this.genotypeConfigurationsCtModel.add(genotypeConfiguration.toString());
			}
			for(int i = 0; i < celltypes.size()-1; ++i){
				StringBuilder genotypeConfiguration = new StringBuilder(String.join("", Collections.nCopies(celltypes.size()-1, "1")));
				genotypeConfiguration.setCharAt(i, '0');
				this.genotypeConfigurationsCtModel.add(genotypeConfiguration.toString());
			}
	
			System.exit(0);

		}else if(genotypeConfigurationType.equals("NONE")){
			// this gets no swapping, so all are 000
			this.genotypeConfigurationsFullModel.add(String.join("", Collections.nCopies(celltypes.size(), "0")));
			this.genotypeConfigurationsCtModel.add(String.join("", Collections.nCopies(celltypes.size()-1, "0")));

		}else{
			throw new RuntimeException("configurationType should be either \"all\" or \"two\", was: "+genotypeConfigurationType);
		}


		for(String genotypeConfiguration : genotypeConfigurationsFullModel){
			genotypeConfigMap.putIfAbsent(genotypeConfiguration, new ArrayList<String>());
			for(int i = 0; i < genotypeConfiguration.length()-1; i++){
				String s = genotypeConfiguration.substring(0, i);
				String s2 = genotypeConfiguration.substring(i+1, genotypeConfiguration.length());
				String newS = s.concat(s2);
				String ctModelName = "ctModel_"+this.getCellCount().getCelltype(i)+"_"+newS;
				genotypeConfigMap.get(genotypeConfiguration).add(ctModelName);
			}
			String newS = genotypeConfiguration.substring(0, genotypeConfiguration.length()-1);
			String ctModelName = "ctModel_"+this.getCellCount().getCelltype(genotypeConfiguration.length()-1)+"_"+newS;
			genotypeConfigMap.get(genotypeConfiguration).add(ctModelName);
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
	 * Add all permutations of genotypes/swappedGenotypes (swappedGenotypes 0=2, 2=0)
	 * 
	 * TODO: Move this to InteractionModel class. Also, merge overlapping code with createObservedValueMatricesCtModel
	 * 
	 * @throws IllegalAccessException	If cell counts file can not be read
	 */
	public void createObservedValueMatricesFullModel(Boolean addGenotypeTerm) 
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
			if(addGenotypeTerm) {
				// add one to number of terms because we are adding the genotype as a separate term
				fullModel = new InteractionModel(numberOfSamples, numberOfTerms+1);
			}
			fullModel.setGenotypeConfiguration(genotypeConfiguration);
			String modelName = String.format("fullModel_%s",genotypeConfiguration);
			fullModel.setModelName(modelName);
			addInteractionModel(fullModel, modelName, true);

			for(int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; ++celltypeIndex){
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

			// number of terms + 1 because for full model all cell types are included
			for (int sampleIndex = 0; sampleIndex <= numberOfSamples-1; ++sampleIndex) {
				for (int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; ++celltypeIndex) {

					double celltypePerc = cellCount.getCellCountPercentages()[sampleIndex][celltypeIndex];
					// if i (cell type index) is the same as m (model index), don't add the interaction term of celltype:GT
					fullModel.addObservedValue(celltypePerc, sampleIndex, celltypeIndex);
					// Have permutation of (2**number of celltypes) as binary ( so 00, 10, 01, 11 ), when 0 do normal genotype, 1 do swapped genotype
					double[] genotypes;
					char genotypeOrderAtCelltype = genotypeConfiguration.charAt(celltypeIndex);
					// Use the binary string permutation to decide if the genotype should be swapped or not
					if(genotypeOrderAtCelltype == '0'){
						genotypes = getGenotypes();
					} else{
						genotypes = getSwappedGenotypes();
					}

					fullModel.addObservedValue(celltypePerc * genotypes[sampleIndex], 
							sampleIndex, numberOfCelltypes + celltypeIndex);


				}
				if(addGenotypeTerm) {
					fullModel.addIndependentVariableName("GT");
					fullModel.addObservedValue(genotypes[sampleIndex], 
							sampleIndex, (numberOfCelltypes*2));
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
			String celltype = ctModelName.split("ctModel_")[1].split("_")[0];
			cellTypeCtModelName.put(celltype, ctModelName);
			ctModelByGenotypeConfiguration.putIfAbsent(bestFullModelgenotypeConfiguration, cellTypeCtModelName);			
		}
	}

	/**
	 * Construct the observed value matrices that are used for calculating the regression
	 * 
	 * TODO: Move this to InteractionModel class. Also, merge overlapping code with createObservedValueMatricesFullModel
	 * 
	 * @throws IllegalAccessException	If cell counts file can not be read
	 */
	public void createObservedValueMatricesCtModels(Boolean addGenotypeTerm) 
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
				if(addGenotypeTerm) {
					// add one to number of terms because we are adding the genotype as a separate term
					ctModel = new InteractionModel(numberOfSamples, numberOfTerms+1);
				}
				ctModel.setGenotypeConfiguration(genotypeConfiguration);
				// calculate p-value and save it, with other information, in a ctModel object. 
				// Then, add it to a list of these models to return as decon results
				String celltypeName = cellCount.getCelltype(modelIndex);
				String modelName = String.format("ctModel_%s_%s", celltypeName, genotypeConfiguration);
				ctModel.setModelName(modelName);
				ctModel.setCelltypeName(celltypeName);
				addInteractionModel(ctModel,ctModel.getModelName(), false);	
				for (int sampleIndex = 0; sampleIndex <= numberOfSamples-1; sampleIndex++) {
					int configurationIndex = 0;
					for (int celltypeIndex = 0; celltypeIndex < numberOfCelltypes; celltypeIndex++) {
						// There is one fullModel including all celltypes add values for celltypePerc and interaction term of
						// celltypePerc * genotypePerc so that you get [[0.3, 0.6], [0.4, 0.8], [0.2, 0.4], [0.1, 0.2]]
						// where numberOfSamples = 1 and numberOfCellTypes = 4 with celltypePerc = 0.3, 0.4, 0.2, and 0.1 and genotype = 2
						// for each cell type is 1 model, celltype% * genotype without 1 celltype.
						// j+1 because j==0 is header
						double celltype_perc = cellCount.getCellCountPercentages()[sampleIndex][celltypeIndex];
						ctModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
						if(sampleIndex == 0){
							// add the celltype name at position i so that it gets in front of the celltype:GT, but once
							try{
								ctModel.addIndependentVariableName(celltypeIndex, celltypeName);
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
								
							}
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

							genotypeCounter++;
						}
						// if i==m there is not celltype:GT interaction term so only one index added to CelltypeVariables
						else if (sampleIndex == 0){
							int[] index = new int[] {celltypeIndex};
							ctModel.addCelltypeVariablesIndex(index);
						}
					}
					if(addGenotypeTerm) {
						ctModel.addIndependentVariableName("GT");
						ctModel.addObservedValue(genotypes[sampleIndex], 
								sampleIndex, (numberOfCelltypes*2)-1);
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
