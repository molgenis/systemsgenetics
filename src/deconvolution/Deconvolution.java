package deconvolution;

import deconvolution.Qtl;
import deconvolution.Validate;
import deconvolution.CellCount;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.time.DurationFormatUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import JSci.maths.statistics.FDistribution;

public class Deconvolution {
	private static CommandLineOptions commandLineOptions = new CommandLineOptions(); 
	private static CellCount cellCounts;
	private static List<String> filteredQTLsOutput = new ArrayList<String>();
	private static int QTLsFiltered = 0;
	private static String outputFolder;
	// factory method for making static variable that can throw an exception


	/**
	 * Deconvolutes a set of QTLs given the expression levels, genotypes,
	 * and cell counts. Calculates the p-values for the deconvoluted QTLs
	 * and writes them to an outfile
	 * 
	 * @param args List of command line arguments
	 * 
	 * TODO: Crashes if output folder does not exist, either add clear error message or make output folder 
	 */
	public static void main(String[] args) throws Exception {
		commandLineOptions.parseCommandLine(args);
		outputFolder = commandLineOptions.getOutfolder();
		cellCounts = new CellCount(commandLineOptions.getCellcountFile());


		if (commandLineOptions.getNumberOfPermutations() > 1) {
			// permutationResults contains a hashmap of all pvalues per celltype gotten through permutation testing
			//PermutationResult permutationResults = permutationTest(expressionFile, genotypeFile, cellcountTable);
			throw new java.lang.UnsupportedOperationException("Not implemented yet.");
		}
		runDeconPerGeneSnpPair();
	}

	/**
	 * For each of the gene-SNP pair in the SnpsToTestFile run deconvolution
	 */
	private static void runDeconPerGeneSnpPair() throws IOException, IllegalAccessException, RuntimeException{
		HashMap<String,ArrayList<String>> geneSnpPairs = Utils.parseSnpPerGeneFile(commandLineOptions.getSnpsToTestFile());
		String expressionFile = commandLineOptions.getExpressionFile();
		DeconvolutionLogger.log.info(String.format("Parse expression data from %s",expressionFile));
		ExpressionData expressionData = new ExpressionData(expressionFile);
		DeconvolutionLogger.log.info("Done");
		String genotypeFile = commandLineOptions.getGenotypeFile();
		DeconvolutionLogger.log.info(String.format("Parse genotype data from %s",genotypeFile));
		GenotypeData genotypeData = new GenotypeData(genotypeFile);
		DeconvolutionLogger.log.info("Done");
		if (!Utils.equalLists(expressionData.getSampleNames(), genotypeData.getSampleNames())){
			throw new RuntimeException("Samplenames not the same in expression and genotype file.\nexpression sampleNames, genotype sampleNames: \n\n"+
					Arrays.toString(expressionData.getSampleNames().toArray())+"\n "+Arrays.toString(genotypeData.getSampleNames().toArray()));
		}
		//file to write all samples in that got filtered out
		Path filteredQTLsFile = Paths.get(outputFolder+"filteredQTLs.csv");
		filteredQTLsOutput.add("QTL\treason");

		int whileIndex = 0;
		long time = System.currentTimeMillis();
		List<DeconvolutionResult> deconvolutionResults = new ArrayList<DeconvolutionResult>();
		int QTLsTotal = 0;
		HashMap<String, double[]> geneExpressionLevels = expressionData.getGeneExpression();
		int skippedGenotypeGeneCombinations = 0;
		for(String gene : geneSnpPairs.keySet()){
			for(String genotype : geneSnpPairs.get(gene)){
				if(commandLineOptions.getTestRun() && whileIndex == 100){
					break;
				}
				if (whileIndex % 100 == 0) {
					long completedIn = System.currentTimeMillis() - time;
					DeconvolutionLogger.log.info(String.format("Processed %d gene-SNP pairs - %s - skipped %d gene-SNP combinations", whileIndex, DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS"), skippedGenotypeGeneCombinations));
				}
				whileIndex++;
				String qtlName = gene+'_'+genotype;
				try{
					QTLsTotal++;
					try{
						deconvolutionResults.add(deconvolution(geneExpressionLevels.get(gene), genotypeData.getGenotypes().get(genotype), qtlName));
					}
					catch(IllegalAccessException e){
						if(commandLineOptions.getSkipGenotypes()){
							skippedGenotypeGeneCombinations++;
							continue;
						}
						else{
							double[] genes = geneExpressionLevels.get(gene);
							double[] genotypes = genotypeData.getGenotypes().get(genotype);
							if (genes == null) {
								DeconvolutionLogger.log.info(String.format("gene %s in SNP-gene pair file but not in expression data file",gene));
							}

							if (genotypes == null) {
								DeconvolutionLogger.log.info(String.format("genotype %s in SNP-gene pair file but not in genotype data file",genotype));
							}
							throw e;
						}
					}
				}
				// If there are not enough samples per genotype, skip this QTL
				catch(NotEnoughGenotypesException e){
					if(!commandLineOptions.getFilterSamples()){
						deconvolutionResults.add(setPvaluesNA(qtlName));
					}
					else{
						QTLsFiltered++;
					}
					filteredQTLsOutput.add(qtlName+"\tNot enough genotypes (e.g. AA and AB but no BB)");
				}
				catch(NotEnoughSamplesPerGenotypeException e){
					if(!commandLineOptions.getFilterSamples()){
						deconvolutionResults.add(setPvaluesNA(qtlName));
					}
					else{
						QTLsFiltered++;
					}
					filteredQTLsOutput.add(qtlName+"\tNot enough samples per genotype");
				}
			}
		}
		writeDeconvolutionResults(deconvolutionResults);
		DeconvolutionLogger.log.info(String.format("Skipped %d gene-SNP combinations (because genotype in SNP-pair file but not in genotype file)",skippedGenotypeGeneCombinations));
		DeconvolutionLogger.log.info(String.format("QTLs passed: %d", QTLsTotal-(QTLsFiltered+skippedGenotypeGeneCombinations)));
		DeconvolutionLogger.log.info(String.format("QTLs filtered: %d", QTLsFiltered));
		DeconvolutionLogger.log.info(String.format("Total: %d",QTLsTotal-skippedGenotypeGeneCombinations));
		Files.write(filteredQTLsFile, filteredQTLsOutput, Charset.forName("UTF-8"));
		String validate = commandLineOptions.getValidationFile();
		if (validate != null){
			DeconvolutionLogger.log.info(String.format("Validating results..."));
			new Validate(deconvolutionResults, validate, outputFolder+"/figures/");
		}
	}

	private static void writeDeconvolutionResults(List<DeconvolutionResult> deconvolutionResults) throws IllegalAccessException, IOException{
		/**
		 * Append the deconvolution pvalues and, if the flag is set, the coefficients to the deconvolution result file
		 * 
		 * @param deconvolutionResult The deconvolutionresult
		 */
		List<String> output = new ArrayList<String>();
		String header = "\t"+Utils.listToTabSeparatedString(cellCounts.getCelltypes(), "_pvalue");
		DeconvolutionLogger.log.info("Getting decon result with full model info for writing the header");
		// celltypes.size()*2 because there are twice as many betas as celltypes (CC% & CC%:GT)
		InteractionModel bestFullModel = deconvolutionResults.get(0).getInteractionModelCollection().getBestFullModel();
		for(int i = 1; i < cellCounts.getNumberOfCelltypes()*2 + 1; i++){
			try{
				header += "\tBeta" + Integer.toString(i) +"_"+bestFullModel.getIndependentVariableNames().get(i-1);
			}
			catch(IndexOutOfBoundsException e){
				DeconvolutionLogger.log.info(String.format("DeconvolutionResult index error with beta %d", i));
				throw e;
			}
		}

		header += "\tgenotypeConfiguration";
		// TODO: Don't know how to get standard error for NNLS so I removed this for now, also we didnt really use the error.
		// 		 but want to leave it here incase I will use it later
		//for(int i = 1; i < cellCounts.getNumberOfCelltypes()*2 + 1; i++){
		//	header += "\tBetaStandardError" + Integer.toString(i); 
		//}

		if(commandLineOptions.getWholeBloodQTL()){
			header += "\tSpearman correlation expression~GT\tSpearman correlation p-value";
		}
		output.add(header);
		for(DeconvolutionResult deconvolutionResult : deconvolutionResults){
			bestFullModel = deconvolutionResult.getInteractionModelCollection().getBestFullModel();
			if(commandLineOptions.getOnlyOutputSignificant() && commandLineOptions.getFilterSamples()){
				if(Collections.min(deconvolutionResult.getPvalues()) > 0.05){
					QTLsFiltered++;
					filteredQTLsOutput.add(deconvolutionResult.getQtlName()+"\tNone of the celltypes had a significant p-value");
					continue;
				}
			}
			String results = "";
			results += deconvolutionResult.getQtlName()+"\t"+Utils.listToTabSeparatedString(deconvolutionResult.getPvalues());
			try{
				results += "\t"+Utils.listToTabSeparatedString(bestFullModel.getEstimateRegressionParameters());
				//results += "\t"+Utils.listToTabSeparatedString(deconvolutionResult.getInteractionModelCollection().getBestFullModel().getEstimateRegressionParametersStandardErrors());
			}catch (java.lang.IllegalAccessException e){
				// if -m is set not all deconvolution resuts will have a full model. If not, set betas to NA
				String str = "\tNA";
				// celltype.size()*4 because for each celltype are 2 betas + 2 standard error of betas
				results += StringUtils.repeat(str, cellCounts.getNumberOfCelltypes()*4);
			}

			results += "\t"+bestFullModel.getGenotypeOrder();
			
			if(commandLineOptions.getWholeBloodQTL()){
				results += "\t"+Double.toString(deconvolutionResult.getWholeBloodQTL());
				results += "\t"+Double.toString(deconvolutionResult.getWholeBloodQTLpvalue());
			}
			output.add(results);	
		}

		Path file = Paths.get(outputFolder+commandLineOptions.getOutfile());
		Files.write(file, output, Charset.forName("UTF-8"));
		DeconvolutionLogger.log.info(String.format("Deconvolution output written to %s", file.toAbsolutePath()));
		DeconvolutionLogger.log.info(String.format("Files with additional info in  %s", outputFolder));
	}

	private static DeconvolutionResult setPvaluesNA(String qtlName) throws IllegalAccessException{
		List<Double> pvalues = new ArrayList<Double>();
		for (int i = 0; i < cellCounts.getNumberOfCelltypes(); i++){
			pvalues.add(333.0);
		}
		InteractionModel dummyModel = new InteractionModel();
		dummyModel.setModelName("dummy");
		dummyModel.setAlltIndependentVariableNames(cellCounts.getCelltypes());
		return(new DeconvolutionResult(cellCounts.getCelltypes(), qtlName, pvalues, dummyModel, 0, 1));
	}

	public static double maximumLikelihoodEstimator() throws IOException, IllegalAccessException {
		/**
		 * Calculate the p-values of the betas by MLE estimation of parameters
		 */

		throw new NotImplementedException("mle not implemented yet");
	}

	private static double anova(double sumOfSquaresModelA, double sumOfSquaresModelB, int degreesOfFreedomA,
			int degreesOfFreedomB, Boolean no_intercept) {
		/**
		 * Compare and return the p-value of two linear models being
		 * significantly different
		 *
		 * From Joris Meys: http://stackoverflow.com/a/35458157/651779 1.
		 * calculate MSE for the largest model by dividing the Residual Sum of
		 * Squares (RSS) by the degrees of freedom (df) 2. calculate the
		 * MSEdifference by substracting the RSS of both models (result is
		 * "Sum of Sq." in the R table), substracting the df for both models
		 * (result is "Df" in the R table), and divide these numbers. 3. Divide
		 * 2 by 1 and you have the F value 4. calculate the p-value using the F
		 * value in 3 and for df the df-difference in the numerator and df of
		 * the largest model in the denominator. For more info:
		 * http://www.bodowinter.com/tutorial/bw_anova_general.pdf
		 * 
		 * @param sumOfSquaresModelA A vector with the genotype of all samples
		 * for *one* eQTL-gene pair
		 * 
		 * @param sumOfSquaresModelB A vector with the expression levels of all
		 * samples for *one* eQTL-gene pair
		 * 
		 * @param degreesOfFreedomA A 2D list with for all samples the different
		 * cell counts
		 * 
		 * @param degreesOfFreedomB A 2D list with for all samples the different
		 * cell counts
		 * 
		 * @return The p-value result from comparing two linear models with the
		 * the Anova test
		 */
		if (no_intercept) {
			// removing the intercept will give another degree of freedom
			degreesOfFreedomA++;
			degreesOfFreedomB++;
		}
		// Within-group Variance
		double meanSquareError = sumOfSquaresModelA / degreesOfFreedomA;

		int degreesOfFreedomDifference = Math.abs(degreesOfFreedomB - degreesOfFreedomA);
		// Between-group Variance
		// 234111286.801326
		double meanSquareErrorDiff = Math.abs((sumOfSquaresModelA - sumOfSquaresModelB) / (degreesOfFreedomDifference));

		/**
		 * F = Between-group Variance / Within-group Variance <- high value if
		 * variance between the models is high, and variance within the models
		 * is low
		 **/
		if(meanSquareError == 0){
			meanSquareError += 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001;
		}
		double Fval = meanSquareErrorDiff / meanSquareError;

		/***
		 * Make an F distribution with degrees of freedom as parameter. If full
		 * model and ctModel have the same number of samples, difference in df
		 * is 1 and degreesOfFreedomB are all the terms of the ctModel (so neut%
		 * + eos% + ... + neut% * GT + eos% * GT With 4 cell types and 1891
		 * samples the dfs are 1883 and 1884, giving the below distribution
		 * http://keisan.casio.com/exec/system/1180573186
		 **/
		FDistribution Fdist = new FDistribution(degreesOfFreedomDifference, degreesOfFreedomB);
		/*** Calculate 1 - the probability of observing a lower Fvalue **/
		double pval = 1 - Fdist.cumulative(Fval);
		return (pval);
	}

	public static DeconvolutionResult deconvolution(Qtl qtl) throws RuntimeException, IllegalAccessException, 
	NotEnoughGenotypesException, IOException, 
	NotEnoughSamplesPerGenotypeException {
		return deconvolution(qtl.getExpressionVector(), qtl.getGenotypeVector(), qtl.getQtlName());
	}

	/**
	 * Make the linear regression models and then do an Anova of the sum of
	 * squares
	 * 
	 * Full model: Exp ~ celltype_1 + celltype_2 + ... + celltype_n +
	 * celltype_1:Gt + celltype_2:Gt + ... + celltype_n:Gt <- without
	 * intercept
	 * 
	 * Compare with anova to Exp ~ celltype_1 + celltype_2 + celtype_n +
	 * celltype_1:Gt + celltype_2:Gt + .. + celltype_n-1 <- without
	 * intercept Exp ~ celltype_1 + celltype_2 + celtype_n + celltype_1:Gt +
	 * .. + celltype_n <- without intercept Exp ~ celltype_1 + celltype_2 +
	 * celtype_n + celltype_2:Gt + .. + celltype_n <- without intercept
	 *
	 * 
	 * @param expression A vector with the expression value per sample
	 * 
	 * @param genotypes A vector with the expression levels of all
	 * samples for *one* eQTL-gene pair. This should include qtl names as in first column, and sample names in first row
	 * 
	 * @param qtlName Name of the QTL (usaully snp name + gene name)
	 * 
	 * @return A list with for each celltype a p-value for the celltype
	 * specific eQTL for one eQTL
	 */
	private static DeconvolutionResult deconvolution(double[] expression, double[] genotypes, String qtlName) throws RuntimeException, NotEnoughGenotypesException, IllegalAccessException, 
	IOException, NotEnoughSamplesPerGenotypeException {


		/** 
		 * If roundDosage option is selected on the command line, round of the dosage to closest integer -> 0.49 = 0, 0.51 = 1, 1.51 = 2. 
		 * If minimumSamplesPerGenotype is selected on the command line, check for current QTL if for each dosage (in case they are not round
		 * the dosages are binned in same way as with roundDosage option) there are at least <minimumSamplesPerGenotype> samples that have it.
		 */

		if (commandLineOptions.getRoundDosage() || commandLineOptions.getMinimumSamplesPerGenotype() > 0 || commandLineOptions.getAllDosages()) {
			int dosage_ref = 0;
			int dosage_heterozygote = 0;
			int dosage_alt = 0;
			for (int i = 0; i < genotypes.length; i++) {
				if (commandLineOptions.getRoundDosage()){
					genotypes[i] = Math.round(genotypes[i]);
				}
				if (commandLineOptions.getMinimumSamplesPerGenotype() > 0 || commandLineOptions.getAllDosages()){
					if(genotypes[i] < 0){
						throw new RuntimeException("Genotype dosage can not be negative, check your dosage input file");
					}
					if(genotypes[i] < 0.5){
						dosage_ref++;
					}
					else if(genotypes[i] < 1.5){
						dosage_heterozygote++;
					}
					else if (genotypes[i] <= 2){
						dosage_alt++;
					}
					else{
						throw new RuntimeException("Genotype dosage can not be larger than 2, check your dosage input file");
					}
				}
			}
			if(commandLineOptions.getAllDosages()){
				// check that all dosages have at least one sample
				if(dosage_ref == 0 || dosage_heterozygote == 0 || dosage_alt == 0){
					throw new NotEnoughGenotypesException("Not all dosages present for this eQTL");
				}
			}
			if(commandLineOptions.getMinimumSamplesPerGenotype() > 0){
				// Check that each genotype has enough samples (AA >= minimum_samples_per_genotype, AB >= minimum_samples_per_genotype, BB >= minimum_samples_per_genotype)
				if(!(dosage_ref >= commandLineOptions.getMinimumSamplesPerGenotype() && dosage_heterozygote >= commandLineOptions.getMinimumSamplesPerGenotype() && dosage_alt >= commandLineOptions.getMinimumSamplesPerGenotype())){
					throw new NotEnoughSamplesPerGenotypeException("Not enough samples for each genotype");
				}}
		}

		InteractionModelCollection interactionModelCollection = new InteractionModelCollection();
		interactionModelCollection.setQtlName(qtlName);
		interactionModelCollection.setCelltypes(cellCounts.getCelltypes());
		interactionModelCollection.setGenotypes(genotypes);
		interactionModelCollection.setExpressionValues(expression);


		/**
		 * For each cell type model, e.g. ctModel 1 -> y = neut% + mono% + neut%:GT; ctModel 2 -> y = neut% + mono% + mono%:GT, one for each cell type, 
		 * where the interaction term (e.g mono%:GT) of the celltype:genotype to test is removed, calculate and save the observations in an observation vector
		 * where the observation vector for the example ctModel 1 is
		 *  
		 * 		celltypeModel = [[sample1_neut%, sample1_mono%, sample1_neut%*sample1_genotype], [sample2_neut%, sample2_mono%, sample2_neut%*sample2_genotype]]
		 *  
		 * with for each sample a cellcount percentage for each cell type and the genotype of the QTL that is being testetd. 
		 * 
		 * Using this observation vector calculate the sum of squares and test with Anova if it is significantly different from the sum of squares of the full model. 
		 * Here the full model includes all interaction terms of the cell type models, e.g. fullModel -> y = neut% + mono% + neut%:GT + mono%:GT so the observation vector
		 * 
		 * 		fullModel = [[sample1_neut%, sample1_mono%, sample1_neut%*sample1_genotype, sample1_mono%*sample1_genotype], [sample2_neut%, ..., etc]]
		 * 
		 */
		createObservedValueMatricesFullModel(interactionModelCollection);
		findBestFullModel(interactionModelCollection);
		createObservedValueMatricesCtModels(interactionModelCollection, interactionModelCollection.getBestFullModel().getGenotypeOrder());

		calculateDeconvolutionPvalue(interactionModelCollection);

		double wholeBloodQTL = 0;
		double wholeBloodQTLpvalue = 0;
		if(commandLineOptions.getWholeBloodQTL()){
			// if true calculate spearman correlation between genotypes and expression values (i.e. whole blood eQTL)
			wholeBloodQTL = new SpearmansCorrelation().correlation(interactionModelCollection.getGenotypes(), interactionModelCollection.getExpessionValues());
			wholeBloodQTLpvalue = Statistics.calculateSpearmanTwoTailedPvalue(wholeBloodQTL, interactionModelCollection.getSampleSize());
		}
		interactionModelCollection.emptyExpressionValues();
		interactionModelCollection.emptyGenotypes();
		DeconvolutionResult deconResult =  new DeconvolutionResult();

		deconResult = new DeconvolutionResult(interactionModelCollection, wholeBloodQTL, wholeBloodQTLpvalue);

		return deconResult;
	}

/*
 * Go through all full models, calculate the regression statistics and select the model with the highest R2 as the new full model
 */
private static void findBestFullModel(InteractionModelCollection interactionModelCollection) throws IllegalAccessException, IOException{
	// set to -1 so that first loop can be initialized
	double sumOfSquares = -1;
	for (String modelName : interactionModelCollection.getFullModelNames()){
		InteractionModel fullModel = interactionModelCollection.getInteractionModel(modelName);
		if(commandLineOptions.getUseNNLS()){
			fullModel.calculateSumOfSquaresNNLS(interactionModelCollection.getExpessionValues());
		}
		else{
			fullModel.calculateSumOfSquaresOLS(interactionModelCollection.getExpessionValues(), true);
		}
		if (sumOfSquares == -1){
			sumOfSquares = fullModel.getSumOfSquares();
			interactionModelCollection.setBestFullModel(fullModel.getModelName());
		}
		if (fullModel.getSumOfSquares() <= sumOfSquares){
			sumOfSquares = fullModel.getSumOfSquares();
			interactionModelCollection.setBestFullModel(fullModel.getModelName());
			fullModel.emptyObservedValues();
		}
		else{
			interactionModelCollection.removeInteractionModel(fullModel.getModelName());
		}
	}
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
	private static void createObservedValueMatricesFullModel(InteractionModelCollection interactionModelCollection) 
			throws IllegalAccessException{
		int numberOfTerms = cellCounts.getNumberOfCelltypes() * 2;
		ArrayList<String> binaryPermutations = new ArrayList<String>();
		if(commandLineOptions.getUseNNLS()){
			// if we use NNLS we need to see what genotype configuration works best, so get configuration permutations for all celltypes
			binaryPermutations = Utils.binaryPermutations("",interactionModelCollection.getCelltypes().size(), new ArrayList<String>());
		}else{
			// if we use OLS we just use default genotype orientation (all 0's)
			binaryPermutations.add(String.join("", Collections.nCopies(interactionModelCollection.getCelltypes().size(), "0")));
		}
		
		// Have to test which genotype combination is the best, so 2**number of celltype loops
		for (String binaryPermutation : binaryPermutations){
			// things neded for fullModel defined outside of loop because every celltype model (ctModel) has to be compared to it
			InteractionModel fullModel = new InteractionModel(cellCounts.getNumberOfSamples(), 
															  numberOfTerms);
			fullModel.setModelName(String.format("fullModel_%s",binaryPermutation));
			interactionModelCollection.addInteractionModel(fullModel, String.format("fullModel_%s",binaryPermutation), true);
			// number of terms + 1 because for full model all cell types are included
			interactionModelCollection.addInteractionModel(fullModel, fullModel.getModelName()); 

			for (int sampleIndex = 0; sampleIndex <= cellCounts.getNumberOfSamples()-1; sampleIndex++) {
				for (int celltypeIndex = 0; celltypeIndex < cellCounts.getNumberOfCelltypes(); celltypeIndex++) {
					double celltype_perc = cellCounts.getCellcountPercentages()[sampleIndex][celltypeIndex];
					// if i (cell type index) is the same as m (model index), don't add the interaction term of celltype:GT
					fullModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
					try {
						if(sampleIndex == 0){
							/** save the index of the variables related to current celltype so that this can be used later to calculate
							 * Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
							 * but for ctModel this is easiest method
							 */
							int[] index = new int[] {celltypeIndex, cellCounts.getNumberOfCelltypes() + celltypeIndex};
							fullModel.addCelltypeVariablesIndex(index);
								// add the celltype name at position i so that it gets in front of the celltype:GT
								fullModel.addIndependentVariableName(celltypeIndex, cellCounts.getCelltypes().get(celltypeIndex));
								fullModel.addIndependentVariableName(cellCounts.getCelltypes().get(celltypeIndex)+":GT");
								
						}
						// Have permutation of (2**number of celltypes) as binary ( so 00, 10, 01, 11 ), when 0 do normal genotype, 1 do swapped genotype
						double[] genotypes;
						fullModel.addGenotypeOrder(interactionModelCollection.getCelltypes().get(celltypeIndex), 
								binaryPermutation.charAt(celltypeIndex));
						// Use the binary string permutation to decide if the genotype should be swapped or not
						if(binaryPermutation.charAt(celltypeIndex) == '0'){
							genotypes = interactionModelCollection.getGenotypes();
						} else{
							genotypes = interactionModelCollection.getSwappedGenotypes();
						}
						fullModel.addObservedValue(celltype_perc * genotypes[sampleIndex], 
												   sampleIndex, cellCounts.getNumberOfCelltypes() + celltypeIndex);					
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
	private static void createObservedValueMatricesCtModels(InteractionModelCollection interactionModelCollection, HashMap<String, Character> genotypeOrder) 
			throws IllegalAccessException{
		int genotypeCounter = cellCounts.getNumberOfCelltypes();
		// -1 because one interaction term is removed
		int numberOfTerms = (cellCounts.getNumberOfCelltypes() * 2) - 1;

		// m = model, there are equally many models as celltypes, the fullModel gets made during the first iteration
		for (int modelIndex = 0; modelIndex < cellCounts.getNumberOfCelltypes(); modelIndex++) {
			InteractionModel ctModel = new InteractionModel(cellCounts.getNumberOfSamples(), numberOfTerms);

			// fullModel will be done in first loop as all models have to be compared to it

			// calculate p-value and save it, with other information, in a ctModel object. Then, add it to a list of these models to return as decon results
			ctModel.setModelName(interactionModelCollection.getCelltypes().get(modelIndex));
			interactionModelCollection.addInteractionModel(ctModel,ctModel.getModelName());

			for (int sampleIndex = 0; sampleIndex <= cellCounts.getNumberOfSamples()-1; sampleIndex++) {
				for (int celltypeIndex = 0; celltypeIndex < cellCounts.getNumberOfCelltypes(); celltypeIndex++) {
					// There is one fullModel including all celltypes add values for celltypePerc and interaction term of
					// celltypePerc * genotypePerc so that you get [[0.3, 0.6], [0.4, 0.8], [0.2, 0.4], [0.1, 0.2]]
					// where numberOfSamples = 1 and numberOfCellTypes = 4 with celltypePerc = 0.3, 0.4, 0.2, and 0.1 and genotype = 2
					// for each cell type is 1 model, celltype% * genotype without 1 celltype.
					// j+1 because j==0 is header
					double celltype_perc = cellCounts.getCellcountPercentages()[sampleIndex][celltypeIndex];
					ctModel.addObservedValue(celltype_perc, sampleIndex, celltypeIndex);
					if(sampleIndex == 0){
						// add the celltype name at position i so that it gets in front of the celltype:GT, but once
						try{
							ctModel.addIndependentVariableName(celltypeIndex, cellCounts.getCelltypes().get(celltypeIndex));
						}
						catch(NullPointerException e){
							DeconvolutionLogger.log.info(String.format("Nullpoint exception with celltype %s", celltypeIndex));
							throw e;
						}
					}

					// if i (cell type index) is the same as m (model index), don't add the interaction term of celltype:GT
					if (celltypeIndex != modelIndex) {
						try {
							// Only add IndependentVariableName once per QTL (j==0)
							if(sampleIndex == 0){

								// Add the interaction term of celltype:genotype
								ctModel.addIndependentVariableName(cellCounts.getCelltypes().get(celltypeIndex)+":GT");
								// save the index of the variables related to current celltype so that this can be used later to calculate
								// Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
								// but for ctModel this is easiest method
								int[] index = new int[] {celltypeIndex, cellCounts.getNumberOfCelltypes()-1+celltypeIndex};
								ctModel.addCelltypeVariablesIndex(index);
								// add the celltype name. This could be done with less code by getting it from IndependentVariableName, but this way 
								// it is explicit. Don't know if better.
							}
							try{
								double genotype = 0;
								char genotypeOrderAtCelltype = genotypeOrder.get(cellCounts.getCelltypes().get(celltypeIndex));
								if(genotypeOrderAtCelltype == '0'){
									genotype = interactionModelCollection.getGenotypes()[sampleIndex];
								}
								else if(genotypeOrderAtCelltype == '1'){
									genotype = interactionModelCollection.getSwappedGenotypes()[sampleIndex];
								}
								else{
									throw new RuntimeException(String.format("Genotype order should be 0 or 1, was: %s", genotypeOrderAtCelltype));
								}
								ctModel.addObservedValue(celltype_perc * genotype, sampleIndex, genotypeCounter);
							}
							catch(NullPointerException e){
								DeconvolutionLogger.log.info(String.format("Nullpoint exception with genotype %s", sampleIndex));
								throw e;
							}
						} catch (ArrayIndexOutOfBoundsException error) {
							throw new RuntimeException(
									"The counts file and expression and/or genotype file do not have equal number of samples or QTLs",
									error);
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
				genotypeCounter = cellCounts.getNumberOfCelltypes();
			}
		ctModel.setModelLength();	
		}
	}
	
	
	/**
	 * get pvalue for each ctmodel
	 * 
	 * @param ctModel InteractionModel object for saving the results
	 * @param m The current model that is being evaluated (for each celltype 1 model)
	 * @param fullModel InteractionModel object that contains information on the fullModel (such as expression values)
	 * @param qtlName Name of the current qtl being calculated
	 */
	private static void calculateDeconvolutionPvalue(InteractionModelCollection interactionModelCollection) 
			throws IllegalAccessException, IOException {
		for (int modelIndex = 0; modelIndex < cellCounts.getNumberOfCelltypes(); modelIndex++) {
			String modelName = interactionModelCollection.getCelltypes().get(modelIndex);
			InteractionModel fullModel = interactionModelCollection.getBestFullModel();
			InteractionModel ctModel = interactionModelCollection.getInteractionModel(modelName);

			if (commandLineOptions.getUseNNLS()){
				ctModel.calculateSumOfSquaresNNLS(interactionModelCollection.getExpessionValues());
			}
			else{
				ctModel.calculateSumOfSquaresOLS(interactionModelCollection.getExpessionValues(), false);
			}
			ctModel.emptyObservedValues();
			
			int expressionLength = interactionModelCollection.getExpessionValues().length;
			if (expressionLength != fullModel.getModelLength()) {
				throw new RuntimeException("expression vector and fullModel have different number of samples.\nexpression: "
						+ expressionLength + "\nfullModel: " + fullModel.getModelLength());
			}

			//double pval = anova(fullModel.getSumOfSquares(), ctModel.getSumOfSquares(), fullModel.getDegreesOfFreedom(), ctModel.getDegreesOfFreedom(), true);
			double pval = anova(fullModel.getSumOfSquares(), ctModel.getSumOfSquares(), fullModel.getDegreesOfFreedom(), ctModel.getDegreesOfFreedom(), true);
			ctModel.setPvalue(pval);
			ctModel.emptyObservedValues();
			interactionModelCollection.setPvalue(pval,modelName);
		}
	}
}
