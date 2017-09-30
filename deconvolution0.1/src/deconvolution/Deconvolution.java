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
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import JSci.maths.statistics.FDistribution;

public class Deconvolution {
	private static CommandLineOptions commandLineOptions = new CommandLineOptions(); 
	private static CellCount cellCounts;
	private static List<String> filteredQTLsOutput = new ArrayList<String>();
	private static int QTLsFiltered = 0;
	// things neded for fullModel defined here because for optimizing pvalue it needs to be available in 2 functions. 
	// this could probably be done a lot neater
	private static double sumOfSquaresFullModel = 0;
	private static int degreesOfFreedomFullModel = 0;
	private static int fullModelLength = 0;
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
				if (whileIndex % 5000 == 0) {
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
				catch(NonNegativeConstraintViolatedException e){
					if(!commandLineOptions.getFilterSamples()){
						deconvolutionResults.add(setPvaluesNA(qtlName));
					}
					else{
						QTLsFiltered++;
						// TODO: Use MLE
					}
					filteredQTLsOutput.add(qtlName+"\tFailed non-negative linear model parameter constraint");
				}
				catch(NotEnoughSamplesPerGenotypeException e){
					if(!commandLineOptions.getFilterSamples()){
						deconvolutionResults.add(setPvaluesNA(qtlName));
					}
					else{
						QTLsFiltered++;
					}
					filteredQTLsOutput.add(qtlName+"\tNot enough samples per genotype");
					// TODO: Use MLE
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
		for(int i = 1; i < cellCounts.getNumberOfCelltypes()*2 + 1; i++){
			try{
				header += "\tBeta" + Integer.toString(i) +"_"+deconvolutionResults.get(0).getFullModel().getIndependentVariableNames().get(i-1);
			}
			catch(IndexOutOfBoundsException e){
				DeconvolutionLogger.log.info(String.format("DeconvolutionResult index error with beta %d", i));
				throw e;
			}
		}

		for(int i = 1; i < cellCounts.getNumberOfCelltypes()*2 + 1; i++){
			header += "\tBetaStandardError" + Integer.toString(i); 
		}

		if(commandLineOptions.getWholeBloodQTL()){
			header += "\tSpearman correlation expression~GT\tSpearman correlation p-value";
		}
		output.add(header);
		for(DeconvolutionResult deconvolutionResult : deconvolutionResults){
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
				results += "\t"+Utils.listToTabSeparatedString(deconvolutionResult.getFullModel().getEstimateRegressionParameters());
				results += "\t"+Utils.listToTabSeparatedString(deconvolutionResult.getFullModel().getEstimateRegressionParametersStandardErrors());
			}catch (java.lang.IllegalAccessException e){
				// if -m is set not all deconvolution resuts will have a full model. If not, set betas to NA
				String str = "\tNA";
				// celltype.size()*4 because for each celltype are 2 betas + 2 standard error of betas
				results += StringUtils.repeat(str, cellCounts.getNumberOfCelltypes()*4);
			}


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
		dummyModel.setCelltypes(cellCounts.getCelltypes());
		dummyModel.setAlltIndependentVariableNames();
		return(new DeconvolutionResult(cellCounts.getCelltypes(), qtlName, pvalues, dummyModel, 0, 1));
	}

	/**
	 * Calculate the sum of squares, using Non-Negative Linear Regression, given a y expression vector with y ~
	 * model. This calculates the Beta parameters three times, for samples with GT = 0, GT = 1, GT = 2
	 * instead of using the interaction terms of the actual genotype.
	 * If no_intercept == true, remove the intercept (equivalent to y * ~ model -1 in R)
	 * 
	 * @return	The parameters per genotype dosage
	 * @param model An InteractionModel object including the y vector expression values and ObservedValues (model)
	 * Such that
	 * test_trait ~ geno_A + lymph% + geno_A:geno_B it can be for one QTL
	 * [[2, 43.4, 86.8], [2, 40.3, 80.6]], for another QTL [[0, 46.7, 0],
	 * [0, 51.5, 0] [0, 48.7, 0]] 
	 */
	public static void calculateSumOfSquaresNNLS(LeastSquareModel model, Boolean plotBetaTimesVariables) throws IOException, IllegalAccessException {
		// OLS = Ordinary Least Squares
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		// if GetIntercept is false, remove the intercept (Beta1) from the linear model
		regression.setNoIntercept(model.getNoIntercept());
		regression.newSampleData(model.getExpessionValuesGt0(), model.getCellCountsGt0());
		regression.estimateRegressionParameters();
		throw new NotImplementedException("NNLS not implemneted");
	}

	/**
	 * Calculate the sum of squares, using Ordinary Linear Regression, given a y expression vector with y ~
	 * model. If no_intercept == true, remove the intercept (equivalent to y
	 * ~ model -1 in R)
	 * 
	 * @param model An InteractionModel object including the y vector expression values and ObservedValues (model)
	 * Such that
	 * test_trait ~ geno_A + lymph% + geno_A:geno_B it can be for one QTL
	 * [[2, 43.4, 86.8], [2, 40.3, 80.6]], for another QTL [[0, 46.7, 0],
	 * [0, 51.5, 0] [0, 48.7, 0]]
	 * 
	 * @return The sum of squares value from running linear regression with
	 * y ~ model
	 */
	private static OLSMultipleLinearRegression multipleLinearRegression(InteractionModel model) throws IOException, IllegalAccessException {
		// OLS = Ordinary Least Squares
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		// if GetIntercept is false, remove the intercept (Beta1) from the linear model
		regression.setNoIntercept(model.getNoIntercept());
		try{
			regression.newSampleData(model.getExpessionValues(), model.getObservedValues());
		}
		catch (DimensionMismatchException e){
			DeconvolutionLogger.log.info(String.format("Length of expression and and genotype data not the same\nexpression length: %d\nobserved values length: %d\n", 
					model.getExpessionValues().length, model.getObservedValues().length));
			throw(e);
		}

		return (regression);
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
	NonNegativeConstraintViolatedException, 
	NotEnoughSamplesPerGenotypeException {
		return deconvolution(qtl.getExpressionVector(), qtl.getGenotypeVector(), qtl.getQtlName());
	}

	//TODO: update documentation on parameters
	private static DeconvolutionResult deconvolution(double[] expression, double[] genotypes, String qtlName) throws RuntimeException, NotEnoughGenotypesException, IllegalAccessException, 
	IOException, NonNegativeConstraintViolatedException, 
	NotEnoughSamplesPerGenotypeException {
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
		 * @param expressionVector A vector with the genotype of all samples for
		 * *one* eQTL-gene pair. This should include qtl names as in first column, and sample names in first row
		 * 
		 * @param genotypeVector A vector with the expression levels of all
		 * samples for *one* eQTL-gene pair. This should include qtl names as in first column, and sample names in first row
		 * 
		 * @param minimum_samples_per_genotype Minimum number of samples per genotype need to count it
		 * 
		 * @return A list with for each celltype a p-value for the celltype
		 * specific eQTL for one eQTL
		 */

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

		List<Double> pvalues = new ArrayList<Double>();
		// things neded for fullModel defined outside of loop because every celltype model (ctModel) has to be compared to it
		InteractionModel fullModel = new InteractionModel();
		fullModel.setModelName("full model");
		fullModel.setQtlName(qtlName);
		fullModel.setCelltypes(cellCounts.getCelltypes());
		List<InteractionModel> ctModels = new ArrayList<InteractionModel>(); 
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
		// m = model, there are equally many models as celltypes, the fullModel gets made during the first iteration
		for (int modelIndex = 0; modelIndex < cellCounts.getNumberOfCelltypes(); modelIndex++) {
			InteractionModel ctModel = new InteractionModel();
			ctModel.setCelltypes(cellCounts.getCelltypes());
			ctModel.setExpressionValues(expression);
			ctModel.setGenotypes(genotypes);
			ctModel.setQtlName(qtlName);
			int numberOfTerms = (cellCounts.getNumberOfCelltypes() * 2) - 1;
			ctModel.InitializeObservedValue(cellCounts.getNumberOfSamples(), numberOfTerms);
			// fullModel will be done in first loop as all models have to be compared to it
			if(modelIndex == 0){
				fullModel.setExpressionValues(expression);
				fullModel.setGenotypes(genotypes);
				// number of terms + 1 because for full model all cell types are included
				fullModel.InitializeObservedValue(cellCounts.getNumberOfSamples(), numberOfTerms+1);
			}
			// calculate p-value and save it, with other information, in a ctModel object. Then, add it to a list of these models to return as decon results
			ctModel = calculateDeconvolutionPvalue(ctModel, modelIndex, fullModel, qtlName);
			pvalues.add(ctModel.getPvalue());
			ctModel.emptyExpressionValues();
			ctModel.emptyGenotypes();
			ctModel.emptyObservedValues();
			ctModels.add(ctModel);
		}
		/**
		 * TODO: implement below function
		 * calculateSumOfSquaresNNLS(leastSquareModel, true);
		 */
		double wholeBloodQTL = 0;
		double wholeBloodQTLpvalue = 0;
		if(commandLineOptions.getWholeBloodQTL()){
			// if true calculate spearman correlation between genotypes and expression values (i.e. whole blood eQTL)
			wholeBloodQTL = new SpearmansCorrelation().correlation(fullModel.getGenotypes(), fullModel.getExpessionValues());
			wholeBloodQTLpvalue = Statistics.calculateSpearmanTwoTailedPvalue(wholeBloodQTL, fullModel.getSampleSize());
		}
		fullModel.emptyExpressionValues();
		fullModel.emptyObservedValues();
		fullModel.emptyGenotypes();
		DeconvolutionResult deconResult =  new DeconvolutionResult();
		deconResult = new DeconvolutionResult(cellCounts.getCelltypes(), qtlName, pvalues, fullModel, ctModels, wholeBloodQTL, wholeBloodQTLpvalue);

		return deconResult;
	}

	private static InteractionModel calculateDeconvolutionPvalue(InteractionModel ctModel, int modelIndex,
			InteractionModel fullModel, String qtlName ) throws IllegalAccessException, IOException, NonNegativeConstraintViolatedException{
		/**
		 * get pvalue for each ctmodel
		 * 
		 * @param ctModel InteractionModel object for saving the results
		 * @param m The current model that is being evaluated (for each celltype 1 model)
		 * @param fullModel InteractionModel object that contains information on the fullModel (such as expression values)
		 * @param qtlName Name of the current qtl being calculated
		 */
		/**
		 * Below nested for loops make the matrices that are necesarry to run the linear model and put them in a model object
		 */
		int genotypeCounter = cellCounts.getNumberOfCelltypes();
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
							ctModel.addObservedValue(celltype_perc * ctModel.getGenotypes()[sampleIndex], sampleIndex, genotypeCounter);
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
				if (modelIndex == 0){
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
						fullModel.addObservedValue(celltype_perc * ctModel.getGenotypes()[sampleIndex], sampleIndex, cellCounts.getNumberOfCelltypes() + celltypeIndex);					
					} catch (ArrayIndexOutOfBoundsException error) {
						throw new RuntimeException(
								"The counts file and expression and/or genotype file do not have equal number of samples or QTLs",
								error);
					}
				}
			}
			// because 1 of numberOfCelltypes + i needs to be skipped,
			// keeping it tracked with separate value is easier
			genotypeCounter = cellCounts.getNumberOfCelltypes();
		}
		Boolean noIntercept = true;
		if(modelIndex == 0){
			// only need to set data of fullModel once, reused every loop of m
			fullModel.setNoIntercept(noIntercept);
			OLSMultipleLinearRegression regression = multipleLinearRegression(fullModel);
			//for (int i = 0; i < estimatedRegressionParameters.length; i++){
			//	DeconvolutionLogger.log.info(String.format("beta: %f\terror: %f\n", estimatedRegressionParameters[i], estimateRegressionParametersStandardErrors[i]));
			//}
			sumOfSquaresFullModel = regression.calculateResidualSumOfSquares();
			degreesOfFreedomFullModel = ctModel.getExpessionValues().length - (fullModel.getObservedValues()[0].length + 1);
			fullModelLength = fullModel.getObservedValues().length;
			double[] estimatedRegressionParameters = regression.estimateRegressionParameters();
			double[] estimateRegressionParametersStandardErrors = regression.estimateRegressionParametersStandardErrors();
			fullModel.setEstimateRegressionParameters(estimatedRegressionParameters);
			fullModel.setEstimateRegressionParametersStandardErrors(estimateRegressionParametersStandardErrors);
		}
		ctModel.setNoIntercept(noIntercept);
		/*** SUM OF SQUARES - CELLTYPE MODEL **/
		ctModel.setModelName("ctModel_"+ Integer.toString(modelIndex));
		OLSMultipleLinearRegression regression = multipleLinearRegression(ctModel);
		double sumOfSquaresCtModel = regression.calculateResidualSumOfSquares();
		int degreesOfFreedomCtModel = ctModel.getExpessionValues().length - (ctModel.getObservedValues()[0].length + 1);

		int expressionLength = ctModel.getExpessionValues().length;

		if (expressionLength != fullModelLength) {
			throw new RuntimeException("expression vector and fullModel have different number of samples.\nexpression: "
					+ expressionLength + "\nfullModel: " + fullModelLength);
		}

		// ANOVA compare full model to celltype model
		/**		if(ctModel.GetQtlName().contains("ENSG00000262539")){
			System.out.println(sumOfSquaresFullModel);
			System.out.println(sumOfSquaresCtModel);
			System.out.println(degreesOfFreedomFullModel);
			System.out.println(degreesOfFreedomCtModel);
		}
		 */			
		double pval = anova(sumOfSquaresFullModel, sumOfSquaresCtModel, degreesOfFreedomFullModel, degreesOfFreedomCtModel, true);
		ctModel.setPvalue(pval);
		ctModel.emptyExpressionValues();
		ctModel.emptyObservedValues();
		return ctModel;
	}
}
