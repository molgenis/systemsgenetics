package deconvolution;

import deconvolution.Qtl;
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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
			Set<String> expressionSamplesSet1 = new HashSet<String>(expressionData.getSampleNames());
			// use 2 times expresionSampleSet because do inplace replacement, so when removing from genotypeSamples set need a new expressionSamples set
			Set<String> expressionSamplesSet2 = new HashSet<String>(expressionData.getSampleNames());
			Set<String> genotypeSamplesSet = new HashSet<String>(genotypeData.getSampleNames());
			expressionSamplesSet1.removeAll(genotypeSamplesSet);
			genotypeSamplesSet.removeAll(expressionSamplesSet2);
			int numberSamplesMissingExpression = expressionSamplesSet1.size();
			int numbeSamplesrMissingGenotype = genotypeSamplesSet.size();
			String expressionSamples = Arrays.toString(expressionSamplesSet1.toArray());
			String genotypeSamples = Arrays.toString(genotypeSamplesSet.toArray());
			throw new RuntimeException(String.format("Samplenames not the same in expression and genotype file.\nexpression samples not in genotypes (%d): %s\ngenotype samples not in expression (%d): %s\n",
					numberSamplesMissingExpression, expressionSamples, numbeSamplesrMissingGenotype, genotypeSamples));
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
				if (whileIndex % 2500 == 0) {
					long completedIn = System.currentTimeMillis() - time;
					DeconvolutionLogger.log.info(String.format("Processed %d gene-SNP pairs - %s - skipped %d gene-SNP combinations", whileIndex, DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS"), skippedGenotypeGeneCombinations));
				}
				++whileIndex;
				String qtlName = gene+'_'+genotype;
				try{
					++QTLsTotal;
					try{
						double[] dosages = genotypeData.getGenotypes().get(genotype);
						if(dosages == null){
							throw new RuntimeException(String.format("SNP %s not in genotype file, is your snpsToTest file correct?", genotype));
						}
						DeconvolutionResult deconResult = deconvolution(geneExpressionLevels.get(gene), 
								dosages, qtlName);
						deconvolutionResults.add(deconResult);
					}
					catch(IllegalAccessException e){
						if(commandLineOptions.getSkipGenotypes()){
							++skippedGenotypeGeneCombinations;
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
						++QTLsFiltered;
					}
					filteredQTLsOutput.add(qtlName+"\tNot enough genotypes (e.g. AA and AB but no BB)");
				}
				catch(NotEnoughSamplesPerGenotypeException e){
					if(!commandLineOptions.getFilterSamples()){
						deconvolutionResults.add(setPvaluesNA(qtlName));
					}
					else{
						++QTLsFiltered;
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
	}

	/**
	 * Write the deconfolution results
	 * 
	 * @param deconvolutionResult The deconvolutionresult
	 */
	private static void writeDeconvolutionResults(List<DeconvolutionResult> deconvolutionResults) throws IllegalAccessException, IOException{
		String header = "\t"+Utils.listToTabSeparatedString(cellCounts.getAllCelltypes(), "_pvalue");
		
		if(!commandLineOptions.getUseBaseModel()){
			header += "\tAIC_fullModel";
		
			for(String celltype : cellCounts.getAllCelltypes()){
				header += "\tAIC_diff_"+celltype;
			}
		}
		DeconvolutionLogger.log.info("Getting decon result with full model info for writing the header");
		// celltypes.size()*2 because there are twice as many betas as celltypes (CC% & CC%:GT)
		InteractionModelCollection firstInteractionModelCollection = deconvolutionResults.get(0).getInteractionModelCollection();
		InteractionModel bestFullModelForHeaderOnly = firstInteractionModelCollection.getBestFullModel();
		if(commandLineOptions.getUseBaseModel()){
			header += "\tBeta1\tBeta2\tBeta3\tBeta4\tBeta5\tBeta6";
		}
		
		else{
			for(int i = 1; i < cellCounts.getNumberOfCelltypes()*2 + 1; ++i){
				header += "\tBeta" + Integer.toString(i) +"_"+bestFullModelForHeaderOnly.getIndependentVariableNames().get(i-1);
			}
		}
		for(String celltype : cellCounts.getAllCelltypes()){
			header += "\teffectDirectionDosage2_"+celltype;
		}
		//header += "\tgenotypeConfiguration";
		//for(String celltype : cellCounts.getAllCelltypes()){
		//	header += "\tgenotypeConfiguration_"+celltype;
		//}\
		
		if(commandLineOptions.getWholeBloodQTL()){
			header += "\tSpearman correlation expression~GT\tSpearman correlation p-value";
		}

		header += "\tStandardError";
		List<String> output = new ArrayList<String>();
		output.add(header);
		for(DeconvolutionResult deconvolutionResult : deconvolutionResults){
			InteractionModelCollection interationModelCollection = deconvolutionResult.getInteractionModelCollection();
			if(commandLineOptions.getOnlyOutputSignificant() && commandLineOptions.getFilterSamples()){
				if(Collections.min(deconvolutionResult.getPvalues()) > 0.05){
					++QTLsFiltered;
					filteredQTLsOutput.add(deconvolutionResult.getQtlName()+"\tNone of the celltypes had a significant p-value");
					continue;
				}
			}
			String results = "";
			results += deconvolutionResult.getQtlName()+"\t"+Utils.listToTabSeparatedString(deconvolutionResult.getPvalues());
			InteractionModel bestFullModel = null;
			if(!commandLineOptions.getUseBaseModel()){
				bestFullModel = interationModelCollection.getBestFullModel();
				double bestFullModelAIC = bestFullModel.getAIC();
				results += "\t"+bestFullModelAIC;
			}
			for(String celltype : cellCounts.getAllCelltypes()){
				//System.out.println(celltype);
				String modelName = deconvolutionResult.getInteractionModelCollection()
													.getCtModelSameGenotypeConfigurationAsBestFullModel(celltype, 
																										commandLineOptions.getUseBaseModel());
				if(!commandLineOptions.getUseBaseModel()){
					results += "\t"+interationModelCollection.getInteractionModel(modelName).getAICdelta();
				}
				if(commandLineOptions.getUseBaseModel()){
					InteractionModel bestFullModelOfCurrentCelltype = interationModelCollection.getBestFullModel(celltype);
					/*
					 *  when using base model, because the model is
					 *  	y ~ cc + (100-cc) + cc:GT + (100-cc):GT
					 *  the beta is always 3rd (index=2) term
					 * 
					 */
					results += "\t"+bestFullModelOfCurrentCelltype.getEstimateRegressionParameters()[2];
				}
			}

			if(!commandLineOptions.getUseBaseModel()){
				results += "\t"+Utils.listToTabSeparatedString(bestFullModel.getEstimateRegressionParameters());
			}

			// check what the genotype configuration is and the beta of the interaction term. 
			// If genotype configuration == 0 and beta == positive, dosage2 effect = positive
			// If genotype configuration == 1 and beta == negative, dosage2 effect = positive
			// else is negative
			int numberOfCelltypes = cellCounts.getNumberOfCelltypes();
			for(int i = 0; i < numberOfCelltypes; ++i){
				char genotypeConfiguration = 0;
				double estimatedRegressionParameter;

				if(commandLineOptions.getUseBaseModel()){
					bestFullModel = interationModelCollection.getBestFullModel(cellCounts.getCelltype(i));
					genotypeConfiguration = bestFullModel.getGenotypeConfiguration().charAt(0);
					estimatedRegressionParameter = bestFullModel.getEstimateRegressionParameters()[2];
				}
				else{
					estimatedRegressionParameter = bestFullModel.getEstimateRegressionParameters()[i+numberOfCelltypes];
					genotypeConfiguration = bestFullModel.getGenotypeConfiguration().charAt(i);
				}
				if (genotypeConfiguration == '0'){
					// add cellCounts.getNumberOfCelltypes() to get the regression parameter for the interaction term (first ones are indepent effect betas)
					if(estimatedRegressionParameter < 0){
						results += "\t-";			
					}
					else{
						results += "\t+";
					}
				}else if(genotypeConfiguration == '1'){
					if(estimatedRegressionParameter < 0){
						results += "\t+";			
					}
					else{
						results += "\t-";
					}
				}
				else{
					throw new RuntimeException(String.format("Genotype configuration should be 0 or 1, not %s", genotypeConfiguration));
				}

			}

			//results += "\t"+bestFullModel.getGenotypeConfiguration();
			//for(String celltype : cellCounts.getAllCelltypes()){
			//	InteractionModel bestCtModel = deconvolutionResult.getInteractionModelCollection().getBestCtModel(celltype); 
			//	results += "\t"+bestCtModel.getGenotypeConfiguration();
			//}
			if(commandLineOptions.getWholeBloodQTL()){
				results += "\t"+deconvolutionResult.getWholeBloodQTL();
				results += "\t"+deconvolutionResult.getWholeBloodQTLpvalue();
			}
			
			results += "\t"+bestFullModel.getEstimatedStandardError();
			output.add(results);	
		}

		Path file = Paths.get(outputFolder+commandLineOptions.getOutfile());
		Files.write(file, output, Charset.forName("UTF-8"));
		
		Boolean writePredictedExpression = commandLineOptions.getOutputPredictedExpression(); 
		if(writePredictedExpression){
			writePredictedExpression(deconvolutionResults);
		}
		DeconvolutionLogger.log.info(String.format("Deconvolution output written to %s", file.toAbsolutePath()));
		DeconvolutionLogger.log.info(String.format("Files with additional info in  %s", outputFolder));
	}

	/*
	 * Write the predicted expression values to a separate file
	 * 
	 * @param deconvolutionResult The deconvolutionresult
	 */
	private static void writePredictedExpression(List<DeconvolutionResult> deconvolutionResults) throws IOException, IllegalAccessException{
		DeconvolutionResult deconResult = deconvolutionResults.get(0);
		String header = "";
		for(String sampleName : deconResult.getInteractionModelCollection().getSampleNames()){
			// counts.get(0) is the sample name
			header += "\t"+sampleName;
			
		}
		
		List<String> output = new ArrayList<String>();
		output.add(header);
		
		for(DeconvolutionResult deconvolutionResult : deconvolutionResults){
			
			String results = "";
			InteractionModel bestFullModel = deconvolutionResult.getInteractionModelCollection().getBestFullModel();
			results += deconvolutionResult.getQtlName()+"\t"+Utils.listToTabSeparatedString(bestFullModel.getPredictedValues());
			output.add(results);	
		}

	
		
		Path file = Paths.get(outputFolder+"predictedExpressionLevels.txt");
		Files.write(file, output, Charset.forName("UTF-8"));
		DeconvolutionLogger.log.info(String.format("predicted expression written to %s", file.toAbsolutePath()));
	}
	
	/*
	 * Incase p-values have to be written as NA (e.g. when they should be filtered)
	 */
	private static DeconvolutionResult setPvaluesNA(String qtlName) throws IllegalAccessException{
		List<Double> pvalues = new ArrayList<Double>();
		for (int i = 0; i < cellCounts.getNumberOfCelltypes(); ++i){
			pvalues.add(333.0);
		}
		InteractionModel dummyModel = new InteractionModel();
		dummyModel.setModelName("dummy");
		dummyModel.setAlltIndependentVariableNames(cellCounts.getAllCelltypes());
		return(new DeconvolutionResult(cellCounts.getAllCelltypes(), qtlName, pvalues, dummyModel, 0, 1));
	}


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
	public static double anova(double sumOfSquaresModelA, double sumOfSquaresModelB, int degreesOfFreedomA,
			int degreesOfFreedomB, Boolean no_intercept) {
		if (no_intercept) {
			// removing the intercept will give another degree of freedom
			++degreesOfFreedomA;
			++degreesOfFreedomB;
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
			for (int i = 0; i < genotypes.length; ++i) {
				if (commandLineOptions.getRoundDosage()){
					genotypes[i] = Math.round(genotypes[i]);
				}
				if (commandLineOptions.getMinimumSamplesPerGenotype() > 0 || commandLineOptions.getAllDosages()){
					if(genotypes[i] < 0){
						throw new RuntimeException("Genotype dosage can not be negative, check your dosage input file");
					}
					if(genotypes[i] < 0.5){
						++dosage_ref;
					}
					else if(genotypes[i] < 1.5){
						++dosage_heterozygote;
					}
					else if (genotypes[i] <= 2){
						++dosage_alt;
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

		InteractionModelCollection interactionModelCollection = new InteractionModelCollection(cellCounts, 
																								commandLineOptions.getGenotypeConfigurationType(),
																								commandLineOptions.getUseBaseModel());
		interactionModelCollection.setQtlName(qtlName);
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
		interactionModelCollection.createObservedValueMatricesFullModel();
		interactionModelCollection.findBestFullModel(commandLineOptions.getUseBaseModel());		
		interactionModelCollection.createObservedValueMatricesCtModels();
		interactionModelCollection.findBestCtModel(commandLineOptions.getUseBaseModel());

		calculateDeconvolutionPvalue(interactionModelCollection);

		double wholeBloodQTL = 0;
		double wholeBloodQTLpvalue = 0;
		if(commandLineOptions.getWholeBloodQTL()){
			// if true calculate spearman correlation between genotypes and expression values (i.e. whole blood eQTL)
			wholeBloodQTL = new SpearmansCorrelation().correlation(interactionModelCollection.getGenotypes(), interactionModelCollection.getExpessionValues());
			wholeBloodQTLpvalue = Statistics.calculateSpearmanTwoTailedPvalue(wholeBloodQTL, cellCounts.getNumberOfSamples());
		}
		DeconvolutionResult deconResult =  new DeconvolutionResult();
		if(!commandLineOptions.getUseBaseModel()){
			interactionModelCollection.setAIC(commandLineOptions.getUseBaseModel());
		}
		interactionModelCollection.cleanUp(!commandLineOptions.getOutputPredictedExpression());
		deconResult = new DeconvolutionResult(interactionModelCollection, wholeBloodQTL, wholeBloodQTLpvalue);
		return deconResult;
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
		for (int modelIndex = 0; modelIndex < cellCounts.getNumberOfCelltypes(); ++modelIndex) {
			String celltypeName = cellCounts.getCelltype(modelIndex);
			InteractionModel fullModel;
			if(commandLineOptions.getUseBaseModel()){
				fullModel = interactionModelCollection.getBestFullModel(celltypeName);
			}
			else{
				fullModel = interactionModelCollection.getBestFullModel();
			}

			int expressionLength = interactionModelCollection.getExpessionValues().length;
			if (expressionLength != fullModel.getModelLength()) {
				throw new RuntimeException("expression vector and fullModel have different number of samples.\nexpression: "
						+ expressionLength + "\nfullModel: " + fullModel.getModelLength());
			}

			InteractionModel ctModel = interactionModelCollection.getBestCtModel(celltypeName);
			double pval = anova(fullModel.getSumOfSquares(), ctModel.getSumOfSquares(), 
								fullModel.getDegreesOfFreedom(),ctModel.getDegreesOfFreedom(), true);
			ctModel.setPvalue(pval);
			interactionModelCollection.setPvalue(pval, ctModel.getCelltypeName());
			interactionModelCollection.setPvalue(pval,interactionModelCollection
														.getBestCtModel(cellCounts.getCelltype(modelIndex)).getModelName());

		}
	}
}
