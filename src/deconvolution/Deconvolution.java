package deconvolution;

import deconvolution.Qtl;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang3.time.DurationFormatUtils;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import JSci.maths.statistics.FDistribution;

public class Deconvolution {
	private static CommandLineOptions commandLineOptions = new CommandLineOptions(); 
	public static void main(String[] args) throws Exception {
		/*
		 * Deconvolutes a set of QTLs given the expression levels, genotypes,
		 * and cell counts. Calculates the p-values for the deconvoluted QTLs
		 * and writes them to an outfile
		 * 
		 * @param args List of command line arguments
		 */

		// begin command line option parsing
		// for some variables need final variable because roundDosage is used in local scope when using parallel threads (http://stackoverflow.com/questions/25894509/problems-with-local-variable-scope-how-to-solve-it)
		commandLineOptions.parseCommandLine(args);

		// End of command line option parsing

		// the cell type names are the first row of cellcount file, extract for
		// later printing
		String cellCountFile = commandLineOptions.getCellcountFile();
		List<List<String>> cellcountTable = readTabDelimitedColumns(cellCountFile);
		String outfilePath = commandLineOptions.getOutfile();

		// output saves all the processed output to write it to a file later
		List<String> output = new ArrayList<String>();

		String expressionFile = commandLineOptions.getExpressionFile();
		String genotypeFile = commandLineOptions.getGenotypeFile();

		System.out.printf("\n");

		// make iterator of files so that we can loop over 2 at the same time
		LineIterator expressionIterator = FileUtils.lineIterator(new File(expressionFile), "UTF-8");
		LineIterator genotypeIterator = FileUtils.lineIterator(new File(genotypeFile), "UTF-8");
		// leastsquare model is for calculating the coefficients that are used for plotting

		// remove the headers (sample names)
		String sampleNames = expressionIterator.next();
		if (!sampleNames.equals(genotypeIterator.next())){
			throw new Exception("Samplenames not the same in expression and genotype file");
		}
		Path file = Paths.get(outfilePath);
		if (commandLineOptions.getNumberOfPermutations() > 0) {
			// permutationResults contains a hasmap of all pvalues per celltype gotten through permutation testing
			PermutationResult permutationResults = permutationTest(expressionFile, genotypeFile, cellcountTable);
		}
		int whileIndex = 0;
		long time = System.currentTimeMillis();
		// do something that takes some time...
		Collection<Qtl> qtlList = new LinkedList<Qtl>();
		List<DeconvolutionResult> deconvolutionResults = new ArrayList<DeconvolutionResult>();
		while (expressionIterator.hasNext() && genotypeIterator.hasNext()) {
			if (whileIndex % 500 == 0) {
				System.out.printf("Processed %d eQTLs - ", whileIndex);
				long completedIn = System.currentTimeMillis() - time;
				System.out.printf("%s\n", DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS"));
			}
			whileIndex++;
			ArrayList<String> expressionStringVector = new ArrayList<String>(
					Arrays.asList(expressionIterator.next().split("\t")));
			ArrayList<String> genotypeStringVector = new ArrayList<String>(
					Arrays.asList(genotypeIterator.next().split("\t")));


			if(commandLineOptions.getNumberOfThreads() <= 1){
				try{
					deconvolutionResults.add(deconvolution(expressionStringVector, genotypeStringVector, cellcountTable));
				}
				// If there are not enough samples per genotype, skip this QTL
				catch(NotEnoughGenotypesException e){}
				// System.out.printf("\n");
			}
			else{qtlList.add(new Qtl(expressionStringVector, genotypeStringVector, cellcountTable));}
		}
		if(commandLineOptions.getNumberOfThreads() >= 2){
			Parallel parallel = new Parallel(commandLineOptions.getNumberOfThreads());
			parallel.For(qtlList, 
					// The operation to perform with each item
					new Parallel.Operation<Qtl>() {
				public void perform(Qtl param) {
					try {
						try{
							deconvolutionResults.add(deconvolution(param));
						}
						// If there are not enough samples per genotype, skip this QTL
						catch(NotEnoughGenotypesException e){}
					} catch (Exception e) {
						e.printStackTrace();
					}
				};
			});
		}
		output.add("\t"+toString(deconvolutionResults.get(0).GetCelltypes()));
		for(DeconvolutionResult deconResult : deconvolutionResults){
			String results = "";
			results += deconResult.GetQtlName()+"\t"+toString(deconResult.GetPvalues());
			output.add(results);
		}

		Files.write(file, output, Charset.forName("UTF-8"));
		System.out.printf("Deconvolution output written to %s\n", file.toAbsolutePath());
	}

	public static PermutationResult permutationTest(String expressionFile, String genotypeFile, List<List<String>> cellcountTable) throws Exception {
		/*
		 * Scramble the expression values for all the samples n number of times
		 * and calculate p-value, for comparison to p-value from unscrambled
		 * p-value
		 * 
		 * @param expressionIterator Iterator containing all rows of expression file 
		 * 
		 * @param genotypeIterator Iterator containing all rows of genotype file, same order as expression iterator 
		 * 
		 * @param cellcountTable A 2D list with for all samples the different cell counts
		 * 
		 * @return PermutationResult, containing hasmap of per cell type all pvalues
		 */
		System.out.println("Starting permutation test...");
		PermutationResult deconResultsPermutation = new PermutationResult();
		long time = System.currentTimeMillis();
		// make iterator of files so that we can loop over 2 at the same time
		for(int i = 0; i < commandLineOptions.getNumberOfPermutations(); i++){
			LineIterator expressionIterator = FileUtils.lineIterator(new File(expressionFile), "UTF-8");
			LineIterator genotypeIterator = FileUtils.lineIterator(new File(genotypeFile), "UTF-8");
			// Skip over the header
			expressionIterator.next();
			genotypeIterator.next();
			// set a seed so that all rows in the expression file can be shuffled identically
			long seed = System.nanoTime();
			// change number after % to get more/less prints of progress of permutations
			if (i % 1 == 0) {
				System.out.printf("Processed %d/%d permutations - ", i, commandLineOptions.getNumberOfPermutations());
				long completedIn = System.currentTimeMillis() - time;
				System.out.printf("%s\n", DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS"));
			}
			Collection<Qtl>  qtlCollection = new LinkedList<Qtl>();  
			while (expressionIterator.hasNext() && genotypeIterator.hasNext()) {

				ArrayList<String> expressionStringVector = new ArrayList<String>(Arrays.asList(expressionIterator.next().split("\t")));

				ArrayList<String> genotypeStringVector = new ArrayList<String>(Arrays.asList(genotypeIterator.next().split("\t")));
				// shuffle on expression or genotype, because at the moment not sure which one should be used for permutation
				// and using this variable is easier than commenting in/out the 3 related lines
				if (commandLineOptions.getPermutationType().equals("expression")){
					String qtlName = expressionStringVector.remove(0);
					// shuffle the expression vector with the previously set seed
					Collections.shuffle(expressionStringVector, new Random(seed));
					// add back the qtl name at the front
					expressionStringVector.add(0, qtlName);
				}
				else if (commandLineOptions.getPermutationType().equals("genotype")){
					String qtlName = genotypeStringVector.remove(0);
					// shuffle the genotype vector with the previously set seed
					Collections.shuffle(genotypeStringVector, new Random(seed));
					// add back the qtl name at the front
					genotypeStringVector.add(0, qtlName);
				}
				else{
					throw new Exception("shuffle should be either expression or genotype, not "+commandLineOptions.getPermutationType());
				}

				if(commandLineOptions.getNumberOfThreads() <= 1){
					// deconvolution() gives for the current SNP-Gene QTL pair for each celltype in cellount table a pvalue. deconResultsPermutation.add() adds
					// the pvalue to a hashmap of pvalues per celltype. After permutation for loop is finished, deconResultsPermutation will contain a hashmap
					// with for each celltype all pvalues found during permutaiton testing
					try{
						deconResultsPermutation.add(deconvolution(expressionStringVector, genotypeStringVector, cellcountTable), commandLineOptions.getNumberOfPermutations());
					}
					catch(NotEnoughGenotypesException e){
						// If there are not enough samples per genotype for this QTL, exclude it from testing
					}
				}
				else{
					Qtl qtl = new Qtl(expressionStringVector, genotypeStringVector,  cellcountTable);
					qtlCollection.add(qtl);
				}
			}
			if(commandLineOptions.getNumberOfThreads() >= 2){
				Parallel parralel = new Parallel(commandLineOptions.getNumberOfThreads());
				parralel.For(qtlCollection, 
						// The operation to perform with each item
						new Parallel.Operation<Qtl>() {
					public void perform(Qtl qtl) {
						try {
							// qtl contians the expression vector, genotype vector, and cell count table for the qtl
							deconResultsPermutation.add(deconvolution(qtl), commandLineOptions.getNumberOfPermutations());
						} catch (Exception e) {
							e.printStackTrace();
						}
					};
				});
			}
		}
		
		
		commandLineOptions.setOutfolder(commandLineOptions.getOutfolder()+"/shuffled_" + commandLineOptions.getPermutationType());
		if(commandLineOptions.getForceNormalExpression()){
			commandLineOptions.setOutfolder(commandLineOptions.getOutfolder()+"_forcedNormalExpression");
			commandLineOptions.setOutfolder(commandLineOptions.getOutfolder()+"_"+commandLineOptions.getNormalizationType());
		}
		if(commandLineOptions.getForceNormalCellcount()){
			commandLineOptions.setOutfolder(commandLineOptions.getOutfolder()+"_forceNormalCellcount");
		}
		new File(commandLineOptions.getOutfolder()).mkdirs();
		Plots.drawHistogram(deconResultsPermutation, commandLineOptions.getOutfolder());
		return(deconResultsPermutation);
	}

	public static void calculateSumOfSquaresNNLS(LeastSquareModel model, Boolean plotBetaTimesVariables) throws IOException, IllegalAccessException {
		/*
		 * Calculate the sum of squares, using Non-Negative Linear Regression, given a y expression vector with y ~
		 * model. This calculates the Beta parameters three times, for samples with GT = 0, GT = 1, GT = 2
		 * instead of using the interaction terms of the actual genotype.
		 * If no_intercept == true, remove the intercept (equivalent to y * ~ model -1 in R)
		 * 
		 * @param model An InteractionModel object including the y vector expression values and ObservedValues (model)
		 * Such that
		 * test_trait ~ geno_A + lymph% + geno_A:geno_B it can be for one QTL
		 * [[2, 43.4, 86.8], [2, 40.3, 80.6]], for another QTL [[0, 46.7, 0],
		 * [0, 51.5, 0] [0, 48.7, 0]]
		 * 
		 * @return The parameters per genotype dosage
		 */
		// OLS = Ordinary Least Squares
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		// if GetIntercept is false, remove the intercept (Beta1) from the linear model
		regression.setNoIntercept(model.GetNoIntercept());
		regression.newSampleData(model.GetExpessionValuesGt0(), model.GetCellCountsGt0());
		regression.estimateRegressionParameters();
	}
	
	public static double calculateSumOfSquaresOLS(InteractionModel model, Boolean plotBetaTimesVariables) throws IOException, IllegalAccessException {
		/*
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
		// OLS = Ordinary Least Squares
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		// if GetIntercept is false, remove the intercept (Beta1) from the linear model
		regression.setNoIntercept(model.GetNoIntercept());
		try{
			regression.newSampleData(model.GetExpessionValues(), model.GetObservedValues());
		}
		catch (DimensionMismatchException e){
			System.out.printf("expression length: %d\nobserved values length: %d\n", model.GetExpessionValues().length, model.GetObservedValues().length);
			throw(e);
		}

		double sumOfSquares = regression.calculateResidualSumOfSquares();

		if(plotBetaTimesVariables){
			String outfolder = "/Users/NPK/UMCG/projects/deconvolution/parameterTimesValue/";
			new File(outfolder).mkdirs();
			if(model.GetQtlName().length() > 50){
				model.SetQtlName(model.GetQtlName().substring(0, 20));
			}
			Plots.boxPlot(regression, model, outfolder+'/'+model.GetQtlName()+"_"+model.GetModelName()+"_betaTimesExplanatoryVariables.PNG");

		}
		//System.exit(0);
		return (sumOfSquares);
	}

	public static double anova(double sumOfSquaresModelA, double sumOfSquaresModelB, int degreesOfFreedomA,
			int degreesOfFreedomB, Boolean no_intercept) {
		/*
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
		double Fval = meanSquareErrorDiff / meanSquareError;

		/**
		 * Make an F distribution with degrees of freedom as parameter. If full
		 * model and ctModel have the same number of samples, difference in df
		 * is 1 and degreesOfFreedomB are all the terms of the ctModel (so neut%
		 * + eos% + ... + neut% * GT + eos% * GT With 4 cell types and 1891
		 * samples the dfs are 1883 and 1884, giving the below distribution
		 * http://keisan.casio.com/exec/system/1180573186
		 **/
		FDistribution Fdist = new FDistribution(degreesOfFreedomDifference, degreesOfFreedomB);
		/** Calculate 1 - the probability of observing a lower Fvalue **/
		double pval = 1 - Fdist.cumulative(Fval);
		return (pval);
	}

	public static DeconvolutionResult deconvolution(Qtl qtl) throws Exception {
		return(deconvolution(qtl.expressionVector, qtl.genotypeVector, qtl.cellcountTable));
	}
	//TODO: update documentation on parameters
	public static DeconvolutionResult deconvolution(ArrayList<String> expressionStringVector, ArrayList<String>genotypeStringVector,
			List<List<String>> cellcountTable) throws Exception {
		/*
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
		 * @param cellcountTable A 2D list with for all samples the different
		 * cell counts. Should include celltype name in first row.
		 * 
		 * @param minimum_samples_per_genotype Minimum number of samples per genotype need to count it
		 * 
		 * @return A list with for each celltype a p-value for the celltype
		 * specific eQTL for one eQTL
		 */

		// remove the eQTL name from the expression and genotype vector (the first column in the expression and genotype file)
		String qtlName = expressionStringVector.remove(0);
		String genotypeQTL = genotypeStringVector.remove(0);
		if (!qtlName.equals(genotypeQTL)) {
			String errorMessage = "eQTL names not the same for expression and genotype\nexpression file: " + qtlName
					+ "\ngenotype file: " + genotypeQTL;
			throw new Exception(errorMessage);
		}

		double[] expressionVector = StringVectorToDoubleVector(expressionStringVector);
		if(commandLineOptions.getForceNormalExpression()){
			Statistics expressionStatistics = new Statistics(expressionVector);
			if(commandLineOptions.getNormalizationType().equals("normalizeAddMean")){
				expressionVector = expressionStatistics.normalizeKeepMean(expressionVector, false);
			}
			else if (commandLineOptions.getNormalizationType().equals("normalizeKeepExponential")){
				expressionVector = expressionStatistics.normalizeKeepExponential(expressionVector, false);
			}
		}
		
		double[] genotypeVector = StringVectorToDoubleVector(genotypeStringVector);
		/* 
		 * If roundDosage option is selected on the command line, round of the dosage to closest integer -> 0.49 = 0, 0.51 = 1, 1.51 = 2. 
		 * If minimumSamplesPerGenotype is selected on the command line, check for current QTL if for each dosage (in case they are not round
		 * the dosages are binned in same way as with roundDosage option) there are at least <minimumSamplesPerGenotype> samples that have it.
		 */

		if (commandLineOptions.getRoundDosage() || commandLineOptions.getMinimumSamplesPerGenotype() > 0 || commandLineOptions.getAllDosages()) {
			int dosage_ref = 0;
			int dosage_heterozygote = 0;
			int dosage_alt = 0;
			for (int i = 0; i < genotypeVector.length; i++) {
				if (commandLineOptions.getRoundDosage()){
					genotypeVector[i] = Math.round(genotypeVector[i]);
				}
				if (commandLineOptions.getMinimumSamplesPerGenotype() > 0 || commandLineOptions.getAllDosages()){
					if(genotypeVector[i] < 0){
						throw new Exception("Genotype dosage can not be negative, check your dosage input file");
					}
					if(genotypeVector[i] < 0.5){
						dosage_ref++;
					}
					else if(genotypeVector[i] < 1.5){
						dosage_heterozygote++;
					}
					else if (genotypeVector[i] <= 2){
						dosage_alt++;
					}
					else{
						throw new Exception("Genotype dosage can not be larger than 2, check your dosage input file");
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
					throw new NotEnoughGenotypesException("Not enough samples for each genotype");
				}}
		}

		List<Double> pvalues = new ArrayList<Double>();
		int numberOfCelltypes = cellcountTable.size();
		List<String>  celltypes = new ArrayList<String>();
		// minus one because the size includes the celltype header
		int numberOfSamples = cellcountTable.get(0).size()-1;	

		// variable for saving which celltypes belong to the specific model
		// things neded for fullModel defined outside of loop because every celltype model (ctModel) has to be compared to it
		double sumOfSquaresFullModel = 0;
		double[][] observedValuesFullModel = null;
		int degreesOfFreedomFullModel = 0;
		int fullModelLength = 0;

		InteractionModel fullModel = new InteractionModel();
		fullModel.SetModelName("full model");
		fullModel.SetQtlName(qtlName);
		LeastSquareModel leastSquareModel = new LeastSquareModel();
		/*
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
		// TODO: For loop looks overly complicated, see if it can be simplified
		for (int m = 0; m < numberOfCelltypes; m++) {
			InteractionModel ctModel = new InteractionModel();
			ctModel.SetExpressionValues(expressionVector);
			int numberOfTerms = (numberOfCelltypes * 2) - 1;
			double[][] observedValues = new double[numberOfSamples][numberOfTerms];
			// fullModel will be done in first loop as all models have to be compared to it
			if(m == 0){
				fullModel.SetExpressionValues(expressionVector);
				// number of terms + 1 because for full model all cell types are included
				observedValuesFullModel = new double[numberOfSamples][numberOfTerms+1];
			}
			int genotypeCounter = numberOfCelltypes;
			for (int j = 0; j <= numberOfSamples-1; j++) {
				for (int i = 0; i < numberOfCelltypes; i++) {
					// There is one fullModel including all celltypes add values for celltypePerc and interaction term of
					// celltypePerc * genotypePerc so that you get [[0.3, 0.6], [0.4, 0.8], [0.2, 0.4], [0.1, 0.2]]
					// where numberOfSamples = 1 and numberOfCellTypes = 4 with celltypePerc = 0.3, 0.4, 0.2, and 0.1 and genotype = 2
					// for each cell type is 1 model, celltype% * genotype without 1 celltype.
					// j+1 because j==0 is header
					double celltype_perc = Double.parseDouble(cellcountTable.get(i).get(j+1));
					observedValues[j][i] = celltype_perc;
					if(j == 0){
						// add the celltype name at position i so that it gets in front of the celltype:GT, but once
						ctModel.AddIndependentVariableName(i, cellcountTable.get(i).get(0));
					}

					// if i (cell type index) is the same as m (model index), don't add the interaction term of celltype:GT
					if (i != m) {
						try {
							// Only add IndependentVariableName once per QTL (j==0)
							if(j == 0){

								// Add the interaction term of celltype:genotype
								ctModel.AddIndependentVariableName(cellcountTable.get(i).get(0)+":GT");
								// save the index of the variables related to current celltype so that this can be used later to calculate
								// Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
								// but for ctModel this is easiest method
								int[] index = new int[] {i, numberOfCelltypes-1+i};
								ctModel.AddCelltypeVariablesIndex(index);
								// add the celltype name. This could be done with less code by getting it from IndependentVariableName, but this way 
								// it is explicit. Don't know if better.
								fullModel.AddCelltype(cellcountTable.get(i).get(0));
							}
							observedValues[j][genotypeCounter] = celltype_perc * genotypeVector[j];
						} catch (ArrayIndexOutOfBoundsException error) {
							throw new Exception(
									"The counts file and expression and/or genotype file do not have equal number of samples or QTLs",
									error);
						}
						genotypeCounter++;
					}
					// if i==m there is not celltype:GT interaction term so only one index added to CelltypeVariables
					else if (j == 0){
						int[] index = new int[] {i};
						ctModel.AddCelltypeVariablesIndex(index);
						// add the celltype name. This could be done with less code by getting it from IndependentVariableName, but this way 
						// it is explicit. Don't know if better.
						fullModel.AddCelltype(cellcountTable.get(i).get(0));
					}
					if (m == 0){
						observedValuesFullModel[j][i] = celltype_perc;
						try {
							if(j == 0){
								/* save the index of the variables related to current celltype so that this can be used later to calculate
								 * Beta1 celltype% + Beta2 * celltype%:GT. For fullModel not so necesarry as it's always <numberOfCelltypes> away,
								 * but for ctModel this is easiest method
								 */
								int[] index = new int[] {i, numberOfCelltypes + i};
								fullModel.AddCelltypeVariablesIndex(index);
								// add the celltype name. This could be done with less code by getting it from IndependentVariableName, but this way 
								// it is explicit. Don't know if better.
								fullModel.AddCelltype(cellcountTable.get(i).get(0));
								// add the celltype name at position i so that it gets in front of the celltype:GT
								fullModel.AddIndependentVariableName(i, cellcountTable.get(i).get(0));
								fullModel.AddIndependentVariableName(cellcountTable.get(i).get(0)+":GT");
								celltypes.add(cellcountTable.get(i).get(0));

							}
							// Add the interaction term of celltype:genotype
							/*if(genotypeVector[j] < 0.5){
								leastSquareModel.AddExpressionValuesGt0(expressionVector[j]);
								leastSquareModel.AddCellCountsGt0(0, celltype_perc);
							}
							else if(genotypeVector[j] < 1.5){
								leastSquareModel.AddExpressionValuesGt1(expressionVector[j]);
								leastSquareModel.AddCellCountsGt1(1, celltype_perc);
							}
							else if(genotypeVector[j] <= 2){
								leastSquareModel.AddExpressionValuesGt2(expressionVector[j]);
								leastSquareModel.AddCellCountsGt2(2, celltype_perc);
							}
							else{
								throw new Exception("Dosage should not be below 0 or above 2");
							}*/
							observedValuesFullModel[j][numberOfCelltypes + i] = celltype_perc * genotypeVector[j];
						} catch (ArrayIndexOutOfBoundsException error) {
							throw new Exception(
									"The counts file and expression and/or genotype file do not have equal number of samples or QTLs",
									error);
						}
					}
				}
				// because 1 of numberOfCelltypes + i needs to be skipped,
				// keeping it tracked with separate value is easier
				genotypeCounter = numberOfCelltypes;
			}
			Boolean noIntercept = true;
			if(m == 0){
				// only need to set data of fullModel once, reused every loop of m
				fullModel.SetObservedValues(observedValuesFullModel);
				fullModel.SetNoIntercept(noIntercept);
				sumOfSquaresFullModel = calculateSumOfSquaresOLS(fullModel, commandLineOptions.getPlotBetas());
				degreesOfFreedomFullModel = expressionVector.length - (observedValuesFullModel[0].length + 1);
				fullModelLength = fullModel.GetObservedValues().length;
			}
			ctModel.SetObservedValues(observedValues);
			ctModel.SetNoIntercept(noIntercept);
			ctModel.SetQtlName(qtlName);
			/** SUM OF SQUARES - CELLTYPE MODEL **/
			ctModel.SetModelName("ctModel_"+ Integer.toString(m));
			double sumOfSquaresCtModel = calculateSumOfSquaresOLS(ctModel, commandLineOptions.getPlotBetas());
			int degreesOfFreedomCtModel = expressionVector.length - (ctModel.GetObservedValues()[0].length + 1);

			int expressionLength = expressionVector.length;

			if (expressionLength != fullModelLength) {
				throw new Exception("expression vector and fullModel have different number of samples.\nexpression: "
						+ expressionLength + "\nfullModel: " + fullModelLength);
			}

			// ANOVA compare full model to celltype model 
			double pval = anova(sumOfSquaresFullModel, sumOfSquaresCtModel, degreesOfFreedomFullModel, degreesOfFreedomCtModel, true);
			pvalues.add(pval);
		}
		//calculateSumOfSquaresNNLS(leastSquareModel, true);
		DeconvolutionResult deconResult = new DeconvolutionResult(celltypes, qtlName, pvalues);
		return deconResult;
	}

	public static List<List<String>> readTabDelimitedColumns(String filepath) throws IOException {
		/*
		 * Reads tab delimited file and returns them as list of list, with [x] =
		 * colummn and [x][y] is value in column. Needed for reading counts
		 * file, as there the rows are the samples, as opposed to expression and
		 * genotype file where the columns are the samples. Needs to be read in
		 * memory, so minimal memory requirement is larger than the size of the
		 * counts file.
		 * 
		 * @param filepath The path to a tab delimited file to read
		 * 
		 * @return A 2D array with each array being one column from filepath
		 */
		List<List<String>> allColumns = new ArrayList<List<String>>();
		// parses file on tabs
		CSVParser parser = new CSVParser(new FileReader(filepath), CSVFormat.newFormat('\t'));
		for (CSVRecord row : parser) {
			// starts at 1 because 1st element of column is the samplename
			for (int i = 1; i < row.size(); i++) {
				// First try to add an element from row to the *i*th list, if
				// the *i*th list does not exist yet
				// catch the IndexOutOfBoundsException, make a new list and it
				// at the *i*th position of the 2D array allColumns
				try {
					allColumns.get(i - 1).add(row.get(i));
				} catch (IndexOutOfBoundsException e) {
					List<String> newColumn = new ArrayList<String>();
					newColumn.add(row.get(i));
					allColumns.add(newColumn);	
				}
			}
		}
		parser.close();
		return allColumns;
	}

	public static double[] StringVectorToDoubleVector(ArrayList<String> vector) {
		/*
		 * Converting a vector of string to a vector of doubles
		 * 
		 * @param vector A vector of strings
		 * 
		 * @return A vector of doubles
		 */
		int vectorLength = vector.size();
		double[] doubles = new double[vectorLength];
		for (int i = 0; i < vectorLength; i++) {
			doubles[i] = Double.parseDouble(vector.get(i));
		}
		return doubles;
	}




	public static <T> String toString(List<T> celltypes)
	{
		/* Turn list into tab separated string*/
		StringBuilder builder = new StringBuilder();
		for(Object o: celltypes)
		{
			builder.append(o+"\t");
		}
		return builder.toString().trim();
	}
}
