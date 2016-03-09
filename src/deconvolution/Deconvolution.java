package deconvolution;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
//import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import JSci.maths.statistics.FDistribution;

/*
 * -o decovolution.txt 
-e /Users/NPK/UMCG/projects/deconvolution/eqtl_exp.tsv 
-g /Users/NPK/UMCG/projects/deconvolution/eqtl_dosage.tsv 
-c /Users/NPK/UMCG/projects/deconvolution/counts1.txt
 */
public class Deconvolution {
	public static void main(String[] args) throws Exception {
		/*
		 * To run: java -jar <expression table> <genotype table> <cell count
		 * table> For all three rows are samples and columns are SNPs
		 * 
		 */
		// begin command line option parsing
		Options options = new Options();
		Option help = new Option("help", "print this message");
		Option eqtl_exp_names_flag = new Option("eqtl_exp_names", "If set, first column of expression file is eQTL names");
		Option eqtl_genotype_names_flag = new Option("eqtl_geno_names", "If set, first column of genotype file is eQTL names");
		Option permute = Option.builder("p").required(false).longOpt("permute").desc("Do permutations").build();
		Option outfile = Option.builder("o").required(true).hasArg().longOpt("outfile").desc("Outfile name")
				.argName("file").build();
		Option expression = Option.builder("e").required(true).hasArg().longOpt("expression")
				.desc("Expression file name").argName("file").build();
		Option genotype = Option.builder("g").required(true).hasArg().longOpt("genotype").desc("Genotype file name")
				.argName("file").build();
		Option cellcount = Option.builder("c").required(true).hasArg().longOpt("cellcount").desc("Cellcount file name")
				.argName("file").build();

		options.addOption(help);
		options.addOption(permute);
		options.addOption(outfile);
		options.addOption(expression);
		options.addOption(genotype);
		options.addOption(cellcount);
		options.addOption(eqtl_exp_names_flag);
		options.addOption(eqtl_genotype_names_flag);
		CommandLineParser cmdLineParser = new DefaultParser();
		CommandLine cmdLine = cmdLineParser.parse(options, args);
		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		if (cmdLine.hasOption("help")) {
			formatter.printHelp("deconvolution", options, true);
		}
		Boolean eqtl_exp_names = false;
		if (cmdLine.hasOption("eqtl_exp_names")) {
			eqtl_exp_names = true;
		}
		Boolean eqtl_genotype_names = false;
		if (cmdLine.hasOption("eqtl_geno_names")) {
			eqtl_genotype_names = true;
		}
		Boolean permutation = false;
		if (cmdLine.hasOption("permute")) {
			permutation = true;
		}
		// End of command line option parsing


		// the cell type names are the first row of cellcount file, extract for
		// later printing
		String cellCountFile =  cmdLine.getOptionValue("cellcount");
		LineIterator cellcountIterator = FileUtils.lineIterator(new File(cellCountFile), "UTF-8");
		ArrayList<String> celltypes = new ArrayList<String>(Arrays.asList(cellcountIterator.next().split("\t")));
		celltypes.removeAll(Arrays.asList("", null));
		cellcountIterator.close();
		List<List<String>> cellcountTable = readTabDelimitedColumns(cellCountFile);
		String outfilePath = cmdLine.getOptionValue("o");


		// output saves all the processed output to write it to a file later
		List<String> output = new ArrayList<String>();
		List<String> output_permutation = new ArrayList<String>();
		String celltypeLine = new String();
		for (int i = 0; i < celltypes.size(); i++) {
			//System.out.printf("%s\t", celltypes.get(i));
			celltypeLine += celltypes.get(i);
			if (i != celltypes.size()-1){
				celltypeLine += "\t";
			}
		}
		if (permutation){
			output_permutation.add(celltypeLine);
		}
		output.add(celltypeLine);
		System.out.printf("\n");
		// make iterator of file so that we can loop over 2 at the same time 
		String expressionFile = cmdLine.getOptionValue("expression");
		String genotypeFile = cmdLine.getOptionValue("genotype");
		LineIterator expressionIterator = FileUtils.lineIterator(new File(expressionFile), "UTF-8");
		LineIterator genotypeIterator = FileUtils.lineIterator(new File(genotypeFile), "UTF-8");
		expressionIterator.next();
		genotypeIterator.next();
		
		while (expressionIterator.hasNext() && genotypeIterator.hasNext()) {
			ArrayList<String> expressionStringVector = new ArrayList<String>(Arrays.asList(expressionIterator.next().split("\t")));
			ArrayList<String> genotypeStringVector = new ArrayList<String>(Arrays.asList(genotypeIterator.next().split("\t")));
			// remove the eQTL name from the vector
			String eQTL = new String();
			if (eqtl_exp_names){
				// if the first column of expression file is eQTL name, 
				// remove it and save it to compare incase genotype file also has eQTL names
				eQTL = expressionStringVector.remove(0);
			}
			if (eqtl_genotype_names){
				// if the first column of genotype file is eQTL name, 
				// remove it. If eQTL name was not in expression file, add it to eQTL string
				String eQTL_geno = genotypeStringVector.remove(0);
				if (eQTL.length() > 0){
					if (!eQTL.equals(eQTL_geno)){
						throw new Exception("eQTLs not the same for expression and genotype");
					}
				}else{
					eQTL = eQTL_geno;
				}
				genotypeStringVector.remove(0);
			}
			double[] expressionVector = StringVectorToDoubleVector(expressionStringVector);
			double[] genotypeVector = StringVectorToDoubleVector(genotypeStringVector);
			List<Double> pvalues_permutation = new ArrayList<Double>();
			if (permutation){
				double[] expressionVectorPermute = expressionVector;
				Collections.shuffle(Arrays.asList(expressionVectorPermute));
				pvalues_permutation = deconvolution(expressionVectorPermute, genotypeVector, cellcountTable);
			}
			List<Double> pvalues = deconvolution(expressionVector, genotypeVector, cellcountTable);
			String results = new String("");
			String results_permutation = new String("");

			if (eQTL.length()>0){
				results += eQTL+"\t";
				if (permutation){
					results_permutation += eQTL+"\t";
				}
				//System.out.printf("%s\t",eQTL);
			}
			
			for (int j = 0; j < pvalues.size(); j++) {
				//System.out.printf("%f\t", pvalues.get(j));
				results += pvalues.get(j) + "\t";
				if (permutation){
					results_permutation += pvalues_permutation.get(j) + "\t";
				}
			}
			output.add(results);
			if (permutation){
				output_permutation.add(results_permutation);
			}
			//System.out.printf("\n");
		}
		Path file = Paths.get(outfilePath);
		Files.write(file, output, Charset.forName("UTF-8"));
		if (permutation){
			String permutation_outpath = outfilePath.replaceAll("(\\.[^.]*)$","_permutation$1");
			file = Paths.get(permutation_outpath);
			Files.write(file, output_permutation, Charset.forName("UTF-8"));
			System.out.printf("Output written to %s", permutation_outpath);
		}
		System.out.printf("Output written to %s", outfilePath);
	}

	public static double calculateSumOfSquares(double[] y, double[][] model, Boolean intercept){
		OLSMultipleLinearRegression regr = new OLSMultipleLinearRegression();
		regr.setNoIntercept(intercept);
		regr.newSampleData(y, model);
		double sumOfSquares = regr.calculateResidualSumOfSquares();
		return (sumOfSquares);
	}
	

	public static double anova(double sumOfSquaresModelA, double sumOfSquaresModelB, int degreesOfFreedomA, int degreesOfFreedomB){
		/** From Joris Meys: http://stackoverflow.com/a/35458157/651779
		 * 1. calculate MSE for the largest model by dividing the Residual Sum of Squares (RSS) by the degrees of freedom (df)
		 * 2. calculate the MSEdifference by substracting the RSS of both models (result is "Sum of Sq." in the R table), substracting the df for both models (result is "Df" in the R table), and divide these numbers.
		 * 3. Divide 2 by 1 and you have the F value
		 * 4.calculate the p-value using the F value in 3 and for df the df-difference in the numerator and df of the largest model in the denominator.
		 * For more info: http://www.bodowinter.com/tutorial/bw_anova_general.pdf
		 **/
		// Within-group Variance
		double meanSquareError = sumOfSquaresModelB / degreesOfFreedomB; 
		int degreesOfFreedomDifference = Math.abs(degreesOfFreedomB - degreesOfFreedomA);
		// Between-group Variance
		double meanSquareErrorDiff = Math.abs((sumOfSquaresModelB - sumOfSquaresModelA) / (degreesOfFreedomDifference));

		/** F = Between-group Variance / Within-group Variance <- high value if variance between the models is high, and
		 * 														  variance within the models is low **/
		double Fval = meanSquareErrorDiff / meanSquareError;
		/**Make an F distribution with degrees of freedom as parameter. If full model and 
		 * ctModel have the same number of samples, difference in df is 1 and degreesOfFreedomB
		 * are all the terms of the ctModel (so neut% + eos% + ... + neut% * GT + eos% * GT
		 * With 4 cell types and 1891 samples the dfs are 1 and 1883, giving the below distribution
		 * http://keisan.casio.com/exec/system/1180573186
		 * **/
		FDistribution Fdist = new FDistribution(degreesOfFreedomDifference, degreesOfFreedomB);
		/**Calculate 1 - the probability of observing a lower Fvalue**/
		double pval = 1 - Fdist.cumulative(Fval);
		return(pval);
	}
	
	public static List<Double> deconvolution(double[] expressionVector, double[] genotypeVector,
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
		 */
		// numberOfCelltypes-1 because sample name is included in the columns of
		// cellcountTable
		int numberOfCelltypes = cellcountTable.size();
		int numberOfSamples = cellcountTable.get(0).size();

		// fullModel [genotype][celltype]
		double[][] fullModel = new double[numberOfSamples][numberOfCelltypes * 2];
		for (int j = 0; j < numberOfSamples; j++) {
			for (int i = 0; i < numberOfCelltypes; i++) {
				// add values for celltypePerc and interaction term of celltypePerc * genotypePerc so that you get
				// [[0.3, 0.6], [0.4, 0.8], [0.2, 0.4], [0.1, 0.2]]
				// where numberOfSamples = 1 and numberOfCellTypes = 4
				// with celltypePerc = 0.3, 0.4, 0.2, and 0.1 and genotype = 2
				double celltype_perc =  Double.parseDouble(cellcountTable.get(i).get(j));
				fullModel[j][i] = celltype_perc;
				//fullModel[j][i] = celltypePercentage[i][j];
				fullModel[j][numberOfCelltypes + i] = celltype_perc * genotypeVector[j];
			}
		}
		/**SUM OF SQUARES - FULL MODEL**/
		double sumOfSquaresA = calculateSumOfSquares(expressionVector, fullModel, true);

		int expressionLength = expressionVector.length;
		int fullModelLength = fullModel.length;
		if (expressionLength != fullModelLength){
			throw new Exception("expression vector and fullModel have different number of samples.\nexpression: "+expressionLength+"\nfullModel: "+fullModelLength);
		}
		int degreesOfFreedomA = expressionVector.length - (fullModel[0].length + 1);
		// df = n - number of coefficients, including intercept
		double[][] ctModel = new double[numberOfSamples][(numberOfCelltypes * 2) - 1];
		List<Double> pvalues = new ArrayList<Double>();
		// should be in same loop as fullModel but couldn't get the
		// genotypeCounter correct
		// m = model, there are equally many models as celltypes
		for (int m = 0; m < numberOfCelltypes; m++) {
			// because 1 of numberOfCelltypes + i needs to be skipped,
			// keeping it tracked with separate value is easier
			int genotypeCounter = numberOfCelltypes;
			for (int j = 0; j < numberOfSamples; j++) {
				for (int i = 0; i < numberOfCelltypes; i++) {
					// for each cell type is 1 model, celltype% * genotype
					// without 1 celltype
					//ctModel[j][i] = celltypePercentage[i][j];
					double celltype_perc = Double.parseDouble(cellcountTable.get(i).get(j));
					ctModel[j][i] = celltype_perc;
					if (i != m) {
						//ctModel[j][genotypeCounter] = celltypePercentage[i][j] * genotypeVector[j];
						ctModel[j][genotypeCounter] = celltype_perc * genotypeVector[j];
						genotypeCounter++;
					}
				}
				// reset genotypeCounter
				genotypeCounter = numberOfCelltypes;
			}
			/**SUM OF SQUARES - CELLTYPE MODEL**/
			double sumOfSquaresB = calculateSumOfSquares(expressionVector, ctModel, true);
			int degreesOfFreedomB = expressionVector.length - (ctModel[0].length + 1);
			/*******ANOVA COMPARE FULL MODEL TO CELLTYPE MODEL************/
			double pval = anova(sumOfSquaresA, sumOfSquaresB, degreesOfFreedomA, degreesOfFreedomB);
			pvalues.add(pval);
		}
		return pvalues;
	}

	public static List<List<String>> readTabDelimitedColumns(String filepath) throws IOException {
		/*
		 * Reads tab delimited file and returns them as list of list, with [x] =
		 * colummn and [x][y] is value in column
		 */
		List<List<String>> allColumns = new ArrayList<List<String>>();
		CSVParser parser = new CSVParser(new FileReader(filepath), CSVFormat.newFormat('\t').withHeader());
		
		for (CSVRecord record : parser) {
			for (int i = 1; i < record.size(); i++){
				try {
					allColumns.get(i-1).add(record.get(i));
				} catch (IndexOutOfBoundsException e) {
					List<String> newColumn = new ArrayList<String>();
					newColumn.add(record.get(i));
					allColumns.add(newColumn);
				}
			}
		}
		parser.close();
		return allColumns;
	}

	public static double[] StringVectorToDoubleVector(ArrayList<String> vector) {
		int vectorLength = vector.size();
		double[] doubles = new double[vectorLength];
		for (int i = 0; i < vectorLength; i++) {
				doubles[i] = Double.parseDouble(vector.get(i));
		}
		return doubles;
	}
}
