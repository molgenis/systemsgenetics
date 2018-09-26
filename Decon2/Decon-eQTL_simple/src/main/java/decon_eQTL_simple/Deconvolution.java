package main.java.decon_eQTL_simple;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.lang3.time.DurationFormatUtils;

import main.java.decon_eQTL_simple.InteractionModel;
import main.java.decon_eQTL_simple.InteractionModelCollection;
import main.java.decon_eQTL_simple.Statistics;
import main.java.decon_eQTL_simple.CellCount;

public class Deconvolution {
	private int QTLsFiltered = 0;
	private String outputFolder;
	public CellCount cellCounts;
	public ExpressionData expressionData;
	public HashMap<String,ArrayList<String>> geneSnpPairs;
	private CommandLineOptions commandLineOptions;
	public GenotypeData genotypeData;
	private cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = null;
	private cern.jet.random.tdouble.StudentT tDistColt = null;
	
	public Deconvolution(CommandLineOptions commandLineOptions) {
		this.commandLineOptions = commandLineOptions;
		randomEngine = new cern.jet.random.tdouble.engine.DRand();
	}

	/*
	 * Read all the input data
	 * 
	 * @throw IllegalAccessException	If file can't be accessed
	 * @throw IOException	If file can not be read
	 */
	public void readInputData() throws IllegalAccessException, IOException {
		outputFolder = commandLineOptions.getOutfolder();
		cellCounts = new CellCount();
		cellCounts.readCellCountData(commandLineOptions.getCellcountFile());

		geneSnpPairs = Utils.parseSnpPerGeneFile(commandLineOptions.getSnpsToTestFile());
		String expressionFile = commandLineOptions.getExpressionFile();
		DeconvolutionLogger.log.info(String.format("Parse expression data from %s",expressionFile));
		expressionData = new ExpressionData(expressionFile);
		DeconvolutionLogger.log.info("Done");
		String genotypeFile = commandLineOptions.getGenotypeFile();
		DeconvolutionLogger.log.info(String.format("Parse genotype data from %s",genotypeFile));
		genotypeData = new GenotypeData(genotypeFile);

		DeconvolutionLogger.log.info("Done");
		if (!expressionData.getSampleNames().equals(genotypeData.getSampleNames())){
			ArrayList<String> differences = Utils.getDifferencesBetweenLists(expressionData.getSampleNames(), genotypeData.getSampleNames());
			throw new RuntimeException(String.format("Samplenames not the same in expression and genotype file, or not in the same order."+
					"\nexpression samples not in genotypes (%d): %s\ngenotype samples not in expression (%d): %s\n",
					differences.get(0).length(), differences.get(0), 
					differences.get(1).length(), differences.get(1)));
		}
	}

	/**
	 * For each of the gene-SNP pair in the SnpsToTestFile run deconvolution
	 * 
	 * @param commandLineOptions	commandLineOptions to run it with 
	 * 
	 * @return Deconvolution results
	 * 
	 * @throws RuntimeException
	 * @throws IllegalAccessException 
	 * @throws IOException 
	 */
	public List<DeconvolutionResult> runDeconPerGeneSnpPair() throws RuntimeException, IllegalAccessException, IOException{
		int whileIndex = 0;
		long time = System.currentTimeMillis();
		List<DeconvolutionResult> deconvolutionResults = new ArrayList<DeconvolutionResult>();
		int QTLsTotal = 0;
		HashMap<String, double[]> geneExpressionLevels = expressionData.getGeneExpression();
		int skippedGenotypeGeneCombinations = 0;
		for(String gene : geneSnpPairs.keySet()){
			for(String genotype : geneSnpPairs.get(gene)){
				if (whileIndex % 500 == 0) {
					long completedIn = System.currentTimeMillis() - time;
					DeconvolutionLogger.log.info(String.format("Processed %d gene-SNP pairs - %s - skipped %d gene-SNP combinations", whileIndex, DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS"), skippedGenotypeGeneCombinations));
				}
				++whileIndex;
				String qtlName = gene+'_'+genotype;
				++QTLsTotal;
				double[] dosages = genotypeData.getGenotypes().get(genotype);
				if(dosages == null){
					DeconvolutionLogger.log.info(String.format("Error: Genotype %s included in gene/snp combinations to test, but not available in the expression file!",genotype));
					throw new RuntimeException(String.format("Error: Genotype %s included in gene/snp combinations to test, but not available in the expression file!",genotype));
				}
				double[] expressionLevels = geneExpressionLevels.get(gene);

				if(expressionLevels == null){
					DeconvolutionLogger.log.info(String.format("Error: Gene %s included in gene/snp combinations to test, but not available in the expression file!",gene));
					throw new RuntimeException(String.format("Gene %s included in gene/snp combinations to test, but not available in the expression file!",gene));

				}
				DeconvolutionResult deconResult = deconvolution(expressionLevels, dosages, qtlName);
				deconvolutionResults.add(deconResult);
			}
		}
		DeconvolutionLogger.log.info(String.format("Skipped %d gene-SNP combinations (because genotype in SNP-pair file but not in genotype file)",skippedGenotypeGeneCombinations));
		DeconvolutionLogger.log.info(String.format("QTLs passed: %d", QTLsTotal-(QTLsFiltered+skippedGenotypeGeneCombinations)));
		DeconvolutionLogger.log.info(String.format("QTLs filtered: %d", QTLsFiltered));
		DeconvolutionLogger.log.info(String.format("Total: %d",QTLsTotal-skippedGenotypeGeneCombinations));
		return deconvolutionResults;
	}

	/**
	 * Write the deconvolution results
	 * 
	 * @param deconvolutionResult The deconvolution result
	 */
	public void writeDeconvolutionResults(List<DeconvolutionResult> deconvolutionResults) throws IllegalAccessException, IOException{
		List<String> celltypes = cellCounts.getAllCelltypes();
		String header = "\t"+Utils.listToTabSeparatedString(celltypes, "_pvalue");

		DeconvolutionLogger.log.info("Getting decon result with full model info for writing the header");
		// celltypes.size()*2 because there are twice as many betas as celltypes (CC% & CC%:GT)


		for(int i = 0; i < cellCounts.getNumberOfCelltypes(); ++i){
			header += "\tBeta" + Integer.toString(i+1) +"_"+cellCounts.getCelltype(i)+":GT";
		}


		//header += "\tStandardError";
		List<String> output = new ArrayList<String>();
		output.add(header);
		for(DeconvolutionResult deconvolutionResult : deconvolutionResults){
			InteractionModelCollection interactionModelCollection = deconvolutionResult.getInteractionModelCollection();
			String results = "";
			results += deconvolutionResult.getQtlName()+"\t"+Utils.listToTabSeparatedString(deconvolutionResult.getPvalues());
			for(String modelName : interactionModelCollection.getModelNames()) {
				InteractionModel model = interactionModelCollection.getInteractionModel(modelName);

				double[] estimateRegressionParameters = model.getEstimateRegressionParameters();

				// because model is y = intercept + snp + cc + snp:cc the interaction beta is at index [3]
				// but have to check what direction the beta of snp has, because if it is negative the interaction
				// term beta has to be flipped.
				if (estimateRegressionParameters[1] < 0) {
					estimateRegressionParameters[3] = estimateRegressionParameters[3]*-1;
				}
				results += "\t"+estimateRegressionParameters[3];
			}
			output.add(results);	
		}
		Path file = Paths.get(outputFolder+"/"+commandLineOptions.getOutfile());
		Files.write(file, output, Charset.forName("UTF-8"));

		DeconvolutionLogger.log.info(String.format("Deconvolution output written to %s", file.toAbsolutePath()));
		DeconvolutionLogger.log.info(String.format("Files with additional info in  %s", outputFolder));
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
	private DeconvolutionResult deconvolution(double[] expression, double[] genotypes, String qtlName) throws RuntimeException, IllegalAccessException, 
	IOException {
		InteractionModelCollection interactionModelCollection = new InteractionModelCollection(cellCounts);
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
		interactionModelCollection.createObservedValueMatrices();
		interactionModelCollection.calculateOlsForAllModels();		
		
		calculateDeconvolutionPvalue(interactionModelCollection);

		DeconvolutionResult deconResult =  new DeconvolutionResult();
		
		interactionModelCollection.cleanUp();
		deconResult = new DeconvolutionResult(interactionModelCollection);
		return deconResult;
	}

	/**
	 * get pvalue for each ctmodel
	 * 
	 * @param interactionModelCollection InteractionModelCollection object that has fullModel and ctModels for ANOVA comparison
	 */
	private void calculateDeconvolutionPvalue(InteractionModelCollection interactionModelCollection) 
			throws IllegalAccessException, IOException {
	
		for (String interactionModelName : interactionModelCollection.getModelNames()) {
			
			InteractionModel interactionModel = interactionModelCollection.getInteractionModel(interactionModelName);
			tDistColt = new cern.jet.random.tdouble.StudentT(interactionModel.getDegreesOfFreedom(), randomEngine);
//			InteractionModel snpModel = interactionModelCollection.getSnpModel(interactionModelName);
			//double pval = Statistics.anova(interactionModel.getSumOfSquares(), snpModel.getSumOfSquares(), 
			//		interactionModel.getDegreesOfFreedom(),snpModel.getDegreesOfFreedom(), 
			//		true);
			double se = interactionModel.getEstimateRegressionParametersStandardErrors()[2];
			double pval = Statistics.convertBetaToP(interactionModel.getEstimateRegressionParameters()[2],se,tDistColt);
			interactionModel.setPvalue(pval);
			interactionModelCollection.setPvalue(pval, interactionModelName);
		}
	}
}
