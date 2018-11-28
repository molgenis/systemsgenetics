package main.java.decon_eQTL;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.lang3.time.DurationFormatUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import main.java.decon_eQTL.CellCount;

public class Deconvolution {
	private int QTLsFiltered = 0;
	private String outputFolder;
	public CellCount cellCounts;
	public ExpressionData expressionData;
	public HashMap<String,ArrayList<String>> geneSnpPairs;
	private CommandLineOptions commandLineOptions;
	public GenotypeData genotypeData;

	public Deconvolution(CommandLineOptions commandLineOptions) {
		this.commandLineOptions = commandLineOptions;
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
				if(commandLineOptions.getTestRun() && whileIndex == 100){
					break;
				}
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
		InteractionModelCollection firstInteractionModelCollection = deconvolutionResults.get(0).getInteractionModelCollection();
		InteractionModel bestFullModelForHeaderOnly = firstInteractionModelCollection.getBestFullModel();

		for(int i = 1; i < cellCounts.getNumberOfCelltypes()*2 + 1; ++i){
			header += "\tBeta" + Integer.toString(i) +"_"+bestFullModelForHeaderOnly.getIndependentVariableNames().get(i-1);
		}
		//for(String celltype : celltypes){
		//	header += "\teffectDirectionDosage2_"+celltype;
		//}
		//header += "\tgenotypeConfiguration";
		//for(String celltype : cellCounts.getAllCelltypes()){
		//	header += "\tgenotypeConfiguration_"+celltype;
		//}\

		if(commandLineOptions.getWholeBloodQTL()){
			header += "\tSpearman correlation expression~GT\tSpearman correlation p-value";
		}

		//header += "\tStandardError";
		List<String> output = new ArrayList<String>();
		output.add(header);
		for(DeconvolutionResult deconvolutionResult : deconvolutionResults){
			InteractionModelCollection interactionModelCollection = deconvolutionResult.getInteractionModelCollection();

			String results = "";
			results += deconvolutionResult.getQtlName()+"\t"+Utils.listToTabSeparatedString(deconvolutionResult.getPvalues());
			InteractionModel bestFullModel = null;

			bestFullModel = interactionModelCollection.getBestFullModel();


			double[] estimateRegressionParameters = bestFullModel.getEstimateRegressionParameters();

			// check what the genotype configuration is and the beta of the interaction term. 
			// If genotype configuration == 0 and beta == positive, dosage2 effect = positive
			// If genotype configuration == 1 and beta == negative, dosage2 effect = positive
			// else is negative
			int numberOfCelltypes = cellCounts.getNumberOfCelltypes();
			// first write out the beta of the cell proportion term
			for(int i = 0; i < numberOfCelltypes; ++i){
				results += "\t"+estimateRegressionParameters[i];
			}
			
			// then write out cell proportion-genotype interaction term with correct sign
			for(int i = 0; i < numberOfCelltypes; ++i){
				char genotypeConfiguration = 0;
				genotypeConfiguration = bestFullModel.getGenotypeConfiguration().charAt(i);
				double interactionTermCurrentCelltype = estimateRegressionParameters[i+numberOfCelltypes];
				if (genotypeConfiguration == '0' || interactionTermCurrentCelltype == 0){
					results += "\t"+interactionTermCurrentCelltype;
				}else if(genotypeConfiguration == '1'){
						results += "\t-"+interactionTermCurrentCelltype;
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

			//results += "\t"+bestFullModel.getEstimatedStandardError();
			output.add(results);	
		}

		Path file = Paths.get(outputFolder+"/"+commandLineOptions.getOutfile());
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
	private void writePredictedExpression(List<DeconvolutionResult> deconvolutionResults) throws IOException, IllegalAccessException{
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
		/** 
		 * If roundDosage option is selected on the command line, round of the dosage to closest integer -> 0.49 = 0, 0.51 = 1, 1.51 = 2. 
		 */
		if (commandLineOptions.getRoundDosage()) {
			for (int i = 0; i < genotypes.length; ++i) {
				if (commandLineOptions.getRoundDosage()){
					genotypes[i] = Math.round(genotypes[i]);
				}
			}
		}

		InteractionModelCollection interactionModelCollection = new InteractionModelCollection(cellCounts, 
				commandLineOptions.getGenotypeConfigurationType());
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
		interactionModelCollection.createObservedValueMatricesFullModel(commandLineOptions.getAddGenotypeTerm());
		interactionModelCollection.findBestFullModel();		
		interactionModelCollection.createObservedValueMatricesCtModels(commandLineOptions.getAddGenotypeTerm());
		interactionModelCollection.findBestCtModel();
		calculateDeconvolutionPvalue(interactionModelCollection);

		double wholeBloodQTL = 0;
		double wholeBloodQTLpvalue = 0;
		if(commandLineOptions.getWholeBloodQTL()){
			// if true calculate spearman correlation between genotypes and expression values (i.e. whole blood eQTL)
			wholeBloodQTL = new SpearmansCorrelation().correlation(interactionModelCollection.getGenotypes(), interactionModelCollection.getExpessionValues());
			wholeBloodQTLpvalue = Statistics.calculateSpearmanTwoTailedPvalue(wholeBloodQTL, cellCounts.getNumberOfSamples());
		}
		DeconvolutionResult deconResult =  new DeconvolutionResult();

		interactionModelCollection.cleanUp(!commandLineOptions.getOutputPredictedExpression());
		deconResult = new DeconvolutionResult(interactionModelCollection, wholeBloodQTL, wholeBloodQTLpvalue);
		return deconResult;
	}


	/**
	 * get pvalue for each ctmodel
	 * 
	 * @param interactionModelCollection InteractionModelCollection object that has fullModel and ctModels for ANOVA comparison
	 */
	private void calculateDeconvolutionPvalue(InteractionModelCollection interactionModelCollection) 
			throws IllegalAccessException, IOException {
		for (int modelIndex = 0; modelIndex < cellCounts.getNumberOfCelltypes(); ++modelIndex) {
			String celltypeName = cellCounts.getCelltype(modelIndex);
			InteractionModel fullModel;

			fullModel = interactionModelCollection.getBestFullModel();

			InteractionModel ctModel = interactionModelCollection.getBestCtModel(celltypeName);
			double pval = Statistics.anova(fullModel.getSumOfSquares(), ctModel.getSumOfSquares(), 
					fullModel.getDegreesOfFreedom(),ctModel.getDegreesOfFreedom(), 
					true);

			ctModel.setPvalue(pval);
			interactionModelCollection.setPvalue(pval, ctModel.getCelltypeName());
			// TODO: why is this method called twice?
			interactionModelCollection.setPvalue(pval,interactionModelCollection
					.getBestCtModel(cellCounts.getCelltype(modelIndex)).getModelName());

		}

	}
}
