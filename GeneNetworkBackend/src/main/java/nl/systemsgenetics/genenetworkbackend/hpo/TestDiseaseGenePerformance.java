/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import static nl.systemsgenetics.genenetworkbackend.ConvertHpoToMatrix.loadHgncToEnsgMap;
import static nl.systemsgenetics.genenetworkbackend.ConvertHpoToMatrix.loadNcbiToEnsgMap;
import static nl.systemsgenetics.genenetworkbackend.div.CalculateGenePredictability.loadSignificantTerms;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.biojava.nbio.ontology.Ontology;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class TestDiseaseGenePerformance {

	/**
	 * @param args the command line arguments
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws Exception {

		final File diseaseGeneHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt");
		final File ncbiToEnsgMapFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgNcbiId.txt");
		final File hgncToEnsgMapFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgnc.txt");
		final File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions.txt.gz");
		final File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_bonSigTerms.txt");
		final double correctedPCutoff = 0.05;
		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");
		final File hpoPredictionInfoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_auc_bonferroni.txt");
		final File hposToExcludeFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoToExclude.txt");
		final File skewnessFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\skewnessSummary.txt");
		final boolean randomize = true;
		final File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		final File backgroundForRandomize = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\Ensembl2Reactome_All_Levels.txt_genesInPathways.txt");
		final boolean randomizeCustomBackground = true;

		final File outputFile;
		final ArrayList<String> backgroundGenes;
		if (randomize) {

			if (randomizeCustomBackground) {
				backgroundGenes = loadBackgroundGenes(backgroundForRandomize);
				outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoDiseaseBenchmarkRandomizedCustomBackground.txt");
			} else {
				backgroundGenes = null;
				outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoDiseaseBenchmarkRandomized.txt");
			}

		} else {
			backgroundGenes = null;
			outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoDiseaseBenchmark.txt");
		}

		final HashMap<String, ArrayList<String>> ncbiToEnsgMap = loadNcbiToEnsgMap(ncbiToEnsgMapFile);
		final HashMap<String, ArrayList<String>> hgncToEnsgMap = loadHgncToEnsgMap(hgncToEnsgMapFile);
		final HashSet<String> exludedHpo = loadHpoExclude(hposToExcludeFile);

		

		final SkewnessInfo skewnessInfo = new SkewnessInfo(skewnessFile);

		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);

		DoubleMatrixDataset<String, String> predictionMatrix = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> predictionMatrixSignificant = predictionMatrix.viewColSelection(significantTerms);

		DoubleMatrixDataset<String, String> predictionMatrixSignificantCorrelationMatrix = predictionMatrixSignificant.calculateCorrelationMatrix();
		
		System.out.println("Done with correlation matrix");
		
		DiseaseGeneHpoData diseaseGeneHpoData = new DiseaseGeneHpoData(diseaseGeneHpoFile, ncbiToEnsgMap, hgncToEnsgMap, exludedHpo);
		if (randomize) {
			diseaseGeneHpoData = diseaseGeneHpoData.getPermutation(1, backgroundGenes, predictionMatrixSignificantCorrelationMatrix, 0.5);
		}
		
		DoubleMatrixDataset<String, String> annotationnMatrix = DoubleMatrixDataset.loadDoubleData(annotationMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> annotationMatrixSignificant = annotationnMatrix.viewColSelection(significantTerms);

		HashMap<String, MeanSd> hpoMeanSds = calculatePathayMeansOfAnnotatedGenes(predictionMatrixSignificant, annotationMatrixSignificant);

		Map<String, PredictionInfo> predictionInfo = HpoFinder.loadPredictionInfo(hpoPredictionInfoFile);

		Ontology hpoOntology = HpoFinder.loadHpoOntology(hpoOboFile);

		HpoFinder hpoFinder = new HpoFinder(hpoOntology, predictionInfo);

		final int totalGenes = predictionMatrixSignificant.rows();
		final double[] geneScores = new double[totalGenes];
		final NaturalRanking naturalRanking = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.MAXIMUM);

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[11];
		int c = 0;
		outputLine[c++] = "Disease";
		outputLine[c++] = "Gene";
		outputLine[c++] = "Rank";
		outputLine[c++] = "Z-score";
		outputLine[c++] = "HPO_skewness";
		outputLine[c++] = "Other_mean_skewness";
		outputLine[c++] = "Other_max_skewness";
		outputLine[c++] = "HPO_phenotypic_match_score";
		outputLine[c++] = "HPO_count";
		outputLine[c++] = "HPO_terms";
		outputLine[c++] = "HPO_terms_match_score";
		writer.writeNext(outputLine);

		for (DiseaseGeneHpoData.DiseaseGene diseaseGene : diseaseGeneHpoData.getDiseaseGeneHpos()) {

			String gene = diseaseGene.getGene();
			String disease = diseaseGene.getDisease();

			if (!predictionMatrixSignificant.containsRow(gene)) {
				continue;
			}

			Set<String> geneHpos = diseaseGeneHpoData.getDiseaseEnsgHpos(diseaseGene);

			LinkedHashSet<String> geneHposPredictable = new LinkedHashSet<>();

			for (String hpo : geneHpos) {
				geneHposPredictable.addAll(hpoFinder.getTermsToNames(hpoFinder.getPredictableTerms(hpo, correctedPCutoff)));
			}

			if (geneHposPredictable.isEmpty()) {
				continue;
			}

			DoubleMatrixDataset<String, String> predictionCaseTerms = predictionMatrixSignificant.viewColSelection(geneHposPredictable);
			DoubleMatrix2D predictionCaseTermsMatrix = predictionCaseTerms.getMatrix();

			double denominator = Math.sqrt(geneHposPredictable.size());

			for (int g = 0; g < totalGenes; ++g) {
				geneScores[g] = predictionCaseTermsMatrix.viewRow(g).zSum() / denominator;
				if (Double.isNaN(geneScores[g])) {
					geneScores[g] = 0;
				}
			}

			double[] geneRanks = naturalRanking.rank(geneScores);

			int diseaseGeneIndex = predictionMatrixSignificant.getRowIndex(gene);

			double zscore = geneScores[diseaseGeneIndex];
			double rank = (totalGenes - geneRanks[diseaseGeneIndex]) + 1;

			double hpoPhenotypicMatchScore = 0;
			StringBuilder individualMatchScore = new StringBuilder();
			boolean notFirst = false;
			int usedHpos = 0;

			for (String hpo : geneHpos) {

				if (!predictionMatrixSignificant.containsCol(hpo)) {
					continue;
				}

				usedHpos++;

				MeanSd hpoMeanSd = hpoMeanSds.get(hpo);

				double hpoPredictionZ = predictionMatrixSignificant.getElement(gene, hpo);

				double hpoPredictionOutlierScore = ((hpoPredictionZ - hpoMeanSd.getMean()) / hpoMeanSd.getSd());

				if (notFirst) {
					notFirst = true;
					individualMatchScore.append(';');
				}
				individualMatchScore.append(hpoPredictionOutlierScore);

				hpoPhenotypicMatchScore += hpoPredictionOutlierScore;

			}

			if (usedHpos == 0) {
				hpoPhenotypicMatchScore = Double.NaN;
			} else {
				hpoPhenotypicMatchScore = hpoPhenotypicMatchScore / usedHpos;
			}

			c = 0;
			outputLine[c++] = disease;
			outputLine[c++] = gene;
			outputLine[c++] = String.valueOf(rank);
			outputLine[c++] = String.valueOf(zscore);
			outputLine[c++] = String.valueOf(skewnessInfo.getHpoSkewness(gene));
			outputLine[c++] = String.valueOf(skewnessInfo.getMeanSkewnessExHpo(gene));
			outputLine[c++] = String.valueOf(skewnessInfo.getMaxSkewnessExHpo(gene));
			outputLine[c++] = String.valueOf(hpoPhenotypicMatchScore);
			outputLine[c++] = String.valueOf(geneHpos.size());
			outputLine[c++] = String.join(";", geneHpos);
			outputLine[c++] = individualMatchScore.toString();
			writer.writeNext(outputLine);

		}

		writer.close();

	}

	public static HashSet<String> loadHpoExclude(File hposToExclude) throws IOException {

		LinkedHashSet<String> hpos = new LinkedHashSet<>();

		BufferedReader reader = new BufferedReader(new FileReader(hposToExclude));

		String line;
		while ((line = reader.readLine()) != null) {
			hpos.add(line);
		}

		return hpos;

	}

	private static HashMap<String, MeanSd> calculatePathayMeansOfAnnotatedGenes(DoubleMatrixDataset<String, String> predictionMatrixSignificant, DoubleMatrixDataset<String, String> annotationMatrixSignificant) {

		HashMap<String, MeanSd> pathwayMeanSdMap = new HashMap<>(predictionMatrixSignificant.columns());

		Mean meanCalculator = new Mean();
		Variance varianceCalculator = new Variance();

		for (String pathway : predictionMatrixSignificant.getColObjects()) {

			meanCalculator.clear();
			varianceCalculator.clear();

			DoubleMatrix1D pathwayPredictions = predictionMatrixSignificant.getCol(pathway);
			DoubleMatrix1D pathwayAnnotations = annotationMatrixSignificant.getCol(pathway);

			for (int g = 0; g < pathwayPredictions.size(); ++g) {
				if (pathwayAnnotations.get(g) != 0) {
					meanCalculator.increment(pathwayPredictions.getQuick(g));
					varianceCalculator.increment(pathwayPredictions.getQuick(g));
				}
			}

			double v = varianceCalculator.getResult();

			pathwayMeanSdMap.put(pathway, new MeanSd(meanCalculator.getResult(), v * v));

		}

		return pathwayMeanSdMap;

	}

	private static class MeanSd {

		private final double mean;
		private final double sd;

		public MeanSd(double mean, double sd) {
			this.mean = mean;
			this.sd = sd;
		}

		public double getMean() {
			return mean;
		}

		public double getSd() {
			return sd;
		}

	}

	public static ArrayList<String> loadBackgroundGenes(File backgroundForRandomize) throws IOException {

		ArrayList<String> genes = new ArrayList<>();

		BufferedReader reader = new BufferedReader(new FileReader(backgroundForRandomize));

		String line;
		while ((line = reader.readLine()) != null) {
			genes.add(line);
		}

		return genes;

	}

}
