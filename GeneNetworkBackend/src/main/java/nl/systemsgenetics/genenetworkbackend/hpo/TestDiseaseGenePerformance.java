/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

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
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoDiseaseBenchmark.txt");

		final HashMap<String, ArrayList<String>> ncbiToEnsgMap = loadNcbiToEnsgMap(ncbiToEnsgMapFile);
		final HashMap<String, ArrayList<String>> hgncToEnsgMap = loadHgncToEnsgMap(hgncToEnsgMapFile);
		final HashSet<String> exludedHpo = loadHpoExclude(hposToExcludeFile);

		final DiseaseGeneHpoData diseaseGeneHpoData = new DiseaseGeneHpoData(diseaseGeneHpoFile, ncbiToEnsgMap, hgncToEnsgMap, exludedHpo);

		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);

		DoubleMatrixDataset<String, String> predictionMatrix = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> predictionMatrixSignificant = predictionMatrix.viewColSelection(significantTerms);

		Map<String, PredictionInfo> predictionInfo = HpoFinder.loadPredictionInfo(hpoPredictionInfoFile);

		Ontology hpoOntology = HpoFinder.loadHpoOntology(hpoOboFile);

		HpoFinder hpoFinder = new HpoFinder(hpoOntology, predictionInfo);

		final int totalGenes = predictionMatrixSignificant.rows();
		final double[] geneScores = new double[totalGenes];
		final NaturalRanking naturalRanking = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.MAXIMUM);

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[5];
		int c = 0;
		outputLine[c++] = "Disease";
		outputLine[c++] = "Gene";
		outputLine[c++] = "Rank";
		outputLine[c++] = "HPO_count";
		outputLine[c++] = "HPO_terms";
		writer.writeNext(outputLine);

		for ( DiseaseGeneHpoData.DiseaseGene diseaseGene : diseaseGeneHpoData.getDiseaseGeneHpos()) {

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

			double rank = (totalGenes - geneRanks[diseaseGeneIndex]) + 1;

			c = 0;
			outputLine[c++] = disease;
			outputLine[c++] = gene;
			outputLine[c++] = String.valueOf(rank);
			outputLine[c++] = String.valueOf(geneHpos.size());
			outputLine[c++] = String.join(";", geneHpos);
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

}
