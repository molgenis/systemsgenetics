/*) 
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.runners;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.depict2.Depict2Options;
import nl.systemsgenetics.depict2.Depict2Step2Results;
import nl.systemsgenetics.depict2.pathway.PathwayDatabase;
import nl.systemsgenetics.depict2.pathway.PathwayEnrichments;
import static nl.systemsgenetics.depict2.runners.Depict2Utilities.loadExistingStep2Restuls;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.MannWhitneyUTest2;

/**
 *
 * @author patri
 */
public class TestCoregulationPerformance {
	
	private static final Logger LOGGER = Logger.getLogger(TestCoregulationPerformance.class);

	public static void testCoreGenePredictionPerformance(Depict2Options options) throws IOException, Exception {
		
		Depict2Step2Results step2 = loadExistingStep2Restuls(options);

		List<PathwayDatabase> pathwayDatabases2 = options.getPathwayDatabases2();
		
		testPredictions(step2.getGenePvalues(), pathwayDatabases2, options, "genePvalues");
		
		for (PathwayEnrichments step2Enrichment : step2.getPathwayEnrichments()) {
			testPredictions( step2Enrichment.getEnrichmentZscores(), pathwayDatabases2, options, step2Enrichment.getPathwayDatabase().getName());
		}
		

	}

	private static void testPredictions(DoubleMatrixDataset<String, String> predictionMatrix, List<PathwayDatabase> pathwayDatabases2, Depict2Options options, String predictionSource) throws IOException {
		ArrayList<String> genesWithPrediciton = predictionMatrix.getRowObjects();
		
		for (PathwayDatabase pathwayDatabase2 : pathwayDatabases2) {
			
			final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase2.getLocation());

			Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();

			final LinkedHashSet<String> sharedGenes = new LinkedHashSet<>();

			for (String gene : genesWithPrediciton) {
				if (pathwayGenes.contains(gene)) {
					sharedGenes.add(gene);
				}
			}

			final int sharedGenesCount = sharedGenes.size();

			final DoubleMatrixDataset<String, String> pathwayMatrix = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes);
			final DoubleMatrixDataset<String, String> gwasCoreGenePredictionsMatched = predictionMatrix.viewRowSelection(sharedGenes);

			final DoubleMatrixDataset<String, String> outputMatrixAuc = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), gwasCoreGenePredictionsMatched.getHashCols());
			final DoubleMatrixDataset<String, String> outputMatrixPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), gwasCoreGenePredictionsMatched.getHashCols());

			final ArrayList<String> pathwayNames = pathwayMatrix.getColObjects();
			
			try (ProgressBar pb = new ProgressBar(predictionSource + "_" + pathwayDatabase2.getName(), pathwayMatrix.columns(), ProgressBarStyle.ASCII)) {

				for (int pathwayI = 0; pathwayI < pathwayMatrix.columns(); ++pathwayI) {

					final DoubleMatrix1D pathwayAnnotation = pathwayMatrix.getCol(pathwayI);
					final int annotatedGenesCount = pathwayAnnotation.cardinality();

					final int pathwayI2 = pathwayI;

					if (annotatedGenesCount < 10) {
						if(LOGGER.isDebugEnabled()){
							LOGGER.debug("Skipping " + pathwayNames.get(pathwayI) + " " + annotatedGenesCount + " genes annotated");
						}
						pb.step();
						continue;
					}
					
					final int pathwayI_2 = pathwayI;

					IntStream.range(0, gwasCoreGenePredictionsMatched.columns()).parallel().forEach(traitI -> {

						final MannWhitneyUTest2 uTest = new MannWhitneyUTest2();

						final DoubleMatrix1D coreGeneScore = gwasCoreGenePredictionsMatched.getCol(traitI);

						final double[] coreGeneScoresAnnotatedGenes = new double[annotatedGenesCount];
						final double[] coreGeneScoresOtherGenes = new double[sharedGenesCount - annotatedGenesCount];

						int x = 0;
						int y = 0;

						for (int g = 0; g < sharedGenesCount; g++) {
							if (pathwayAnnotation.getQuick(g) > 0) {
								coreGeneScoresAnnotatedGenes[x++] = coreGeneScore.get(g);
							} else {
								coreGeneScoresOtherGenes[y++] = coreGeneScore.get(g);
							}
						}

						uTest.setData(coreGeneScoresAnnotatedGenes, coreGeneScoresOtherGenes);

						final double auc = uTest.getAuc();
						final double pval = uTest.getP();

						if(LOGGER.isDebugEnabled() && Double.isNaN(auc)){
							LOGGER.debug(pathwayNames.get(pathwayI_2) + " NaN AUC " + Arrays.toString(coreGeneScoresAnnotatedGenes));
						}

						outputMatrixAuc.setElementQuick(pathwayI2, traitI, auc);
						outputMatrixPvalues.setElementQuick(pathwayI2, traitI, pval);

					});
					
					pb.step();

				}

				outputMatrixAuc.save(options.getOutputBasePath() + "_" + predictionSource + "_auc_" + pathwayDatabase2.getName() + ".txt");
				outputMatrixPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_aucPvals_" + pathwayDatabase2.getName() + ".txt");
			}
		}
	}

}
