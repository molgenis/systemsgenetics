/*) 
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.MannWhitneyUTest2;

/**
 *
 * @author patri
 */
public class testCoregulationPerformance {

	public static void testCoreGenePredictionPerformance(Depict2Options options) throws IOException {

		DoubleMatrixDataset<String, String> gwasCoreGenePredictions = DoubleMatrixDataset.loadDoubleBinaryData(options.getGwasZscoreMatrixPath());
		List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();

		ArrayList<String> genesWithPrediciton = gwasCoreGenePredictions.getRowObjects();

		for (PathwayDatabase pathwayDatabase : pathwayDatabases) {

			final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase.getLocation());

			Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();

			final LinkedHashSet<String> sharedGenes = new LinkedHashSet<>();

			for (String gene : genesWithPrediciton) {
				if (pathwayGenes.contains(gene)) {
					sharedGenes.add(gene);
				}
			}

			final int sharedGenesCount = sharedGenes.size();

			final DoubleMatrixDataset<String, String> pathwayMatrix = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes);
			final DoubleMatrixDataset<String, String> gwasCoreGenePredictionsMatched = gwasCoreGenePredictions.viewRowSelection(sharedGenes);

			final DoubleMatrixDataset<String, String> outputMatrix = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), gwasCoreGenePredictionsMatched.getHashCols());

			try (ProgressBar pb = new ProgressBar(pathwayDatabase.getName(), pathwayMatrix.columns(), ProgressBarStyle.ASCII)) {

				for (int pathwayI = 0; pathwayI < pathwayMatrix.columns(); ++pathwayI) {

					final DoubleMatrix1D pathwayAnnotation = pathwayMatrix.getCol(pathwayI);
					final int annotatedGenesCount = pathwayAnnotation.cardinality();

					final int pathwayI2 = pathwayI;

					if (annotatedGenesCount < 10) {
						continue;
					}

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

						outputMatrix.setElementQuick(pathwayI2, traitI, uTest.getAuc());

					});
					
					pb.step();

				}

				outputMatrix.save(options.getOutputBasePath() + "_" + pathwayDatabase.getName() + ".txt");

			}
		}

	}

}
