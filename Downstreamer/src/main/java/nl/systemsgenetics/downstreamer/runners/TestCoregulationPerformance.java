/*) 
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.io.CoregeneEnrichmentExcelWriter;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import static nl.systemsgenetics.downstreamer.runners.DownstreamerUtilities.loadExistingStep2Results;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.MannWhitneyUTest2;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class TestCoregulationPerformance {

	private static final Logger LOGGER = Logger.getLogger(TestCoregulationPerformance.class);

	public static void testCoreGenePredictionPerformance(DownstreamerOptions options) throws IOException, Exception {

		DownstreamerStep2Results step2 = loadExistingStep2Results(options);

		List<PathwayDatabase> pathwayDatabases2 = options.getPathwayDatabases2();
		testPredictionsGenePvalues(step2.getGenePvalues(), pathwayDatabases2, options, "gwasGenePvalues");

		for (PathwayEnrichments step2Enrichment : step2.getPathwayEnrichments()) {
			testPredictions(step2Enrichment, pathwayDatabases2, options, step2Enrichment.getPathwayDatabase().getName());
		}

	}

	private static void testPredictions(PathwayEnrichments step2Enrichment, List<PathwayDatabase> pathwayDatabases2, DownstreamerOptions options, String predictionSource) throws IOException {

		DoubleMatrixDataset<String, String> predictionZscores = step2Enrichment.getEnrichmentZscores();
		DoubleMatrixDataset<String, String> predictionPvalues = step2Enrichment.getpValues();
		DoubleMatrixDataset<String, String> predictionQvalues = step2Enrichment.getqValues();

		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2Auc = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2Utest = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfOdds = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfFisherP = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2FdrOdds = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2FdrFisherP = new HashMap<>();

		ArrayList<String> genesWithPrediciton = predictionZscores.getRowObjects();

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
			final double bonfSigThreshold = 0.05d / sharedGenesCount;
			final double fdrSigThreshold = 0.05;

			final DoubleMatrixDataset<String, String> pathwayMatrix;

			{
				final DoubleMatrixDataset<String, String> pathwayMatrixUnfiltered = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes);
				//Filter on pathways with atleast 10 genes
				ArrayList<String> allPathways = pathwayMatrixUnfiltered.getColObjects();
				List<String> includedPathways = Collections.synchronizedList(new ArrayList<String>());

				IntStream.range(0, pathwayMatrixUnfiltered.columns()).parallel().forEach(pathwayI -> {
					if (pathwayMatrixUnfiltered.getCol(pathwayI).cardinality() >= 10) {
						includedPathways.add(allPathways.get(pathwayI));
					}
				});
				pathwayMatrix = pathwayMatrixUnfiltered.viewColSelection(includedPathways);

			}

			final DoubleMatrixDataset<String, String> predictionZscoresMatched = predictionZscores.viewRowSelection(sharedGenes);
			final DoubleMatrixDataset<String, String> predictionPvaluesMatched = predictionPvalues.viewRowSelection(sharedGenes);
			final DoubleMatrixDataset<String, String> predictionQvaluesMatched = predictionQvalues.viewRowSelection(sharedGenes);

			final DoubleMatrixDataset<String, String> outputMatrixAuc = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> outputMatrixPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());

			final DoubleMatrixDataset<String, String> bonfExactPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> fdrExactPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> bonfOdds = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> fdrOdds = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());

			final ArrayList<String> pathwayNames = pathwayMatrix.getColObjects();

			try (ProgressBar pb = new ProgressBar(predictionSource + "_" + pathwayDatabase2.getName(), pathwayMatrix.columns(), ProgressBarStyle.ASCII)) {

				IntStream.range(0, pathwayMatrix.columns()).parallel().forEach(pathwayI -> {

					final DoubleMatrix1D pathwayAnnotation = pathwayMatrix.getCol(pathwayI);
					final int annotatedGenesCount = pathwayAnnotation.cardinality();

					for (int traitI = 0; traitI < predictionZscoresMatched.columns(); ++traitI) {

						final MannWhitneyUTest2 uTest = new MannWhitneyUTest2();

						final DoubleMatrix1D coreGeneZScores = predictionZscoresMatched.getCol(traitI);
						final DoubleMatrix1D coreGenePvalues = predictionPvaluesMatched.getCol(traitI);
						final DoubleMatrix1D coreGeneQvalues = predictionQvaluesMatched.getCol(traitI);

						final double[] coreGeneScoresAnnotatedGenes = new double[annotatedGenesCount];
						final double[] coreGeneScoresOtherGenes = new double[sharedGenesCount - annotatedGenesCount];

						int x = 0;
						int y = 0;

						int inPathwayBonfSig = 0;
						int inPathwayNotBonfSig = 0;
						int notPathwayBonfSig = 0;
						int notPathwayNotBonfSig = 0;

						int inPathwayFdrSig = 0;
						int inPathwayNotFdrSig = 0;
						int notPathwayFdrSig = 0;
						int notPathwayNotFdrSig = 0;

						for (int g = 0; g < sharedGenesCount; g++) {
							if (pathwayAnnotation.getQuick(g) > 0) {
								//Gene is annotated to pathway
								coreGeneScoresAnnotatedGenes[x++] = coreGeneZScores.get(g);

								if (coreGeneZScores.get(g) >= 0 && coreGenePvalues.get(g) <= bonfSigThreshold) {
									inPathwayBonfSig++;
								} else {
									inPathwayNotBonfSig++;
								}

								if (coreGeneZScores.get(g) >= 0 && coreGeneQvalues.get(g) <= fdrSigThreshold) {
									inPathwayFdrSig++;
								} else {
									inPathwayNotFdrSig++;
								}

							} else {
								//Gene is not annoated to pathway
								coreGeneScoresOtherGenes[y++] = coreGeneZScores.get(g);

								if (coreGeneZScores.get(g) >= 0 && coreGenePvalues.get(g) <= bonfSigThreshold) {
									notPathwayBonfSig++;
								} else {
									notPathwayNotBonfSig++;
								}

								if (coreGeneZScores.get(g) >= 0 && coreGeneQvalues.get(g) <= fdrSigThreshold) {
									notPathwayFdrSig++;
								} else {
									notPathwayNotFdrSig++;
								}
							}
						}

						uTest.setData(coreGeneScoresAnnotatedGenes, coreGeneScoresOtherGenes);

						final double auc = uTest.getAuc();
						final double pval = uTest.getP();

						if (LOGGER.isDebugEnabled() && Double.isNaN(auc)) {
							LOGGER.debug(pathwayNames.get(pathwayI) + " NaN AUC " + Arrays.toString(coreGeneScoresAnnotatedGenes));
						}

						FisherExactTest ft = new FisherExactTest();
						//first do this before single sided has a value;
						ft.getFisherPValue(inPathwayBonfSig, inPathwayNotBonfSig, notPathwayBonfSig, notPathwayNotBonfSig);
						final double bonFp = ft.getFisherRightTail();
						final double bonOr = (double) (inPathwayBonfSig * notPathwayNotBonfSig) / (double) (inPathwayNotBonfSig * notPathwayBonfSig);

						ft.getFisherPValue(inPathwayFdrSig, inPathwayNotFdrSig, notPathwayFdrSig, notPathwayNotFdrSig);
						final double fdrFp = ft.getFisherRightTail();
						final double fdrOr = (double) (inPathwayFdrSig * notPathwayNotFdrSig) / (double) (inPathwayNotFdrSig * notPathwayFdrSig);

						bonfOdds.setElementQuick(pathwayI, traitI, bonOr);
						bonfExactPvalues.setElementQuick(pathwayI, traitI, bonFp);
						fdrOdds.setElementQuick(pathwayI, traitI, fdrOr);
						fdrExactPvalues.setElementQuick(pathwayI, traitI, fdrFp);
						outputMatrixAuc.setElementQuick(pathwayI, traitI, auc);
						outputMatrixPvalues.setElementQuick(pathwayI, traitI, pval);

					}

					pb.step();

				});

				pathwayDatabase2BonfFisherP.put(pathwayDatabase2.getName(), bonfExactPvalues);
				pathwayDatabase2BonfOdds.put(pathwayDatabase2.getName(), bonfOdds);
				pathwayDatabase2FdrFisherP.put(pathwayDatabase2.getName(), fdrExactPvalues);
				pathwayDatabase2FdrOdds.put(pathwayDatabase2.getName(), fdrOdds);
				pathwayDatabase2Auc.put(pathwayDatabase2.getName(), outputMatrixAuc);
				pathwayDatabase2Utest.put(pathwayDatabase2.getName(), outputMatrixPvalues);

				bonfExactPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactPvals_" + pathwayDatabase2.getName() + ".txt");
				fdrExactPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_fdrFexactPvals_" + pathwayDatabase2.getName() + ".txt");
				bonfOdds.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactOr_" + pathwayDatabase2.getName() + ".txt");
				fdrOdds.save(options.getOutputBasePath() + "_" + predictionSource + "_fdrFexactOr_" + pathwayDatabase2.getName() + ".txt");
				outputMatrixAuc.save(options.getOutputBasePath() + "_" + predictionSource + "_auc_" + pathwayDatabase2.getName() + ".txt");
				outputMatrixPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_utestPvals_" + pathwayDatabase2.getName() + ".txt");
			}
		}

		CoregeneEnrichmentExcelWriter.write(options, pathwayDatabase2Auc, pathwayDatabase2Utest, pathwayDatabase2BonfOdds, pathwayDatabase2BonfFisherP, pathwayDatabase2FdrOdds, pathwayDatabase2FdrFisherP, predictionPvalues.getColObjects(), pathwayDatabases2);

	}

	private static void testPredictionsGenePvalues(DoubleMatrixDataset<String, String> predictionMatrix, List<PathwayDatabase> pathwayDatabases2, DownstreamerOptions options, String predictionSource) throws IOException {
		ArrayList<String> genesWithPrediciton = predictionMatrix.getRowObjects();

		for (PathwayDatabase pathwayDatabase2 : pathwayDatabases2) {

			final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase2.getLocation());

			Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();

			final LinkedHashSet<String> sharedGenes = new LinkedHashSet<>();
			final double bonfSigThreshold = 0.05d / sharedGenes.size();

			for (String gene : genesWithPrediciton) {
				if (pathwayGenes.contains(gene)) {
					sharedGenes.add(gene);
				}
			}

			final int sharedGenesCount = sharedGenes.size();

			final DoubleMatrixDataset<String, String> pathwayMatrix;

			{
				final DoubleMatrixDataset<String, String> pathwayMatrixUnfiltered = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes);
				//Filter on pathways with atleast 10 genes
				ArrayList<String> allPathways = pathwayMatrixUnfiltered.getColObjects();
				List<String> includedPathways = Collections.synchronizedList(new ArrayList<String>());

				IntStream.range(0, pathwayMatrixUnfiltered.columns()).parallel().forEach(pathwayI -> {
					if (pathwayMatrixUnfiltered.getCol(pathwayI).cardinality() >= 10) {
						includedPathways.add(allPathways.get(pathwayI));
					}
				});
				pathwayMatrix = pathwayMatrixUnfiltered.viewColSelection(includedPathways);

			}

			final DoubleMatrixDataset<String, String> gwasCoreGenePredictionsMatched = predictionMatrix.viewRowSelection(sharedGenes);

			final DoubleMatrixDataset<String, String> outputMatrixAuc = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), gwasCoreGenePredictionsMatched.getHashCols());
			final DoubleMatrixDataset<String, String> outputMatrixPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), gwasCoreGenePredictionsMatched.getHashCols());

			final DoubleMatrixDataset<String, String> bonfExactPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), gwasCoreGenePredictionsMatched.getHashCols());
			final DoubleMatrixDataset<String, String> bonfOdds = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), gwasCoreGenePredictionsMatched.getHashCols());

			final ArrayList<String> pathwayNames = pathwayMatrix.getColObjects();

			try (ProgressBar pb = new ProgressBar(predictionSource + "_" + pathwayDatabase2.getName(), pathwayMatrix.columns(), ProgressBarStyle.ASCII)) {

				for (int pathwayI = 0; pathwayI < pathwayMatrix.columns(); ++pathwayI) {

					final DoubleMatrix1D pathwayAnnotation = pathwayMatrix.getCol(pathwayI);
					final int annotatedGenesCount = pathwayAnnotation.cardinality();

					final int pathwayI2 = pathwayI;

//					if (annotatedGenesCount < 10) {
//						if (LOGGER.isDebugEnabled()) {
//							LOGGER.debug("Skipping " + pathwayNames.get(pathwayI) + " " + annotatedGenesCount + " genes annotated");
//						}
//						pb.step();
//						continue;
//					}
					final int pathwayI_2 = pathwayI;

					IntStream.range(0, gwasCoreGenePredictionsMatched.columns()).parallel().forEach(traitI -> {

						final MannWhitneyUTest2 uTest = new MannWhitneyUTest2();

						final DoubleMatrix1D gwasGeneZscore = gwasCoreGenePredictionsMatched.getCol(traitI);

						final double[] coreGeneScoresAnnotatedGenes = new double[annotatedGenesCount];
						final double[] coreGeneScoresOtherGenes = new double[sharedGenesCount - annotatedGenesCount];

						int x = 0;
						int y = 0;

						int inPathwayBonfSig = 0;
						int inPathwayNotBonfSig = 0;
						int notPathwayBonfSig = 0;
						int notPathwayNotBonfSig = 0;

						for (int g = 0; g < sharedGenesCount; g++) {
							if (pathwayAnnotation.getQuick(g) > 0) {

								if (gwasGeneZscore.get(g) <= bonfSigThreshold) {
									inPathwayBonfSig++;
								} else {
									inPathwayNotBonfSig++;
								}

								coreGeneScoresAnnotatedGenes[x++] = gwasGeneZscore.get(g);
							} else {

								if (gwasGeneZscore.get(g) <= bonfSigThreshold) {
									notPathwayBonfSig++;
								} else {
									notPathwayNotBonfSig++;
								}

								coreGeneScoresOtherGenes[y++] = gwasGeneZscore.get(g);
							}
						}

						uTest.setData(coreGeneScoresAnnotatedGenes, coreGeneScoresOtherGenes);

						final double auc = uTest.getAuc();
						final double pval = uTest.getP();

						if (LOGGER.isDebugEnabled() && Double.isNaN(auc)) {
							LOGGER.debug(pathwayNames.get(pathwayI_2) + " NaN AUC " + Arrays.toString(coreGeneScoresAnnotatedGenes));
						}

						FisherExactTest ft = new FisherExactTest();
						//first do this before single sided has a value;
						ft.getFisherPValue(inPathwayBonfSig, inPathwayNotBonfSig, notPathwayBonfSig, notPathwayNotBonfSig);
						final double bonFp = ft.getFisherRightTail();
						final double bonOr = (double) (inPathwayBonfSig * notPathwayNotBonfSig) / (double) (inPathwayNotBonfSig * notPathwayBonfSig);

						bonfOdds.setElementQuick(pathwayI2, traitI, bonOr);
						bonfExactPvalues.setElementQuick(pathwayI2, traitI, bonFp);

						outputMatrixAuc.setElementQuick(pathwayI2, traitI, auc);
						outputMatrixPvalues.setElementQuick(pathwayI2, traitI, pval);

					});

					pb.step();

				}

				bonfExactPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactPvals_" + pathwayDatabase2.getName() + ".txt");
				bonfOdds.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactOr_" + pathwayDatabase2.getName() + ".txt");
				outputMatrixAuc.save(options.getOutputBasePath() + "_" + predictionSource + "_auc_" + pathwayDatabase2.getName() + ".txt");
				outputMatrixPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_utestPvals_" + pathwayDatabase2.getName() + ".txt");
			}
		}
	}

}
