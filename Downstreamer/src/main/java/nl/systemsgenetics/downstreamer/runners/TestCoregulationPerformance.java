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
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.io.CoregeneEnrichmentExcelWriter;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.MannWhitneyUTest2;

/**
 * @author patri
 */
@Deprecated
public class TestCoregulationPerformance {

	private static final Logger LOGGER = LogManager.getLogger(TestCoregulationPerformance.class);

	public static void testCoreGenePredictionPerformance(DownstreamerOptionsDeprecated options) throws IOException, Exception {

		DownstreamerStep2Results step2 = IoUtils.loadExistingStep2Results(options);

		List<PathwayDatabase> pathwayDatabases2 = options.getPathwayDatabases2();
		testPredictionsGenePvalues(step2.getGenePvalues(), pathwayDatabases2, options, "gwasGenePvalues");

		for (PathwayEnrichments step2Enrichment : step2.getPathwayEnrichments()) {
			testPredictions(step2Enrichment, pathwayDatabases2, options, step2Enrichment.getPathwayDatabase().getName());
		}

	}

	private static void testPredictions(PathwayEnrichments step2Enrichment, List<PathwayDatabase> pathwayDatabases2, DownstreamerOptionsDeprecated options, String predictionSource) throws IOException {

		DoubleMatrixDataset<String, String> predictionZscores = step2Enrichment.getEnrichmentZscores();
		DoubleMatrixDataset<String, String> predictionPvalues = step2Enrichment.getpValues();
		DoubleMatrixDataset<String, String> predictionQvalues = step2Enrichment.getqValues();

		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2Auc = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2Utest = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfOverlap = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfOdds = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfFisherP = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfCisOverlap = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfCisOdds = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfCisFisherP = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfTransOverlap = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfTransOdds = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2BonfTransFisherP = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2FdrOdds = new HashMap<>();
		HashMap<String, DoubleMatrixDataset<String, String>> pathwayDatabase2FdrFisherP = new HashMap<>();

		ArrayList<String> genesWithPrediciton = predictionZscores.getRowObjects();

		HashMap<String, HashMap<String, DownstreamerUtilities.NearestVariant>> distanceGeneToTopCisSnpPerTrait = DownstreamerUtilities.getDistanceGeneToTopCisSnpPerTrait(options);

		for (PathwayDatabase pathwayDatabase2 : pathwayDatabases2) {

			final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase2.getLocation());

			Set<String> pathwayGenes = pathwayMatrixLoader.getAllRowIdentifiers();

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
			final DoubleMatrixDataset<String, String> bonfExactOverlap = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> bonfOdds = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());

			final DoubleMatrixDataset<String, String> bonfCisExactPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> bonfCisExactOverlap = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> bonfCisOdds = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());

			final DoubleMatrixDataset<String, String> bonfTransExactPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> bonfTransExactOverlap = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> bonfTransOdds = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());


			final DoubleMatrixDataset<String, String> fdrExactPvalues = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());
			final DoubleMatrixDataset<String, String> fdrOdds = new DoubleMatrixDataset<>(pathwayMatrix.getHashCols(), predictionZscoresMatched.getHashCols());

			final ArrayList<String> pathwayNames = pathwayMatrix.getColObjects();
			final ArrayList<String> traits = predictionZscoresMatched.getColObjects();
			final ArrayList<String> geneOrder = pathwayMatrix.getRowObjects();

			try (ProgressBar pb = new ProgressBar(predictionSource + "_" + pathwayDatabase2.getName(), pathwayMatrix.columns(), ProgressBarStyle.ASCII)) {

				IntStream.range(0, pathwayMatrix.columns()).parallel().forEach(pathwayI -> {

					final DoubleMatrix1D pathwayAnnotation = pathwayMatrix.getCol(pathwayI);
					final int annotatedGenesCount = pathwayAnnotation.cardinality();

					for (int traitI = 0; traitI < predictionZscoresMatched.columns(); ++traitI) {

						String trait = traits.get(traitI);

						HashMap<String, DownstreamerUtilities.NearestVariant> distanceGeneToTopCisSnpPer = distanceGeneToTopCisSnpPerTrait.get(trait);

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

						int inPathwayCisBonfSig = 0;
						int inPathwayCisNotBonfSig = 0;
						int notPathwayCisBonfSig = 0;
						int notPathwayCisNotBonfSig = 0;

						int inPathwayTransBonfSig = 0;
						int inPathwayTransNotBonfSig = 0;
						int notPathwayTransBonfSig = 0;
						int notPathwayTransNotBonfSig = 0;

						int inPathwayFdrSig = 0;
						int inPathwayNotFdrSig = 0;
						int notPathwayFdrSig = 0;
						int notPathwayNotFdrSig = 0;

						for (int g = 0; g < sharedGenesCount; g++) {

							String gene = geneOrder.get(g);

							boolean cisGene = distanceGeneToTopCisSnpPer.containsKey(gene) && distanceGeneToTopCisSnpPer.get(gene).getDistance() >= 0;

							if (pathwayAnnotation.getQuick(g) > 0) {
								//Gene is annotated to pathway
								coreGeneScoresAnnotatedGenes[x++] = coreGeneZScores.get(g);

								if (coreGeneZScores.get(g) >= 0 && coreGenePvalues.get(g) <= bonfSigThreshold) {
									inPathwayBonfSig++;
									if (cisGene) {
										inPathwayCisBonfSig++;
									} else {
										inPathwayTransBonfSig++;
									}
								} else {
									inPathwayNotBonfSig++;
									if (cisGene) {
										inPathwayCisNotBonfSig++;
									} else {
										inPathwayTransNotBonfSig++;
									}
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
									if (cisGene) {
										notPathwayCisBonfSig++;
									} else {
										notPathwayTransBonfSig++;
									}
								} else {
									notPathwayNotBonfSig++;
									if (cisGene) {
										notPathwayCisNotBonfSig++;
									} else {
										notPathwayTransNotBonfSig++;
									}
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

						ft.getFisherPValue(inPathwayCisBonfSig, inPathwayCisNotBonfSig, notPathwayCisBonfSig, notPathwayCisNotBonfSig);
						final double bonCisFp = ft.getFisherRightTail();
						final double bonCisOr = (double) (inPathwayCisBonfSig * notPathwayCisNotBonfSig) / (double) (inPathwayCisNotBonfSig * notPathwayCisBonfSig);

						ft.getFisherPValue(inPathwayTransBonfSig, inPathwayTransNotBonfSig, notPathwayTransBonfSig, notPathwayTransNotBonfSig);
						final double bonTransFp = ft.getFisherRightTail();
						final double bonTransOr = (double) (inPathwayTransBonfSig * notPathwayTransNotBonfSig) / (double) (inPathwayTransNotBonfSig * notPathwayTransBonfSig);

						ft.getFisherPValue(inPathwayFdrSig, inPathwayNotFdrSig, notPathwayFdrSig, notPathwayNotFdrSig);
						final double fdrFp = ft.getFisherRightTail();
						final double fdrOr = (double) (inPathwayFdrSig * notPathwayNotFdrSig) / (double) (inPathwayNotFdrSig * notPathwayFdrSig);

						bonfOdds.setElementQuick(pathwayI, traitI, bonOr);
						bonfExactPvalues.setElementQuick(pathwayI, traitI, bonFp);
						bonfExactOverlap.setElementQuick(pathwayI, traitI, inPathwayBonfSig);

						bonfCisOdds.setElementQuick(pathwayI, traitI, bonCisOr);
						bonfCisExactPvalues.setElementQuick(pathwayI, traitI, bonCisFp);
						bonfCisExactOverlap.setElementQuick(pathwayI, traitI, inPathwayCisBonfSig);

						bonfTransOdds.setElementQuick(pathwayI, traitI, bonTransOr);
						bonfTransExactPvalues.setElementQuick(pathwayI, traitI, bonTransFp);
						bonfTransExactOverlap.setElementQuick(pathwayI, traitI, inPathwayTransBonfSig);

						fdrOdds.setElementQuick(pathwayI, traitI, fdrOr);
						fdrExactPvalues.setElementQuick(pathwayI, traitI, fdrFp);
						outputMatrixAuc.setElementQuick(pathwayI, traitI, auc);
						outputMatrixPvalues.setElementQuick(pathwayI, traitI, pval);

					}

					pb.step();

				});

				pathwayDatabase2BonfFisherP.put(pathwayDatabase2.getName(), bonfExactPvalues);
				pathwayDatabase2BonfOverlap.put(pathwayDatabase2.getName(), bonfExactOverlap);
				pathwayDatabase2BonfOdds.put(pathwayDatabase2.getName(), bonfOdds);

				pathwayDatabase2BonfCisFisherP.put(pathwayDatabase2.getName(), bonfCisExactPvalues);
				pathwayDatabase2BonfCisOverlap.put(pathwayDatabase2.getName(), bonfCisExactOverlap);
				pathwayDatabase2BonfCisOdds.put(pathwayDatabase2.getName(), bonfCisOdds);

				pathwayDatabase2BonfTransFisherP.put(pathwayDatabase2.getName(), bonfTransExactPvalues);
				pathwayDatabase2BonfTransOverlap.put(pathwayDatabase2.getName(), bonfTransExactOverlap);
				pathwayDatabase2BonfTransOdds.put(pathwayDatabase2.getName(), bonfTransOdds);

				pathwayDatabase2FdrFisherP.put(pathwayDatabase2.getName(), fdrExactPvalues);
				pathwayDatabase2FdrOdds.put(pathwayDatabase2.getName(), fdrOdds);
				pathwayDatabase2Auc.put(pathwayDatabase2.getName(), outputMatrixAuc);
				pathwayDatabase2Utest.put(pathwayDatabase2.getName(), outputMatrixPvalues);

				bonfExactPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactPvals_" + pathwayDatabase2.getName() + ".txt.gz");
				fdrExactPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_fdrFexactPvals_" + pathwayDatabase2.getName() + ".txt.gz");
				bonfOdds.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactOr_" + pathwayDatabase2.getName() + ".txt.gz");
				fdrOdds.save(options.getOutputBasePath() + "_" + predictionSource + "_fdrFexactOr_" + pathwayDatabase2.getName() + ".txt.gz");
				outputMatrixAuc.save(options.getOutputBasePath() + "_" + predictionSource + "_auc_" + pathwayDatabase2.getName() + ".txt.gz");
				outputMatrixPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_utestPvals_" + pathwayDatabase2.getName() + ".txt.gz");
			}
		}

		CoregeneEnrichmentExcelWriter.write(options, pathwayDatabase2Auc, pathwayDatabase2Utest,
				pathwayDatabase2BonfOverlap, pathwayDatabase2BonfOdds, pathwayDatabase2BonfFisherP,
				pathwayDatabase2BonfCisOverlap, pathwayDatabase2BonfCisOdds, pathwayDatabase2BonfCisFisherP,
				pathwayDatabase2BonfTransOverlap, pathwayDatabase2BonfTransOdds, pathwayDatabase2BonfTransFisherP,
				pathwayDatabase2FdrOdds, pathwayDatabase2FdrFisherP, predictionPvalues.getColObjects(), pathwayDatabases2, predictionSource);

	}

	private static void testPredictionsGenePvalues(DoubleMatrixDataset<String, String> predictionMatrix, List<PathwayDatabase> pathwayDatabases2, DownstreamerOptionsDeprecated options, String predictionSource) throws IOException {
		ArrayList<String> genesWithPrediciton = predictionMatrix.getRowObjects();

		for (PathwayDatabase pathwayDatabase2 : pathwayDatabases2) {

			final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase2.getLocation());

			Set<String> pathwayGenes = pathwayMatrixLoader.getAllRowIdentifiers();

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

				bonfExactPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactPvals_" + pathwayDatabase2.getName() + ".txt.gz");
				bonfOdds.save(options.getOutputBasePath() + "_" + predictionSource + "_bonFexactOr_" + pathwayDatabase2.getName() + ".txt.gz");
				outputMatrixAuc.save(options.getOutputBasePath() + "_" + predictionSource + "_auc_" + pathwayDatabase2.getName() + ".txt.gz");
				outputMatrixPvalues.save(options.getOutputBasePath() + "_" + predictionSource + "_utestPvals_" + pathwayDatabase2.getName() + ".txt.gz");
			}
		}
	}

}
