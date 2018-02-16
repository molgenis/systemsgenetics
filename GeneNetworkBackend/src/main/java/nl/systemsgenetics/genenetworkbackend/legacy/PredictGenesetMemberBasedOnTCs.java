package nl.systemsgenetics.genenetworkbackend.legacy;

import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Probability;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import nl.systemsgenetics.genenetworkbackend.PredictGenesetMembersOptions;
import org.apache.commons.cli.ParseException;

public class PredictGenesetMemberBasedOnTCs {

	public static void main(String[] args) throws IOException {

		PredictGenesetMembersOptions options;

		if (args.length == 0) {
			PredictGenesetMembersOptions.printHelp();
			return;
		}

		try {
			options = new PredictGenesetMembersOptions(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			PredictGenesetMembersOptions.printHelp();
			return;
		}

		options.printOptions();

		final String eigenVectorPath = options.getEigenVectorFile().getAbsolutePath();
		final String pathwayMatrixPath = options.getPathwayMatrixFile().getAbsolutePath();
		
		final HashMap genesToInclude = options.getBackgroudGenesFile() == null ? null : readFileLinesToHashMap(options.getBackgroudGenesFile());

		//Load the matrix with the eigenvectors (rows = genes, columns = components, missing values are coded as NaN, components from multiple species can be combined)
		String filenameEigenvector = eigenVectorPath;
		DoubleRandomEngine randomEngine = new DRand();
		ExpressionDataset eigenVectorDataset = new ExpressionDataset(filenameEigenvector, "\t", genesToInclude);
		ExpressionDataset eigenVectorDatasetTransposed = new ExpressionDataset(eigenVectorDataset.fileName, "\t", genesToInclude);
		eigenVectorDatasetTransposed.transposeDataset();
		
		if(genesToInclude != null && eigenVectorDataset.nrProbes != genesToInclude.size()){
			System.err.println("Not all genes in background file found in eigenvector matrix" );
		}

		//Identify strata in the dataset (mix of different species, that cause different missing values per species):
		int nrStrata = 0;
		int[] nrTCsPerStrata = new int[eigenVectorDataset.nrSamples];
		nrTCsPerStrata[nrStrata]++;
		for (int s = 1; s < eigenVectorDataset.nrSamples; s++) {
			for (int p = 0; p < eigenVectorDataset.nrProbes; p++) {
				if (Double.isNaN(eigenVectorDataset.rawData[p][s]) != Double.isNaN(eigenVectorDataset.rawData[p][s - 1])) {
					nrStrata++;
					break;
				}
			}
			nrTCsPerStrata[nrStrata]++;
		}
		nrStrata++;
		int[] strataStartColumn = new int[nrStrata];
		if (nrStrata > 1) {
			for (int s = 0; s < nrStrata; s++) {
				if (s > 0) {
					strataStartColumn[s] = strataStartColumn[s - 1] + nrTCsPerStrata[s - 1];
				}
				System.out.println("Number of species:\t" + s + "\t" + nrTCsPerStrata[s] + "\t" + strataStartColumn[s]);
			}
		}

		//Load the matrix with genesets (rows = genes, columns = genesets, rows should be identical to rows filenameEigenvector, coding = 1 or 0):
		ExpressionDatasetBool datasetGeneset = new ExpressionDatasetBool(pathwayMatrixPath, "\t", genesToInclude);
		ExpressionDatasetBool datasetGenesetT = new ExpressionDatasetBool(pathwayMatrixPath, "\t", genesToInclude);
		datasetGenesetT.transposeDataset();

		if(genesToInclude != null && datasetGeneset.nrProbes != genesToInclude.size()){
			System.err.println("Not all genes in background file found in pathway matrix" );
		}
		
		ExpressionDataset result = new ExpressionDataset(pathwayMatrixPath);

		//Load the Ensembl to HGNC annotation:
//        HashMap hashEnsemblToHGNC = new HashMap();
//        HashMap hashHGNCToEnsembl = new HashMap();
//        try {
//           java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File("EnsemblToHGNCAnnotation.txt")));
//           String str = in.readLine();
//           while ((str = in.readLine()) != null) {
//               String[] data = str.split("\t");
//                hashEnsemblToHGNC.put(data[0],data[1]);
//                hashHGNCToEnsembl.put(data[1],data[0]);
//            }
//        } catch (Exception e) {
//            System.out.println("Error:\t" + e.getMessage());
//            e.printStackTrace();
//        }        
		//Now process each of the genesets:
		String[] queryGenesets = datasetGeneset.sampleNames;
		for (int q = 0; q < queryGenesets.length; q++) {
			int geneset = ((Integer) datasetGeneset.hashSamples.get(queryGenesets[q])).intValue();
			System.out.println("Inventorizing geneset:\t" + queryGenesets[q] + "\t" + geneset);

			//Check whether there are at least 10 genes for each TC with data for this geneset:
			int minNrGenesPerStrata = eigenVectorDataset.nrProbes;
			int[] nrGenesInGenesetPerStrata = new int[nrStrata];
			int[] nrGenesNotInGenesetPerStrata = new int[nrStrata];
			for (int s = 0; s < nrStrata; s++) {
				int nrGenesAvailable = 0;
				for (int gene = 0; gene < eigenVectorDataset.nrProbes; gene++) {
					if (!Double.isNaN(eigenVectorDatasetTransposed.rawData[strataStartColumn[s]][gene])) {
						if (datasetGeneset.rawData[gene][geneset]) {
							nrGenesAvailable++;
							nrGenesInGenesetPerStrata[s]++;
						} else {
							nrGenesNotInGenesetPerStrata[s]++;
						}
					}
				}
				minNrGenesPerStrata = Math.min(nrGenesAvailable, minNrGenesPerStrata);
			}

			System.out.println("Inventorizing geneset:\t" + queryGenesets[q] + "\t" + geneset + "\t" + minNrGenesPerStrata);

			if (minNrGenesPerStrata >= 10) {

				System.out.println("Processing geneset:\t" + queryGenesets[q] + "\t" + geneset + "\t" + minNrGenesPerStrata);
				System.out.println("TC\tTCIndex\tNrGenesInGeneset\tT-Statistic\tZScore\tpValue\tdf");

				double[] zScoreProfile = new double[eigenVectorDataset.nrSamples];
				for (int s = 0; s < nrStrata; s++) {

					double[] vals1 = new double[nrGenesInGenesetPerStrata[s]];
					double[] vals2 = new double[nrGenesNotInGenesetPerStrata[s]];
					for (int tc = 0; tc < nrTCsPerStrata[s]; tc++) {
						int vals1Itr = 0;
						int vals2Itr = 0;
						for (int gene = 0; gene < eigenVectorDataset.nrProbes; gene++) {
							if (!Double.isNaN(eigenVectorDatasetTransposed.rawData[tc + strataStartColumn[s]][gene])) {
								if (datasetGeneset.rawData[gene][geneset]) {
									vals1[vals1Itr] = eigenVectorDatasetTransposed.rawData[tc + strataStartColumn[s]][gene];
									vals1Itr++;
								} else {
									vals2[vals2Itr] = eigenVectorDatasetTransposed.rawData[tc + strataStartColumn[s]][gene];
									vals2Itr++;
								}
							}
						}
						double n1 = vals1.length;
						double n2 = vals2.length;
						double mean1 = mean(vals1);
						double mean2 = mean(vals2);
						double var1 = variance(vals1, mean1);
						double var2 = variance(vals2, mean2);
						double t = (mean1 - mean2) / Math.sqrt(var1 / n1 + var2 / n2);
						//Assume unequal variance (Welch T-Test):
						double df = ((var1 / n1 + var2 / n2) * (var1 / n1 + var2 / n2)) / (((var1 / n1) * (var1 / n1)) / (n1 - 1) + ((var2 / n2) * (var2 / n2)) / (n2 - 1));

						StudentT tDist = new StudentT(df, randomEngine);
						double pValue;
						double zScore;
						if (t < 0) {
							pValue = tDist.cdf(t);
							if (pValue < 2.0E-323) {
								pValue = 2.0E-323;
							}
							zScore = Probability.normalInverse(pValue);
						} else {
							pValue = tDist.cdf(-t);
							if (pValue < 2.0E-323) {
								pValue = 2.0E-323;
							}
							zScore = -Probability.normalInverse(pValue);
						}
						pValue *= 2;
						zScoreProfile[tc + strataStartColumn[s]] = zScore;

						System.out.println(eigenVectorDataset.sampleNames[tc + strataStartColumn[s]] + "\t" + (tc + strataStartColumn[s]) + "\t" + nrGenesInGenesetPerStrata[s] + "\t" + t + "\t" + zScore + "\t" + pValue + "\t" + df);
					}
				}

				System.out.println("\n\n");
				//System.out.println("ProbeIndex\tEnsemblGene\tCorrelation\tP-Value\tZ-Score\tGeneInGenesetInput");
				double[] zScoresIndividualGenes = new double[eigenVectorDataset.nrProbes];

				int countPreventedOverwfitting = 0;
				long timeInVector = 0;
				long timeInCor = 0;
				long timeinT = 0;
				long timeinPreventOverfitA = 0;
				long timeinPreventOverfitB = 0;
				long timeinPreventOverfitBb = 0;
				long timeinPreventOverfitC = 0;
				long timeinPreventOverfitCb = 0;
				long timeinPreventOverfitCc = 0;
				long timeinPreventOverfitD = 0;

				long timeStart;
				long timeStop;

				for (int p = 0; p < eigenVectorDataset.nrProbes; p++) {

//					if (p % 1000 == 0) {
//						System.out.println("Processed: " + p + " genes of which " + countPreventedOverwfitting + " used overfitting protection");
//						System.out.println(" - Time in prefent overfit A: " + timeinPreventOverfitA);
//						System.out.println(" - Time in prefent overfit B: " + timeinPreventOverfitB);
//						System.out.println(" - Time in prefent overfit Bb: " + timeinPreventOverfitBb);
//						System.out.println(" - Time in prefent overfit C: " + timeinPreventOverfitC);
//						System.out.println(" - Time in prefent overfit Cb: " + timeinPreventOverfitCb);
//						System.out.println(" - Time in prefent overfit Cc: " + timeinPreventOverfitCc);
//						System.out.println(" - Time in prefent overfit D: " + timeinPreventOverfitD);
//						System.out.println(" - Time in vector: " + timeInVector);
//						System.out.println(" - Time in Cor: " + timeInCor);
//						System.out.println(" - Time in T: " + timeinT);
//					}
					double[] zScoreProfileToUse;

					//Prevent overfitting:
					if (datasetGeneset.rawData[p][geneset]) {

						timeStart = System.currentTimeMillis();

						++countPreventedOverwfitting;
						zScoreProfileToUse = new double[eigenVectorDataset.nrSamples];
						int[] nrGenesInGenesetPerStrataToUse = new int[nrStrata];
						int[] nrGenesNotInGenesetPerStrataToUse = new int[nrStrata];
						for (int s = 0; s < nrStrata; s++) {
							for (int gene = 0; gene < eigenVectorDataset.nrProbes; gene++) {
								if (p != gene) {
									if (!Double.isNaN(eigenVectorDatasetTransposed.rawData[strataStartColumn[s]][gene])) {
										if (datasetGeneset.rawData[gene][geneset]) {
											nrGenesInGenesetPerStrataToUse[s]++;
										} else {
											nrGenesNotInGenesetPerStrataToUse[s]++;
										}
									}
								}
							}
						}

						timeStop = System.currentTimeMillis();
						timeinPreventOverfitA += (timeStop - timeStart);

						for (int s = 0; s < nrStrata; s++) {
							double[] vals1 = new double[nrGenesInGenesetPerStrataToUse[s]];
							double[] vals2 = new double[nrGenesNotInGenesetPerStrataToUse[s]];
							for (int tc = 0; tc < nrTCsPerStrata[s]; tc++) {

								timeStart = System.currentTimeMillis();

								double[] x = eigenVectorDatasetTransposed.rawData[tc + strataStartColumn[s]];
								boolean[] y = datasetGenesetT.rawData[geneset];

								timeStop = System.currentTimeMillis();
								timeinPreventOverfitB += (timeStop - timeStart);

								timeStart = System.currentTimeMillis();

								//double dummy;
								int vals1Itr = 0;
								int vals2Itr = 0;
								for (int gene = 0; gene < eigenVectorDataset.nrProbes; gene++) {
									if (p != gene) {
										//if (!Double.isNaN(x[gene])) {
										if (y[gene]) {
											vals1[vals1Itr++] = x[gene];
										} else {
											vals2[vals2Itr++] = x[gene];
										}
										//}
									}
								}

								timeStop = System.currentTimeMillis();
								timeinPreventOverfitBb += (timeStop - timeStart);

								timeStart = System.currentTimeMillis();

								double n1 = vals1.length;
								double n2 = vals2.length;
								double mean1 = mean(vals1);
								double mean2 = mean(vals2);
								double var1 = variance(vals1, mean1);
								double var2 = variance(vals2, mean2);

								timeStop = System.currentTimeMillis();
								timeinPreventOverfitC += (timeStop - timeStart);

								timeStart = System.currentTimeMillis();

								double t = (mean1 - mean2) / Math.sqrt(var1 / n1 + var2 / n2);

								timeStop = System.currentTimeMillis();
								timeinPreventOverfitCb += (timeStop - timeStart);

								timeStart = System.currentTimeMillis();

								double df = ((var1 / n1 + var2 / n2) * (var1 / n1 + var2 / n2)) / (((var1 / n1) * (var1 / n1)) / (n1 - 1) + ((var2 / n2) * (var2 / n2)) / (n2 - 1));

								timeStop = System.currentTimeMillis();
								timeinPreventOverfitCc += (timeStop - timeStart);

								timeStart = System.currentTimeMillis();

								//double df = n1 + n2 - 2;
								StudentT tDist = new StudentT(df, randomEngine);
								double pValue;
								double zScore;
								if (t < 0) {
									pValue = tDist.cdf(t);
									if (pValue < 2.0E-323) {
										pValue = 2.0E-323;
									}
									zScore = Probability.normalInverse(pValue);
								} else {
									pValue = tDist.cdf(-t);
									if (pValue < 2.0E-323) {
										pValue = 2.0E-323;
									}
									zScore = -Probability.normalInverse(pValue);
								}
								zScoreProfileToUse[tc + strataStartColumn[s]] = zScore;

								timeStop = System.currentTimeMillis();
								timeinPreventOverfitD += (timeStop - timeStart);

							}

						}

					} else {
						zScoreProfileToUse = zScoreProfile;
					}

					timeStart = System.currentTimeMillis();

					TDoubleList vec1 = new TDoubleArrayList();
					TDoubleList vec2 = new TDoubleArrayList();
					for (int tc = 0; tc < eigenVectorDataset.nrSamples; tc++) {
						if (!Double.isNaN(eigenVectorDataset.rawData[p][tc])) {
							vec1.add(eigenVectorDataset.rawData[p][tc]);
							vec2.add(zScoreProfileToUse[tc]);
						}
					}
					double[] vals1 = new double[vec1.size()];
					double[] vals2 = new double[vec2.size()];
					for (int v = 0; v < vec1.size(); v++) {
						vals1[v] = vec1.get(v);
						vals2[v] = vec2.get(v);
						//vals1[v * 2 + 1] = -vec1.get(v);
						//vals2[v * 2 + 1] = -vec2.get(v);
					}

					timeStop = System.currentTimeMillis();
					timeInVector += (timeStop - timeStart);

					timeStart = System.currentTimeMillis();

					double correlation = JSci.maths.ArrayMath.correlation(vals1, vals2);

					timeStop = System.currentTimeMillis();
					timeInCor += (timeStop - timeStart);

					timeStart = System.currentTimeMillis();

					StudentT tDistColt = new StudentT(vals1.length / 2 - 2, randomEngine);
					double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (vals1.length / 2 - 2)));
					double pValue;
					double zScore;
					if (t < 0) {
						pValue = tDistColt.cdf(t);
						if (pValue < 2.0E-323) {
							pValue = 2.0E-323;
						}
						zScore = Probability.normalInverse(pValue);
					} else {
						pValue = tDistColt.cdf(-t);
						if (pValue < 2.0E-323) {
							pValue = 2.0E-323;
						}
						zScore = -Probability.normalInverse(pValue);
					}
					pValue *= 2;

					timeStop = System.currentTimeMillis();
					timeinT += (timeStop - timeStart);

					zScoresIndividualGenes[p] = zScore;
					// System.out.println(p + "\t" + eigenVectorDataset.probeNames[p] + "\t" + correlation + "\t" + pValue + "\t" + zScore + "\t" + datasetGeneset.rawData[p][geneset]);
				}

				System.out.println("Done with genes");

				System.out.println("\n\n");
				TDoubleList vecInGeneset = new TDoubleArrayList();
				TDoubleList vecNotInGeneset = new TDoubleArrayList();
				for (int p = 0; p < eigenVectorDataset.nrProbes; p++) {
					if (datasetGeneset.rawData[p][geneset]) {
						vecInGeneset.add(zScoresIndividualGenes[p]);
					} else {
						vecNotInGeneset.add(zScoresIndividualGenes[p]);
					}
				}
				double[] vals1 = new double[vecInGeneset.size()];
				for (int v = 0; v < vals1.length; v++) {
					vals1[v] = vecInGeneset.get(v);
				}
				double[] vals2 = new double[vecNotInGeneset.size()];
				for (int v = 0; v < vals2.length; v++) {
					vals2[v] = vecNotInGeneset.get(v);
				}

				System.out.println("Starting wilcox test");

				WilcoxonMannWhitney wmw = new WilcoxonMannWhitney();
				double pValue = wmw.returnWilcoxonMannWhitneyPValue(vals1, vals2);
				System.out.println("PredictionPerformanceGeneset:\t" + queryGenesets[q] + "\tGenesetNr:\t" + geneset + "\t\tNrGenesInGeneset\t" + vals1.length + "\tNrGenesOutsideGeneset:\t" + vals2.length + "\tP-Value:\t" + pValue + "\tAUC:\t" + wmw.getAUC());
				System.out.println();
				//Replace the 1 or 0 values with the Z-Scores:
				for (int p = 0; p < eigenVectorDataset.nrProbes; p++) {
					result.rawData[p][geneset] = zScoresIndividualGenes[p];
				}

			}

		}

		result.save(options.getPredictionsFile().getAbsolutePath());

	}

	private static double mean(double[] v) {
		double sum = 0;
		for (int k = 0; k < v.length; k++) {
			sum += v[k];
		}
		return (sum / (double) v.length);
	}

	private static double variance(double[] v, double mean) {
		double ans = 0.0;
		for (int i = 0; i < v.length; i++) {
			ans += (v[i] - mean) * (v[i] - mean);
		}
		return ans / (v.length - 1);
	}

	/**
	 * HashMap construction to be compatible with old Lude matrix
	 *
	 * @param file
	 * @return
	 */
	private static HashMap readFileLinesToHashMap(File file) throws FileNotFoundException, IOException {

		HashMap fileContent = new HashMap();

		BufferedReader reader = new BufferedReader(new FileReader(file));

		String line;
		while ((line = reader.readLine()) != null) {
			fileContent.put(line, line);
		}

		return fileContent;

	}

}
