/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.originalLude;

import java.io.File;
import java.util.HashMap;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author lude
 */
public class PerformDEPICT2Analysis {

	public PerformDEPICT2Analysis() {

		performDEPICTAnalysis();

		System.exit(0);

	}

	public void performDEPICTAnalysis() {

		//This is the code that we used for studying the trans-eQTLs. For the trans-eQTL analysis we had studied a set of 10,000 SNPs, 
		//and ascertained whether these SNPs were affecting the expression of 20,000 genes that were expressed in blood.
		//As such the nomenclature is that we are dealing with SNPs and genes.
		//The directory holding the data:
		String fileDir = "/Users/lude/Documents/Genetica/eQTLGen/ZScoreMarices-20180125/";

		//Load eQTL matrix, first the Z-scores, subsequently the number of samples used for ascertaining the relationship between the SNP and genes:
		//First load the dataset that holds the Z-score matrix, SNPs are rows, phenotypes are columns
		//IMPORTANT: The file that is loaded here is a file that contains permuted results! The real Z-score matrix, based on actual GWAS results is loaded at a later point!
		ExpressionDataset dataset = new ExpressionDataset(fileDir + "/ZScoreMatrix-Permutation.txt.binary", "\t", null, null);
		//Now load the dataset that holds the NrSamples matrix, SNPs are rows, phenotypes are columns
		ExpressionDataset datasetNrSamples = new ExpressionDataset(fileDir + "/ZScoreMatrixNrSamples.txt.binary");

		//Now load the genomic positions of each of the SNPs. IMPORTANT: This file needs to be sorted, first per chromosome, and subsequently per chromosome position!!!
		HashMap hashColumns = new HashMap();
		hashColumns.put("CHR_HG19", null);
		hashColumns.put("POS_HG19", null);
		ExpressionDataset datasetSNPPos = new ExpressionDataset("/Users/lude/Documents/Genetica/eQTLGen/ZScoreMatrices/help_files/genetic_risk_factors_cleaned_traits_standardized_sorted_chromosomal_position.txt", "\t", dataset.hashProbes, hashColumns);
		hashColumns = new HashMap();
		hashColumns.put("Chr", null);
		hashColumns.put("ChrStart", null);
		hashColumns.put("ChrEnd", null);

		//Load the genomic positions of each of the genes tested:
		ExpressionDataset datasetGenePos = new ExpressionDataset("/Users/lude/Documents/Genetica/eQTLGen/BIOSData/HelpFiles/ProbeAnnotation_STARv2.3.0e_Ensembl71.txt", "\t", dataset.hashSamples, hashColumns);

		//Determine list of genes that are not co-expressed, based on the 10 permuted Z-score matrices:
		//ProbeCovariance-Perms1To10.binary is a symmetric matrix holding the correlations between each of the genes. 
		//This file is used to prune these genes, and arrive at a set of genes that are not correlated:
		HashMap hashPrunedGenes = new HashMap();
		ExpressionDataset datasetCoexpression = new ExpressionDataset(fileDir + "/GeneCorrelation/ProbeCovariance-Perms1To10.binary");
		boolean[] geneExcluded = new boolean[dataset.nrSamples];
		for (int p = 0; p < dataset.nrSamples; p++) {
			if (!geneExcluded[p]) {
				for (int q = 0; q < dataset.nrSamples; q++) {
					if (p != q && !geneExcluded[q]) {
						if (Math.abs(datasetCoexpression.rawData[p][q]) > 0.05) {
							geneExcluded[q] = true;
						}
					}
				}
				hashPrunedGenes.put(dataset.sampleNames[p], null);
			}
		}
		System.out.println("Number of pruned genes:\t" + hashPrunedGenes.size());
		datasetCoexpression.rawData = null;
		datasetCoexpression = null;

		File file = new File(dataset.fileName + ".SNPSetXGenePValues-100kbDistance.binary");
		if (!file.exists()) {
			System.out.println("Generating the matrix with the SNPset x Gene P-Values, bassed on 100kb distance:");
			//Generate a matrix that permits accurate determination of LD between SNPs. We use the first 5 permuted Z-score matrices, and only use genes that are not co-expressed:
			int nrPrunedGenes = hashPrunedGenes.size();
			ExpressionDataset datasetNull = new ExpressionDataset(dataset.nrProbes, nrPrunedGenes * 5);
			for (int perm = 0; perm < 5; perm++) {
				ExpressionDataset datasetPerm = new ExpressionDataset(fileDir + "/ZScoreMatrix-Permutation" + (perm + 1) + ".txt.binary", "\t", null, hashPrunedGenes);
				datasetNull.probeNames = datasetPerm.probeNames;
				for (int p = 0; p < datasetPerm.nrProbes; p++) {
					for (int s = 0; s < datasetPerm.nrSamples; s++) {
						datasetNull.rawData[p][s + perm * nrPrunedGenes] = datasetPerm.rawData[p][s];
					}
				}
			}

			//Determine per gene which SNPs are in close physical proximity:
			int[][] geneToSNP = new int[dataset.nrSamples][dataset.nrProbes];//samples is genes and probes is snps
			int[] geneNrSNPsAssigned = new int[dataset.nrSamples];//number of snps per gene
			int nrPerms = 100000000;
			double[] geneChi2SumNull = new double[nrPerms];//For the current gene the permutation chi2 values
			int[] pValueDistribution = new int[21];//used to create histogram 

			ExpressionDataset datasetSNPSetXGenePValues = new ExpressionDataset(datasetGenePos.nrProbes, dataset.nrSamples);// probes genomic region (aggregate over SNPs around a gene, samples is number of phenotypes)
			datasetSNPSetXGenePValues.probeNames = datasetGenePos.probeNames;
			datasetSNPSetXGenePValues.sampleNames = dataset.sampleNames;
			datasetSNPSetXGenePValues.recalculateHashMaps();
			for (int gene = 0; gene < datasetGenePos.nrProbes; gene++) { // This loop is mostly done when loading genotype data but make sure to prune out near perfect LD variants
				for (int s = 0; s < dataset.nrSamples; s++) {
					datasetSNPSetXGenePValues.rawData[gene][s] = Double.NaN;
				}
				for (int snp = 0; snp < dataset.nrProbes; snp++) {
					geneToSNP[gene][snp] = -1;
				}
				int geneChr = (int) datasetGenePos.rawData[gene][0];
				int geneChrPosStart = (int) datasetGenePos.rawData[gene][1];
				int geneChrPosEnd = (int) datasetGenePos.rawData[gene][2];
				for (int snp = 0; snp < datasetSNPPos.nrProbes; snp++) {
					if (datasetSNPPos.rawData[snp][0] == geneChr) {
						int snpPos = (int) datasetSNPPos.rawData[snp][1];
						if (snpPos >= geneChrPosStart - 100000 && snpPos <= geneChrPosEnd + 100000) {
							boolean includeSNP = true;
							int snpIndex = ((Integer) dataset.hashProbes.get(datasetSNPPos.probeNames[snp])).intValue();;
							for (int q = 0; q < geneNrSNPsAssigned[gene]; q++) {//check back over already included SNPs and only include this SNP if not in strong LD
								double corr = JSci.maths.ArrayMath.correlation(datasetNull.rawData[snpIndex], datasetNull.rawData[geneToSNP[gene][q]]);
								if (Math.abs(corr) > 0.95) {
									includeSNP = false;
									break;
								}
							}
							if (includeSNP) {
								geneToSNP[gene][geneNrSNPsAssigned[gene]] = snpIndex;
								geneNrSNPsAssigned[gene]++;
							}
						}
					}
				}
				System.out.print(gene + "\t" + datasetGenePos.probeNames[gene] + "\t" + geneNrSNPsAssigned[gene]);
//				for (int p = 0; p < geneNrSNPsAssigned[gene]; p++) {
//					int snpIndexP = geneToSNP[gene][p];
//					System.out.print("\t" + datasetNull.probeNames[snpIndexP]);
//				}
				System.out.println("");
				if (geneNrSNPsAssigned[gene] > 1) {
					double[][] cov = new double[geneNrSNPsAssigned[gene]][geneNrSNPsAssigned[gene]];
					for (int p = 0; p < geneNrSNPsAssigned[gene]; p++) {
						int snpIndexP = geneToSNP[gene][p];
						cov[p][p] = 1;
						for (int q = p + 1; q < geneNrSNPsAssigned[gene]; q++) {
							int snpIndexQ = geneToSNP[gene][q];
							cov[p][q] = JSci.maths.ArrayMath.correlation(datasetNull.rawData[snpIndexP], datasetNull.rawData[snpIndexQ]);
							cov[q][p] = cov[p][q];
						}
					}
					Jama.EigenvalueDecomposition eig = eigenValueDecomposition(cov);
					double[] eigenValues = eig.getRealEigenvalues();

					int nrThreads = Runtime.getRuntime().availableProcessors();
					java.util.concurrent.ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
					CompletionService<DoubleArrayIntegerObject> pool = new ExecutorCompletionService<DoubleArrayIntegerObject>(threadPool);
					int nrTasks = 0;
					for (int taskNr = 0; taskNr < 10000; taskNr++) {
						EstimateChi2SumDistUsingCorrelatedVariablesThread task = new EstimateChi2SumDistUsingCorrelatedVariablesThread(eigenValues, 10000, taskNr);
						//pool.submit(task);
						nrTasks++;
					}
					try {
						for (int task = 0; task < nrTasks; task++) {
							try {
								DoubleArrayIntegerObject result = pool.take().get();
								int taskNr = result.intValue;
								double[] chi2SumDist = result.doubleArray;
								for (int w = 0; w < chi2SumDist.length; w++) {
									geneChi2SumNull[taskNr * 10000 + w] = chi2SumDist[w];
								}
							} catch (ExecutionException ex) {
								ex.printStackTrace();
								System.out.println(ex.getMessage());
							}
						}
						threadPool.shutdown();
					} catch (Exception e) {
						e.printStackTrace();
						System.out.println(e.getMessage());
					}

				}
				for (int phenotypeIndex = 0; phenotypeIndex < dataset.nrSamples; phenotypeIndex++) {
					//for (int phenotypeIndex = 0; phenotypeIndex<100; phenotypeIndex++) {
					if (geneNrSNPsAssigned[gene] > 1) {
						double sumChi2 = 0;
						boolean missingData = false;
						for (int p = 0; p < geneNrSNPsAssigned[gene]; p++) {
							if (datasetNrSamples.rawData[geneToSNP[gene][p]][phenotypeIndex] == 0) {
								missingData = true;
								break;
							}
							double z = dataset.rawData[geneToSNP[gene][p]][phenotypeIndex];
							sumChi2 += z * z;
						}
						if (!missingData) {
							double p = 0.5;
							int nrPermsHere = 1000; 
							for (int perm = 0; perm < nrPermsHere; perm++) {
								if (geneChi2SumNull[perm] >= sumChi2) {
									p += 1;
								}
								if (perm == nrPermsHere - 1 && p < 5 && nrPermsHere < nrPerms) {
									nrPermsHere *= 10;
								}
							}
							p /= (double) nrPermsHere + 1;
							pValueDistribution[(int) (20d * p)]++;
							datasetSNPSetXGenePValues.rawData[gene][phenotypeIndex] = p;
							if (p < 1e-8) {
								System.out.println(gene + "\t" + datasetGenePos.probeNames[gene] + "\t" + geneNrSNPsAssigned[gene] + "\t" + sumChi2 + "\t" + p);
							}
						}
					} else if (geneNrSNPsAssigned[gene] == 1) {
						if (datasetNrSamples.rawData[geneToSNP[gene][0]][phenotypeIndex] != 0) {
							double z = dataset.rawData[geneToSNP[gene][0]][phenotypeIndex];
							//double p = 2d * cern.jet.stat.Probability.normal(-Math.abs(z));
							double p = ZScores.zToP(-Math.abs(z));//Patrick used differnt method to get to p-value because of old lib
							if (p == 1) {
								p = 0.99999d;
							}
							if (p < 1e-300d) {
								p = 1e-300d;
							}
							pValueDistribution[(int) (20d * p)]++;
							datasetSNPSetXGenePValues.rawData[gene][phenotypeIndex] = p;
							//System.out.println(gene + "\t" + datasetGenePos.probeNames[gene] + "\t" + geneNrSNPsAssigned[gene] + "\t" + z*z + "\t" + p);
						}

					}
				}
				if (gene % 10 == 0) {
					System.out.println("\n\n");
					for (int a = 0; a <= 20; a++) {
						System.out.println(gene + "\t" + a / 20d + "\t" + pValueDistribution[a]);
					}
				}
			}
			datasetSNPSetXGenePValues.save(dataset.fileName + ".SNPSetXGenePValues-100kbDistance.binary");
			System.exit(0);
		}

		ExpressionDataset datasetGenesets = null;

		//We now have a matrix where rows = SNPs collapsed into SNPsets (sorted by gene), cols = gene expression phenotypes.
		//We can now calculate to what extend these SNPsets are correlated, enabling us to correct for correlated SNPsets.
		//Reload the matrix, but only the pruned genes:
		ExpressionDataset datasetSNPSetXGenePValues = new ExpressionDataset(dataset.fileName + ".SNPSetXGenePValues-100kbDistance.binary", "\t", datasetGenesets.hashProbes, hashPrunedGenes);
		//Only include SNP genesets that have been tested, but exclude the HLA region:
		HashMap hashValidSNPSets = new HashMap();
		for (int p = 0; p < datasetSNPSetXGenePValues.nrProbes; p++) {
			int chr = (int) datasetGenePos.rawData[((Integer) datasetGenePos.hashProbes.get(datasetSNPSetXGenePValues.probeNames[p])).intValue()][0];
			int chrPosStart = (int) datasetGenePos.rawData[((Integer) datasetGenePos.hashProbes.get(datasetSNPSetXGenePValues.probeNames[p])).intValue()][1];
			int chrPosEnd = (int) datasetGenePos.rawData[((Integer) datasetGenePos.hashProbes.get(datasetSNPSetXGenePValues.probeNames[p])).intValue()][1];
			if (chr == 6 && chrPosEnd > 20000000 && chrPosStart <= 40000000) {
				//Exclude the HLA. We did not check whether this is needed statistically, but for many traits the HLA shows strong association, that might inflate results.
			} else {
				for (int s = 0; s < datasetSNPSetXGenePValues.nrSamples; s++) {
					if (!Double.isNaN(datasetSNPSetXGenePValues.rawData[p][s])) {
						hashValidSNPSets.put(datasetSNPSetXGenePValues.probeNames[p], null);
						break;
					}
				}
			}
		}
		datasetSNPSetXGenePValues = new ExpressionDataset(datasetSNPSetXGenePValues.fileName, "\t", hashValidSNPSets, hashPrunedGenes);
		for (int p = 0; p < datasetSNPSetXGenePValues.nrProbes; p++) {
			for (int s = 0; s < datasetSNPSetXGenePValues.nrSamples; s++) {
				if (Double.isNaN(datasetSNPSetXGenePValues.rawData[p][s])) {
					//datasetSNPSetXGenePValues.rawData[p][s] = 0;
				} else {
					datasetSNPSetXGenePValues.rawData[p][s] = cern.jet.stat.Probability.normalInverse(datasetSNPSetXGenePValues.rawData[p][s]);
				}
			}
		}

		ExpressionDataset datasetSNPSetXGenePValuesReal = new ExpressionDataset(fileDir + "/ZScoreMatrix.txt.binary" + ".SNPSetXGenePValues-100kbDistance.binary", "\t", hashValidSNPSets, null);
		for (int p = 0; p < datasetSNPSetXGenePValuesReal.nrProbes; p++) {
			for (int s = 0; s < datasetSNPSetXGenePValuesReal.nrSamples; s++) {
				if (Double.isNaN(datasetSNPSetXGenePValuesReal.rawData[p][s])) {
					//datasetSNPSetXGenePValuesReal.rawData[p][s] = 0;
				} else {
					if (datasetSNPSetXGenePValuesReal.rawData[p][s] < 5e-9d) {
						datasetSNPSetXGenePValuesReal.rawData[p][s] = 5e-9d;
					}
					double zScore = -cern.jet.stat.Probability.normalInverse(datasetSNPSetXGenePValuesReal.rawData[p][s] / 2);
					datasetSNPSetXGenePValuesReal.rawData[p][s] = zScore * zScore;
				}
			}
		}

		//Calculate the weight of the individual SNPsets in order to correct for genes in close proximity that yield highly similar, and thus correlated results:
		double[] weightSNPSet = new double[datasetSNPSetXGenePValues.nrProbes];
		if (1 == 1) {
			System.out.println("Calculating the weights of the individual genes, necessary in order to account for LD between genes");
			int[] centromerePos = {125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000, 49000000, 40200000, 53700000, 35800000, 17900000, 17600000, 19000000, 36600000, 24000000, 17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 12500000};
			//Use PCA to decompose per chromosome the gene correlation matrix
			//Subsequently we can infer the weights of each gene by summing the squared factorloadings weighted for the inverse of the eigenvalue, if the eigenvalue > 1;
			for (int chr = 1; chr <= 22; chr++) {
				for (int shortLong = 0; shortLong < 2; shortLong++) {
					int counter = 0;
					int[] geneIndex = new int[datasetSNPSetXGenePValues.nrProbes];
					for (int gene = 0; gene < datasetSNPSetXGenePValues.nrProbes; gene++) {
						int geneChr = (int) datasetGenePos.rawData[((Integer) datasetGenePos.hashProbes.get(datasetSNPSetXGenePValues.probeNames[gene])).intValue()][0];
						int geneChrPos = (int) datasetGenePos.rawData[((Integer) datasetGenePos.hashProbes.get(datasetSNPSetXGenePValues.probeNames[gene])).intValue()][1];
						if (geneChr == chr) {
							if (shortLong == 0 && geneChrPos < centromerePos[chr - 1]) {
								geneIndex[counter] = gene;
								counter++;
							}
							if (shortLong == 1 && geneChrPos >= centromerePos[chr - 1]) {
								geneIndex[counter] = gene;
								counter++;
							}
						}
					}
					System.out.print("Processing chromosome:\t" + chr + "\t" + shortLong + "\tNr of genes with nearby SNPs:\t" + counter);
					if (counter > 0) {
						double[][] corrMatrix = new double[counter][counter];
						for (int x = 0; x < counter; x++) {
							int p = geneIndex[x];
							for (int y = x; y < counter; y++) {
								int q = geneIndex[y];
								double correlation = 0;
								int nr = 0;
								for (int s = 0; s < datasetSNPSetXGenePValues.nrSamples; s++) {
									if (!Double.isNaN(datasetSNPSetXGenePValues.rawData[p][s]) && !Double.isNaN(datasetSNPSetXGenePValues.rawData[q][s])) {
										correlation += datasetSNPSetXGenePValues.rawData[p][s] * datasetSNPSetXGenePValues.rawData[q][s];
										nr++;
									}
								}
								correlation /= (double) (nr - 1);
								corrMatrix[x][y] = correlation;
								corrMatrix[y][x] = correlation;
							}
						}
						//Calculate the factorloadings:
						Jama.EigenvalueDecomposition eig = eigenValueDecomposition(corrMatrix);
						double[] eigenValues = eig.getRealEigenvalues();
						double[][] factorLoadings = new double[counter][counter];
						for (int comp = 0; comp < counter; comp++) {
							double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
							if (eigenvalue < 0) {
								eigenvalue = 0;
							}
							double sqrtEigenvalue = Math.sqrt(eigenvalue);
							factorLoadings[comp] = getEigenVector(eig, comp);
							for (int a = 0; a < counter; a++) {
								factorLoadings[comp][a] *= sqrtEigenvalue;
							}
						}
						//Calculate the weights of the individual SNPs:
						double[] weights = new double[counter];
						for (int p = 0; p < counter; p++) {
							double weight = 0;
							for (int comp = 0; comp < counter; comp++) {
								double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
								if (eigenvalue < 1) {
									eigenvalue = 1;
								}
								weight += factorLoadings[comp][p] * factorLoadings[comp][p] / eigenvalue;
							}
							weightSNPSet[geneIndex[p]] = weight;
							weights[p] = weight;
						}
						System.out.println("\tAverage weight:\t" + JSci.maths.ArrayMath.mean(weights) + "\tMax weight:\t" + JSci.maths.ArrayMath.max(weights));
					} else {
						System.out.println("");
					}
				}
			}
		}
		double weightsAverage = JSci.maths.ArrayMath.mean(weightSNPSet);
		System.out.println("Average weight of the genes (i.e. effective number of independent, unlinked genes):\t" + weightsAverage);
		System.out.println("\n\n\n\n");

		System.out.println("Now performing pathway enrichment analysis, using pathway information and gene p-values");

		if (1 == 1) {
			double[][] genesetZScores = new double[datasetGenesets.nrSamples][datasetSNPSetXGenePValues.nrProbes];
			for (int p = 0; p < datasetSNPSetXGenePValues.nrProbes; p++) {
				int index = ((Integer) datasetGenesets.hashProbes.get(datasetSNPSetXGenePValues.probeNames[p])).intValue();
				for (int geneset = 0; geneset < datasetGenesets.nrSamples; geneset++) {
					genesetZScores[geneset][p] = datasetGenesets.rawData[index][geneset];
				}
			}
			if (1 == 2) {
				//Do we want do a parametric or non-parametric enrichment test?
				//If we want to do a non-parametric enrichment test we can enfore a normal distributiopn on the geneset Z-scores. This will guard us against strong outliers, but might yield redcuced power.

				//Needs differnt rank algorithm
//				for (int geneset=0; geneset<datasetGenesets.nrSamples; geneset++) {
//                    //Enforce normal distribution on the genesets:
//                    double[] vals = new double[datasetSNPSetXGenePValues.nrProbes];
//                    for (int s=0; s<genesetZScores[geneset].length; s++) {
//                        vals[s] = genesetZScores[geneset][s];
//                    }
//                    jsc.util.Rank rank = new jsc.util.Rank(vals, 0d);
//                    double[] rankedValues = rank.getRanks();
//                    for (int s=0; s<genesetZScores[geneset].length; s++) {
//                        double pValue = (0.5d + rankedValues[s] - 1d) / (double) (rankedValues.length);
//                        double zScore = cern.jet.stat.Probability.normalInverse(pValue);  
//                        genesetZScores[geneset][s] = zScore;
//                    }
//                }
			}

			datasetSNPSetXGenePValuesReal.transposeDataset();
			for (int phenotype = 0; phenotype < datasetSNPSetXGenePValuesReal.nrProbes; phenotype++) {
				double mean = 0;
				int nrSNPSets = 0;
				for (int s = 0; s < datasetSNPSetXGenePValuesReal.nrSamples; s++) {
					if (!Double.isNaN(datasetSNPSetXGenePValuesReal.rawData[phenotype][s])) {
						mean += datasetSNPSetXGenePValuesReal.rawData[phenotype][s];
						nrSNPSets++;
					}
				}
				mean /= (double) nrSNPSets;
				double stdev = 0;
				for (int s = 0; s < datasetSNPSetXGenePValuesReal.nrSamples; s++) {
					if (!Double.isNaN(datasetSNPSetXGenePValuesReal.rawData[phenotype][s])) {
						datasetSNPSetXGenePValuesReal.rawData[phenotype][s] -= mean;
						stdev += datasetSNPSetXGenePValuesReal.rawData[phenotype][s] * datasetSNPSetXGenePValuesReal.rawData[phenotype][s];
					}
				}
				stdev = Math.sqrt(stdev / (double) (nrSNPSets - 1));
				for (int s = 0; s < datasetSNPSetXGenePValuesReal.nrSamples; s++) {
					if (!Double.isNaN(datasetSNPSetXGenePValuesReal.rawData[phenotype][s])) {
						datasetSNPSetXGenePValuesReal.rawData[phenotype][s] /= stdev;
					}
				}
			}
			double[] betas = new double[datasetGenesets.nrSamples];
			for (int geneset = 0; geneset < datasetGenesets.nrSamples; geneset++) {
				double r2Sum = 0;

				for (int phenotype = 0; phenotype < datasetSNPSetXGenePValuesReal.nrProbes; phenotype++) {

					//Calculate the weighted mean for the trans-eQTL phenotype of interest:
					double phenotypeWeightedMean = 0;
					double genesetWeightedMean = 0;
					double sumWeights = 0;
					for (int s = 0; s < datasetSNPSetXGenePValues.nrProbes; s++) {
						if (!Double.isNaN(datasetSNPSetXGenePValuesReal.rawData[phenotype][s])) {
							phenotypeWeightedMean += weightSNPSet[s] * datasetSNPSetXGenePValuesReal.rawData[phenotype][s];
							genesetWeightedMean += weightSNPSet[s] * genesetZScores[geneset][s];
							sumWeights += weightSNPSet[s];
						}
					}
					phenotypeWeightedMean /= sumWeights;
					genesetWeightedMean /= sumWeights;
					double covXX = 0;
					double covYY = 0;
					double covXY = 0;
					for (int s = 0; s < datasetSNPSetXGenePValues.nrProbes; s++) {
						if (!Double.isNaN(datasetSNPSetXGenePValuesReal.rawData[phenotype][s])) {
							covXX += weightSNPSet[s] * (datasetSNPSetXGenePValuesReal.rawData[phenotype][s] - phenotypeWeightedMean) * (datasetSNPSetXGenePValuesReal.rawData[phenotype][s] - phenotypeWeightedMean);
							covYY += weightSNPSet[s] * (genesetZScores[geneset][s] - genesetWeightedMean) * (genesetZScores[geneset][s] - genesetWeightedMean);
							covXY += weightSNPSet[s] * (datasetSNPSetXGenePValuesReal.rawData[phenotype][s] - phenotypeWeightedMean) * (genesetZScores[geneset][s] - genesetWeightedMean);
						}
					}
					covXX /= sumWeights;
					covYY /= sumWeights;
					covXY /= sumWeights;
					double corr = covXY / Math.sqrt(covXX * covYY);

					betas[geneset] = corr;
					double r2 = corr * corr;
					r2Sum += r2;
					if (Math.abs(corr) > 0.055) {
						//Significant correlation:
						double zScorePhenotypeGene = 0;
						String phenotypeGene = datasetSNPSetXGenePValuesReal.probeNames[phenotype];
						if (datasetGenesets.hashProbes.containsKey(phenotypeGene)) {
							zScorePhenotypeGene = datasetGenesets.rawData[((Integer) datasetGenesets.hashProbes.get(phenotypeGene)).intValue()][geneset];
						}
						System.out.println(phenotype + "\t" + phenotypeWeightedMean + "\t" + geneset + "\t" + corr + "\t" + r2 + "\t" + sumWeights + "\t" + datasetSNPSetXGenePValuesReal.probeNames[phenotype] + "\t" + datasetGenesets.sampleNames[geneset] + "\t" + zScorePhenotypeGene);
					}
				}
			}
		}

		System.exit(0);
	}

	public void orthogonalizeDataset(String inputFile) {
		orthogonalizeDataset(inputFile, -1, true);
	}

	public void orthogonalizeDataset(String inputFile, int maxNrComps, boolean standardNormalize) {

		ExpressionDataset dataset = new ExpressionDataset(inputFile);
		dataset.transposeDataset();
		if (standardNormalize) {
			dataset.standardNormalizeData();
		}
		int nrVars = dataset.nrProbes;
		int nrSamples = dataset.nrSamples;
		int nrPCsToOutput = nrVars;
		if (maxNrComps != -1) {
			if (maxNrComps < nrPCsToOutput) {
				nrPCsToOutput = maxNrComps;
			}
		}

		double[][] matrix = new double[nrVars][nrSamples];
		for (int s = 0; s < nrVars; s++) {
			for (int sample = 0; sample < nrSamples; sample++) {
				matrix[s][sample] = dataset.rawData[s][sample];
			}
		}
		double[][] correlationMatrix = new double[nrVars][nrVars];
		for (int p = 0; p < nrVars; p++) {
			correlationMatrix[p][p] = 1d;
			for (int q = p + 1; q < nrVars; q++) {
				double covariance = 0;
				for (int sample = 0; sample < nrSamples; sample++) {
					covariance += matrix[p][sample] * matrix[q][sample];
				}
				covariance /= (double) (nrSamples - 1);
				correlationMatrix[p][q] = covariance;
				correlationMatrix[q][p] = covariance;
			}
		}
		Jama.EigenvalueDecomposition eig = eigenValueDecomposition(correlationMatrix);
		double[] eigenValues = eig.getRealEigenvalues();

		double[][] eigenVectors = new double[correlationMatrix.length][correlationMatrix.length];
		ExpressionDataset datasetEigenvectors = new ExpressionDataset(nrPCsToOutput, correlationMatrix.length);
		ExpressionDataset datasetEigenvalues = new ExpressionDataset(correlationMatrix.length, 3);
		for (int pca = 0; pca < correlationMatrix.length; pca++) {
			if (pca < nrPCsToOutput) {
				datasetEigenvectors.probeNames[pca] = "Comp" + (pca + 1);
			}
			datasetEigenvalues.probeNames[pca] = "Comp" + (pca + 1);
			datasetEigenvectors.sampleNames[pca] = dataset.probeNames[pca];
		}
		datasetEigenvalues.sampleNames[0] = "Eigenvalues";
		datasetEigenvalues.sampleNames[1] = "ExplainedVariance";
		datasetEigenvalues.sampleNames[2] = "CumExplainedVariance";
		for (int pca = 0; pca < correlationMatrix.length; pca++) {
			if (pca < nrPCsToOutput) {
				datasetEigenvectors.rawData[pca] = getEigenVector(eig, pca);
			}
			datasetEigenvalues.rawData[pca][0] = eigenValues[eigenValues.length - 1 - pca];
			datasetEigenvalues.rawData[pca][1] = getEigenValueVar(eigenValues, pca);
			datasetEigenvalues.rawData[pca][2] = datasetEigenvalues.rawData[pca][1];
			if (pca > 0) {
				datasetEigenvalues.rawData[pca][2] += datasetEigenvalues.rawData[pca - 1][2];
			}
			System.out.println(pca + "\tExplainedVariance:\t" + getEigenValueVar(eigenValues, pca) + "\tEigenvalue:\t" + eigenValues[eigenValues.length - 1 - pca]);
		}
		datasetEigenvectors.transposeDataset();
		datasetEigenvectors.save(inputFile + ".Eigenvectors.txt");
		datasetEigenvalues.save(inputFile + ".Eigenvalues.txt");

		//Calculate principal components:
		ExpressionDataset datasetPCs = new ExpressionDataset(dataset.nrSamples, nrPCsToOutput);
		for (int pca = 0; pca < nrPCsToOutput; pca++) {
			datasetPCs.sampleNames[pca] = "Comp" + (pca + 1);
		}
		for (int p = 0; p < datasetPCs.nrProbes; p++) {
			datasetPCs.probeNames[p] = dataset.sampleNames[p];
		}
		for (int pca = 0; pca < nrPCsToOutput; pca++) {
			for (int p = 0; p < dataset.nrProbes; p++) {
				for (int s = 0; s < dataset.nrSamples; s++) {
					datasetPCs.rawData[s][pca] += datasetEigenvectors.rawData[p][pca] * dataset.rawData[p][s];
				}
			}
		}
		datasetPCs.save(dataset.fileName + ".PrincipalComponents.txt");

		ExpressionDataset datasetFactorloadings = new ExpressionDataset(correlationMatrix.length, nrPCsToOutput);
		datasetPCs.transposeDataset();
		for (int p = 0; p < dataset.nrProbes; p++) {
			datasetFactorloadings.probeNames[p] = dataset.probeNames[p];
		}
		for (int pca = 0; pca < nrPCsToOutput; pca++) {
			datasetFactorloadings.sampleNames[pca] = "Comp" + (pca + 1);
			for (int p = 0; p < dataset.nrProbes; p++) {
				datasetFactorloadings.rawData[p][pca] = JSci.maths.ArrayMath.correlation(datasetPCs.rawData[pca], dataset.rawData[p]);
			}
		}
		datasetFactorloadings.save(dataset.fileName + ".Factorloadings.txt");

	}

	public ExpressionDataset orthogonalizeMatrix(ExpressionDataset dataset) {

		dataset.standardNormalizeData();
		int nrVars = dataset.nrProbes;
		int nrSamples = dataset.nrSamples;
		double[][] matrix = new double[nrVars][nrSamples];
		for (int s = 0; s < nrVars; s++) {
			for (int sample = 0; sample < nrSamples; sample++) {
				matrix[s][sample] = dataset.rawData[s][sample];
			}
		}
		double[][] correlationMatrix = new double[nrVars][nrVars];
		for (int p = 0; p < nrVars; p++) {
			correlationMatrix[p][p] = 1d;
			for (int q = p + 1; q < nrVars; q++) {
				double covariance = 0;
				for (int sample = 0; sample < nrSamples; sample++) {
					covariance += matrix[p][sample] * matrix[q][sample];
				}
				covariance /= (double) (nrSamples - 1);
				if (covariance > 1) {
					covariance = 1d;
				}
				if (covariance < -1) {
					covariance = -1d;
				}
				correlationMatrix[p][q] = covariance;
				correlationMatrix[q][p] = covariance;
			}
		}
		Jama.EigenvalueDecomposition eig = eigenValueDecomposition(correlationMatrix);
		double[] eigenValues = eig.getRealEigenvalues();
		int nrCompsWithPositiveEigenvalues = 0;
		for (int e = 0; e < eigenValues.length; e++) {
			//System.out.println(e + "\t" + eigenValues[e]);
			if (eigenValues[e] > 1e-7) {
				nrCompsWithPositiveEigenvalues++;
			}
		}

		ExpressionDataset datasetEigenvectors = new ExpressionDataset(correlationMatrix.length, correlationMatrix.length);
		for (int pca = 0; pca < correlationMatrix.length; pca++) {
			datasetEigenvectors.rawData[pca] = getEigenVector(eig, pca);
		}
		datasetEigenvectors.transposeDataset();

		//Calculate principal components:
		ExpressionDataset datasetPCs = new ExpressionDataset(dataset.nrSamples, nrCompsWithPositiveEigenvalues);
		for (int pca = 0; pca < nrCompsWithPositiveEigenvalues; pca++) {
			datasetPCs.sampleNames[pca] = "Comp" + (pca + 1);
		}
		for (int p = 0; p < datasetPCs.nrProbes; p++) {
			datasetPCs.probeNames[p] = dataset.sampleNames[p];
		}
		for (int pca = 0; pca < nrCompsWithPositiveEigenvalues; pca++) {
			for (int p = 0; p < dataset.nrProbes; p++) {
				for (int s = 0; s < dataset.nrSamples; s++) {
					datasetPCs.rawData[s][pca] += datasetEigenvectors.rawData[p][pca] * dataset.rawData[p][s];
				}
			}
		}
		datasetPCs.transposeDataset();
		return datasetPCs;

	}

	private Jama.EigenvalueDecomposition eigenValueDecomposition(double[][] data) {
		Jama.Matrix m = new Jama.Matrix(data);
		Jama.EigenvalueDecomposition eig = m.eig();
		return eig;
	}

	private double[] getEigenVector(Jama.EigenvalueDecomposition eig, double[] eigenValues, int pca) {
		Jama.Matrix eigenValueMatrix = eig.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
		}
		return eigenVector;
	}

	private double[] getEigenVector(Jama.EigenvalueDecomposition eig, int pca) {
		Jama.Matrix eigenValueMatrix = eig.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
		}
		return eigenVector;
	}

	private double getEigenValueVar(double[] eigenValues, int pca) {
		double sumEigenvalues = 0.0;
		for (Double d : eigenValues) {
			sumEigenvalues += Math.abs(d);
		}
		double result = eigenValues[eigenValues.length - 1 - pca] / sumEigenvalues;
		return result;
	}

	private double[] getEigenVectorSVD(Jama.SingularValueDecomposition svd, double[] singularValues, int pca) {
		Jama.Matrix eigenValueMatrix = svd.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][pca] * Math.sqrt(singularValues[pca]);
		}
		return eigenVector;
	}

	public double[] getLinearRegressionCoefficients(double[] xVal, double[] yVal) {
		double n = (double) xVal.length;
		double sumX = 0;
		double sumXX = 0;
		double sumY = 0;
		double sumXY = 0;
		for (int x = 0; x < xVal.length; x++) {
			sumX += xVal[x];
			sumXX += xVal[x] * xVal[x];
			sumY += yVal[x];
			sumXY += xVal[x] * yVal[x];
		}
		double sXX = sumXX - sumX * sumX / n;
		double sXY = sumXY - sumX * sumY / n;
		double a = sXY / sXX;
		double b = (sumY - a * sumX) / n;
		double[] regressionCoefficients = new double[2];
		regressionCoefficients[0] = a;
		regressionCoefficients[1] = b;
		return regressionCoefficients;
	}

	private double mean(double[] v) {
		double sum = 0;
		for (int k = 0; k < v.length; k++) {
			sum += v[k];
		}
		return (sum / (double) v.length);
	}

	private double variance(double[] v, double mean) {
		double ans = 0.0;
		for (int i = 0; i < v.length; i++) {
			ans += (v[i] - mean) * (v[i] - mean);
		}
		return ans / (v.length - 1);
	}

	public double weightedCorrelation(double[] x, double[] y, double[] weights) {
		double wmX = weightedMean(x, weights);
		double wmY = weightedMean(y, weights);
		double sumWeights = JSci.maths.ArrayMath.mass(weights);
		double covXX = 0;
		double covXY = 0;
		double covYY = 0;
		for (int s = 0; s < x.length; s++) {
			covXX += weights[s] * (x[s] - wmX) * (x[s] - wmX);
			covXY += weights[s] * (x[s] - wmX) * (y[s] - wmY);
			covYY += weights[s] * (y[s] - wmY) * (y[s] - wmY);
		}
		covXX /= sumWeights;
		covXY /= sumWeights;
		covYY /= sumWeights;
		double corr = covXY / (Math.sqrt(covXX * covYY));
		return corr;
	}

	public double weightedMean(double[] x, double[] weights) {
		double m = 0;
		double sumWeights = 0;
		for (int s = 0; s < x.length; s++) {
			m += x[s] * weights[s];
			sumWeights += weights[s];
		}
		return m / sumWeights;
	}

}
