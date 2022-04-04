package betaqtl.junk;//package betaqtl.junk;
//
//import betaqtl.*;
//import betaqtl.vcf.VCFTabix;
//import betaqtl.vcf.VCFVariant;
//import umcg.genetica.console.ProgressBar;
//import umcg.genetica.containers.Triple;
//import umcg.genetica.enums.Chromosome;
//import umcg.genetica.features.Feature;
//import umcg.genetica.io.text.TextFile;
//import umcg.genetica.math.stats.Correlation;
//import umcg.genetica.math.stats.ZScores;
//import umcg.genetica.util.RankArray;
//import umontreal.iro.lecuyer.probdist.BetaDist;
//
//import java.io.IOException;
//import java.lang.reflect.Array;
//import java.util.*;
//import java.util.concurrent.atomic.AtomicBoolean;
//import java.util.stream.IntStream;
//
//public class BetaQTL2Partest2 extends QTLAnalysis {
//
//	private int nrPermutations = 1000;
//	private int cisWindow = 1000000;
//	private long randomSeed = 123456789;
//	private boolean rankData = true;
//	private boolean outputAll = false;
//	private boolean replaceMissingGenotypes = false;
//	private VCFTabix tabix;
//
//
//	public BetaQTL2Partest2(String vcfFile, int chromosome, String linkfile, String snpLimitFile, String geneLimitFile, String snpGeneLimitFile, String geneExpressionDataFile, String geneAnnotationFile, String outfile) throws IOException {
//		super(vcfFile, chromosome, linkfile, snpLimitFile, geneLimitFile, snpGeneLimitFile, geneExpressionDataFile, geneAnnotationFile, outfile);
//	}
//
//	public void setNrPermutations(int nrPermutations) {
//		this.nrPermutations = nrPermutations;
//	}
//
//	public void setCisWindow(int cisWindow) {
//		this.cisWindow = cisWindow;
//	}
//
//	public void setRandomSeed(long randomSeed) {
//		this.randomSeed = randomSeed;
//	}
//
//	public void setRankData(boolean rankData) {
//		this.rankData = rankData;
//	}
//
//	public void setOutputAll(boolean outputAll) {
//		this.outputAll = outputAll;
//	}
//
//	public void setReplaceMissingGenotypes(boolean replaceMissingGenotypes) {
//		this.replaceMissingGenotypes = replaceMissingGenotypes;
//	}
//
//
//	public void run() throws IOException {
//        /*
//        TODO:
//        - function to test specific SNPs
//        - output dataset Z-scores and sample sizes
//        - plotting?
//        - check whether code runs with no permutations
//        - a way to run on a VCF containing all chromosomes
//        - check proper way to shuffle when there is missing data in the genotype as well as the phenotype data
//        */
//
//		System.out.println();
//		System.out.println("-----------------------------------------------------");
//		System.out.println("QTL analysis with Beta distribution approximated null");
//		System.out.println("-----------------------------------------------------");
//		System.out.println("MAF:\t" + mafthreshold);
//		System.out.println("Call-rate:\t" + callratethreshold);
//		System.out.println("HWE-P:\t" + hwepthreshold);
//		System.out.println("Min nr datasets:\t" + minNumberOfDatasets);
//		System.out.println("Nr Permutations:\t" + nrPermutations);
//		System.out.println("Cis window:\t" + cisWindow);
//		System.out.println("Random seed:\t" + randomSeed);
//		System.out.println("Ranking data:\t" + rankData);
////        System.out.println("Re-ranking data:\t");
//		System.out.println("Replacing missing genotypes: " + replaceMissingGenotypes);
//		System.out.println("Min observations: " + minObservations);
//		System.out.println("Writing all snp/feature pairs: " + outputAll);
//		System.out.println();
//
//		Chromosome chromosomeObj = Chromosome.parseChr("" + chromosome);
//		System.out.println("Processing: " + vcfFile);
//
//
//		// initialize permutation seeds
//		long[] seed = new long[nrPermutations];
//		Random rand = new Random(randomSeed);
//		for (int i = 0; i < seed.length; i++) {
//			seed[i] = rand.nextLong();
//		}
//
//		// initialize output
//		ProgressBar pb = new ProgressBar(expressionData.genes.length, "Processing " + expressionData.genes.length + " genes...");
//		TextFile outTopFx = new TextFile(outputPrefix + "-TopEffects.txt", TextFile.W);
//		String header = "Gene\t" +
//				"GeneChr\t" +
//				"GenePos\t" +
//				"GeneSymbol\t" +
//				"NrTestedSNPs\t" +
//				"BestSNP\t" +
//				"BestSNPChr\t" +
//				"BestSNPPos\t" +
//				"BestSNPAlleles\t" +
//				"BestSNPEffectAllele\t" +
//				"LowestMetaP\t" +
//				"LowestMetaPN\t" +
//				"LowestMetaPZ\t" +
//				"LowestMetaBeta\t" +
//				"LowestMetaBetaSE\t" +
//				"LowestMetaPNrDatasets\t" +
//				"ProportionBetterPermPvals\t" +
//				"BetaDistAlpha\t" +
//				"BetaDistBeta\t" +
//				"BetaAdjustedMetaP";
//		outTopFx.writeln(header);
//		TextFile outAll = null;
//		if (outputAll) {
//			outAll = new TextFile(outputPrefix + "-AllEffects.txt.gz", TextFile.W);
////            String headerAll = "Gene\tGeneSymbol\tSNP\tSNPAlleles\tSNPEffectAllele\tMetaP\tMetaPN\tMetaPZ\tMetaBeta\tMetaSE\tNrDatasets\tProportionBetterPermPvals\tBetaAdjustedMetaP";
//			String headerAll = "Gene\t" +
//					"GeneChr\t" +
//					"GenePos\t" +
//					"GeneSymbol\t" +
//					"SNP\t" +
//					"SNPAlleles\t" +
//					"SNPEffectAllele\t" +
//					"MetaP\t" +
//					"MetaPN\t" +
//					"MetaPZ\t" +
//					"MetaBeta\t" +
//					"MetaSE\t" +
//					"NrDatasets";
//			outAll.writeln(headerAll);
//		}
//
//		// iterate genes
//		TextFile finalOutAll = outAll;
//		tabix = new VCFTabix(vcfFile);
//		IntStream.range(0, expressionData.genes.length).parallel().forEach(g -> {
//			try {
//				String gene = expressionData.genes[g];
//
//				Integer geneAnnotationId = geneAnnotation.getGeneId(gene);
//				if (geneAnnotationId != null) {
//					double[] expData = expressionData.data[g];
//
//					// define CIS window
//					int pos = geneAnnotation.getPos(geneAnnotationId);
//					String geneSymbol = geneAnnotation.getSymbol(geneAnnotationId);
//					int start = pos - cisWindow;
//					if (start < 0) {
//						start = 0;
//					}
//					int stop = pos + cisWindow;
//					Feature cisRegion = new Feature(chromosomeObj, start, stop);
//
//					// iterate SNPs
//					double[] permutationPvals = new double[nrPermutations];
//					Arrays.fill(permutationPvals, 1);
//
//					int nrTestedSNPs = 0;
//
//					// split expression data per dataset
//					double[][] expressionPerDataset = new double[datasets.length][];
//					IntStream.range(0, datasets.length).forEach(d -> {
//						Dataset thisDataset = datasets[d];
//						double[] datasetExp = thisDataset.select(expData, thisDataset.expressionIds);
//						double[] datasetExpRanked = datasetExp;
//						if (rankData) {
//							RankArray ranker = new RankArray();
//							datasetExpRanked = ranker.rank(datasetExp, true); // does this work with NaNs? answer: no
//						}
//						expressionPerDataset[d] = datasetExpRanked;
//					});
//
//					final UnpermutedResult topUnpermutedResult = new UnpermutedResult(); // this is safe, because values are only changed once per SNP, when permutation == -1
//
//					ArrayList<VCFVariant> variants = getVariants(cisRegion, genotypeSamplesToInclude, snpLimitSet, gene);
//
//					for (VCFVariant variant : variants) {
////						VCFVariant variant = snpIterator.next();
//						if (variant != null || (snpGeneLimitSet != null && snpGeneLimitSet.get(gene) != null && snpGeneLimitSet.get(gene).contains(variant.getId()))) {
//
//							final double[] genotypes = getGenotype(variant.getGenotypesAsByteVector());
//							final double[] dosages = getDosage(variant.getDosage());
//
//							// split genotype data per dataset, perform QC
//							double[][] genotypesPerDataset = new double[datasets.length][];
//							double[][] dosagesPerDataset = new double[datasets.length][];
//							VariantQCObj[] qcobjs = new VariantQCObj[datasets.length];
//							IntStream.range(0, datasets.length).forEach(d -> {
//								Dataset thisDataset = datasets[d];
//								dosagesPerDataset[d] = thisDataset.select(dosages, thisDataset.genotypeIds); // select required dosages
//								genotypesPerDataset[d] = thisDataset.select(genotypes, thisDataset.genotypeIds); // select required genotype IDs
//
//								VariantQCObj qcobj = checkVariant(genotypesPerDataset[d]);
//								if (qcobj.passqc) {
//									if (replaceMissingGenotypes) {
//										// only replace missing genotypes on variants that pass the qc thresholds
//										double meanDosage = JSci.maths.ArrayMath.mean(dosagesPerDataset[d]);
//										double meanGenotype = Math.round(JSci.maths.ArrayMath.mean(genotypesPerDataset[d]));
//										for (int i = 0; i < dosagesPerDataset[d].length; i++) {
//											if (genotypesPerDataset[d][i] == -1) {
//												genotypesPerDataset[d][i] = meanGenotype;
//												dosagesPerDataset[d][i] = meanDosage;
//											}
//										}
//									}
//
//									// prune the data here once, to check if there are enough values to go ahead with this snp/gene combo
//									// but only if the variant is passing the QC in the first place for the samples selected in this dataset
//									Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(genotypesPerDataset[d],
//											dosagesPerDataset[d],
//											expressionPerDataset[d]);
//
//									// check the variant again, taking into account missingness in the expression data
//									qcobj = checkVariant(prunedDatasetData.getLeft());
//
//									// require minimum number of observations, otherwise kick out dataset from analysis
//									if (prunedDatasetData.getLeft().length < minObservations) {
//										qcobj.passqc = false;
//									}
//								}
//								qcobjs[d] = qcobj;
//							});
//
//							// run permutations, and non-permuted result (permutation == -1)
//							AtomicBoolean tested = new AtomicBoolean(false);
////							TextFile finalOutAll = outAll;
//							IntStream.range(-1, nrPermutations).forEach(permutation -> {
//								double[] zscores = new double[datasets.length];
//								int[] samplesizes = new int[datasets.length];
//
//								// iterate datasets
//								int nrAltAlleles = 0;
//								int nrTotalAlleles = 0;
//								for (int d = 0; d < datasets.length; d++) {
//									Dataset thisDataset = datasets[d];
//									double[] datasetGt = genotypesPerDataset[d]; // thisDataset.select(genotypes, thisDataset.genotypeIds); // select required genotype IDs
//									VariantQCObj qcobj = qcobjs[d]; // check maf, hwep, call-rate, number of genotypes per genotype group
//									if (qcobj.passqc) {
//
//										double[] datasetExp = expressionPerDataset[d];
//										double[] datasetExpCopy = new double[datasetExp.length];
//										System.arraycopy(datasetExp, 0, datasetExpCopy, 0, datasetExpCopy.length);
//
//										double[] datasetDs = dosagesPerDataset[d];
//
//										// if this is a permutation, shuffle the data
//										if (permutation != -1) {
//											Util.shuffleArray(datasetExpCopy, seed[permutation]);
//										}
//
//										// prune the data (remove missing values)
//										// can't prune the data earlier (would save a lot of compute time) because shuffling is performed over all available samples for this dataset
//										// this is because the order of permuted samples should be equal across all SNPs
//										Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(datasetGt,
//												datasetDs,
//												datasetExpCopy);
//
//										// re-rank data here? original EMP does not, but it is the right thing to do...
//										double[] datasetExpPruned = prunedDatasetData.getRight();
////                                    if (rankData) {
////                                        RankArray ranker = new RankArray();
////                                        datasetExpPruned = ranker.rank(datasetExpPruned, true); // does this work with NaNs? answer: no
////                                    }
//										datasetExpPruned = Util.centerScale(datasetExpPruned);
//										double[] datasetDsPruned = Util.centerScale(prunedDatasetData.getMiddle());
//										double[] datasetGtPruned = Util.centerScale(prunedDatasetData.getLeft());
//
//
//										// count the number of alleles, used later to estimate Beta and SE from MetaZ
//										if (permutation == -1) {
//											for (int i = 0; i < datasetGtPruned.length; i++) {
//												if (datasetGtPruned[i] >= 0.5 && datasetGtPruned[i] < 1.5) {
//													nrAltAlleles += 1;
//												} else if (datasetGtPruned[i] > 1.5) {
//													nrAltAlleles += 2;
//												}
//											}
//											nrTotalAlleles += datasetGtPruned.length * 2;
//										}
//
//										// perform correlation
//										double r = Correlation.correlate(datasetDsPruned, datasetExpPruned);
//										double p = PVal.getPvalue(r, datasetExpPruned.length - 2);
//										double z = ZScores.pToZTwoTailed(p); // p value is already two-tailed, so need to use this other p-value conversion method... :/; returns negative z-scores by default
//										if (r > 0) {
//											z *= -1; // flip z-score if correlation is positive because p-value conversion returns only negative z-scores
//										}
//										zscores[d] = z;
//										samplesizes[d] = datasetExpPruned.length;
//									} else {
//										zscores[d] = Double.NaN;
//										samplesizes[d] = 0;
//									}
//								} // ENDIF: test every dataset
//
//								// determine number of datasets with data
//								int nDatasets = 0;
//								int n = 0;
//								for (int d = 0; d < zscores.length; d++) {
//									n += samplesizes[d];
//									if (!Double.isNaN(zscores[d])) {
//										nDatasets++;
//									}
//								}
//								if (nDatasets >= minNumberOfDatasets) {
//									// meta-analyze, weight by sample size
//									double metaZ = ZScores.getWeightedZ(zscores, samplesizes);
//									double metaP = ZScores.zToP(metaZ);
//									if (permutation != -1) { // this is a permuted result
//										if (metaP < permutationPvals[permutation]) {
//											permutationPvals[permutation] = metaP;
//										}
//									} else { // this is a non-permuted result
//										tested.getAndSet(true);
//
//										// calculate overall MAF
//										double overallMaf = (double) nrAltAlleles / nrTotalAlleles;
//										double[] betaAndSEEstimate = ZScores.zToBeta(metaZ, overallMaf, n);
//
//										// non-permuted p-value
//										if (metaP <= topUnpermutedResult.metaP) {
//											boolean replace = true;
//											// if the SNP is in perfect LD (has equal pvalue), select the closest one to the gene
//											if (metaP == topUnpermutedResult.metaP && topUnpermutedResult.snpID != null) {
//												// if the Z-score is sufficiently large (>40) we exceed the range of the normal distribution, returning a p-value of ~2x10-232
//												// in that case, compare the absolute Z-scores to determine the top effect for this gene
//												if (Math.abs(metaZ) < Math.abs(topUnpermutedResult.metaPZ)) {
//													replace = false;
//												} else {
//													int genePos = geneAnnotation.getPos(geneAnnotationId);
//													int tssDist = Math.abs(genePos - variant.getPos());
//													int tssDist2 = Math.abs(genePos - topUnpermutedResult.snpPos);
//													if (tssDist > tssDist2) {
//														replace = false;
//													}
//												}
//											}
//											if (replace) {
//												topUnpermutedResult.metaP = metaP;
//												topUnpermutedResult.metaPN = n;
//												topUnpermutedResult.metaPZ = metaZ;
//												topUnpermutedResult.metaPD = nDatasets;
//												topUnpermutedResult.metaBeta = betaAndSEEstimate[0];
//												topUnpermutedResult.metaBetaSE = betaAndSEEstimate[1];
//												topUnpermutedResult.snpID = variant.getId();
//												topUnpermutedResult.snpPos = variant.getPos();
//												topUnpermutedResult.snpAlleles = variant.getAlleles()[0] + "/" + variant.getAlleles()[1];
//												topUnpermutedResult.snpEffectAllele = variant.getAlleles()[1];
//											}
//										}
//										if (outputAll) { // this code only runs when in the 'not-permuted' iteration
//											String snpAlleles = variant.getAlleles()[0] + "/" + variant.getAlleles()[1];
//											String snpEffectAllele = variant.getAlleles()[1];
//
//											String outln = gene + "\t" + chromosome + "\t" + pos + "\t" + geneSymbol + "\t" + variant.getId() + "\t" + chromosome + "\t" + variant.getPos() + "\t" + snpAlleles + "\t" + snpEffectAllele
//													+ "\t" + metaP + "\t" + n + "\t" + metaZ + "\t" + betaAndSEEstimate[0] + "\t" + betaAndSEEstimate[1] + "\t" + nDatasets;
//
//											try {
//												finalOutAll.writelnsynced(outln);
//											} catch (IOException e) {
//												e.printStackTrace();
//											}
//
//										}
//									}
//								}
//							}); // ENDIF: for each permutation in parallel
//
//							// determine if SNP was tested somehow...
//							if (tested.get()) {
//								nrTestedSNPs++;
//							}
//						} // ENDIF: if variant != null
//					} // ENDIF: while snpiterator has next
//
//					if (nrTestedSNPs == 0) {
//						System.out.println(gene + " has zero SNPs tested?");
//					} else {
//						// determine beta distribution etc
//						double propBetterPvals = 0;
//						double betaAdjPval = 1;
//						double[] shape = new double[]{Double.NaN, Double.NaN};
//						BetaDist betaDistribution = null;
//						if (nrPermutations > 1) { // permutations are required for the following step
//							for (int p = 0; p < permutationPvals.length; p++) {
//								if (permutationPvals[p] < topUnpermutedResult.metaP) {
//									propBetterPvals++;
//								}
//							}
//							propBetterPvals /= permutationPvals.length;
//							BetaDistributionMLE mle = new BetaDistributionMLE();
//							shape = mle.fit(permutationPvals);
//							betaDistribution = new BetaDist(shape[0], shape[1]);
//							betaAdjPval = betaDistribution.cdf(topUnpermutedResult.metaP);
//							if (betaAdjPval < 2.0E-323D) {
//								betaAdjPval = 2.0E-323D;
//							}
//						} else {
//							propBetterPvals = 1;
//						}
//
//						String outln = gene + "\t" + chromosome + "\t" + pos + "\t" + geneSymbol + "\t" + nrTestedSNPs
//								+ "\t" + topUnpermutedResult.snpID + "\t" + chromosome + "\t" + topUnpermutedResult.snpPos + "\t" + topUnpermutedResult.snpAlleles + "\t" + topUnpermutedResult.snpEffectAllele
//								+ "\t" + topUnpermutedResult.metaP + "\t" + topUnpermutedResult.metaPN + "\t" + topUnpermutedResult.metaPZ
//								+ "\t" + topUnpermutedResult.metaBeta + "\t" + topUnpermutedResult.metaBetaSE + "\t" + topUnpermutedResult.metaPD
//								+ "\t" + propBetterPvals + "\t" + shape[0] + "\t" + shape[1] + "\t" + betaAdjPval;
//
//						outTopFx.writelnsynced(outln);
//
//					}
//				} // ENDIF: if the gene has an annotation
//
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
////			pb.set(g + 1);
//			pb.iterateSynched();
//			pb.print();
//
//		}); // ENDIF: iterate genes
//		tabix.close();
//		outTopFx.close();
//		if (outputAll) {
//			outAll.close();
//		}
//		pb.close();
//	}
//
//	int getVariantCtr = 1;
//
//	private synchronized ArrayList<VCFVariant> getVariants(Feature cisRegion, boolean[] genotypeSamplesToInclude, Set<String> snpLimitSet, String gene) throws IOException {
//		System.out.println();
//		System.out.println(getVariantCtr + " - Getting variants: " + cisRegion);
//		Iterator<VCFVariant> snpIterator = tabix.getVariants(cisRegion, genotypeSamplesToInclude, snpLimitSet);
//		ArrayList<VCFVariant> variants = new ArrayList<>();
//		while (snpIterator.hasNext()) {
//			VCFVariant variant = snpIterator.next();
//			if (variant != null || (snpGeneLimitSet != null && snpGeneLimitSet.get(gene) != null && snpGeneLimitSet.get(gene).contains(variant.getId()))) {
//				variants.add(variant);
//			}
//		}
//		System.out.println(getVariantCtr + " - Getting variants: " + cisRegion + " --> " + variants.size());
//		System.out.println();
//		getVariantCtr++;
//		return variants;
//	}
//
//	private class UnpermutedResult {
//		public double metaBeta;
//		public double metaBetaSE;
//		public String snpID;
//		public String snpAlleles;
//		public String snpEffectAllele;
//		public int snpPos;
//		double metaP = 1;
//		double metaPN = 0;
//		double metaPZ = 0;
//		double metaPD = 0;
//
//		@Override
//		public String toString() {
//			return "UnpermutedResult{" +
//					"metaP=" + metaP +
//					", metaPN=" + metaPN +
//					", metaPZ=" + metaPZ +
//					", metaPD=" + metaPD +
//					", snpID=" + snpID +
//					'}';
//		}
//
//	}
//
//}
