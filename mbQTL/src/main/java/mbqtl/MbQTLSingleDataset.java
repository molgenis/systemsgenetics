package mbqtl;

import mbqtl.datastructures.Dataset;
import mbqtl.stat.BetaDistributionMLE;
import mbqtl.stat.PVal;
import mbqtl.vcf.VCFTabix;
import mbqtl.vcf.VCFVariant;

import umcg.genetica.containers.Triple;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.util.RankArray;
import umcg.genetica.util.RunTimer;
import umontreal.iro.lecuyer.probdist.BetaDist;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.stream.IntStream;

public class MbQTLSingleDataset extends QTLAnalysis {

	private int nrPermutations = 1000;
	private int cisWindow = 1000000;
	private long randomSeed = 123456789;

	public MbQTLSingleDataset(String vcfFile, int chromosome, String linkfile, String geneLimitFile, String geneExpressionDataFile, String geneAnnotationFile, String outfile) throws IOException {
		super(vcfFile, chromosome, linkfile, null, geneLimitFile, null, geneExpressionDataFile, geneAnnotationFile, outfile);
		System.out.println("Running beta-analysis for single dataset");
		System.out.println("Permutations: " + nrPermutations);
		System.out.println("Ciswindow: " + cisWindow);

	}

	public void fastqtlclone2() throws IOException {

		// iterate genes
		Chromosome chromosomeObj = Chromosome.parseChr("" + chromosome);
		System.out.println("Processing: " + vcfFile);
		VCFTabix tabix = new VCFTabix(vcfFile);

		TextFile output = new TextFile(outputPrefix + "TopFxPerGene.txt", TextFile.W);
		TextFile outputAll = new TextFile(outputPrefix + "AllAssociations.txt", TextFile.W);
		outputAll.writeln("Gene\tSNP\tN\tR\tP");

		String header = "Gene\tNrGenotypes\tBestSNP\tBestSNPPval\tBestPvalN\tBestPvalR\tBdist-Alpha\tBdist-Beta\tProportionBetterPermPvals\tBetaApproxPval";
		output.writeln(header);
		long[] seed = new long[nrPermutations];
		Random rand = new Random(randomSeed);
		for (int i = 0; i < seed.length; i++) {
			seed[i] = rand.nextLong();
		}

		// iterate genes
		RunTimer timer = new RunTimer();
		timer.start();

		for (int g = 0; g < expressionData.genes.length; g++) {
			String gene = expressionData.genes[g];

//            Integer geneAnnotationId = geneAnnotation.getGeneId(gene);
			Gene geneAnnotationObj = geneAnnotation.getGene(gene);
			if (geneAnnotationObj != null) {
				double[] expData = expressionData.data[g];

				// rank once, then center+scale
				int d = 0;
				Dataset thisDataset = datasets[d];
				double[] datasetExpressionData = thisDataset.select(expData, thisDataset.getExpressionIds());
				RankArray ranker = new RankArray();
				double[] expDataPerDatasetRanked = ranker.rank(datasetExpressionData, true);

				// get variants, 1mb up and downstream
				int pos = geneAnnotationObj.getStart(); //.getStartPos(geneAnnotationId);
				int start = pos - cisWindow;
				if (start < 0) {
					start = 0;
				}
				int stop = pos + cisWindow;

				int vctr = 0;
				int tctr = 0;

				double bestPval = 1;
				double bestPvalR = 0;
				double bestPvalN = 0;
				String bestPvalSNP = null;

				double[] bestPermPVals = new double[nrPermutations];
				Arrays.fill(bestPermPVals, 1);
				double[] bestPermPValsN = new double[nrPermutations];

				// load genotype data
				ArrayList<double[]> genotypes = new ArrayList<>();
				ArrayList<double[]> genotypesDosages = new ArrayList<>();
				ArrayList<String> genotypesIds = new ArrayList<>();
				Feature region = new Feature(chromosomeObj, start, stop);
				Iterator<VCFVariant> snpIterator = tabix.getVariants(region, genotypeSamplesToInclude);
				while (snpIterator.hasNext()) {
					VCFVariant variant = snpIterator.next();
					if (variant != null) {
						String snpid = variant.getId();
						// collect data per dataset
						double[] datasetGenotypeData = thisDataset.select(variant.getGenotypesAsByteVector(), thisDataset.getGenotypeIds());

						// do some QC checks: test MAF, Callrate, HWE-P
						VariantQCObj qcobj = checkVariant(datasetGenotypeData);
						double[] datasetDosageValues = thisDataset.select(variant.getDosage(), thisDataset.getGenotypeIds());

						if (qcobj.passqc) {
							genotypes.add(datasetGenotypeData);
							genotypesDosages.add(datasetDosageValues);
							genotypesIds.add(snpid);
						}
					}
				}

				// nominal pass
				// final int df = datasets[0].expressionIds.length - 2;
				for (int gt = 0; gt < genotypes.size(); gt++) {
					double[] datasetGenotypeData = genotypes.get(gt);
					double[] datasetGenotypeDosageData = genotypesDosages.get(gt);
					Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(datasetGenotypeData,
							datasetGenotypeDosageData,
							expDataPerDatasetRanked);
					double[] expReRanked = prunedDatasetData.getRight();
					double[] expDataCenterScale = Util.centerScale(expReRanked);
					double[] gtDataCenterScale = Util.centerScale(prunedDatasetData.getMiddle());

					// perform correlation
					double r = Correlation.correlate(expDataCenterScale, gtDataCenterScale);
					double p = PVal.getPvalue(r, expReRanked.length - 2);
					if (p < bestPval) {
						bestPval = p;
						bestPvalR = r;
						bestPvalN = expData.length;
						bestPvalSNP = genotypesIds.get(gt);
					}
					outputAll.writeln(gene + "\t" + genotypesIds.get(gt) + "\t" + expData.length + "\t" + r + "\t" + p);
				}

				// permutation passes
				double[] finalExpDataPerDatasetRanked = expDataPerDatasetRanked;
				IntStream.range(0, nrPermutations).parallel().forEach(perm -> {
					double[] exp = new double[finalExpDataPerDatasetRanked.length];
					System.arraycopy(finalExpDataPerDatasetRanked, 0, exp, 0, exp.length);
					Util.shuffleArray(exp, seed[perm]);

					for (int gt = 0; gt < genotypes.size(); gt++) {
						double[] datasetsGenotypeData = genotypes.get(gt);
						double[] datasetsGenotypeDosageData = genotypesDosages.get(gt);
						Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(datasetsGenotypeData,
								datasetsGenotypeDosageData,
								exp);
						double[] expReRanked = prunedDatasetData.getRight();
						double[] expDataCenterScale = Util.centerScale(expReRanked);
						double[] gtDataCenterScale = Util.centerScale(prunedDatasetData.getMiddle());
						double rperm = Correlation.correlate(expDataCenterScale, gtDataCenterScale);
						double pperm = PVal.getPvalue(rperm, exp.length - 2);

						if (pperm < bestPermPVals[perm]) {
							bestPermPVals[perm] = pperm;
							bestPermPValsN[perm] = exp.length;
						}
					}
				});

				System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + genotypes.size() + " variants tested. T: " + timer.getTimeDesc() + "\n");
				if (genotypes.size() > 0) {

					TextFile geneoutput = new TextFile(outputPrefix + "GenePerm-" + gene + ".txt", TextFile.W);
					geneoutput.writeln("Perm\tBestPval\tBestPValN");
					for (int p = 0; p < nrPermutations; p++) {
						geneoutput.writeln(p + "\t" + bestPermPVals[p] + "\t" + bestPermPValsN[p]);
					}
					geneoutput.close();

					// calculate permutation pvalues with adjusted DF
//					double trueDF = thisDataset.expressionIds.length - 2;
//					double varPermCor = JSci.maths.ArrayMath.variance(bestPermCorrelations);
//					DFMLE dfmle = new DFMLE();
//					if (varPermCor > 0) {
//						// learn actual DF
//						trueDF = dfmle.learnDF(bestPermCorrelations, trueDF);
//					}
//					System.out.println("DF: " + trueDF + "\t" + (thisDataset.expressionIds.length - 2));
//					double[] permP = new double[bestPermCorrelations.length];
//					double propBetterCorrel = 0;
//					for (int p = 0; p < bestPermCorrelations.length; p++) {
//						permP[p] = PVal.getPvalue(bestPermCorrelations[p], trueDF);
//						if (Math.abs(bestPermCorrelations[p]) >= Math.abs(bestCorrelation)) {
//							propBetterCorrel++;
//						}
//					}
//					propBetterCorrel /= nrPermutations;

//					double meanPermP = JSci.maths.ArrayMath.mean(permP);
//					double varPermP = JSci.maths.ArrayMath.variance(permP);
//					double beta_shape1 = meanPermP * (meanPermP * (1 - meanPermP) / varPermP - 1);
//					double beta_shape2 = beta_shape1 * (1 / meanPermP - 1);
//					System.out.println("Shape1 " + beta_shape1);
//					System.out.println("Shape2 " + beta_shape2);
//					BetaDist bd = BetaDist.getInstanceFromMLE(permP, permP.length);
//					BetaDistributionMLE bdmle = new BetaDistributionMLE();
//					double[] betaParams = bdmle.fit(permP, beta_shape1, beta_shape2);
//					BetaDist bdmledist = new BetaDist(betaParams[0], betaParams[1]);
//
//					double pval_fdo = PVal.getPvalue(bestCorrelation, trueDF);
//					double pval_nom = PVal.getPvalue(bestCorrelation, thisDataset.expressionIds.length - 2);
//					double betaApproxPvalAdj = bd.cdf(pval_fdo);
//					double betaMLEApproxPvalAdj = bdmledist.cdf(pval_fdo);
//
//					double proportionBetterPvalsAdj = 0;
//					for (int q = 0; q < permP.length; q++) {
//						if (permP[q] < pval_fdo) {
//							proportionBetterPvalsAdj++;
//						}
//					}
//					proportionBetterPvalsAdj /= permP.length;

					// the same with unadjusted DF
					double proportionBetterPvals = 0;
					for (int q = 0; q < bestPermPVals.length; q++) {
						if (bestPermPVals[q] <= bestPval) {
							proportionBetterPvals++;
						}
					}
					proportionBetterPvals /= bestPermPVals.length;

					double meanPermP2 = JSci.maths.ArrayMath.mean(bestPermPVals);
					double varPermP2 = JSci.maths.ArrayMath.variance(bestPermPVals);
					double beta_shape12 = meanPermP2 * (meanPermP2 * (1 - meanPermP2) / varPermP2 - 1);
					double beta_shape22 = beta_shape12 * (1 / meanPermP2 - 1);
					BetaDistributionMLE bdmle = new BetaDistributionMLE();
					double[] betamleparams2 = bdmle.fit(bestPermPVals, beta_shape12, beta_shape22);
					BetaDist bd22 = new BetaDist(betamleparams2[0], betamleparams2[1]);
					double betaApproxPval2 = bd22.cdf(bestPval);

					String outln = gene + "\t" + genotypes.size()
							+ "\t" + bestPvalSNP + "\t" + bestPval + "\t" + bestPvalN + "\t" + bestPvalR
							+ "\t" + bd22.getAlpha() + "\t" + bd22.getBeta() + "\t" + proportionBetterPvals + "\t" + betaApproxPval2;

					output.writeln(outln);
					// output somewhere
					// // "Gene\tSNP\tAlleles\tEffectallele\tMetaN\tMetaZ\tMetaP\tBDistAlpha\tBDistBeta\tAdjP"
					System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + tctr + " variants tested. T: " + timer.getTimeDesc() + "\n");

					output.flush();
//                    System.exit(0);
				} else {
					System.err.println("Error no tested SNPs for gene " + gene);
				}


			}
		}
		outputAll.close();
		output.close();
		tabix.close();
	}

	public void setNrPermutations(int nrPermutations) {
		this.nrPermutations = nrPermutations;
	}

	public void setCisWindow(int cisWindow) {
		this.cisWindow = cisWindow;
	}

	public void setRandomSeed(long randomSeed) {
		this.randomSeed = randomSeed;
	}


}
