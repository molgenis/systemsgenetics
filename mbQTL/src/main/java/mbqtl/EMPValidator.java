package mbqtl;

import mbqtl.datastructures.Dataset;
import mbqtl.vcf.VCFTabix;
import mbqtl.vcf.VCFVariant;
import umcg.genetica.containers.Triple;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.RankArray;
import umcg.genetica.util.RunTimer;

import java.io.IOException;
import java.util.*;

public class EMPValidator extends QTLAnalysis {


	public EMPValidator(String vcfFile, int chromosome, String linkfile, String geneLimitFile, String geneExpressionDataFile, String geneAnnotationFile, String outfile) throws IOException {
		super(vcfFile, chromosome, linkfile, null, geneLimitFile, null, geneExpressionDataFile, geneAnnotationFile, outfile);
	}

	public void validate(String vcfFile,
						 int chromosome,
						 String linkfile,
						 String geneLimitFile,
						 String geneExpressionDataFile,
						 String geneAnnotationFile,
						 String outfile) throws IOException {



		// iterate genes
		Chromosome chromosomeObj = Chromosome.parseChr("" + chromosome);
		System.out.println("Processing: " + vcfFile);
		VCFTabix tabix = new VCFTabix(vcfFile);

		TextFile output = new TextFile(outfile, TextFile.W);
		TextFile outputsnplog = new TextFile(outfile + "SNPLog.txt.gz", TextFile.W);
		output.writeln("Gene\tSNP\tAlleles\tEffectallele\t" +
				"MetaN\tMetaZ\tMetaP\t" +
				"MetaZPreRanked\tMetaPMPreRanked\t" +
				"JointN1\tJointZ1\tJointP1\t" +
				"JointN2\tJointZ2\tJointP2");

		String snplogheader = "SNP\tDatasetsPassQC\tNTotal\tNJoint1\tMAFJoint1\tCallRateJoint1\tHWEPJoint1\tNJoint1\tMAFJoint1\tCallRateJoint1\tHWEPJoint1";
		HashSet<String> seenSNP = new HashSet<>();
		for (int d = 0; d < datasets.length; d++) {
			String name = datasets[d].getName();
			snplogheader += "\t" + name + "-MAF";
			snplogheader += "\t" + name + "-CallRate";
			snplogheader += "\t" + name + "-HWEP";
		}
		outputsnplog.writeln(snplogheader);

		// iterate genes
		RunTimer timer = new RunTimer();
		timer.start();
		RankArray ranker = new RankArray();
		for (int g = 0; g < expressionData.genes.length; g++) {
			String gene = expressionData.genes[g];

			Gene geneAnnotationObj = geneAnnotation.getGene(gene);
			if (geneAnnotationObj != null) {
				double[] expData = expressionData.data[g];

				// rank once, then center+scale
				double[][] expDataPerDatasetRanked = new double[datasets.length][];
				double[][] expDataPerDatasetRankedCenterScaled = new double[datasets.length][];

				for (int d = 0; d < datasets.length; d++) {
					Dataset thisDataset = datasets[d];
					double[] datasetExpressionData = thisDataset.select(expData, thisDataset.getExpressionIds());
					expDataPerDatasetRanked[d] = ranker.rank(datasetExpressionData, true);
					expDataPerDatasetRankedCenterScaled[d] = Util.centerScale(expDataPerDatasetRanked[d]);
				}

				// get variants, 1mb up and downstream
				int pos = geneAnnotationObj.getStart();
				int start = pos - 1000001;
				if (start < 0) {
					start = 0;
				}
				int stop = pos + 1000001;
				Feature region = new Feature(chromosomeObj, start, stop);
				Iterator<VCFVariant> snpIterator = tabix.getVariants(region, genotypeSamplesToInclude);
				int vctr = 0;
				int tctr = 0;
				boolean print = false;
				while (snpIterator.hasNext()) {
					VCFVariant variant = snpIterator.next();

					if (variant != null) {
						String snpid = variant.getId();
						if (snpid.equals("22:18908190:rs424512:T_C")) {
							System.out.println();
							System.out.println("-------");
							print = true;
						} else {
							print = false;
						}
						// collect data per dataset
						double[][] datasetsGenotypeData = new double[datasets.length][];
						double[][] datasetsGenotypeDosageData = new double[datasets.length][];
						VariantQCObj[] datasetsQCObjs = new VariantQCObj[datasets.length];
						int datasetsPassQC = 0;

//                        ArrayList<Double> genotypeDataForJointAnalysis = new ArrayList<>();
//                        ArrayList<Double> genotypeDosageDataForJointAnalysis = new ArrayList<>();
//                        ArrayList<Double> expressionDataForJointAnalysis = new ArrayList<>();
//
//                        ArrayList<Double> genotypeDataForJointAnalysisDatasetsPassingQC = new ArrayList<>();
//                        ArrayList<Double> genotypeDosageDataForJointAnalysisDatasetsPassingQC = new ArrayList<>();
//                        ArrayList<Double> expressionDataForJointAnalysisDatasetsPassingQC = new ArrayList<>();

						// System.out.println(variant.getId() + "\t" + gene + "\t" + Strings.concat(variant.getAlleles(), Strings.backwardslash));
						for (int d = 0; d < datasets.length; d++) {
							Dataset thisDataset = datasets[d];

							double[] datasetGenotypeData = thisDataset.select(variant.getGenotypesAsByteVector(), thisDataset.getGenotypeIds());

							// do some QC checks? test MAF, Callrate, HWE-P
							VariantQCObj qcobj = checkVariant(datasetGenotypeData);
							datasetsQCObjs[d] = qcobj;
							datasetsGenotypeData[d] = datasetGenotypeData;
//                            double[] datasetDosageValues = thisDataset.select(variant.getGenotypeDosage(), thisDataset.genotypeIds);
							double[] datasetDosageValues = thisDataset.select(variant.getDosage(), thisDataset.getGenotypeIds());
							datasetsGenotypeDosageData[d] = datasetDosageValues;

//                            System.out.println(variant.getGenotypesAsByteVector().length);
							if (print) {
								System.out.println(thisDataset.getName() + "\t" + datasetGenotypeData.length + "\t" + qcobj.toString());
							}
//                            for (int q = 0; q < datasetGenotypeData.length; q++) {
//                                System.out.println(datasetGenotypeData[q]);
//                            }
//                            System.exit(-1);

							// System.out.println(datasets[d].name + "\t" + qcobj.toString());
							for (int v = 0; v < datasetGenotypeData.length; v++) {
//                                genotypeDataForJointAnalysis.add(datasetGenotypeData[v]);
//                                genotypeDosageDataForJointAnalysis.add(datasetDosageValues[v]);
//                                expressionDataForJointAnalysis.add(expDataPerDatasetRankedCenterScaled[d][v]);
								if (qcobj.passqc) {
//                                    genotypeDataForJointAnalysisDatasetsPassingQC.add(datasetGenotypeData[v]);
//                                    genotypeDosageDataForJointAnalysisDatasetsPassingQC.add(datasetDosageValues[v]);
//                                    expressionDataForJointAnalysisDatasetsPassingQC.add(expDataPerDatasetRankedCenterScaled[d][v]);
								}
							}
							if (qcobj.passqc) {
								datasetsPassQC++;
							}
						}

						if (datasetsPassQC > 1) {
							int totalN = 0;
							// perform different types of analysis:

							////////////////////////////////
							// 1. EMP style meta-analysis //
							////////////////////////////////
							double[] zScores = new double[datasets.length];
							int[] datasetSampleSizes = new int[datasets.length];

							int dctr = 0;
							for (int d = 0; d < datasets.length; d++) {
								if (datasetsQCObjs[d].passqc) {
									// prune the dataset; remove missing values

									Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(datasetsGenotypeData[d], datasetsGenotypeDosageData[d], expDataPerDatasetRanked[d]);
									if (prunedDatasetData.getRight().length == 0 || prunedDatasetData.getLeft().length == 0) {
										System.out.println(datasets[d].getName() + " has 0 values, but " + datasets[d].getExpressionIds().length + " samples");
										for (int v = 0; v < datasetsGenotypeData[d].length; v++) {
											System.out.println(datasetsGenotypeData[d][v] + "\t" + expDataPerDatasetRanked[d][v]);
										}
										System.out.println(datasetsQCObjs[d].toString());
										System.exit(-1);
									}
									double[] expDataCenterScale = Util.centerScale(prunedDatasetData.getRight());
									double[] gtDataCenterScale = Util.centerScale(prunedDatasetData.getMiddle());


									// perform correlation
									double r = Correlation.correlate(expDataCenterScale, gtDataCenterScale);
									// convert correlation to z-score
									int n = gtDataCenterScale.length;
									datasetSampleSizes[d] = n;
									totalN += n;
									double z = Correlation.convertCorrelationToZScore(n, r);
									zScores[d] = z;
//                                    System.out.println(datasets[d].name + "\t" + n + "\t" + r + "\t" + z + "\t" + datasetsQCObjs[d]);
									dctr++;
								} else {
									zScores[d] = Double.NaN;
//                                    System.out.println(datasets[d].name + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + datasetsQCObjs[d]);
								}
								if (print) {
									System.out.println(datasets[d].getName() + "\t" + zScores[d]);
								}
							}


							// 1.2 meta-analyze z-scores
							double metaZ = ZScores.getWeightedZ(zScores, datasetSampleSizes);
							double metaP = ZScores.zToP(metaZ);
							if (print) {
								System.out.println(metaZ + "\t" + metaP);
								System.out.println("-------");
//                                System.exit(0);
							}
//                            System.out.println(metaZ);
//                            System.out.println(metaP);
//                            System.exit(0);

							String[] alleles = variant.getAlleles();
							String alleleStr = alleles[0] + "/" + alleles[1];
							tctr++;

							////////////////////////////////////////////////////////////////
							// 2. joint analysis including all datasets, regardless of QC //
							////////////////////////////////////////////////////////////////

							// prune the data
//                            Pair<double[], double[]> prunedData = pruneMissing(Primitives.toPrimitiveArr(genotypeDataForJointAnalysis),
//                                    Primitives.toPrimitiveArr(genotypeDosageDataForJointAnalysis),
//                                    Primitives.toPrimitiveArr(expressionDataForJointAnalysis));
//                            // perform correlation
//                            QCObj combinedQCObj1 = checkVariant(Primitives.toPrimitiveArr(genotypeDataForJointAnalysis));
							double zJoint1 = Double.NaN;
							double pJoint1 = Double.NaN;
							int nJoint1 = 0;
//                            //if (combinedQCObj1.passqc) {
//                            double[] genotypeDataCenterScale = centerScale(prunedData.getLeft());
//                            double rJoint1 = Correlation.correlate(prunedData.getRight(), genotypeDataCenterScale);
//                            nJoint1 = genotypeDataCenterScale.length;
//                            zJoint1 = Correlation.convertCorrelationToZScore(nJoint1, rJoint1);
//                            pJoint1 = ZScores.zToP(zJoint1);
							//}

							//////////////////////////////////////////////////////////////
							// 3. joint analysis including only datasets passing SNP QC //
							//////////////////////////////////////////////////////////////

							// prune the data
//                            Pair<double[], double[]> prunedData2 = pruneMissing(Primitives.toPrimitiveArr(genotypeDataForJointAnalysisDatasetsPassingQC),
//                                    Primitives.toPrimitiveArr(genotypeDosageDataForJointAnalysisDatasetsPassingQC),
//                                    Primitives.toPrimitiveArr(expressionDataForJointAnalysisDatasetsPassingQC));
//                            // perform correlation
//                            QCObj combinedQCObj2 = checkVariant(Primitives.toPrimitiveArr(genotypeDataForJointAnalysisDatasetsPassingQC));
							double zJoint2 = Double.NaN;
							double pJoint2 = Double.NaN;
							int nJoint2 = 0;
//
//                            genotypeDataCenterScale = centerScale(prunedData2.getLeft());
//                            double rJoint2 = Correlation.correlate(prunedData2.getRight(), genotypeDataCenterScale);
//                            nJoint2 = genotypeDataCenterScale.length;
//                            zJoint2 = Correlation.convertCorrelationToZScore(nJoint2, rJoint2);
//                            pJoint2 = ZScores.zToP(zJoint2);


							// "Gene\tSNP\tAlleles\tEffectallele\tMetaN\tMetaZ\tMetaP\tRankedOverAllSamplesZ\tRankedOverAllSamplesP\tRankedPerDatasetZ\tRankedPerDatasetP"
							output.writeln(gene + "\t" + snpid + "\t" + alleleStr + "\t" + alleles[1] + "\t" +
									totalN + "\t" + metaZ + "\t" + metaP + "\t" +
									nJoint1 + "\t" + zJoint1 + "\t" + pJoint1 + "\t" +
									nJoint2 + "\t" + zJoint2 + "\t" + pJoint2);

							StringBuilder snplogstr = null;
							if (!seenSNP.contains(snpid)) {
								// SNP	DatasetsPassQC	NTotal	NJoint MAFJoint	CallRateJoint	HWEPJoint
//                                snplogstr = new StringBuilder().append(snpid)
//                                        .append("\t").append(datasetsPassQC)
//                                        .append("\t").append(totalN)
//                                        .append("\t").append(nJoint1)
//                                        .append("\t").append(combinedQCObj1.maf)
//                                        .append("\t").append(combinedQCObj1.cr)
//                                        .append("\t").append(combinedQCObj1.hwep)
//                                        .append("\t").append(nJoint2)
//                                        .append("\t").append(combinedQCObj2.maf)
//                                        .append("\t").append(combinedQCObj2.cr)
//                                        .append("\t").append(combinedQCObj2.hwep);
								snplogstr = new StringBuilder().append(snpid)
										.append("\t").append(datasetsPassQC)
										.append("\t").append(totalN)
										.append("\t").append(nJoint1)
										.append("\t").append(0)
										.append("\t").append(0)
										.append("\t").append(1)
										.append("\t").append(nJoint2)
										.append("\t").append(0)
										.append("\t").append(0)
										.append("\t").append(1);
								for (int d = 0; d < datasetsQCObjs.length; d++) {
									snplogstr.append("\t").append(datasetsQCObjs[d].maf)
											.append("\t").append(datasetsQCObjs[d].cr)
											.append("\t").append(datasetsQCObjs[d].hwep);
								}
								outputsnplog.writeln(snplogstr.toString());
								seenSNP.add(snpid);
							}
						}


					}
					vctr++;
					if (vctr % 1000 == 0) {
						System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + tctr + " variants tested. T: " + timer.getTimeDesc() + "\r");
					}
				}
				System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + tctr + " variants tested. T: " + timer.getTimeDesc() + "\n");
			}
		}
		outputsnplog.close();
		output.close();
		tabix.close();

	}


}
