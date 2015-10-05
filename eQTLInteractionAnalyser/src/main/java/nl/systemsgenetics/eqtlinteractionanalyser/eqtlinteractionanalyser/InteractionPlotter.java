package nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import static nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser.TestEQTLDatasetForInteractions.getEqtls;
import static nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser.TestEQTLDatasetForInteractions.getLinearRegressionCoefficients;
import static nl.systemsgenetics.eqtlinteractionanalyser.eqtlinteractionanalyser.TestEQTLDatasetForInteractions.orthogonalizeDataset;
import org.apache.commons.math3.stat.ranking.NaturalRanking;

/**
 *
 * @author Patrick Deelen
 */
public class InteractionPlotter {

	static String inputDir = null;
	static String outputDir = null;

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		//makeInteractionPlot("D:\\tmp\\test.png", new double[]{0,0,0,0.2,1,1,1,1,2,2,2}, new double[]{5,4,3,0.2,8,12,6,7,23,5,7}, new double[]{3,2,1,0.2,2,6,4,6,20,2,5});

		inputDir = args[0];
		outputDir = args[1];
		String eQTLfileName = args[2];
		String covariate = args[3];
		File genesFile = new File(args[4]);



		System.out.println("Input dir: " + inputDir);
		System.out.println("Output dir: " + outputDir);
		System.out.println("eQTL file: " + eQTLfileName);
		System.out.println("covariate: " + covariate);
		System.out.println("genes file: " + genesFile.getAbsolutePath());

		String[] covsToCorrect = {"gender", "GC", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "CEU", "GBR", "FIN", "TSI", "YRI"};
		//String[] covsToCorrect = {"age", "gender", "GC", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "LLdeep", "LLS", "RS", "CODAM"};
		HashMap hashEQTLs = getEqtls(eQTLfileName);

		HashMap hashSamples = new HashMap();

		if (1 == 1) {

			System.out.println("Removing outlier samples!!!");
			HashMap hashCovariates = new HashMap();
			hashCovariates.put("MEDIAN_5PRIME_BIAS", null);
			hashCovariates.put("MEDIAN_3PRIME_BIAS", null);
			ExpressionDataset datasetCovariates = new ExpressionDataset(inputDir + "/covariateTableLude.txt.Covariates.binary", '\t', hashCovariates, null);
			hashSamples = new HashMap();
			for (int s = 0; s < datasetCovariates.nrSamples; s++) {
				if (datasetCovariates.rawData[0][s] != 0) {
					hashSamples.put(datasetCovariates.sampleNames[s], null);
				}
			}
			datasetCovariates = new ExpressionDataset(inputDir + "/covariateTableLude.txt.Covariates.binary", '\t', hashCovariates, hashSamples);
			HashMap hashSamplesToExclude = new HashMap();
			if (1 == 1) {
				int index = ((Integer) datasetCovariates.hashProbes.get("MEDIAN_5PRIME_BIAS")).intValue();
				double mean = JSci.maths.ArrayMath.mean(datasetCovariates.rawData[index]);
				double stdev = JSci.maths.ArrayMath.standardDeviation(datasetCovariates.rawData[index]);
				for (int s = 0; s < datasetCovariates.nrSamples; s++) {
					double z = (datasetCovariates.rawData[index][s] - mean) / stdev;
					if (Math.abs(z) > 3) {
						hashSamplesToExclude.put(datasetCovariates.sampleNames[s], null);
					}
				}
			}
			if (1 == 1) {
				int index = ((Integer) datasetCovariates.hashProbes.get("MEDIAN_3PRIME_BIAS")).intValue();
				double mean = JSci.maths.ArrayMath.mean(datasetCovariates.rawData[index]);
				double stdev = JSci.maths.ArrayMath.standardDeviation(datasetCovariates.rawData[index]);
				for (int s = 0; s < datasetCovariates.nrSamples; s++) {
					double z = (datasetCovariates.rawData[index][s] - mean) / stdev;
					if (Math.abs(z) > 3) {
						hashSamplesToExclude.put(datasetCovariates.sampleNames[s], null);
					}
				}
			}
			hashSamples = new HashMap();
			for (int s = 0; s < datasetCovariates.nrSamples; s++) {
				if (!hashSamplesToExclude.containsKey(datasetCovariates.sampleNames[s])) {
					hashSamples.put(datasetCovariates.sampleNames[s], null);
					hashSamples.put(datasetCovariates.sampleNames[s] + "_exp", null);
					hashSamples.put(datasetCovariates.sampleNames[s] + "_dosage", null);
				}
			}
		}

		ExpressionDataset datasetGenotypes = new ExpressionDataset(inputDir + "/bigTableLude.txt.Genotypes.binary", '\t', hashEQTLs, hashSamples);
		ExpressionDataset datasetExpression = new ExpressionDataset(inputDir + "/bigTableLude.txt.Expression.binary", '\t', hashEQTLs, hashSamples);
		ExpressionDataset datasetCovariates = new ExpressionDataset(inputDir + "/covariateTableLude.txt.Covariates.binary", '\t', null, hashSamples);

		org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression regression = new org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression();
		int nrSamples = datasetGenotypes.nrSamples;


		if (1 == 1) {
			//Define a set of covariates that we want to use as correction:
			System.out.println("Correcting gene expression data for cohort specific effects and top 25 components");
			//String[] cohorts = {"LLDeep", "LLS", "RS", "CODAM"};
			int nrCompsToCorrectFor = 25;
			ExpressionDataset datasetCovariatesToCorrectFor = new ExpressionDataset(nrCompsToCorrectFor, datasetGenotypes.nrSamples);
			datasetCovariatesToCorrectFor.sampleNames = datasetGenotypes.sampleNames;
//			for (int p = 0; p < cohorts.length; p++) {
//				for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
//					if (datasetGenotypes.sampleNames[s].startsWith(cohorts[p])) {
//						datasetCovariatesToCorrectFor.rawData[p][s] = 1;
//					}
//				}
//			}
			if (nrCompsToCorrectFor > 0) {
				for (int comp = 0; comp < nrCompsToCorrectFor; comp++) {
					for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
						datasetCovariatesToCorrectFor.rawData[comp][s] = datasetCovariates.rawData[datasetCovariates.nrProbes - 51 + comp][s];
					}
				}
			}

			datasetCovariatesToCorrectFor.transposeDataset();

			datasetCovariatesToCorrectFor.save(inputDir + "/CovariatesToCorrectFor.txt");
			orthogonalizeDataset(inputDir + "/CovariatesToCorrectFor.txt");
			datasetCovariatesToCorrectFor = new ExpressionDataset(inputDir + "/CovariatesToCorrectFor.txt.PrincipalComponents.txt");
			datasetCovariatesToCorrectFor.transposeDataset();
			ExpressionDataset datasetCovariatesToCorrectForEigenvalues = new ExpressionDataset(inputDir + "/CovariatesToCorrectFor.txt.Eigenvalues.txt");
			for (int snp = 0; snp < datasetExpression.nrProbes; snp++) {
				for (int cov = 0; cov < datasetCovariatesToCorrectFor.nrProbes; cov++) {
					if (datasetCovariatesToCorrectForEigenvalues.rawData[cov][0] > 1E-5) {
						double[] rc = getLinearRegressionCoefficients(datasetCovariatesToCorrectFor.rawData[cov], datasetExpression.rawData[snp]);
						for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
							datasetExpression.rawData[snp][s] -= rc[0] * datasetCovariatesToCorrectFor.rawData[cov][s];
						}
					}
				}
			}


		}


		double[] mainEQTLCorr = new double[datasetGenotypes.nrProbes];
		if (1 == 1) {
			System.out.println("Enforcing for every eQTL that the genotype dosage positively correlated with gene expression levels:");
			for (int snp = 0; snp < datasetGenotypes.nrProbes; snp++) {
				double corr = JSci.maths.ArrayMath.correlation(datasetGenotypes.rawData[snp], datasetExpression.rawData[snp]);
				//System.out.println(datasetExpression.probeNames[snp] + "\t" + snp + "\t" + corr);

				if (corr < 0) {
					corr = -corr;
					for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
						datasetGenotypes.rawData[snp][s] = 2 - datasetGenotypes.rawData[snp][s];
					}
				}

				mainEQTLCorr[snp] = corr;
			}
		}

		if (1 == 1) {

			if (1 == 1) {
				System.out.println("Correcting covariate data for cohort specific effects:");
//                String[] cohorts = {"LLDeep","LLS","RS","CODAM"};
				ExpressionDataset datasetCovariatesToCorrectFor = new ExpressionDataset(covsToCorrect.length, datasetGenotypes.nrSamples);
				datasetCovariatesToCorrectFor.sampleNames = datasetGenotypes.sampleNames;
//                for (int p=0; p<cohorts.length; p++) {
//                    for (int s=0; s<datasetGenotypes.nrSamples; s++) {
//                        if (datasetGenotypes.sampleNames[s].startsWith(cohorts[p])) {
//                            datasetCovariatesToCorrectFor.rawData[p][s]=1;
//                        }
//                    }
//                }
				HashMap hashCovsToCorrect = new HashMap();
				int[] covsToCorrectIndex = new int[covsToCorrect.length];
				for (int c = 0; c < covsToCorrect.length; c++) {
					hashCovsToCorrect.put(covsToCorrect[c], null);
					covsToCorrectIndex[c] = ((Integer) datasetCovariates.hashProbes.get(covsToCorrect[c])).intValue();
					for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
						datasetCovariatesToCorrectFor.rawData[c][s] = datasetCovariates.rawData[covsToCorrectIndex[c]][s];
					}
				}

				datasetCovariatesToCorrectFor.transposeDataset();

				datasetCovariatesToCorrectFor.save(inputDir + "/CovariatesToCorrectFor.txt");
				orthogonalizeDataset(inputDir + "/CovariatesToCorrectFor.txt");
				datasetCovariatesToCorrectFor = new ExpressionDataset(inputDir + "/CovariatesToCorrectFor.txt.PrincipalComponents.txt");
				datasetCovariatesToCorrectFor.transposeDataset();
				ExpressionDataset datasetCovariatesToCorrectForEigenvalues = new ExpressionDataset(inputDir + "/CovariatesToCorrectFor.txt.Eigenvalues.txt");

				for (int p = 0; p < datasetCovariates.nrProbes; p++) {
					if (!hashCovsToCorrect.containsKey(datasetCovariates.probeNames[p])) {
						for (int cov = 0; cov < datasetCovariatesToCorrectFor.nrProbes; cov++) {
							if (datasetCovariatesToCorrectForEigenvalues.rawData[cov][0] > 1E-5) {
								double[] rc = getLinearRegressionCoefficients(datasetCovariatesToCorrectFor.rawData[cov], datasetCovariates.rawData[p]);
								for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
									datasetCovariates.rawData[p][s] -= rc[0] * datasetCovariatesToCorrectFor.rawData[cov][s];
								}
							}
						}
						double stdev = JSci.maths.ArrayMath.standardDeviation(datasetCovariates.rawData[p]);
						double mean = JSci.maths.ArrayMath.mean(datasetCovariates.rawData[p]);
						if (stdev < 1E-5) {
							for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
								datasetCovariates.rawData[p][s] = mean;
							}
						}
					}
				}


			}

			if (1 == 1) {
				System.out.println("Correcting covariate data for cis-eQTL effects:");
				for (int p = 0; p < datasetCovariates.nrProbes; p++) {
					if (datasetExpression.hashProbes.containsKey(datasetCovariates.probeNames[p])) {
						int index = ((Integer) datasetExpression.hashProbes.get(datasetCovariates.probeNames[p])).intValue();
						double[] rc = getLinearRegressionCoefficients(datasetGenotypes.rawData[index], datasetCovariates.rawData[p]);
						for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
							datasetCovariates.rawData[p][s] -= rc[0] * datasetGenotypes.rawData[index][s];
						}
					}
				}
			}

			if (1 == 2) {
				datasetCovariates.save(inputDir + "/CovariatesCorrected.txt");
				HashMap hashProbesToFilter = new HashMap();
				for (int p = 0; p < datasetCovariates.nrProbes; p++) {
					if (datasetCovariates.probeNames[p].startsWith("ENSG")) {
						hashProbesToFilter.put(datasetCovariates.probeNames[p], null);
					}
				}
				ExpressionDataset datasetCovariatesCorrected = new ExpressionDataset(inputDir + "/CovariatesCorrected.txt", '\t', hashProbesToFilter, null);
				datasetCovariatesCorrected.transposeDataset();
				datasetCovariatesCorrected.save(inputDir + "/CovariatesCorrected.txt");
				System.exit(0);
			}

			if (1 == 2) {
				ExpressionDataset datasetICA = new ExpressionDataset("/Users/lude/Documents/ICA/mixingmatrix.txt");
				//ExpressionDataset datasetICA = new ExpressionDataset("/Users/lude/Documents/ICA/signals.txt");
				datasetICA.transposeDataset();
				for (int p = 0; p < datasetICA.nrProbes; p++) {
					datasetCovariates.rawData[p] = datasetICA.rawData[p];
					datasetCovariates.probeNames[p] = datasetICA.probeNames[p];
					if (p == 7) {
						for (int q = 0; q < datasetCovariates.nrProbes; q++) {
							double corr = JSci.maths.ArrayMath.correlation(datasetICA.rawData[p], datasetCovariates.rawData[q]);
							System.out.println(p + "\t" + datasetICA.probeNames[p] + "\t" + q + "\t" + datasetCovariates.probeNames[q] + "\t" + corr + "\t" + corr * corr);
						}
					}
				}

				orthogonalizeDataset("/Users/lude/Documents/ICA/mixingmatrix.txt");
				//System.exit(0);
			}

			System.out.println("Enforcing normal distribution on covariates");

			NaturalRanking ranker = new NaturalRanking();

			for (int p = 0; p < datasetCovariates.nrProbes; p++) {
				//Rank order the expression values:
				double[] values = new double[datasetCovariates.nrSamples];
				for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
					values[s] = datasetCovariates.rawData[p][s];
				}
				double[] rankedValues = ranker.rank(values);
				//Replace the original expression value with the standard distribution enforce:
				for (int s = 0; s < datasetGenotypes.nrSamples; s++) {
					//Convert the rank to a proportion, with range <0, 1>
					double pValue = (0.5d + rankedValues[s] - 1d) / (double) (rankedValues.length);
					//Convert the pValue to a Z-Score:
					double zScore = cern.jet.stat.tdouble.Probability.normalInverse(pValue);
					datasetCovariates.rawData[p][s] = zScore; //Replace original expression value with the Z-Score
				}
			}

		}

		cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();

		ExpressionDataset datasetExpressionBeforeEQTLCorrection = new ExpressionDataset(datasetExpression.nrProbes, datasetExpression.nrSamples);
		for (int p = 0; p < datasetExpression.nrProbes; p++) {
			for (int s = 0; s < datasetExpression.nrSamples; s++) {
				datasetExpressionBeforeEQTLCorrection.rawData[p][s] = datasetExpression.rawData[p][s];
			}
		}

		if (1 == 1) {
			System.out.println("Correcting expression data for predefined gene environment interaction effects (GC content, Gender, 5'Median Bias, 3'Median Bias):");
			int[] covsToCorrectIndex = new int[covsToCorrect.length];
			for (int c = 0; c < covsToCorrect.length; c++) {
				covsToCorrectIndex[c] = ((Integer) datasetCovariates.hashProbes.get(covsToCorrect[c])).intValue();
			}
			for (int snp = 0; snp < datasetGenotypes.nrProbes; snp++) {
				double[][] valsX = new double[nrSamples][1 + covsToCorrect.length * 2]; //store genotypes, covariates, interactions
				for (int s = 0; s < nrSamples; s++) {
					valsX[s][0] = datasetGenotypes.rawData[snp][s]; //genotypes
				}
				for (int c = 0; c < covsToCorrect.length; c++) {
					for (int s = 0; s < nrSamples; s++) {
						valsX[s][c * 2 + 1] = datasetCovariates.rawData[covsToCorrectIndex[c]][s]; //covariate
						valsX[s][c * 2 + 2] = valsX[s][0] * valsX[s][c * 2 + 1]; //interction
					}
				}
				double[] valsY = datasetExpression.rawData[snp];
				regression.newSampleData(valsY, valsX);
				datasetExpression.rawData[snp] = regression.estimateResiduals();
			}
		}


		if (1 == 1) {
			System.out.println("Enforcing normal distribution on expression data:");

			NaturalRanking ranker = new NaturalRanking();

			for (int p = 0; p < datasetExpression.nrProbes; p++) {
				//Rank order the expression values:
				double[] values = new double[datasetExpression.nrSamples];
				for (int s = 0; s < datasetExpression.nrSamples; s++) {
					values[s] = datasetExpression.rawData[p][s];
				}

				double[] rankedValues = ranker.rank(values);
				//Replace the original expression value with the standard distribution enforce:
				for (int s = 0; s < datasetExpression.nrSamples; s++) {
					//Convert the rank to a proportion, with range <0, 1>
					double pValue = (0.5d + rankedValues[s] - 1d) / (double) (rankedValues.length);
					//Convert the pValue to a Z-Score:
					double zScore = cern.jet.stat.tdouble.Probability.normalInverse(pValue);
					datasetExpression.rawData[p][s] = zScore; //Replace original expression value with the Z-Score
				}
			}

			System.out.println("Expression data now force normal");

		}


		CSVReader reader = new CSVReader(new FileReader(genesFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);
		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			String eQtlGene = nextLine[0];
			
			System.out.println(eQtlGene);

			Integer eQtlGeneI = datasetExpression.hashProbes.get(eQtlGene);
			Integer covariateI = datasetCovariates.hashProbes.get(covariate);
			Integer snpI = eQtlGeneI;

			makeInteractionPlot(outputDir + "/" + covariate + "-" + eQtlGene + ".png", datasetGenotypes.rawData[snpI], datasetExpression.rawData[eQtlGeneI], datasetCovariates.rawData[covariateI]);

		}

	}

	public static void makeInteractionPlot(String fileName, double[] genotype, double[] expression, double[] covariate) {

		int nrSamples = genotype.length;

//		int[] cohortIndex = new int[4];
//		String[] cohorts = {"LLDeep", "LLS", "RS", "CODAM"};
//		for (int cohort = 0; cohort < cohorts.length; cohort++) {
//			for (int s = 0; s < nrSamples; s++) {
//				if (sampleNames[s].startsWith(cohorts[cohort])) {
//					cohortIndex[cohort] = s;
//					break;
//				}
//			}
//		}

		int marginLeft = 100;
		int marginRight = 200;
		int marginTop = 100;
		int marginBottom = 100;
		int innerHeight = 500;
		int innerWidth = 500;
		int docWidth = marginLeft + marginRight + innerWidth;
		int docHeight = marginTop + marginBottom + innerHeight;

		BufferedImage bimage = new BufferedImage(docWidth, docHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D g2d = bimage.createGraphics();

		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setColor(new Color(255, 255, 255));
		g2d.fillRect(0, 0, docWidth, docHeight);
		java.awt.AlphaComposite alphaComposite10 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.10f);
		java.awt.AlphaComposite alphaComposite25 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f);
		java.awt.AlphaComposite alphaComposite50 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f);
		java.awt.AlphaComposite alphaComposite100 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC, 1.00f);

		float fontSize = 12f;
		java.awt.Font font = new java.awt.Font("Gill Sans MT", java.awt.Font.PLAIN, (int) fontSize);
		java.awt.Font fontBold = new java.awt.Font("Gill Sans MT", java.awt.Font.BOLD, (int) fontSize);
		java.awt.Font fontSmall = new java.awt.Font("Gill Sans MT", java.awt.Font.PLAIN, 8);
		java.awt.Font fontBoldSmall = new java.awt.Font("Gill Sans MT", java.awt.Font.BOLD, 8);

		java.awt.Color dataColor[] = new Color[10];
		dataColor[0] = new java.awt.Color(167, 72, 20);
		dataColor[1] = new java.awt.Color(62, 138, 20);
		dataColor[2] = new java.awt.Color(228, 171, 0);
		dataColor[3] = new java.awt.Color(0, 148, 183);
		dataColor[4] = new java.awt.Color(119, 80, 152);
		dataColor[5] = new java.awt.Color(106, 106, 106);
		dataColor[6] = new java.awt.Color(212, 215, 10);
		dataColor[7] = new java.awt.Color(210, 111, 0);
		dataColor[8] = new java.awt.Color(0, 0, 141);
		dataColor[9] = new java.awt.Color(190, 190, 190);

		g2d.setComposite(alphaComposite50);
		g2d.setColor(new Color(0, 0, 0));
		g2d.drawLine(marginLeft, marginTop, marginLeft, marginTop + innerHeight);
		g2d.drawLine(marginLeft, marginTop + innerHeight, marginLeft + innerWidth, marginTop + innerHeight);

		double minX = JSci.maths.ArrayMath.min(covariate);
		double maxX = JSci.maths.ArrayMath.max(covariate);
		double minY = JSci.maths.ArrayMath.min(expression);
		double maxY = JSci.maths.ArrayMath.max(expression);

		g2d.setComposite(alphaComposite10);
		for (int rep = 0; rep >= 0; rep--) {
			for (int s = 0; s < nrSamples; s++) {
				int posY = marginTop + innerHeight - (int) ((expression[s] - minY) / (maxY - minY) * innerHeight);
				int posX = marginLeft + (int) ((covariate[s] - minX) / (maxX - minX) * innerWidth);
				if (genotype[s] < 0.5) {
					g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_ATOP, 0.30f - (float) rep / 10f));
					g2d.setColor(new Color(204, 86, 78));
				} else {
					if (genotype[s] > 1.5) {
						g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_ATOP, 0.30f - (float) rep / 10f));
						g2d.setColor(new Color(171, 178, 114));
					} else {
						g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_ATOP, 0.30f - (float) rep / 10f));
						g2d.setColor(new Color(98, 175, 255));
					}
				}

				g2d.fillOval(posX - 5 - rep * 4, posY - 5 - rep * 4, 7 + rep * 8, 7 + rep * 8);

			}
		}

		//Draw the four independent cohorts seperately:
		//int[] cohortIndex = {0,626,1280,1933};
//		for (int rep = 2; rep >= 0; rep--) {
//			for (int s = 0; s < nrSamples; s++) {
//				int cohort = 0;
//				for (int c = 0; c < cohortIndex.length; c++) {
//					if (s >= cohortIndex[c]) {
//						cohort = c;
//					}
//				}
//
//				int posY = marginTop + 100 + cohort * 125 - (int) ((expression[s] - minY) / (maxY - minY) * 100);
//				int posX = marginLeft + innerWidth + 50 + (int) ((covariate[s] - minX) / (maxX - minX) * 100);
//				if (genotype[s] < 0.5) {
//					g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_ATOP, 0.30f - (float) rep / 10f));
//					g2d.setColor(new Color(204, 86, 78));
//				} else {
//					if (genotype[s] > 1.5) {
//						g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_ATOP, 0.30f - (float) rep / 10f));
//						g2d.setColor(new Color(171, 178, 114));
//					} else {
//						g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_ATOP, 0.30f - (float) rep / 10f));
//						g2d.setColor(new Color(98, 175, 255));
//					}
//				}
//				g2d.fillOval(posX - 1 - rep * 2, posY - 1 - rep * 2, 3 + rep * 4, 3 + rep * 4);
//
//			}
//		}


		g2d.setComposite(alphaComposite50);
		double[][] valsX = new double[nrSamples][3];
		for (int s = 0; s < nrSamples; s++) {
			valsX[s][0] = genotype[s];
			valsX[s][1] = covariate[s];
			valsX[s][2] = valsX[s][0] * valsX[s][1];
		}
		double[] valsY = expression;
		org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression regression = new org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression();
		regression.newSampleData(valsY, valsX);
		double[] betas = regression.estimateRegressionParameters();
		double betaInteraction = betas[3];
		double seInteraction = regression.estimateRegressionParametersStandardErrors()[3];
		double tInteraction = betaInteraction / seInteraction;
		double pValueInteraction = 1;
		double zScoreInteraction = 0;
		cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
		cern.jet.random.tdouble.StudentT tDistColt = new cern.jet.random.tdouble.StudentT(genotype.length - 4, randomEngine);
		if (tInteraction < 0) {
			pValueInteraction = tDistColt.cdf(tInteraction);
			if (pValueInteraction < 2.0E-323) {
				pValueInteraction = 2.0E-323;
			}
			zScoreInteraction = cern.jet.stat.tdouble.Probability.normalInverse(pValueInteraction);
		} else {
			pValueInteraction = tDistColt.cdf(-tInteraction);
			if (pValueInteraction < 2.0E-323) {
				pValueInteraction = 2.0E-323;
			}
			zScoreInteraction = -cern.jet.stat.tdouble.Probability.normalInverse(pValueInteraction);
		}
		pValueInteraction *= 2;

		String pValueString = (new java.text.DecimalFormat("0.##E0", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(pValueInteraction);
		if (pValueInteraction > 0.001) {
			pValueString = (new java.text.DecimalFormat("##.###;-##.###", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(pValueInteraction);
		}
		g2d.setFont(new java.awt.Font("Arial", java.awt.Font.BOLD, 14));
		g2d.setColor(new Color(0, 0, 0));
		int posX = marginLeft;
		int posY = marginTop + innerHeight + 20;
		g2d.drawString("Interaction P-Value: " + pValueString, posX, posY);


		for (int g = 0; g <= 2; g++) {

			double valMin = betas[0] + betas[1] * g + minX * betas[2] + betas[3] * g * minX;
			double valMax = betas[0] + betas[1] * g + maxX * betas[2] + betas[3] * g * maxX;
			int posXMin = marginLeft + (int) ((minX - minX) / (maxX - minX) * innerWidth);
			int posYMin = marginTop + innerHeight - (int) ((valMin - minY) / (maxY - minY) * innerHeight);
			int posXMax = marginLeft + (int) ((maxX - minX) / (maxX - minX) * innerWidth);
			int posYMax = marginTop + innerHeight - (int) ((valMax - minY) / (maxY - minY) * innerHeight);

			g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.8f));
			g2d.setColor(new Color(255, 255, 255));
			g2d.setStroke(new java.awt.BasicStroke(5.0f, java.awt.BasicStroke.CAP_ROUND, java.awt.BasicStroke.JOIN_ROUND));
			g2d.drawLine(posXMin, posYMin, posXMax, posYMax);
			if (g < 0.5) {
				g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.30f));
				g2d.setColor(new Color(204, 86, 78));
			} else {
				if (g > 1.5) {
					g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f));
					g2d.setColor(new Color(171, 178, 114));
				} else {
					g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f));
					g2d.setColor(new Color(98, 175, 255));
				}
			}
			g2d.setStroke(new java.awt.BasicStroke(3.0f, java.awt.BasicStroke.CAP_ROUND, java.awt.BasicStroke.JOIN_ROUND));
			g2d.drawLine(posXMin, posYMin, posXMax, posYMax);

		}

		try {
			javax.imageio.ImageIO.write(bimage, "png", new File(fileName));
		} catch (IOException e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
		}


	}
}
