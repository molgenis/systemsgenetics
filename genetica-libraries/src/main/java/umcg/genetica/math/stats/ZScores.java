/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import JSci.maths.statistics.NormalDistribution;
import cern.jet.random.tdouble.StudentT;
//import cern.jet.random.tdouble.engine.DRand;
//import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Probability;

import java.util.ArrayList;

/**
 * @author juha
 */
public class ZScores {
	
	public static double zToR(double z, int n) {
		double r = z * Math.pow(10, 0.00332727459306132 + -0.500803673176124 * Math.log10(n)) - Math.pow(z, 3) * Math.pow(10, -0.600314180613649 + -1.50039529106925 * Math.log10(n)) + Math.pow(z, 5) * Math.pow(10, -1.3415016045602 + -2.48579979106749 * Math.log10(n)) - Math.pow(z, 7) * Math.pow(10, -2.53051421275022 + -3.39996348904233 * Math.log10(n));
		return r;
	}
	
	private static NormalDistribution normDist;
	
	/**
	 * Calculates a weighted Z-score according to Whitlock's paper:
	 * http://www.ncbi.nlm.nih.gov/pubmed/16135132 Square root of the sample
	 * size is used as the weight for each test.
	 *
	 * @param zScores     Z-scores from individual tests
	 * @param sampleSizes sample sizes of these tests
	 * @return
	 */
	public static double getWeightedZ(float[] zScores, int[] sampleSizes) {
		if (zScores.length != sampleSizes.length) {
			throw new IllegalArgumentException("Zscores and sample sizes should have same length!");
		}
		double weightedZ = 0;
		double sampleSizeSum = 0;
//        int nrNans = 0;
		for (int j = 0; j < zScores.length; j++) {
			if (!Float.isNaN(zScores[j])) {
//                nrNans++;
				weightedZ += Math.sqrt(sampleSizes[j]) * zScores[j];
				sampleSizeSum += sampleSizes[j];
			}
		}
		
		weightedZ /= Math.sqrt(sampleSizeSum);
		return weightedZ;
	}
	
	/**
	 * Calculates a weighted Z-score according to Whitlock's paper:
	 * http://www.ncbi.nlm.nih.gov/pubmed/16135132 Square root of the sample
	 * size is used as the weight for each test.
	 *
	 * @param zScores     Z-scores from individual tests
	 * @param sampleSizes sample sizes of these tests
	 * @return
	 */
	public static double getWeightedZ(double[] zScores, int[] sampleSizes) {
		if (zScores.length != sampleSizes.length) {
			throw new IllegalArgumentException("Zscores and sample sizes should have same length!");
		}
		double weightedZ = 0;
		double sampleSizeSum = 0;
		int nrNans = 0;
		for (int j = 0; j < zScores.length; j++) {
			if (!Double.isNaN(zScores[j])) {
				nrNans++;
				weightedZ += Math.sqrt(sampleSizes[j]) * zScores[j];
				sampleSizeSum += sampleSizes[j];
			}
		}
		
		weightedZ /= Math.sqrt(sampleSizeSum);
		return weightedZ;
	}
	
	/**
	 * Calculates a weighted Z-score according to Whitlock's paper:
	 * http://www.ncbi.nlm.nih.gov/pubmed/16135132 Square root of the sample
	 * size is used as the weight for each test.
	 *
	 * @param zScores     Z-scores from individual tests
	 * @param sampleSizes sample sizes of these tests
	 * @return
	 */
	public static double getWeightedZ(double[] zScores, int[] sampleSizes, double[] weights) {
		if (zScores.length != sampleSizes.length) {
			throw new IllegalArgumentException("Zscores and sample sizes should have same length!");
		}
		double weightedZ = 0;
		double sampleSizeSum = 0;
		for (int j = 0; j < zScores.length; j++) {
			if (!Double.isNaN(zScores[j])) {
				weightedZ += Math.sqrt(sampleSizes[j]) * weights[j] * zScores[j];
				sampleSizeSum += (sampleSizes[j] * weights[j]);
			}
		}
		weightedZ /= Math.sqrt(sampleSizeSum);
		return weightedZ;
	}
	
	/**
	 * Returns the p-value corresponding to the given Z-score assuming a two
	 * tailed-test.
	 *
	 * @param z A Z-score to convert
	 * @return The p-value corresponding to the Z-score
	 */
	public static double zToP(double z) {
		
		if (normDist == null) {
			normDist = new NormalDistribution();
			System.out.println("Creating new Normal Dist");
		}
		double p;
		if (z > 0) {
			p = normDist.cumulative(-z);
		} else {
			p = normDist.cumulative(z);
		}
		if (p > 0.5) {
			p = 1 - p;
		}
		p *= 2.0d;
		
		return p;
	}
	
	/**
	 * Returns the absolute Z-score for a given p-value using a normal
	 * distribution.
	 *
	 * @param p p-value
	 * @return absolute Z-score
	 */
	public static double pToZ(double p) {
		if (p == 1d) {
			return 0;
		} else if (p < 0 || p > 1d) {
			throw new IllegalStateException("P-value should be between 0 and 1.");
		} else if (p == 0d) {
			return 40d;
		}
		
		return Probability.normalInverse(p);
	}
	
	/**
	 * Returns the absolute Z-score for a given p-value using a normal
	 * distribution.
	 *
	 * @param p p-value
	 * @return absolute Z-score
	 */
	public static double pToZTwoTailed(double p) {
		
		p = p / 2;
		return pToZ(p);
	}
	
	/**
	 * Returns the Z-score for a given correlation coefficient from Student's t
	 * distribution.
	 *
	 * @param correlation correlation coefficient (r)
	 * @param nrSamples   number of samples used
	 * @return Z-score
	 */
	public static double correlationToZ(double correlation, int nrSamples, StudentT tDist) {
		
		double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (nrSamples - 2)));
		double pValue = 0;
		double zScore = 0;
		if (t < 0) {
			pValue = tDist.cdf(t);
			if (pValue < 2.0E-323) {
				pValue = 2.0E-323;
			}
			zScore = Probability.normalInverse(pValue);
		} else {
			pValue = tDist.cdf(-t); //Take two sided P-Value
			if (pValue < 2.0E-323) {
				pValue = 2.0E-323;
			}
			zScore = -Probability.normalInverse(pValue);
		}
		return zScore;
	}
	
	/**
	 * Returns the two-sided p-value for a given correlation coefficient from
	 * Student's t distribution.
	 *
	 * @param correlation correlation coefficient (r)
	 * @param nrSamples   number of samples used
	 * @return Z-score
	 */
	public static double correlationToP(double correlation, int nrSamples, StudentT tDist) {
		
		double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (nrSamples - 2)));
		double pValue = 0;
		if (t < 0) {
			pValue = tDist.cdf(t);
			if (pValue < 2.0E-323) {
				pValue = 2.0E-323;
			}
		} else {
			pValue = tDist.cdf(-t); //Take two sided P-Value
			if (pValue < 2.0E-323) {
				pValue = 2.0E-323;
			}
		}
		return pValue;
	}
	
	//    private static void zScoreToCorrelation(int nrSamples) {
//        if (m_zScoreToCorrelation == null || m_zScoreToCorrelation.length != nrSamples) {
//            //Fast look-up service to determine P-Value for individual Spearman correlation coefficient, given sample size:
//            cern.jet.random.StudentT tDistColt = new cern.jet.random.StudentT(nrSamples - 2, (new cern.jet.random.engine.DRand()));
//
//            double[] zScoreToCorrelation = new double[2000001];
//            for (int corrInt = 0; corrInt < 2000001; corrInt++) {
//                double correlation = (double) (corrInt - 1000000) / 1000000d;
//                double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (nrSamples - 2)));
//                double pValue = 0;
//                double zScore = 0;
//                if (t < 0) {
//                    pValue = tDistColt.cdf(t);
//                    if (pValue < 2.0E-323) {
//                        pValue = 2.0E-323;
//                    }
//                    zScore = cern.jet.stat.Probability.normalInverse(pValue);
//                } else {
//                    pValue = tDistColt.cdf(-t); //Take two sided P-Value
//                    if (pValue < 2.0E-323) {
//                        pValue = 2.0E-323;
//                    }
//                    zScore = -cern.jet.stat.Probability.normalInverse(pValue);
//                }
//                int zScoreInt = (int) Math.round(zScore * 10000d + 1000000d);
//                if (zScoreInt < 0) {
//                    zScoreInt = 0;
//                }
//                if (zScoreInt > 2000000) {
//                    zScoreInt = 2000000;
//                }
//                zScoreToCorrelation[zScoreInt] = correlation;
//            }
//            for (int i = 0; i < zScoreToCorrelation.length; i++) {
//
//                double zForBin = (i - 1000000d) / 10000;
//                double p = zToP(zForBin);
//
//                double r = 1;
//                if (p > 0.2 && p < 1) {
//                    System.out.println(p+"\t"+nrSamples);
//                    double t = cern.jet.stat.Probability.studentTInverse(p, nrSamples);
//                    // solve the eqn:
//                    int nrSamplesDf = nrSamples - 2;
//                    double tsquared = t * t;
//                    double tsquaredDF = nrSamplesDf * tsquared;
//                    r = -1 / Math.sqrt(1 + (((double) nrSamples / t)*((double) nrSamples / t) ));
//                    System.out.println(i + "\t" + zForBin + "\t" + p + "\t" + r + "\t" + zScoreToCorrelation[i]);
//                }
//                
//                // System.out.println(zScoreToCorrelation[i] + "\t" + r);
//                // t = r / (Math.sqrt(1-r2/nrSamplesDf));
//                // t * Math.sqrt(1-r2/nrSamplesDf) = r
//                // t2 * ((1-r2)/nrSamplesDf) = r2;
//                // t2 * nrSamplesDf * (1-r2) = r2 * nrSamplesDf;
//                // tsquaredDF - tsquaredDF * r2 = r2 * nrSamplesDf;
//                // tsquaredDF - (tsquaredDF-nrSamplesDf)*r2 = 0;
//                // tsquaredDF / (tsquaredDF-nrSamplesDf) - r2 = 0;
//                // sqrt(tsquaredDF / (tsquaredDF-nrSamplesDf))) = r;
//
//            }
//
//            m_zScoreToCorrelation = zScoreToCorrelation;
//        }
//
//
//    }
	private static JSci.maths.statistics.TDistribution tDist = null;
	private static final double FINDROOT_ACCURACY = 1.0e-15;
	private static final int FINDROOT_MAX_ITERATIONS = 15000;
	
	public static double zScoreToCorrelation(double obsZ, int nrSamples) {
		if (tDist == null || tDist.getDegreesOfFreedom() != nrSamples - 2) {
			tDist = new JSci.maths.statistics.TDistribution(nrSamples - 2);
		}
		
		
		//For a given Z-Score calculate the correlation:
		double obsPValue = zToP(obsZ);
		
		if (obsPValue == 0) {
			obsPValue = Double.MIN_VALUE;
		} else if (obsPValue == 0.5) {
			return 0;
		} else {
			obsPValue /= 2; // cannot divide 0 or close to 0 by 2
		}

//        System.out.println(obsZ + "\t" + obsPValue);
		
		
		double obsT = inverseT(obsPValue);
		if (obsZ < 0) {
			obsT *= -1;
		}
		double corr = Math.sqrt(obsT * obsT / (obsT * obsT + nrSamples));
		if (obsT > 0) {
			corr = -corr;
		}
//        System.out.println("p " + obsPValue + "\tz " + obsZ + "\tc " + corr + "\tt" + obsT);
		if (Double.isNaN(corr)) {
			return 0;
		}
		return corr;
	}
	
	public static double extrapolateZScore(int originalSampleSize, int newSampleSize, double originalZ) {
		double correlation = zScoreToCorrelation(originalZ, originalSampleSize);
		Correlation.correlationToZScore(newSampleSize);
		double newZ = Correlation.convertCorrelationToZScore(newSampleSize, correlation);
//        System.out.println(newZ);
		return newZ;
	}
	
	private static double inverseT(double probability) {
		if (probability < 0 || probability > 1) {
			throw new IllegalArgumentException("Error: probability out of normal range (should be 0 >= x <= 1): " + probability);
		}
		if (probability == 0.0) {
			return -Double.MAX_VALUE;
		}
		if (probability == 1.0) {
			return Double.MAX_VALUE;
		}
		if (probability == 0.5) {
			return 0.0;
		}
		return findRoot(probability, 0.0, -0.5 * Double.MAX_VALUE, 0.5 * Double.MAX_VALUE);
	}
	
	private static double findRoot(double prob, double guess, double xLo, double xHi) {
		double x = guess, xNew = guess;
		double error, pdf, dx = 1.0;
		int i = 0;
		while (Math.abs(dx) > FINDROOT_ACCURACY && i++ < FINDROOT_MAX_ITERATIONS) {
			// Apply Newton-Raphson step
			error = tDist.cumulative(x) - prob;
			if (error < 0.0) {
				xLo = x;
			} else {
				xHi = x;
			}
			pdf = tDist.probability(x);
			if (pdf != 0.0) { // Avoid division by zero
				dx = error / pdf;
				xNew = x - dx;
			}
			// If the Newton-Raphson fails to converge (which for example may be the
			// case if the initial guess is to rough) we apply a bisection
			// step to determine a more narrow interval around the root.
			if (xNew < xLo || xNew > xHi || pdf == 0.0) {
				xNew = (xLo + xHi) / 2.0;
				dx = xNew - x;
			}
			x = xNew;
		}
		return x;
	}
	
	public static double betaToZ(double b, double se, double n) {
		return b / (se * Math.sqrt(n));
	}
	
	
}
