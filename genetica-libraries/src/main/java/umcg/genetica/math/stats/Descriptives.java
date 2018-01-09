/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import umcg.genetica.containers.Pair;

/**
 * @author harmjan
 */
public class Descriptives {
	
	public static double[] m_zScoreToPValue = null;
	public static double[] m_sqrtSample = null;
	
	/**
	 * Calculate the mean for an array of floats
	 *
	 * @param v The array of floats to calculate the mean for
	 * @return The mean of the float array
	 */
	public static double mean(float[] v) {
		double sum = 0;
		for (int k = 0; k < v.length; k++) {
			sum += v[k];
		}
		return (sum / (double) v.length);
	}
	
	/**
	 * Calculate the variance for an array of floats
	 *
	 * @param v    The array of floats for which the variance should be calculated
	 * @param mean The mean of the array of floats
	 * @return The variance for the array of floats
	 */
	public static double variance(float[] v, double mean) {
		double ans = 0.0;
		for (int i = 0; i < v.length; i++) {
			ans += (v[i] - mean) * (v[i] - mean);
		}
		return ans / (v.length - 1);
	}
	
	public static void lookupSqrt(int totalNrSamples) {
		//Fast look-up service for squared root of number of samples:
		double[] sqrtSample = new double[totalNrSamples + 1];
		for (int nrSamples = 0; nrSamples <= totalNrSamples; nrSamples++) {
			sqrtSample[nrSamples] = Math.sqrt((double) nrSamples);
		}
		m_sqrtSample = sqrtSample;
	}
	
	public static void initializeZScoreToPValue() {
		//Fast look-up service for normal distribution:
		JSci.maths.statistics.NormalDistribution normDist = new JSci.maths.statistics.NormalDistribution();
		m_zScoreToPValue = new double[376501];
		
		for (int z = 0; z <= 376500; z++) {
			double zScore = ((double) z - 188250d) / 5000d;
			double pValue;
			//Avoid complementation:
			if (zScore > 0) {
				pValue = normDist.cumulative(-zScore);
			} else {
				pValue = normDist.cumulative(zScore);
			}
			//Two sided test:
			if (pValue > 0.5) {
				pValue = 1 - pValue;
			}
			pValue *= 2.0d;
			//Eventual P Value is stored:
			m_zScoreToPValue[z] = pValue;
		}
	}
	
	public static double convertZscoreToPvalue(double zScore) {
		if (Double.isNaN(zScore)) {
			return 1;
		}
		
		// REMOVED: this is not thread safe!
		if (m_zScoreToPValue == null) {
			throw new IllegalArgumentException("ZScore to P-value lookup table is not initialized.");
		}
		return m_zScoreToPValue[getZScorePvalueIndex(zScore)];
	}
	
	public static int getZScorePvalueIndex(double zScore) {
		int zScoreIndex = (int) ((zScore * 5000.0d) + 188250);
		if (zScoreIndex < 0) {
			zScoreIndex = 0;
		}
		if (zScoreIndex > 376500) {
			zScoreIndex = 376500;
		}
		return zScoreIndex;
	}
	
	
	public static double zScore(double value, double mean, double variance) {
		if (variance > 0.0 && mean > 0.0) {
			double sd = Math.sqrt(variance);
			double val = value - mean;
			double zscore = val / sd;
			zscore = Math.sqrt(Math.pow(zscore, 2));
			return zscore;
		} else {
			return Double.MAX_VALUE;
		}
	}
	
	public static double mean(double[] v) {
		double sum = 0;
		for (int k = 0; k < v.length; k++) {
			sum += v[k];
		}
		return (sum / (double) v.length);
	}
	
	/**
	 * Simultainiusly calculates mean for two double[] of the same lengths
	 * Warning: Does not check if the length is identical!
	 *
	 * @param v
	 * @param w
	 * @return
	 */
	public static Pair<Double, Double> mean(double[] v, double[] w) {
		double sumV = 0;
		double sumW = 0;
		for (int k = 0; k < v.length; k++) {
			sumV += v[k];
			sumW += w[k];
		}
		return (new Pair<Double, Double>((sumV / (double) v.length), (sumW / (double) v.length)));
	}
	
	public static double variance(double[] v) {
		double mean = mean(v);
		return variance(v, mean);
	}
	
	public static double variance(double[] v, double mean) {
		double ans = 0.0;
		for (int i = 0; i < v.length; i++) {
			ans += (v[i] - mean) * (v[i] - mean);
		}
		return ans / (v.length - 1);
	}
	
	public static double getSqrt(int nrTotalSamples) {
		if (m_sqrtSample == null) {
			System.out.println("ERROR: square-root table not correctly initialized.");
			System.exit(-1);
		}
		return m_sqrtSample[nrTotalSamples];
	}
	
	public static double sum(double[] v) {
		double sum = 0;
		for (double d : v) {
			sum += d;
		}
		return sum;
	}
	
	public static double absSum(double[] v) {
		double sum = 0;
		for (double d : v) {
			sum += Math.abs(d);
		}
		return sum;
	}
	
	
	public static double cityBlockDistance(double[] v, double[] w) {
		double summedDistance = 0;
		for (int i = 0; i < v.length; ++i) {
			summedDistance += Math.abs(v[i] - w[i]);
		}
		return summedDistance;
	}
	
	public static double BrayCurtisDistance(double[] v, double[] w) {
		double summedDistance = 0;
		double sum = 0;
		for (int i = 0; i < v.length; ++i) {
			summedDistance += Math.abs(v[i] - w[i]);
			sum += v[i];
			sum += w[i];
		}
		return (summedDistance / sum);
	}
}
