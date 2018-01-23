/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import cern.jet.stat.tdouble.Probability;
import umcg.genetica.containers.Pair;
import umcg.genetica.util.RankDoubleArray;


/**
 * @author harmjan
 */
public class Correlation {
	
	public static double[][] m_correlationToZScore;
	
	public static double rankCorrelate(double[] x, double[] y) {
		
		RankDoubleArray rda = new RankDoubleArray();
		if (x == null || y == null) {
			throw new IllegalArgumentException("Error calculating correlation: x or y is null!");
		} else if (y.length != x.length) {
			throw new IllegalArgumentException("Error calculating correlation: x and y should have the same length!");
		}
		double[] xNew = new double[x.length];
		double[] yNew = new double[y.length];
		for (int a = 0; a < x.length; a++) {
			xNew[a] = x[a];
			yNew[a] = y[a];
		}
		
		double[] xranked = rda.rank(xNew);
		double[] yranked = rda.rank(yNew);
		
		return correlate(xranked, yranked);
	}
	
	public static void correlationToZScore(int maxNrSamples, int granularity) {
		//Fast look-up service to determine P-Value for individual Spearman correlation coefficient, given sample size:
		if (m_correlationToZScore == null || m_correlationToZScore.length < maxNrSamples + 1) {
//            System.out.println("Creating new LookupTable: "+maxNrSamples);
//            if(m_correlationToZScore != null){
//                System.out.println("Length: "+m_correlationToZScore.length+"\treq: "+maxNrSamples);
//            }
			double[][] correlationToZScore = new double[maxNrSamples + 1][granularity];
			cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
			
			for (int nrSamples = 0; nrSamples <= maxNrSamples; nrSamples++) {
				if (nrSamples < 3) {
					for (int sc = 0; sc <= 2000; sc++) {
						correlationToZScore[nrSamples][sc] = 0;
					}
				} else {
					cern.jet.random.tdouble.StudentT tDistColt = new cern.jet.random.tdouble.StudentT(nrSamples - 2, randomEngine);
					//JSci.maths.statistics.TDistribution tDist = new JSci.maths.statistics.TDistribution(nrSamples - 2);
					for (int s = 0; s <= 2000; s++) {
						
						//Determine Spearman correlation R:
						double spearman = (double) (s - 1000) / 1000d;
						if (Math.abs(spearman - -1.0d) < .0000001) {
							spearman = -0.9999;
						}
						if (Math.abs(spearman - 1.0d) < .0000001) {
							spearman = +0.9999;
						}
						
						//Calculate T score:
						double t = spearman / (Math.sqrt((1 - spearman * spearman) / (nrSamples - 2)));
						
						//Look up P value, avoid complementation, due to Double inaccuracies.
						//And lookup Z score, and avoid complementation here as well.
						double pValueColt = 1;
						double zScoreColt = 0;
						if (t < 0) {
							pValueColt = tDistColt.cdf(t);
							if (pValueColt < 2.0E-323) {
								pValueColt = 2.0E-323;
							}
							zScoreColt = Probability.normalInverse(pValueColt);
							//zScoreColt = normDist.inverse(pValueColt);
						} else {
							pValueColt = tDistColt.cdf(-t);
							if (pValueColt < 2.0E-323) {
								pValueColt = 2.0E-323;
							}
							zScoreColt = -Probability.normalInverse(pValueColt);
							//zScoreColt = -normDist.inverse(pValueColt);
						}
						
						
						//Store Z score:
						correlationToZScore[nrSamples][s] = zScoreColt;
					}
				}
			}
			
			m_correlationToZScore = correlationToZScore;
		}
		
	}
	
	public static void correlationToZScore(int maxNrSamples) {
		correlationToZScore(maxNrSamples, 2001);
	}
	
	//Requires mean of 0
	public static double correlateMeanCenteredData(double[] x, double[] y, double varX, double varY) {
		double covarianceInterim = 0;
		for (int i = 0; i < x.length; i++) {
			covarianceInterim += x[i] * y[i];
		}
		return (covarianceInterim / (x.length - 1)) / Math.sqrt(varX * varY);
	}
	
	//Requires mean of 0
	public static double correlateMeanCenteredData(double[] x, double[] y, double sdXsdY) {
		double covarianceInterim = 0;
		for (int i = 0; i < x.length; i++) {
			covarianceInterim += x[i] * y[i];
		}
		return (covarianceInterim / (x.length - 1)) / (sdXsdY);
	}

//    public static double correlateMeanCenteredData(double[] x, double[] y) {
//        
//        double[] xNew = new double[x.length];
//        double[] yNew = new double[y.length];
//        for (int a=0; a<x.length; a++) {
//            xNew[a] = x[a];
//            yNew[a] = y[a];
//        }
//        
//        double meanX = Descriptives.mean(xNew);
//
//        double varX = Descriptives.variance(xNew, meanX);
//
//        double meanY = Descriptives.mean(yNew);
//        double varY = Descriptives.variance(yNew, meanY);
//        for (int i = 0; i < xNew.length; i++) {
//            xNew[i] -= meanX;
//            yNew[i] -= meanY;
//        }
//        return correlateMeanCenteredData(xNew, yNew, varX, varY);
//    }
	
	/**
	 * Fast correlation of two double[]
	 *
	 * @param x
	 * @param y
	 * @return
	 */
	public static double correlate(double[] x, double[] y) {
		if (x.length == y.length) {
			Pair<Double, Double> tmpMean = Descriptives.mean(x, y);
			double meanX = tmpMean.getLeft();
			double meanY = tmpMean.getRight();
			
			double varX = 0;
			double varY = 0;
			
			double covarianceInterim = 0;
			
			for (int a = 0; a < x.length; a++) {
				double varXT = (x[a] - meanX);
				double varYT = (y[a] - meanY);
				
				covarianceInterim += varXT * varYT;
				
				varX += varXT * varXT;
				varY += varYT * varYT;
			}
			
			varY = varY / (y.length - 1);
			varX = varX / (x.length - 1);
			
			double denominator = Math.sqrt(varX * varY);
			double covariance = covarianceInterim / (x.length - 1);
			double correlation = covariance / denominator;
			return correlation;
		} else {
			System.out.println("Warning two arrays of non identical length are put in for correlation.");
			System.out.println("Returning NaN");
			return (Double.NaN);
		}
		
	}
	
	/**
	 * Fast correlation of two double[] with their mean values
	 *
	 * @param x
	 * @param y
	 * @param meanX
	 * @param meanY
	 * @return
	 */
	public static double correlate(double meanX, double meanY, double[] x, double[] y) {
		if (x.length == y.length) {
			double varX = 0;
			double varY = 0;
			
			double covarianceInterim = 0;
			
			for (int a = 0; a < x.length; a++) {
				double varXT = (x[a] - meanX);
				double varYT = (y[a] - meanY);
				
				covarianceInterim += varXT * varYT;
				
				varX += varXT * varXT;
				varY += varYT * varYT;
			}
			
			varY = varY / (y.length - 1);
			varX = varX / (x.length - 1);
			
			double denominator = Math.sqrt(varX * varY);
			double covariance = covarianceInterim / (x.length - 1);
			double correlation = covariance / denominator;
			return correlation;
		} else {
			System.out.println("Warning two arrays of non identical length are put in for correlation.");
			System.out.println("Returning NaN");
			return (Double.NaN);
		}
		
	}
	
	/**
	 * Fast covariation of two double[] ToDo change correlateMeanCenteredData code to covariate
	 * code! ToDo add test
	 *
	 * @param x
	 * @param y
	 * @return
	 */
	public static double covariate(double[] x, double[] y) {
		if (x.length == y.length) {
			Pair<Double, Double> tmpMean = Descriptives.mean(x, y);
			double meanX = tmpMean.getLeft();
			double meanY = tmpMean.getRight();
			
			double covarianceInterim = 0;
			
			for (int a = 0; a < x.length; a++) {
				double varXT = (x[a] - meanX);
				double varYT = (y[a] - meanY);
				
				covarianceInterim += varXT * varYT;
			}
			
			double covariance = covarianceInterim / (x.length - 1);
			return covariance;
		} else {
			System.out.println("Warning two arrays of non identical length are put in for correlation.");
			System.out.println("Returning NaN");
			return (Double.NaN);
		}
		
	}
	
	/**
	 * Fast covariation of two double[] with their mean values ToDo change
	 * correlateMeanCenteredData code to covariate code! ToDo add test
	 *
	 * @param x
	 * @param y
	 * @param meanX
	 * @param meanY
	 * @return
	 */
	public static double covariate(double meanX, double meanY, double[] x, double[] y) {
		if (x.length == y.length) {
			double covarianceInterim = 0;
			
			for (int a = 0; a < x.length; a++) {
				double varXT = (x[a] - meanX);
				double varYT = (y[a] - meanY);
				
				covarianceInterim += varXT * varYT;
			}
			
			double covariance = covarianceInterim / (x.length - 1);
			return covariance;
		} else {
			System.out.println("Warning two arrays of non identical length are put in for correlation.");
			System.out.println("Returning NaN");
			return (Double.NaN);
		}
		
	}
	
	public static double convertCorrelationToZScore(int length, double correlation) {
		int correlationIndex = (int) Math.round(((correlation + 1.0d) * 1000d));
		
		if (m_correlationToZScore == null) {
			throw new IllegalArgumentException("Correlation to ZScore lookup table is not initialized.");
		} else {
			if (m_correlationToZScore.length < length) {
				
				throw new IllegalArgumentException("Correlation to ZScore lookup table is not initialized properly. Expected: " + length + ". Actual length: " + m_correlationToZScore.length);
				
			} else if (correlationIndex > m_correlationToZScore[length].length - 1) {
				throw new IllegalArgumentException("ERROR! correlation: " + correlation + " does not fit Z-score table for " + m_correlationToZScore[length] + " samples (length is: " + m_correlationToZScore[length].length + ")");
			}
			return m_correlationToZScore[length][correlationIndex];
		}
	}
	
	public static double correlate(double[] a1, double[] a2, double mean1, double mean2, double var1, double var2) {
		double denom = Math.sqrt(var1 * var2);
		if (denom != 0) {
			if (a1.length != a2.length) {
				throw new IllegalArgumentException("Arrays must have the same length : " + a1.length + ", " + a2.length);
			}
			double ans = 0.0;
			for (int i = 0; i < a1.length; i++) {
				ans += (a1[i] - mean1) * (a2[i] - mean2);
			}
			return ans / (a1.length - 1) / denom;
		} else {
			if (var1 == 0 && var2 == 0) {
				return (1.0);
			} else {
				return (0.0); // impossible to correlateMeanCenteredData a null signal with another
			}
		}
	}
}
