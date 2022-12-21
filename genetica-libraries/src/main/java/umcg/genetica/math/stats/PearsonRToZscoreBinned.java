/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Probability;
import smile.stat.distribution.TDistribution;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;


/**
 *
 * @author patri
 */
public class PearsonRToZscoreBinned {

	private final double[] zscoreLookupTable;
	private final double[] rLookupTable;

	private final int totalNumberOfBins;
	private final int binsPerSide;
	private final double halfStep;
	private final int maxBin;

	private final double maxZscore;

	/**
	 * It is recommended to name instances of this class r2Zscore ;)
	 * Max z-score returned is ~38.
	 *
	 * @param numberOfBins determines precision.
	 * @param samplesUsedForCor
	 * @throws java.lang.Exception
	 */
	public PearsonRToZscoreBinned(final int numberOfBins, final int samplesUsedForCor) throws Exception {

		final int df = samplesUsedForCor - 2;
		DoubleRandomEngine randomEngine = new DRand();

		this.binsPerSide = numberOfBins;
		this.totalNumberOfBins = (numberOfBins * 2) + 1;
		this.maxBin = totalNumberOfBins - 1;
		zscoreLookupTable = new double[totalNumberOfBins];
		rLookupTable = new double[totalNumberOfBins];

		final double stepSize = 2d / totalNumberOfBins;
		halfStep = stepSize / 2;
		StudentT tDistColt = new StudentT(df, randomEngine);

		for (int bin = 0; bin < totalNumberOfBins; ++bin) {
			if (bin == this.binsPerSide) {
				//Enforce r=0 results in z=0
				zscoreLookupTable[bin] = 0;
			} else {
				final double corBinCenter = -1 + (stepSize * bin) + halfStep;
				double t = corBinCenter / (Math.sqrt((1 - corBinCenter * corBinCenter) / (double) (df)));
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

				zscoreLookupTable[bin] = zScore;
				//System.out.println("Bin: " + bin + "\tcenter r: " + corBinCenter + "\tZscore: " + zScore);
			}
		}


		maxZscore = zscoreLookupTable[zscoreLookupTable.length-1];

		// Important to use this T distribution as the others gave massive issues with numerical preciscion
		// either above p < 1e-16 for colt, or arround p < 1e-60? for JSci
		// The values this produces match well with R
		TDistribution tDist = new TDistribution(df);

		for (int bin = 0; bin < totalNumberOfBins; ++bin) {
			if (bin == this.binsPerSide) {
				//Enforce r=0 results in z=0
				rLookupTable[bin] = 0;
			} else {
				final double corBinCenter = -1 + (stepSize * bin) + halfStep;
				double curZ = corBinCenter * maxZscore;
				rLookupTable[bin] = exactZscoreToR(curZ, df, tDist);

				// Hack to deal with slight differences in numerical precision.
				// The inverse conversion caps out arround ~37.5 instead of 38.4
				// In these cases return an R of 1.
				if(Double.isNaN(rLookupTable[bin])) {
					rLookupTable[bin] = 1;
				}

			}
		}

	}

	/**
	 * 
	 * 
	 * @param r pearson r. An r > 1 or r small -1 is accepted for imprecise r calculations
	 * @return  Z-score
	 */
	public double lookupZscoreForR(double r) {

		final long bin;

		if (r >= 1) {
			//this is needed because due to rounding a r of 1 will not fall into any bin. 
			bin = maxBin;
		} else if (r < -1) {
			//-1 exactly will always round properly, no need to test this.
			bin = 0;
		} else {
			bin = Math.round((r + 1) * binsPerSide - halfStep);
		}
//
//		System.out.println("r: " + r);
//		System.out.println("bin: " + bin);
//		System.out.println("z: " + zscoreLookupTable[(int) bin]);
//		System.out.println("---");
		return zscoreLookupTable[(int) bin];

	}

	/**
	 * Lookup the pearson r matching a z-score.
	 * @param z
	 * @return
	 */
	public double lookupRforZscore(double z) {

		final long bin;

		if (z >= maxZscore) {
			bin = maxBin;
		} else if (z < -maxZscore) {
			bin = 0;
		} else {
			bin = Math.round(((z/maxZscore) + 1) * binsPerSide -halfStep);
		}

		return rLookupTable[(int) bin];
	}

	/**
	 * Inplace replacement of Z-scores to pearson r
	 *
	 * @param dataset
	 */
	public void inplaceRToZ(DoubleMatrixDataset dataset) {

		DoubleMatrix2D matrix = dataset.getMatrix();

		final int rows = matrix.rows();
		final int cols = matrix.columns();

		for (int r = 0; r < rows; ++r) {
			for (int c = 0; c < cols; ++c) {
				matrix.setQuick(r, c, lookupZscoreForR(matrix.getQuick(r, c)));
			}
		}

	}

	public void inplaceZToR(DoubleMatrixDataset dataset) {

		DoubleMatrix2D matrix = dataset.getMatrix();

		final int rows = matrix.rows();
		final int cols = matrix.columns();

		for (int r = 0; r < rows; ++r) {
			for (int c = 0; c < cols; ++c) {
				matrix.setQuick(r, c, lookupRforZscore(matrix.getQuick(r, c)));
			}
		}

	}


	public static double exactZscoreToR(double z, double df, TDistribution tDist) {
		double p = Probability.normal(-Math.abs(z));
		double t = tDist.quantile(p);
		t = t*t;
		double cor = Math.sqrt(t / (t+df));

		if (z < 0) {
			cor = -cor;
		}

		return cor;
	}

}
