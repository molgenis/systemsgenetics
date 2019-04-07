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
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class PearsonRToZscoreBinned {

	private final double[] zscoreLookupTable;
	private final int numberOfBins;
	private final double halfStep;
	private final int maxBin;

	/**
	 * It is recommended to name instances of this class r2Zscore ;)
	 *
	 * @param numberOfBins determines precision
	 * @param samplesUsedForCor
	 * @throws java.lang.Exception
	 */
	public PearsonRToZscoreBinned(final int numberOfBins, final int samplesUsedForCor) throws Exception {

//		double power = Math.log10(numberOfBins);
//
//		if (power != Math.round(power)) {
//			throw new Exception("Number of bins must be power of 10");
//		}

		final int df = samplesUsedForCor - 2;
		DoubleRandomEngine randomEngine = new DRand();

		this.numberOfBins = numberOfBins;
		this.maxBin = numberOfBins - 1;
		zscoreLookupTable = new double[numberOfBins];

		final double stepSize = 1d / numberOfBins;
		halfStep = stepSize / 2;

		for (int bin = 0; bin < numberOfBins; ++bin) {

			final double corBinCenter = stepSize * bin + halfStep;

			StudentT tDistColt = new StudentT(df, randomEngine);
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
			
//			System.out.println("Bin: " + bin + "\tcenter r: " + corBinCenter + "\tZscore: " + zScore);

		}
		
		

	}

	public double lookupZscoreForR(double r) {

		
		
		long bin = Math.round(Math.abs(r) * numberOfBins - halfStep);
		
//		System.out.println("r:" + r);
//		System.out.println("bin:" + bin);

		//this is needed because due to rounding a r of 1 will not fall into any bin. Also a r > 1 is accepted for imprecise r calculations
		if (bin > maxBin) {
			bin = maxBin;
		}

		return zscoreLookupTable[(int) bin];

	}

	/**
	 * Inplace replacement of pearson r to Z-scores
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

}
