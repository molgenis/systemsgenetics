package umcg.genetica.math.stats;

import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.util.Primitives;

import java.util.ArrayList;

/**
 * This class was freely adapted from Abecasis' METAL package for meta-analysis
 * (http://www.sph.umich.edu/csg/abecasis/Metal/)
 *
 * @author harmjan
 */
public class Heterogeneity {
	
	public static Triple<Double, Double,Integer> getISq(ArrayList<Double> isqd, ArrayList<Integer> isqn) {
		
		double[] z = Primitives.toPrimitiveArr(isqd.toArray(new Double[0]));
		int[] n = Primitives.toPrimitiveArr(isqn.toArray(new Integer[0]));
		
		return getISq(z, n);
	}
	
	public static Triple<Double, Double,Integer> getISq(double[] datasetZ, int[] datasetWeights) {
		
		double weightedZ = 0;
		int totalSample = 0;
		for (int d = 0; d < datasetZ.length; d++) {
			if (!Double.isNaN(datasetZ[d])) {
				weightedZ += Math.sqrt(datasetWeights[d]) * datasetZ[d];
				totalSample += datasetWeights[d];
			}
		}
		
		double hetSum = 0;
		int hetDf = 0;
		for (int d = 0; d < datasetZ.length; d++) {
			if (!Double.isNaN(datasetZ[d])) {
				double expectedZ = Math.sqrt(datasetWeights[d]) * weightedZ / totalSample;
				
				hetSum += (datasetZ[d] - expectedZ) * (datasetZ[d] - expectedZ);
				hetDf++;
			}
		}
		
		double p = 1d;
		double i = 0d;
		
		double iExp = ((hetSum - hetDf + 1) / hetSum) * 100d;
		if (hetDf <= 1 || hetSum < 1E-7) {
			p = 1;
		} else {
			p = ChiSquare.getP(hetDf - 1, hetSum);
		}
		
		
		if (hetSum <= (hetDf - 1) || hetDf <= 1) {
			i = 0;
		} else {
			i = iExp;
		}
		
		if (i > 0) {
			i /= 100;
		}
		return new Triple<Double, Double,Integer>(i, p, totalSample);
	}
	
	
}


/*
// sum Zscores
statistics[marker] += sqrt(w) * z;
weights[marker] += w;

// for each dataset:
ez = sqrt(dsSampleSize) * statistics[marker] / weights[marker]
hetstatistic += (z-ez) * (z-ez);
hetDegreesOfFreedom[marker]++;

// final ISq / p
double p =
               (hetStatistic[marker] < 1e-7 || hetDegreesOfFreedom[marker] <= 1) ?
               1.0 : chidist(hetStatistic[marker], hetDegreesOfFreedom[marker] - 1);
double I2 =
               (hetStatistic[marker] <= hetDegreesOfFreedom[marker] - 1) || hetDegreesOfFreedom[marker] <= 1 ?
               0.0 : (hetStatistic[marker] - hetDegreesOfFreedom[marker] + 1) / hetStatistic[marker] * 100.;
 */
