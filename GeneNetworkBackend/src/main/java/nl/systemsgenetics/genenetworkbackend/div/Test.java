/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Probability;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

/**
 *
 * @author patri
 */
public class Test {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {

		DoubleRandomEngine randomEngine = new DRand();

		double[] vals1 = {0d, 10d, 15d,50d, 5, 3};
		double[] vals2 = {3,11,12,10, 6, 4};

		double correlation = JSci.maths.ArrayMath.correlation(vals1, vals2);

		StudentT tDistColt = new StudentT(vals1.length / 2 - 2, randomEngine);
		double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (vals1.length / 2 - 2)));
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
		pValue *= 2;
		
		System.out.println("cor: " + correlation);
		System.out.println("T: " + t);
		System.out.println("p-value: " + pValue);
		System.out.println("zScore: " + zScore);
		

	}

}
