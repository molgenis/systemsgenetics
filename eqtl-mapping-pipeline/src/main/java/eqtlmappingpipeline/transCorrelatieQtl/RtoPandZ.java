/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.transCorrelatieQtl;

import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Probability;

/**
 *
 * @author patri
 */
public class RtoPandZ {

	private static final DoubleRandomEngine randomEngine = new DRand();

	public static PvalueZscore calculatePandZforCorrelationR(double r, long n) {

		StudentT tDistColt = new StudentT(n - 2, randomEngine);
		double t = r / (Math.sqrt((1 - r * r) / (double) (n - 2)));
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
		
		return new PvalueZscore(pValue, zScore);

	}

}
