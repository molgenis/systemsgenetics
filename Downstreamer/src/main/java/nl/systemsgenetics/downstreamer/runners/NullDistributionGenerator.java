/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import org.apache.commons.math3.stat.descriptive.moment.Kurtosis;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;

/**
 *
 * @author patri
 */
public class NullDistributionGenerator {

	private int n;
	private double mean;
	private double sumOfSquaredDiffs;
	private Kurtosis k;
	private Skewness s; 

	public NullDistributionGenerator() {
		n = 0;
		mean = 0.0;
		sumOfSquaredDiffs = 0.0;
		k = new Kurtosis();
		s = new Skewness();
	}

	public synchronized void addPermutation(final double x) {
		double delta = x - mean;
		mean += delta / ++n;//here we first increase n
		sumOfSquaredDiffs += delta * (x - mean);
		k.increment(x);
		s.increment(x);
	}

	public double getStandardDeviation() {
		if (n < 2) {
			return Double.NaN;
		} else {
			return Math.sqrt(sumOfSquaredDiffs / (n - 1));
		}
	}

	public double getMean() {
		return mean;
	}

	public double getKurtosis(){
		return k.getResult();
	}
	
	public double getSwewness(){
		return s.getResult();
	}

}
