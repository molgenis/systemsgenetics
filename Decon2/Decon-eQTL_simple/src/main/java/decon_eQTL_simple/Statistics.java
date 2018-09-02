package main.java.decon_eQTL_simple;

import org.apache.commons.math3.distribution.NormalDistribution;
import JSci.maths.statistics.FDistribution;
import java.lang.Math;
import cern.jet.random.tdouble.StudentT;

/**
 * Several functions for statistical computations
 */
public class Statistics 
{
	/**
	 * based on https://www.researchgate.net/post/How_do_you_calculate_a_p_value_for_spearmans_rank_correlation
	 * and https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient#Determining_significance
	 * 
	 * this is an approximation, so won't give exactly same results as R cor.test(), but with many samples and
	 * when there is a correlation they approximate eachother
	 * 
	 * @param spearmanCorrelation Spearman correlation from SpearmansCorrelation()
	 * 
	 * @param sampleSize Number of samples used to calculate correlation
	 * 
	 * @return Spearman two-tailed p-value from correlation
	 */
    public static double calculateSpearmanTwoTailedPvalue(double spearmanCorrelation, int sampleSize){
    	double z = Math.sqrt((sampleSize-3)/1.06) * atanh(spearmanCorrelation);
    	NormalDistribution normalDistribution = new NormalDistribution();
    	double p = 2*normalDistribution.cumulativeProbability(-Math.abs(z));

    	if (Double.isNaN(z)){
    		p = 0;
    	}
    	if (Double.isNaN(spearmanCorrelation)){
    		p = 1;
    	}
    	return p;
    }
    
	/**
	 * Calculate hyperbolic Tangent of value (from https://github.com/maths/dragmath/blob/master/lib/jep/src/org/nfunk/jep/function/ArcTanH.java)
	 * 
	 * @param x Value to calculate atanh for
	 */
    private static double atanh(double x){
    	return Math.log((1+x)/(1-x))/2;
    }

	
	/**
	 * Compare and return the p-value of two linear models being
	 * significantly different
	 *
	 * From Joris Meys: http://stackoverflow.com/a/35458157/651779 1.
	 * calculate MSE for the largest model by dividing the Residual Sum of
	 * Squares (RSS) by the degrees of freedom (df) 2. calculate the
	 * MSEdifference by substracting the RSS of both models (result is
	 * "Sum of Sq." in the R table), substracting the df for both models
	 * (result is "Df" in the R table), and divide these numbers. 3. Divide
	 * 2 by 1 and you have the F value 4. calculate the p-value using the F
	 * value in 3 and for df the df-difference in the numerator and df of
	 * the largest model in the denominator. For more info:
	 * http://www.bodowinter.com/tutorial/bw_anova_general.pdf
	 * 
	 * @param sumOfSquaresModelA A vector with the genotype of all samples
	 * for *one* eQTL-gene pair
	 * 
	 * @param sumOfSquaresModelB A vector with the expression levels of all
	 * samples for *one* eQTL-gene pair
	 * 
	 * @param degreesOfFreedomA A 2D list with for all samples the different
	 * cell counts
	 * 
	 * @param degreesOfFreedomB A 2D list with for all samples the different
	 * cell counts
	 * 
	 * @param no_intercept	If intercept was removed to calculate the sum of squares
	 * 
	 * @return The p-value result from comparing two linear models with the
	 * the Anova test
	 */
	public static double anova(double sumOfSquaresModelA, double sumOfSquaresModelB, int degreesOfFreedomA,
			int degreesOfFreedomB, Boolean no_intercept) {
		if (no_intercept) {
			// removing the intercept will give another degree of freedom
			++degreesOfFreedomA;
			++degreesOfFreedomB;
		}
		// Within-group Variance
		double meanSquareError = sumOfSquaresModelA / degreesOfFreedomA;

		int degreesOfFreedomDifference = Math.abs(degreesOfFreedomB - degreesOfFreedomA);
		// Between-group Variance
		// 234111286.801326
		double meanSquareErrorDiff = Math.abs((sumOfSquaresModelA - sumOfSquaresModelB) / (degreesOfFreedomDifference));

		/**
		 * F = Between-group Variance / Within-group Variance <- high value if
		 * variance between the models is high, and variance within the models
		 * is low
		 **/
		if(meanSquareError == 0){
			throw new RuntimeException("meanSquareError should not be 0, no variance in the data?");
		}
		double Fval = meanSquareErrorDiff / meanSquareError;
		/***
		 * Make an F distribution with degrees of freedom as parameter. If full
		 * model and ctModel have the same number of samples, difference in df
		 * is 1 and degreesOfFreedomB are all the terms of the ctModel (so neut%
		 * + eos% + ... + neut% * GT + eos% * GT With 4 cell types and 1891
		 * samples the dfs are 1883 and 1884, giving the below distribution
		 * http://keisan.casio.com/exec/system/1180573186
		 **/
		FDistribution Fdist = new FDistribution(degreesOfFreedomDifference, degreesOfFreedomB);
		/*** Calculate 1 - the probability of observing a lower Fvalue **/
		double pval = 1 - Fdist.cumulative(Fval);
		return pval;
	}
	
	/**
	 * Convert a beta to zscore, and the zscore to a p-value. Adjusted from https://github.com/molgenis/systemsgenetics/blob/master/eqtl-mapping-pipeline/src/main/java/eqtlmappingpipeline/interactionanalysis/InteractionAnalysisTask.java
	 * 
	 * @param beta The beta to convert
	 * @param se The standard error
	 * @param tDistColt t-distribution
	 * @return p-value
	 */
	public static Double convertBetaToP(double beta, double se, StudentT tDistColt) {

		double t = beta / se;
		double p = 1;
		if (t < 0) {
			p = tDistColt.cdf(t);
			if (p < 2.0E-323) {
				p = 2.0E-323;

			}
		} else {
			p = tDistColt.cdf(-t);
			if (p < 2.0E-323) {
				p = 2.0E-323;

			}
		}
		return p;
	}
}

