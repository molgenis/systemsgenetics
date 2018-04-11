package decon_eQTL;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.lang.Math;

/**
 * Several functions for statistical computations
 */
public class Statistics 
{
    public Statistics(){}   
	
    /**
	 * Takes the log10 of all values in vector data and returns the vector of logged values
	 * 
	 * @param data Vector to calculate log10 for
	 * 
	 * @return log10 of input vector
	 */    
    public double[] log10(double[] data){
    	for (int i = 0; i < data.length; i++){
    		data[i] = Math.log10(data[i]);
    	}
    	return data;
    }

	/** log modulus transformation, so that negative numbers are preserved
	 * Take log of absolute value, put back sign
	 * 
	 * @param value value to take log modulus of
	 * 
	 * @return Log modulus of value
	 */
    static public double logmodulus(double value){
		Boolean negative = false;
		if(value < 0){
			negative = true;
		}
		double logValue = Math.log10(Math.abs(value)+1);
		if (negative){
			// if original value was negative, put sign b ack
			logValue *= -1;
		}
		return logValue;
    }
    
	/**
	 * Normalize a vector
	 * 
	 * @param data vector to normalize
	 * 
	 * @return Normalized vector
	 */
    static public double[] normalize(double[] data){
    	double std = getStdDev(data);
    	double mean = getMean(data);
    	for(int i = 0; i < data.length; i++){
    		data[i] = (data[i]-mean)/std;
    	}
    	return data;
    }
    
	/**
	 * Normalize an arraylist
	 * 
	 * @param data Arraylist top normalize
	 * 
	 * @return Normalized ArrayList
	 */
    static public ArrayList<Double> normalize(ArrayList<Double> data){
    	double std = getStdDev(data);
    	double mean = getMean(data);
    	for(int i = 0; i < data.size(); i++){
    		data.set(i, (data.get(i)-mean)/std);
    	}
    	return data;
    }
    
	/**
	 * Normalize vector, reinsert the mean
	 * 
	 * @param expression	Vector of expression values
	 * 
	 * @return Normalized vector with mean reinserted
	 */
    public double[] normalizeKeepMean(double[] expression){
    	double mean = getMean(expression);
    	double[] normalizedData = normalize(expression);
    	// put back the power
    	for(int i = 0; i < normalizedData.length; i++){
    		normalizedData[i] = normalizedData[i] + mean;
    	}
    	return expression;
    }
    
    public double[] normalizeKeepExponential(double[] expression){
    	double[] normalizedData = normalize(expression);
    	// put back the power
    	for(int i = 0; i < normalizedData.length; i++){
    		normalizedData[i] = Math.pow(normalizedData[i], 2);
    	}
    	return expression;
    }
    
    static double getMean(double[] data)
    {
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return sum/data.length;
    }

    static double getMean(List<Double> data)
    {
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return sum/data.size();
    }

    
    private static double getVariance(double[] data)
    {
        double mean = getMean(data);
        double temp = 0;
        for(double a :data)
            temp += (mean-a)*(mean-a);
        return temp/data.length;
    }
    
    private static double getVariance(List<Double> data)
    {
        double mean = getMean(data);
        double temp = 0;
        for(double a :data)
            temp += (mean-a)*(mean-a);
        return temp/data.size();
    }

    static double getStdDev(double[] data)
    {
        return Math.sqrt(getVariance(data));
    }
    
    static double getStdDev(List<Double> data){
        return Math.sqrt(getVariance(data));
    }


	/**
	 * Calculate median of a vector
	 * 
	 * @param data Vector to calculate median for
	 * 
	 * @return median of vector
	 */
    public double median(double[] data){
	    Arrays.sort(data);
	
	    if (data.length % 2 == 0){
	       return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
	    } 
	    else{
	       return data[data.length / 2];
	    }
    }
    
	/**
	 * Transpose a matrix
	 * 
	 * @param cellcountTable Matrix to transpose
	 * 
	 * @return Transposed Matrix
	 */
    public static List<List<String>> transposeStringMatrix(List<List<String>> cellcountTable){
        List<List<String>> temp = new ArrayList<List<String>>();
        for (int i = 0; i < cellcountTable.size(); i++)
            for (int j = 0; j < cellcountTable.get(0).size(); j++)
                temp.get(i).set(j, cellcountTable.get(i).get(j));
        return temp;
    }
    
	/**
	 * Transpose a matrix
	 * 
	 * @param cellcountTable Matrix to transpose
	 * 
	 * @return Transposed matrix
	 */
    public static List<List<Double>> transposeDoubleMatrix(List<List<Double>> cellcountTable){
        List<List<Double>> temp = new ArrayList<List<Double>>();
        for (int i = 0; i < cellcountTable.size(); i++)
            for (int j = 0; j < cellcountTable.get(0).size(); j++)
                temp.get(i).set(j, cellcountTable.get(i).get(j));
        return temp;
    }
   
	/**
	 * based on https://www.researchgate.net/post/How_do_you_calculate_a_p_value_for_spearmans_rank_correlation
	 * and https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient#Determining_significance
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

    	if (Double.isNaN(p)){
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
    
    
	/*
	 * Calculate log-likelihood from linear models residuals (same as logLik in R)
	 * 
	 * @param residuals residuals of linear model
	 */
	private static double logLik(double[] residuals){
		double logSummedResiduals = 0;
		for(double residual : residuals){
			logSummedResiduals += Math.pow(residual, 2);
		}
		logSummedResiduals = Math.log(logSummedResiduals);
			
		double logLikelihood = 0.5 * (0 - residuals.length * (Math.log(2 * Math.PI) + 1 - Math.log(residuals.length)+logSummedResiduals));
		return logLikelihood;
	}
	
	/*
	 * Calculate Akaike's "An Information Criterion" with k = 2 (same is AIC in R)
	 * 
	 * @param residuals residuals of linear model
	 * @param npar represents the number of parameters in the fitted model
	 */
	public static double AIC(double[] residuals, int npar){
		int k=2;
		return -2*logLik(residuals) + k*npar;
	}
}

