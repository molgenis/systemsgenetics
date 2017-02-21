package deconvolution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.lang.Math;

public class Statistics 
{
	/*
	 * Several functions for statistical computations
	 */
    public Statistics(){}   
    
    public double[] log10(double[] data){
    	/*
    	 * Takes the log10 of all values in vector data and returns the vector of logged values
    	 * 
    	 * @param data Vector to calculate log10 for
    	 */
    	for (int i = 0; i < data.length; i++){
    		data[i] = Math.log10(data[i]);
    	}
    	return(data);
    }

    static public double logmodulus(double value){
    	/* log modulus transformation, so that negative numbers are preserved
		 * Take log of absolute value, put back sign
		 * 
		 * @param value value to take log modulus of
		 */
		
		Boolean negative = false;
		if(value < 0){
			negative = true;
		}
		double logValue = Math.log10(Math.abs(value)+1);
		if (negative){
			// if original value was negative, put sign b ack
			logValue *= -1;
		}
		return(logValue);
    }
    
    static public double[] normalize(double[] data){
    	/*
    	 * Normalize a vector
    	 * 
    	 * @param data vector to normalize
    	 * 
    	 */
    	double std = getStdDev(data);
    	double mean = getMean(data);
    	for(int i = 0; i < data.length; i++){
    		data[i] = (data[i]-mean)/std;
    	}
    	return(data);
    }
    
    static public ArrayList<Double> normalize(ArrayList<Double> data){
    	/*
    	 * Normalize an arraylist
    	 * 
    	 * @param data Arraylist top normalize
    	 */
    	double std = getStdDev(data);
    	double mean = getMean(data);
    	for(int i = 0; i < data.size(); i++){
    		data.set(i, (data.get(i)-mean)/std);
    	}
    	return(data);
    }
    
    public double[] normalizeKeepMean(double[] expression, Boolean testNormality){
    	/*
    	 * Normalize vector, reinsert the mean
    	 */
    	double mean = getMean(expression);
    	double[] normalizedData = normalize(expression);
    	if(testNormality){
    		//TODO: Test normality
    	}
    	// put back the power
    	for(int i = 0; i < normalizedData.length; i++){
    		normalizedData[i] = normalizedData[i] + mean;
    	}
    	return(expression);
    }
    
    public double[] normalizeKeepExponential(double[] expression, Boolean testNormality){
    	double[] normalizedData = normalize(expression);
    	if(testNormality){
    		//TODO: Test normality
    	}
    	// put back the power
    	for(int i = 0; i < normalizedData.length; i++){
    		normalizedData[i] = Math.pow(normalizedData[i], 2);
    	}
    	return(expression);
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


    public double median(double[] data){
    	/*
    	 * Calculate median of a vector
    	 * 
    	 * @param data Vector to calculate median for
    	 */
	    Arrays.sort(data);
	
	    if (data.length % 2 == 0){
	       return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
	    } 
	    else{
	       return data[data.length / 2];
	    }
    }
    
    public static List<List<String>> transposeStringMatrix(List<List<String>> cellcountTable){
    	/*
    	 * Transpose a matrix
    	 * 
    	 * @param cellcountTable Matrix to transpose
    	 */
        List<List<String>> temp = new ArrayList<List<String>>();
        for (int i = 0; i < cellcountTable.size(); i++)
            for (int j = 0; j < cellcountTable.get(0).size(); j++)
                temp.get(i).set(j, cellcountTable.get(i).get(j));
        return temp;
    }
    
    public static List<List<Double>> transposeDoubleMatrix(List<List<Double>> cellcountTable){
    	/*
    	 * Transpose a matrix
    	 * 
    	 * @param cellcountTable Matrix to transpose
    	 */
        List<List<Double>> temp = new ArrayList<List<Double>>();
        for (int i = 0; i < cellcountTable.size(); i++)
            for (int j = 0; j < cellcountTable.get(0).size(); j++)
                temp.get(i).set(j, cellcountTable.get(i).get(j));
        return temp;
    }
   
    public static double calculateSpearmanTwoTailedPvalue(double spearmanCorrelation, int sampleSize){
    	/*
    	 * based on https://www.researchgate.net/post/How_do_you_calculate_a_p_value_for_spearmans_rank_correlation
    	 * and https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient#Determining_significance
    	 * 
    	 * @param spearmanCorrelation Spearman correlation from SpearmansCorrelation()
    	 * 
    	 * @param sampleSize Number of samples used to calculate correlation
    	 */
    	double z = Math.sqrt((sampleSize-3)/1.06) * atanh(spearmanCorrelation);
    	NormalDistribution normalDistribution = new NormalDistribution();
    	double p = 2*normalDistribution.cumulativeProbability(-Math.abs(z));

    	if (Double.isNaN(p)){
    		p = 1;
    	}
    	return(p);
    }
    
    private static double atanh(double x){
    	/*
    	 * Calculate hyperbolic Tangent of value (from https://github.com/maths/dragmath/blob/master/lib/jep/src/org/nfunk/jep/function/ArcTanH.java)
    	 * 
    	 * @param x Value to calculate atanh for
    	 */
    	return Math.log((1+x)/(1-x))/2;
    }
    
}

