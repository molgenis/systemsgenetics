package deconvolution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.lang.Math;

public class Statistics 
{

    public Statistics(double[] data) 
    {
    }   
    
    public double[] log10(double[] data){
    	/*
    	 * Takes the log10 of all values in vector data and returns the vector of logged values 
    	 */
    	for (int i = 0; i < data.length; i++){
    		data[i] = Math.log10(data[i]);
    	}
    	return(data);
    }

    static public double logmodulus(double value){
		Boolean negative = false;
		/* log modulus transformation, so that negative numbers are preserved
		 * Take log of absolute value, put back sign
		 */
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
    	double std = getStdDev(data);
    	double mean = getMean(data);
    	for(int i = 0; i < data.length; i++){
    		data[i] = (data[i]-mean)/std;
    	}
    	return(data);
    }
    
    static public List<Double>normalize(List<Double> data){
    	double std = getStdDev(data);
    	double mean = getMean(data);
    	for(int i = 0; i < data.size(); i++){
    		data.set(i, (data.get(i)-mean)/std);
    	}
    	return(data);
    }
    
    public double[] normalizeKeepMean(double[] data, Boolean testNormality){
    	double mean = getMean(data);
    	double[] normalizedData = normalize(data);
    	if(testNormality){
    		//TODO: Test normality
    	}
    	// put back the power
    	for(int i = 0; i < normalizedData.length; i++){
    		normalizedData[i] += mean;
    	}
    	return(data);
    }
    
    public double[] normalizeKeepExponential(double[] data, Boolean testNormality){
    	double[] normalizedData = normalize(data);
    	if(testNormality){
    		//TODO: Test normality
    	}
    	// put back the power
    	for(int i = 0; i < normalizedData.length; i++){
    		normalizedData[i] = Math.pow(normalizedData[i], 2);
    	}
    	return(data);
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
    
    static double getStdDev(List<Double> data)
    {
        return Math.sqrt(getVariance(data));
    }


    public double median(double[] data) 
    {
       
       Arrays.sort(data);

       if (data.length % 2 == 0) 
       {
          return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
       } 
       else 
       {
          return data[data.length / 2];
       }
    }
    
    public static List<List<String>> transposeStringMatrix(List<List<String>> cellcountTable){
        List<List<String>> temp = new ArrayList<List<String>>();
        for (int i = 0; i < cellcountTable.size(); i++)
            for (int j = 0; j < cellcountTable.get(0).size(); j++)
                temp.get(i).set(j, cellcountTable.get(i).get(j));
        return temp;
    }
    
    public static List<List<Double>> transposeDoubleMatrix(List<List<Double>> cellcountTable){
        List<List<Double>> temp = new ArrayList<List<Double>>();
        for (int i = 0; i < cellcountTable.size(); i++)
            for (int j = 0; j < cellcountTable.get(0).size(); j++)
                temp.get(i).set(j, cellcountTable.get(i).get(j));
        return temp;
    }
    
}

