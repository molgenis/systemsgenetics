package deconvolution;

import java.util.Arrays;
import java.lang.Math;
// modified from http://stackoverflow.com/a/7988556/651779
public class Statistics 
{
    private double[] data;
    private int size;   

    public Statistics(double[] data) 
    {
    }   
    
    private double[] log10(double[] data){
    	/*
    	 * Takes the log10 of all values in vector data and returns the vector of logged values 
    	 */
    	for (int i = 0; i < data.length; i++){
    		data[i] = Math.log10(data[i]);
    	}
    	return(data);
    }

    private double[] normalize(double data[]){
    	double std = getStdDev(data);
    	double mean = getMean();
    	for(int i = 0; i < data.length; i++){
    		data[i] = (data[i]-mean)/std;
    	}
    	return(data);
    }
    
    public double[] normalizeKeepMean(double[] data, Boolean testNormality){
    	double mean = getMean();
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
    
    private double getMean()
    {
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return sum/size;
    }

    private double getVariance()
    {
        double mean = getMean();
        double temp = 0;
        for(double a :data)
            temp += (mean-a)*(mean-a);
        return temp/size;
    }

    private double getStdDev(double[] data)
    {
    	this.data = data;
    	this.size = data.length;
        return Math.sqrt(getVariance());
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
}

