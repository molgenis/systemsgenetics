/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.methylation;

/**
 *
 * @author Marc Jan
 */
public class DeepCopy {
	
	public static double[][] CopyDoubleArrayArray(double[][] data){
		double[][] newData = new double[data.length][data[1].length];
		
		for(int i=0; i< data.length;++i){
			newData[i]=data[i].clone();
		}
		
		return newData;
	}
	
}
