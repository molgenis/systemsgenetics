/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.math;

/**
 *
 * @author harmjan
 */
public class SVD {
    private static Jama.SingularValueDecomposition svd;
    
    public void PCA(){
        
    }
    
    public static void eigenValueDecomposition(double[][] data) {  
//        System.out.println("Matrix is: "+data.length+"x"+data[data.length-1].length);
        Jama.Matrix m = new Jama.Matrix(data);                                                                                                                                                                                                                                                                                                                                                                                        
//        System.out.println("Performing decomposition");
        svd = m.svd();
    }       
    
    public static double[][] getU(){
        Jama.Matrix m = svd.getU();
        
        return m.getArray();
    }
  
}
