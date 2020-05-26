/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.mixupmapper.stat;

/**
 *
 * @author harm-jan
 */
public class PCA {

    Jama.EigenvalueDecomposition eig;
    
    public void eigenValueDecomposition(double[][] data) {
        Jama.Matrix m = new Jama.Matrix(data);
        eig = m.eig();
    }

    public double[] getRealEigenvalues(){
        return eig.getRealEigenvalues();
    }


    public double[] getEigenVector(double[] eigenValues) {
        Jama.Matrix eigenValueMatrix = eig.getV();
        double[][] eigenValueMat = eigenValueMatrix.getArray();
        double[] eigenVector = new double[eigenValueMat.length];
        for (int i = 0; i < eigenValueMat.length; i++) {
            eigenVector[i] = Math.abs(eigenValueMat[i][eigenValueMat.length - 1] * Math.sqrt(eigenValues[eigenValues.length - 1]));
        }
        return eigenVector;
    }

    public double getEigenValueVar(double[] eigenValues) {
        double sumEigenvalues = 0.0;
        for (Double d : eigenValues) {
            sumEigenvalues += Math.abs(d);
        }
        double result = eigenValues[eigenValues.length - 1] / sumEigenvalues;
        return result;
    }

    

}
