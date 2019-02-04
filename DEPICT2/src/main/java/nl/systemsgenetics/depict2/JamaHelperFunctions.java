/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

/**
 *
 * @author patri
 */
public class JamaHelperFunctions {

	/**
	 *
	 * @param data must be a square matrix, no check
	 * @return
	 */
	public static Jama.EigenvalueDecomposition eigenValueDecomposition(double[][] data) {
		Jama.Matrix m = new Jama.Matrix(data, data.length, data.length);
		Jama.EigenvalueDecomposition eig = m.eig();
		return eig;
	}

	public static double[] getEigenVector(Jama.EigenvalueDecomposition eig, double[] eigenValues, int pca) {
		Jama.Matrix eigenValueMatrix = eig.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
		}
		return eigenVector;
	}

	public static double[] getEigenVector(Jama.EigenvalueDecomposition eig, int pca) {
		Jama.Matrix eigenValueMatrix = eig.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
		}
		return eigenVector;
	}
	
	public static double getEigenValueVar(double[] eigenValues, int pca) {
		double sumEigenvalues = 0.0;
		for (double d : eigenValues) {
			sumEigenvalues += Math.abs(d);
		}
		return eigenValues[eigenValues.length - 1 - pca] / sumEigenvalues;
	}


}
