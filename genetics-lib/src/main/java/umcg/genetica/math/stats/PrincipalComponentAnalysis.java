/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package umcg.genetica.math.stats;

/**
 *
 * @author harmjan
 */
public class PrincipalComponentAnalysis {

    private Jama.EigenvalueDecomposition eig;
    private double[][] eigenVectors;

    public void eigenValueDecomposition(double[][] data, int numPCAs) {
        Jama.Matrix m = new Jama.Matrix(data);
        eig = m.eig();
	//Performe eigenvalue decomposition:

	eigenVectors = new double[data.length][data.length];
	for (int pca=0; pca<numPCAs; pca++) {
	    eigenVectors[pca] = getEigenVector(eig, pca);
	}
    }

    private double[] getEigenVector(Jama.EigenvalueDecomposition eig, int pca) {
        Jama.Matrix eigenValueMatrix = eig.getV();
        double[][] eigenValueMat = eigenValueMatrix.getArray();
        double[] eigenVector = new double[eigenValueMat.length];
        for (int i = 0; i < eigenValueMat.length; i++) {
            eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca];
        }
        return eigenVector;
    }

    public double getEigenValueVar(double[] eigenValues, int pca) {
        double sumEigenvalues = 0.0;
        for (Double d : eigenValues) {
            sumEigenvalues += Math.abs(d);
        }
        double result = eigenValues[eigenValues.length - 1 - pca] / sumEigenvalues;
        return result;
    }

    public double[][] getDataMatrixPCScores(double[][] dataMatrix, int size, int sampleCount) {
	double[][] dataMatrixPCScores = new double[size][sampleCount];
	for (int sample = 0; sample < sampleCount; sample++) {
	    for (int p=0; p<size; p++) {
		for (int snp=0; snp<size; snp++) {
		    double probeCoefficient = eigenVectors[p][snp];
		    dataMatrixPCScores[p][sample]+= dataMatrix[snp][sample] * probeCoefficient;
		}
	    }
	}
	return dataMatrixPCScores;
    }

    public double[] getRealEigenvalues() {
	return eig.getRealEigenvalues();
    }
}
