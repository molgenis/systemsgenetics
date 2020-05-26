/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math;

/**
 *
 * @author harmjan
 */
public class PCA {

    private Jama.EigenvalueDecomposition eig;
    private double[][] eigenVectors;

    public void PCA() {
    }

    public static Jama.EigenvalueDecomposition eigenValueDecomposition(double[][] data) {
//        System.out.println("Matrix is: "+data.length+"x"+data[data.length-1].length);
        Jama.Matrix m = new Jama.Matrix(data);
//        System.out.println("Performing decomposition");
        Jama.EigenvalueDecomposition eig = m.eig();
        return eig;
    }

    public static double[] getRealEigenvalues(Jama.EigenvalueDecomposition eig) {
        return eig.getRealEigenvalues();
    }

    public static double[] getEigenVector(Jama.EigenvalueDecomposition eig, double[] eigenValues, int pca) {
//        Jama.Matrix eigenValueMatrix = eig.getV();
//        double[][] eigenValueMat = eigenValueMatrix.getArray();
//        double[] eigenVector = new double[eigenValueMat.length];
//        for (int i = 0; i < eigenValueMat.length; i++) {
//            eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
//        }
//        return eigenVector;
        return getEigenVector(eig,pca);
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

    public static double  getEigenValueVar(double[] eigenValues, int pca) {
        double sumEigenvalues = 0.0;
        for (Double d : eigenValues) {
            sumEigenvalues += Math.abs(d);
        }
        double result = eigenValues[eigenValues.length - 1 - pca] / sumEigenvalues;
        return result;
    }

    public static double[] getEigenVectorSVD(Jama.SingularValueDecomposition svd, double[] singularValues, int pca) {
        Jama.Matrix eigenValueMatrix = svd.getV();
        double[][] eigenValueMat = eigenValueMatrix.getArray();
        double[] eigenVector = new double[eigenValueMat.length];
        for (int i = 0; i < eigenValueMat.length; i++) {
            eigenVector[i] = eigenValueMat[i][pca] * Math.sqrt(singularValues[pca]);
        }
        return eigenVector;
    }

    public void eigenValueDecomposition(double[][] data, int numPCAs) {
        Jama.Matrix m = new Jama.Matrix(data);
        eig = m.eig();
        //Performe eigenvalue decomposition:

        eigenVectors = new double[data.length][data.length];
        for (int pca = 0; pca < numPCAs; pca++) {
            eigenVectors[pca] = getEigenVector(eig, pca);
        }
    }

    public double[][] getDataMatrixPCScores(double[][] dataMatrix, int size, int sampleCount) {
        double[][] dataMatrixPCScores = new double[size][sampleCount];
        for (int sample = 0; sample < sampleCount; sample++) {
            for (int p = 0; p < size; p++) {
                for (int snp = 0; snp < size; snp++) {
                    double probeCoefficient = eigenVectors[p][snp];
                    dataMatrixPCScores[p][sample] += dataMatrix[snp][sample] * probeCoefficient;
                }
            }
        }
        return dataMatrixPCScores;
    }

    public double[] getRealEigenvalues() {
        return eig.getRealEigenvalues();
    }
}
