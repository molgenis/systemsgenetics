package umcg.genetica.math;

import Jama.EigenvalueDecomposition;
import com.sun.j3d.utils.geometry.Primitive;
import org.ojalgo.matrix.PrimitiveMatrix;
import org.ojalgo.matrix.decomposition.Eigenvalue;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.PrimitiveDenseStore;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.util.RunTimer;

public class PCAojAlgo {


    private double[] eigenvalues;
    private Eigenvalue<Double> eig;

//    public static void main(String[] args) {
//
//        try {
//            DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData("D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM-first50-5x5.txt");
//            ConcurrentCorrelation c = new ConcurrentCorrelation();
//            DoubleMatrixDataset<String, String> cormat = c.pairwiseCorrelation(ds);
//            double[][] cormatarr = cormat.getMatrix().toArray();
//
//
//            boolean performance = false;
//            if (performance) {
//                RunTimer rt = new RunTimer();
//                int nriter = 5;
//                rt = new RunTimer();
//
//                System.out.println("OjAlgo");
//                for (int i = 0; i < nriter; i++) {
//                    rt.start();
//                    Eigenvalue<Double> pcaoj = PCAojAlgo.eigenValueDecomposition(cormatarr);
//                    System.out.println(i + "\t" + rt.getTimeDesc(rt.getTimeDiff()));
//                }
//                System.out.println();
//                System.out.println("JAMA");
//                for (int i = 0; i < nriter; i++) {
//                    rt.start();
//                    EigenvalueDecomposition pcajama = PCA.eigenValueDecomposition(cormatarr);
//                    System.out.println(i + "\t" + rt.getTimeDesc(rt.getTimeDiff()));
//                }
//
//                System.out.println("OjAlgo");
//                for (int i = 0; i < nriter; i++) {
//                    rt.start();
//                    Eigenvalue<Double> pcaoj = PCAojAlgo.eigenValueDecomposition(cormatarr);
//                    System.out.println(i + "\t" + rt.getTimeDesc(rt.getTimeDiff()));
//                }
//                System.out.println();
//
//            }
//
//            System.out.println("JAMA");
//            EigenvalueDecomposition pcajama = PCA.eigenValueDecomposition(cormatarr);
//
//            System.out.println("OJ");
//            Eigenvalue<Double> pcaoj = PCAojAlgo.eigenValueDecomposition(cormatarr);
//            double[] eigvaroj = PCAojAlgo.getRealEigenValues();
//            double[] eigvarjama = PCA.getRealEigenvalues(pcajama);
//
//            System.out.println("EigenValues");
//            System.out.println(eigvarjama.length + "\t" + eigvaroj.length);
//            for (int i = 0; i < eigvaroj.length; i++) {
//                System.out.println(i + "\t" + eigvarjama[i] + "\t" + eigvaroj[i]);
//            }
//
//            System.out.println("EigenVectors");
//
//
//            for (int i = 0; i < eigvaroj.length; i++) {
//                double[] eigenvectorJama = PCA.getEigenVector(pcajama, i);
//                double[] eigenvectorOJ = PCAojAlgo.getEigenVector(i);
//
//                for (int j = 0; j < eigenvectorJama.length; j++) {
//                    System.out.println(i + "\t" + j + "\t" + eigenvectorJama[j] + "\t" + eigenvectorOJ[j]);
//                }
//            }
//
//            System.out.println();
//            System.out.println("variance");
//            for (int i = 0; i < eigvaroj.length; i++) {
//                double varoj = PCAojAlgo.getEigenValueVar(i);
//                double varjama = PCAojAlgo.getEigenValueVar(i);
//
//                System.out.println(i + "\t" + varjama + "\t" + varoj);
//            }
//
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//
//    }


    public void PCAojAlgo() {

    }


    public Eigenvalue<Double> eigenValueDecomposition(double[][] data) {

        PrimitiveMatrix matrix = PrimitiveMatrix.FACTORY.rows(data);

        eig = Eigenvalue.make(matrix, true);
        if (!eig.decompose(matrix)) {
            throw new RuntimeException("Decomposition failed");
        }
        return eig;
    }

    public double[] getRealEigenValues() {
        if (eig == null) {
            throw new RuntimeException("Eigenvalues requested but no decomposition performed");
        }
        if (eigenvalues == null) {
            double[] eigenvaluestmp = eig.getEigenvalues().toRawCopy1D();
            eigenvalues = new double[eigenvaluestmp.length];
            for (int d = 0; d < eigenvalues.length; d++) {
                eigenvalues[d] = eigenvaluestmp[eigenvaluestmp.length - 1 - d];
            }
        }
        return eigenvalues;
    }

    public double[] getEigenVector(int pca) {
        if (eig == null) {
            throw new RuntimeException("Eigenvector requested but no decomposition performed");
        }
        MatrixStore<Double> eigenValueMatrix = eig.getV();
        int nrrows = (int) eigenValueMatrix.countRows();
        double[] eigenVector = new double[nrrows];
        for (int i = 0; i < nrrows; i++) {
            eigenVector[i] = eigenValueMatrix.get(i, pca); //eigenValueMat[i][pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
        }
        return eigenVector;
    }

    public double getEigenValueVar(int pca) {

        if (eigenvalues == null) {
            getRealEigenValues();
        }
        double sumEigenvalues = 0.0;
        for (Double d : eigenvalues) {
            sumEigenvalues += Math.abs(d);
        }
        double result = eigenvalues[pca] / sumEigenvalues;
        return result;
    }


}
