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

    public static void main(String[] args) {

        try {
            DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData("D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM-first50-5x5.txt");
            ConcurrentCorrelation c = new ConcurrentCorrelation();
            DoubleMatrixDataset<String, String> cormat = c.pairwiseCorrelation(ds);
            double[][] cormatarr = cormat.getMatrix().toArray();

            System.out.println("JAMA");
            EigenvalueDecomposition pcajama = PCA.eigenValueDecomposition(cormatarr);

            PCAojAlgo ojalgo = new PCAojAlgo();
            System.out.println("OJ");
            Eigenvalue<Double> pcaoj = ojalgo.eigenValueDecomposition(cormatarr);
            double[] eigvaroj = ojalgo.getRealEigenValues();
            double[] eigvarjama = PCA.getRealEigenvalues(pcajama);

            System.out.println("EigenValues");
            System.out.println(eigvarjama.length + "\t" + eigvaroj.length);
            for (int i = 0; i < eigvaroj.length; i++) {
                System.out.println(i + "\t" + eigvarjama[i] + "\t" + eigvaroj[i]);
            }

            System.out.println("EigenVectors");

            for (int i = 0; i < eigvaroj.length; i++) {
                double[] eigenvectorJama = PCA.getEigenVector(pcajama, i);
                double[] eigenvectorOJ = ojalgo.getEigenVector(i);

                for (int j = 0; j < eigenvectorJama.length; j++) {
                    System.out.println(i + "\t" + j + "\t" + eigenvectorJama[j] + "\t" + eigenvectorOJ[j]);
                }
            }

            System.out.println();
            System.out.println("variance");
            for (int i = 0; i < eigvaroj.length; i++) {
                double varoj = ojalgo.getEigenValueVar(i);
                double varjama = PCA.getEigenValueVar(pcajama.getRealEigenvalues(), i);

                System.out.println(i + "\t" + varjama + "\t" + varoj);
            }
//

            double cumExpVarPCAOj = 0;
            double cumExpVarPCAJa = 0;
            double[] eigenValuesOj = ojalgo.getRealEigenValues();
            double[] eigenValuesJa = PCA.getRealEigenvalues(pcajama);
            for (int pca = 0; pca < eigvaroj.length; pca++) {
                double expVarPCAOj = ojalgo.getEigenValueVar(pca);
                double expVarPCAJa = PCA.getEigenValueVar(eigenValuesJa, pca);
//                double[] pca1ExpEigenVector = ojalgo.getEigenVector(pca);

                int pcaNr = pca + 1;
                cumExpVarPCAOj += expVarPCAOj;
                cumExpVarPCAJa += expVarPCAJa;

                System.out.println("PCA:\t" + pcaNr + "\t" + eigenValuesOj[eigenValuesOj.length - 1 - pca] + "\t" + eigenValuesJa[eigenValuesJa.length - 1 - pca] + "\t\t" + expVarPCAOj + "\t" + cumExpVarPCAOj + "\t\t" + expVarPCAJa + "\t" + cumExpVarPCAJa);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }


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
        double result = eigenvalues[eigenvalues.length - 1 - pca] / sumEigenvalues;
        return result;
    }


}
