package umcg.genetica.math;

import Jama.EigenvalueDecomposition;
import com.sun.j3d.utils.geometry.Primitive;
import org.ojalgo.matrix.PrimitiveMatrix;
import org.ojalgo.matrix.decomposition.Eigenvalue;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.PrimitiveDenseStore;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.text.Strings;
import umcg.genetica.util.RunTimer;

public class PCAojAlgo {


	private double[] eigenvalues;
	private Eigenvalue<Double> eig;
	private MatrixStore<Double> eigenValueMatrix;

	public static void main(String[] args) {

		try {

//
			DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData("D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM-first50-5x5.txt");
//            ds.normalizeRows();

			ConcurrentCorrelation c = new ConcurrentCorrelation();
			DoubleMatrixDataset<String, String> cormat = c.pairwiseCorrelation(ds.viewDice());
			double[][] cormatarr = cormat.getMatrix().toArray();

			System.out.println("JAMA");
			EigenvalueDecomposition pcajama = PCA.eigenValueDecomposition(cormatarr);

			PCAojAlgo ojalgo = new PCAojAlgo();
			System.out.println("OJ");
			ojalgo.eigenValueDecomposition(cormatarr);
			double[] eigvaroj = ojalgo.getRealEigenValues();
			double[] eigvarjama = PCA.getRealEigenvalues(pcajama);

			System.out.println("EigenValues");
			System.out.println(eigvarjama.length + "\t" + eigvaroj.length);
			for (int i = 0; i < eigvaroj.length; i++) {
				System.out.println(i + "\t" + eigvarjama[i] + "\t" + eigvaroj[i]);
			}

			System.out.println("EigenVectors OJ");
			for (int i = 0; i < eigvaroj.length; i++) {
				double[] eigenvectorJama = PCA.getEigenVector(pcajama, i);
				double[] eigenvectorOJ = ojalgo.getEigenVector(i);

//				for (int j = 0; j < eigenvectorJama.length; j++) {
//					System.out.println(i + "\t" + j + "\t" + eigenvectorJama[j] + "\t" + eigenvectorOJ[j]);
//				}
				System.out.println(Strings.concat(eigenvectorOJ, Strings.tab));
			}

			System.out.println();
			System.out.println("EigenVectors JAMA");
			for (int i = 0; i < eigvaroj.length; i++) {
				double[] eigenvectorJama = PCA.getEigenVector(pcajama, i);
				double[] eigenvectorOJ = ojalgo.getEigenVector(i);

//				for (int j = 0; j < eigenvectorJama.length; j++) {
//					System.out.println(i + "\t" + j + "\t" + eigenvectorJama[j] + "\t" + eigenvectorOJ[j]);
//				}
				System.out.println(Strings.concat(eigenvectorJama, Strings.tab));
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

	boolean ascendingorder = false;

	public void eigenValueDecomposition(double[][] data) {
		System.out.println("Performing eigenvector decomposition on " + data.length + " x " + data[data.length - 1].length + " matrix ");
		PrimitiveMatrix matrix = PrimitiveMatrix.FACTORY.rows(data);


		eig = Eigenvalue.make(matrix, true);

		if (!eig.decompose(matrix)) {
			throw new RuntimeException("Decomposition failed");
		}

		// check if eigenvalues are in descending order
		eigenvalues = eig.getEigenvalues().toRawCopy1D();

		for (int i = 0; i < eigenvalues.length - 1; i++) {
			if (eigenvalues[i] < eigenvalues[eigenvalues.length - 1]) {
				ascendingorder = true;
			}
		}


		if (ascendingorder) {
			System.out.println("WARNING: eigenvalues are in ascending order. Will flip them for you.");
			// invert eigenvalues
			double[] tmpeig = new double[eigenvalues.length];
			for (int d = 0; d < eigenvalues.length; d++) {
				tmpeig[d] = eigenvalues[eigenvalues.length - 1 - d];
			}
			eigenvalues = tmpeig;
		} else {
			System.out.println("Eigenvalues are in descending order. ");
		}

		eigenValueMatrix = eig.getV();
	}

	public double[] getRealEigenValues() {
		if (eig == null) {
			throw new RuntimeException("Eigenvalues requested but no decomposition performed");
		}
		return eigenvalues;
	}

	public double[] getEigenVector(int pca) {
		if (eig == null) {
			throw new RuntimeException("Eigenvector requested but no decomposition performed");
		}

		int nrcols = (int) eigenValueMatrix.countColumns();
		double[] eigenVector = null;
		if (ascendingorder) {
			return eigenValueMatrix.sliceColumn(nrcols - 1 - pca).toRawCopy1D();
		} else {
			return eigenValueMatrix.sliceColumn(pca).toRawCopy1D();
		}

	}

	public double getEigenValueVar(int pca) {
		if (eig == null) {
			throw new RuntimeException("Eigenvalue variance requested but no decomposition performed");
		}

		double sumEigenvalues = 0.0;
		for (Double d : eigenvalues) {
			sumEigenvalues += Math.abs(d);
		}
		double result = eigenvalues[pca] / sumEigenvalues;
		return result;
	}


}
