/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DoubleStatistic;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import java.util.LinkedHashMap;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class PcaColt {

	private final DoubleMatrixDataset<String, String> eigenvectors;
	private final DoubleMatrixDataset<String, String> eigenValues;
	private final DoubleMatrixDataset<String, String> pcs;

	public PcaColt(DoubleMatrixDataset<String, String> dataset, final boolean center) {

		if (center) {
			dataset = dataset.duplicate();
			dataset.centerColumns();
		}

		final DoubleMatrix2D covarianceMatrix = DoubleStatistic.covariance(dataset.getMatrix());

		final DenseDoubleEigenvalueDecomposition eigDecom = new DenseDoubleEigenvalueDecomposition(covarianceMatrix);

		DoubleMatrix2D eigenVectorsMatrix = eigDecom.getV();
		DoubleMatrix1D eigenValuesVector = eigDecom.getRealEigenvalues();

		LinkedHashMap<String, Integer> pcNames = new LinkedHashMap<>(dataset.columns());
		for (int i = 0; i < dataset.columns(); ++i) {
			pcNames.put("Comp" + (i + 1), i);
		}
		
		LinkedHashMap<String, Integer> eigenValueColName = new LinkedHashMap<>(1);
		eigenValueColName.put("EigenValue", 0);
		
		this.eigenValues = new DoubleMatrixDataset<>(pcNames, eigenValueColName);
		this.eigenValues.getCol(0).assign(eigenValuesVector.viewFlip());

		this.eigenvectors = new DoubleMatrixDataset<>(eigenVectorsMatrix.viewColumnFlip(), dataset.getHashCols(), pcNames);
		
		DoubleMatrix2D pcsMatrix = eigenVectorsMatrix.viewDice().zMult(dataset.getMatrix().viewDice(), null);
		
		this.pcs = new DoubleMatrixDataset<>(pcsMatrix.viewDice().viewColumnFlip(), dataset.getHashRows(), pcNames);
		
	}

	public DoubleMatrixDataset<String, String> getEigenvectors() {
		return eigenvectors;
	}

	public DoubleMatrixDataset<String, String> getEigenValues() {
		return eigenValues;
	}

	public DoubleMatrixDataset<String, String> getPcs() {
		return pcs;
	}
	
}
