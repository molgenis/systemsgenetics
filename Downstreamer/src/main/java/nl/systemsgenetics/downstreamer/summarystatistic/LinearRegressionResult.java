package nl.systemsgenetics.downstreamer.summarystatistic;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import org.apache.commons.lang3.NotImplementedException;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.Collection;
import nl.systemsgenetics.downstreamer.runners.DownstreamerRegressionEngine;
import org.apache.commons.math3.distribution.TDistribution;

public class LinearRegressionResult {

	private final String name;
	private final DoubleMatrixDataset<String, String> beta;
	private final DoubleMatrixDataset<String, String> standardError;
	private final DoubleMatrixDataset<String, String> explainedVariance;
	private final int degreesOfFreedom;

	private final DoubleMatrixDataset<String, String> tStatisticCache;
	private final DoubleMatrixDataset<String, String> pValueCache;

	public LinearRegressionResult(Collection<String> rownames, Collection<String> colnames, int degreesOfFreedom, String name) {
		this.beta = new DoubleMatrixDataset<>(rownames, colnames);
		this.standardError = new DoubleMatrixDataset<>(rownames, colnames);
		this.explainedVariance = null;
		this.degreesOfFreedom = degreesOfFreedom;
		this.name = name;
		this.tStatisticCache = null;
		this.pValueCache = null;

	}

	public DoubleMatrixDataset<String, String> getBeta() {
		return beta;
	}

	public DoubleMatrixDataset<String, String> getStandardError() {
		return standardError;
	}

	public DoubleMatrixDataset<String, String> getExplainedVariance() {
		throw new NotImplementedException("Not yet implemented");
	}

	public int getDegreesOfFreedom() {
		return degreesOfFreedom;
	}

	public String getName() {
		return name;
	}

	public void appendBetas(int rowNumber, double[] betas) {
		//int rowNumber = beta.getRowIndex(row);
		for (int i = 0; i < betas.length; i++) {
			beta.setElementQuick(rowNumber, i, betas[i]);
		}
	}

	public void appendSe(int rowNumber, double[] se) {
		//int rowNumber = standardError.getRowIndex(row);
		for (int i = 0; i < se.length; i++) {
			standardError.setElementQuick(rowNumber, i, se[i]);
		}
	}

	public DoubleMatrixDataset<String, String> getTstats() {
		throw new NotImplementedException("Not yet implemented");
	}

	public DoubleMatrixDataset<String, String> getPvalues() {
		throw new NotImplementedException("Not yet implemented");
	}

	public DoubleMatrix1D getBetaForMainEffect() {
		return beta.getCol(DownstreamerRegressionEngine.MAIN_EFFECT_COL_NAME);	
	}
	
	public DoubleMatrix1D getPvalueForMainEffect() {

		final DoubleMatrix1D mainBeta = beta.getCol(DownstreamerRegressionEngine.MAIN_EFFECT_COL_NAME);
		final DoubleMatrix1D mainSe = standardError.getCol(DownstreamerRegressionEngine.MAIN_EFFECT_COL_NAME);
		final DoubleMatrix1D mainPvalues = mainBeta.like();
		
		final int numberPathways = beta.rows();

		for (int i = 0; i < numberPathways; ++i) {
			mainPvalues.setQuick(i, new TDistribution(degreesOfFreedom).cumulativeProbability(-Math.abs(mainBeta.getQuick(i) / mainSe.getQuick(i))) * 2);
		}

		return mainPvalues;
		
	}

	public void save(String basePath, boolean isBinary) throws Exception {

		String finalPath = basePath + "_df_" + degreesOfFreedom;

		if (isBinary) {
			beta.saveBinary(finalPath + "_betas");
			standardError.saveBinary(finalPath + "_se");
		} else {
			beta.save(finalPath + "_betas");
			standardError.save(finalPath + "_se");
		}
	}
}
