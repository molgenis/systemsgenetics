/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.util.Arrays;
import static nl.systemsgenetics.depict2.JamaHelperFunctions.eigenValueDecomposition;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import static nl.systemsgenetics.depict2.GenePvalueCalculator.createHashColsFromList;

/**
 *
 * @author patri
 */
public class GeneWeightCalculator {

	/**
	 * Code from lude to get weights for genes
	 * 
	 * @param genePvaluesNullGwas
	 * @return 
	 */
	public static DoubleMatrixDataset<String, String> calculateGeneWeights(final DoubleMatrixDataset<String, String> genePvaluesNullGwas) {

		final int numberOfGenes = genePvaluesNullGwas.rows();
		final DoubleMatrixDataset<String, String> genePvaluesNullGwasT = genePvaluesNullGwas.viewDice();

		DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneCorrelation = genePvaluesNullGwasT.calculateCorrelationMatrix();

		final Jama.EigenvalueDecomposition eig = eigenValueDecomposition(genePvaluesNullGwasGeneCorrelation.getMatrixAs2dDoubleArray());
		final double[] eigenValues = eig.getRealEigenvalues();
		double[][] factorLoadings = new double[numberOfGenes][numberOfGenes];

		//Calculate the factorloadings:
		for (int comp = 0; comp < numberOfGenes; comp++) {
			double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
			if (eigenvalue < 0) {
				eigenvalue = 0;
			}
			double sqrtEigenvalue = Math.sqrt(eigenvalue);
			factorLoadings[comp] = JamaHelperFunctions.getEigenVector(eig, comp);
			for (int a = 0; a < numberOfGenes; a++) {
				factorLoadings[comp][a] *= sqrtEigenvalue;
			}
		}

		//Calculate the weights of the individual genes, to be used for the weighed correlation, to account for co-expression between genes:
		DoubleMatrixDataset<String, String> weights = new DoubleMatrixDataset<>(genePvaluesNullGwasGeneCorrelation.getHashColsCopy(), createHashColsFromList(Arrays.asList(new String[]{"weight"})));
		for (int p = 0; p < numberOfGenes; p++) {
			double weight = 0;
			for (int comp = 0; comp < numberOfGenes; comp++) {
				double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
				if (eigenvalue < 1) {
					eigenvalue = 1;
				}
				weight += factorLoadings[comp][p] * factorLoadings[comp][p] / eigenvalue;
			}

			weights.setElementQuick(p, 0, weight);

		}
		
		return weights;

	}

}
