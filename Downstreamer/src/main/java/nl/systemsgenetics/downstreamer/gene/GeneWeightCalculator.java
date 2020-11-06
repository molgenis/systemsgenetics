/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.gene;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import static nl.systemsgenetics.downstreamer.gene.JamaHelperFunctions.eigenValueDecomposition;

import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import static nl.systemsgenetics.downstreamer.gene.GenePvalueCalculator.createHashColsFromList;
import org.apache.log4j.Logger;

/**
 *
 * @author patri
 */
public class GeneWeightCalculator {

	private static final Logger LOGGER = Logger.getLogger(Downstreamer.class);

	/**
	 * Code from lude to get weights for genes
	 *
	 * @param genePvaluesNullGwas
	 * @param genes
	 * @return
	 */
	public static DoubleMatrixDataset<String, String> calculateGeneWeights(final DoubleMatrixDataset<String, String> genePvaluesNullGwas, List<Gene> genes, DownstreamerOptions options) {

		Map<String, ArrayList<String>> chrArmToGeneMapping = createChrArmGeneMapping(genes, genePvaluesNullGwas.getHashRows());

		final DoubleMatrixDataset<String, String> weights;

		try (ProgressBar pb = new ProgressBar("Gene weight calculation", chrArmToGeneMapping.size(), ProgressBarStyle.ASCII)) {

			weights = new DoubleMatrixDataset<>(genePvaluesNullGwas.getHashRowsCopy(), createHashColsFromList(Arrays.asList(new String[]{"weight"})));
			chrArmToGeneMapping.keySet().parallelStream().forEach((String chrArm) -> {

				ArrayList<String> armGenes = chrArmToGeneMapping.get(chrArm);

				DoubleMatrixDataset<String, String> genePvaluesNullGwasArm = genePvaluesNullGwas.viewRowSelection(armGenes);
				final int numberOfGenesInArm = genePvaluesNullGwasArm.rows();

				final DoubleMatrixDataset<String, String> genePvaluesNullGwasArmT = genePvaluesNullGwasArm.viewDice();

				DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelation = genePvaluesNullGwasArmT.calculateCorrelationMatrix();

				if (LOGGER.isDebugEnabled()) {
					try {
						genePvaluesNullGwasGeneArmCorrelation.save(new File(options.getOutputBasePath() + "_" + chrArm + "_GeneCorMatrix.txt"));
					} catch (IOException ex) {
						throw new RuntimeException(ex);
					}
				}

				final Jama.EigenvalueDecomposition eig = eigenValueDecomposition(genePvaluesNullGwasGeneArmCorrelation.getMatrixAs2dDoubleArray());
				final double[] eigenValues = eig.getRealEigenvalues();
				double[][] factorLoadings = new double[numberOfGenesInArm][numberOfGenesInArm];

				//Calculate the factorloadings:
				for (int comp = 0; comp < numberOfGenesInArm; comp++) {
					double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
					if (eigenvalue < 0) {
						eigenvalue = 0;
					}
					double sqrtEigenvalue = Math.sqrt(eigenvalue);
					factorLoadings[comp] = JamaHelperFunctions.getEigenVector(eig, comp);
					for (int a = 0; a < numberOfGenesInArm; a++) {
						factorLoadings[comp][a] *= sqrtEigenvalue;
					}
				}

				//Calculate the weights of the individual genes, to be used for the weighed correlation, to account for co-expression between genes:
				for (int p = 0; p < numberOfGenesInArm; p++) {
					double weight = 0;
					for (int comp = 0; comp < numberOfGenesInArm; comp++) {
						double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
						if (eigenvalue < 1) {
							eigenvalue = 1;
						}
						weight += factorLoadings[comp][p] * factorLoadings[comp][p] / eigenvalue;
					}

					weights.setElementQuick(weights.getRowIndex(armGenes.get(p)), 0, weight);

				}

				pb.step();

			});
		}

		return weights;

	}

	private static Map<String, ArrayList<String>> createChrArmGeneMapping(List<Gene> genes, LinkedHashMap<String, Integer> hashRows) {
		Map<String, ArrayList<String>> chrArmToGeneMapping = new HashMap<>(25);
		for (Gene gene : genes) {

			if (hashRows.containsKey(gene.getGene())) {

				String chrArm = gene.getChrAndArm();

				ArrayList<String> armGenes = chrArmToGeneMapping.get(chrArm);
				if (armGenes == null) {
					armGenes = new ArrayList<>();
					chrArmToGeneMapping.put(chrArm, armGenes);
				}

				armGenes.add(gene.getGene());

			}

		}
		return chrArmToGeneMapping;
	}

}
