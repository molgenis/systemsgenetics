/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.gene;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class CalculateGeneInvCorMatrix {

	private static final Logger LOGGER = Logger.getLogger(Downstreamer.class);

	/**
	 *
	 * @param genePvaluesNullGwas
	 * @param genes
	 * @param options
	 * @return inv cor matrix per chr arm
	 */
	public static Map<String, DoubleMatrixDataset<String, String>> CalculateGeneInvCorMatrix(final DoubleMatrixDataset<String, String> genePvaluesNullGwas, List<Gene> genes, DownstreamerOptions options) {

		final Map<String, ArrayList<String>> chrArmToGeneMapping = createChrArmGeneMapping(genes, genePvaluesNullGwas.getHashRows());

		final Map<String, DoubleMatrixDataset<String, String>> invCorMatrixPerChrArm = Collections.synchronizedMap(new HashMap<>(chrArmToGeneMapping.size()));

		//LinkedHashMap<String, Integer> geneHash = genePvaluesNullGwas.getHashRowsCopy();
		//final DoubleMatrixDataset<String, String> invCorMatrix = new DoubleMatrixDataset<>(geneHash, geneHash);
		try (ProgressBar pb = new ProgressBar("Gene inv cor calculation", chrArmToGeneMapping.size(), ProgressBarStyle.ASCII)) {

			chrArmToGeneMapping.keySet().parallelStream().forEach((String chrArm) -> {

				final ArrayList<String> armGenes = chrArmToGeneMapping.get(chrArm);

				final DoubleMatrixDataset<String, String> genePvaluesNullGwasArm = genePvaluesNullGwas.viewRowSelection(armGenes);

				//final DoubleMatrixDataset<String, String> invCorMatrixArmGenes = invCorMatrix.viewSelection(armGenes, armGenes);
				final DoubleMatrixDataset<String, String> genePvaluesNullGwasArmT = genePvaluesNullGwasArm.viewDice();

				final DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelation = genePvaluesNullGwasArmT.calculateCorrelationMatrix();

				if (LOGGER.isDebugEnabled()) {
					try {
						genePvaluesNullGwasGeneArmCorrelation.save(new File(options.getOutputBasePath() + "_" + chrArm + "_GeneCorMatrix.txt"));
					} catch (IOException ex) {
						throw new RuntimeException(ex);
					}
				}

				//We need to take the inverse of the correlation matrix. To do that the correlation between genes can't be correlated
				//Simply removing highly correlated genes did not always work, therefor:
				//(1) create correlation matrix of correlations
				//(2) identifie genes that have correlated correlation
				//(3) prune gene correlation matrix
				DoubleMatrixDataset<String, String> correlationOfCorrelations = genePvaluesNullGwasGeneArmCorrelation.calculateCorrelationMatrix();

				ArrayList<String> variantNames = correlationOfCorrelations.getRowObjects();
				LinkedHashSet<String> includedGenes = new LinkedHashSet<>(correlationOfCorrelations.rows());

				rows:
				for (int r = 0; r < correlationOfCorrelations.rows(); ++r) {
					cols:
					for (int c = 0; c < r; ++c) {
						if (Math.abs(correlationOfCorrelations.getElementQuick(r, c)) >= 0.95 && includedGenes.contains(variantNames.get(c))) {
							continue rows;
						}
					}
					includedGenes.add(variantNames.get(r));
				}

				final DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelationPruned = genePvaluesNullGwasGeneArmCorrelation.viewSelection(includedGenes, includedGenes);

				final DoubleMatrix2D genePvaluesNullGwasGeneArmCorrelationInverseMatrix;
				try {
					genePvaluesNullGwasGeneArmCorrelationInverseMatrix = new DenseDoubleAlgebra().inverse(genePvaluesNullGwasGeneArmCorrelationPruned.getMatrix());
				} catch (RuntimeException ex) {

					LOGGER.fatal("Error during matrix inverse of: " + chrArm);
					LOGGER.fatal("Marix before prune: " + genePvaluesNullGwasGeneArmCorrelation.getMatrix().toStringShort());
					LOGGER.fatal("Matrix info: " + genePvaluesNullGwasGeneArmCorrelationPruned.getMatrix().toStringShort());

//					try {
//						genePvaluesNullGwasGeneArmCorrelationPruned.save(options.getOutputBasePath() + "singularMatrix.txt");
//					} catch (IOException ex1) {
//						throw ex;
//					}
					throw ex;
				}
				final DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelationInverse = new DoubleMatrixDataset(genePvaluesNullGwasGeneArmCorrelationInverseMatrix, genePvaluesNullGwasGeneArmCorrelationPruned.getHashRows(), genePvaluesNullGwasGeneArmCorrelationPruned.getHashCols());

				try {
					genePvaluesNullGwasGeneArmCorrelationInverse.save(options.getOutputBasePath() + "_geneInvCor_" + chrArm + ".txt");
				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}

				invCorMatrixPerChrArm.put(chrArm, genePvaluesNullGwasGeneArmCorrelationInverse);

				pb.step();

			});
		}

		return invCorMatrixPerChrArm;

	}

	protected static Map<String, ArrayList<String>> createChrArmGeneMapping(List<Gene> genes, LinkedHashMap<String, Integer> hashRows) {
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
