/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class CalculateGeneInvCorMatrix {

	private static final Logger LOGGER = Logger.getLogger(Depict2.class);

	/**
	 * 
	 * @param genePvaluesNullGwas
	 * @param genes
	 * @param options
	 * @return 
	 */
	public static DoubleMatrixDataset<String, String> CalculateGeneInvCorMatrix(final DoubleMatrixDataset<String, String> genePvaluesNullGwas, List<Gene> genes, Depict2Options options) {

		Map<String, ArrayList<String>> chrArmToGeneMapping = createChrArmGeneMapping(genes, genePvaluesNullGwas.getHashRows());

		LinkedHashMap<String, Integer> geneHash = genePvaluesNullGwas.getHashRowsCopy();

		final DoubleMatrixDataset<String, String> invCorMatrix = new DoubleMatrixDataset<>(geneHash, geneHash);

		try (ProgressBar pb = new ProgressBar("Gene inv cor calculation", chrArmToGeneMapping.size(), ProgressBarStyle.ASCII)) {

			chrArmToGeneMapping.keySet().parallelStream().forEach((String chrArm) -> {

				final ArrayList<String> armGenes = chrArmToGeneMapping.get(chrArm);

				final DoubleMatrixDataset<String, String> genePvaluesNullGwasArm = genePvaluesNullGwas.viewRowSelection(armGenes);
				
				final DoubleMatrixDataset<String, String> invCorMatrixArmGenes = invCorMatrix.viewSelection(armGenes, armGenes);
			
				final DoubleMatrixDataset<String, String> genePvaluesNullGwasArmT = genePvaluesNullGwasArm.viewDice();

				final DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelation = genePvaluesNullGwasArmT.calculateCorrelationMatrix();

				if (LOGGER.isDebugEnabled()) {
					try {
						genePvaluesNullGwasGeneArmCorrelation.save(new File(options.getOutputBasePath() + "_" + chrArm + "_GeneCorMatrix.txt"));
					} catch (IOException ex) {
						throw new RuntimeException(ex);
					}
				}

				DoubleMatrix2D genePvaluesNullGwasGeneArmCorrelationInverse = new DenseDoubleAlgebra().inverse(genePvaluesNullGwasGeneArmCorrelation.getMatrix());
				
				invCorMatrixArmGenes.getMatrix().assign(genePvaluesNullGwasGeneArmCorrelationInverse);
				
				pb.step();

			});
		}

		return invCorMatrix;

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
