/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import static nl.systemsgenetics.depict2.Depict2.readMatrixAnnotations;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.WeightedCorrelations;

/**
 *
 * @author patri
 */
public class PathwayEnrichments {

	public static void doEnrichments(final DoubleMatrixDataset<String, String> genePvalues, final DoubleMatrixDataset<String, String> geneWeights, final List<PathwayDatabase> pathwayDatabases, final String outputBasePath) {

		try (ProgressBar pb = new ProgressBar("Pathway enrichtment analysis", pathwayDatabases.size(), ProgressBarStyle.ASCII)) {

			pathwayDatabases.parallelStream().forEach((PathwayDatabase pathwayDatabase) -> {

				try {

					final List<String> genesInPathwayMatrix = readMatrixAnnotations(new File(pathwayDatabase.getLocation() + ".rows.txt"));

					Iterator<String> pathwayGeneIterator = genesInPathwayMatrix.iterator();
					String pathwayGene;
					while (pathwayGeneIterator.hasNext()) {
						pathwayGene = pathwayGeneIterator.next();
						if (!genePvalues.containsRow(pathwayGene)) {
							pathwayGeneIterator.remove();
						}
					}
					//Now genesInPathwayMatrix will only contain genes that are also in the gene p-value matrix
					
					final DoubleMatrixDataset<String, String> genePvaluesSubset = genePvalues.viewRowSelection(genesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> geneWeightsSubset = geneWeights.viewRowSelection(genesInPathwayMatrix);

					final DoubleMatrixDataset<String, String> pathwayMatrix = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(pathwayDatabase.getLocation(), genesInPathwayMatrix);
					
					//All three matrices should now contain the same genes in the same order
					DoubleMatrixDataset<String, String> enrichment = WeightedCorrelations.weightedCorrelationColumnsOf2Datasets(genePvaluesSubset, pathwayMatrix, geneWeightsSubset);
					
					enrichment.save(outputBasePath + "_" + pathwayDatabase.getName() + "Enrichment.txt" );

				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}

				pb.step();

			});

		}

	}

}
