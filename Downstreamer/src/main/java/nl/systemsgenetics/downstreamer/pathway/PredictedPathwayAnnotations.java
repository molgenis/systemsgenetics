/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * Used to expand 0/1 pathway matrix with gene annotations to include strongly
 * predicted genes
 *
 * @author patri
 */
public class PredictedPathwayAnnotations {

	private static final Logger LOGGER = LogManager.getLogger(PredictedPathwayAnnotations.class);

	public static void expandAnnotations(DownstreamerOptionsDeprecated options) throws Exception {

		List<PathwayDatabase> pathwayAnnotations = options.getPathwayDatabases();
		List<PathwayDatabase> gnPredictions = options.getPathwayDatabases2();
		HashMap<String, PathwayDatabase> gnPredictionsMap = new HashMap<>();
		for (PathwayDatabase gnp : gnPredictions) {
			gnPredictionsMap.put(gnp.getName(), gnp);
		}

		LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
		LOGGER.info("Loaded " + genes.size() + " genes");

		for (PathwayDatabase pd : pathwayAnnotations) {
			final PathwayDatabase gnp = gnPredictionsMap.get(pd.getName());
			if (gnp == null) {
				throw new Exception("No predictions specified for: " + pd.getName());
			}

			LOGGER.info(pd.getName());

			DoubleMatrixDatasetFastSubsetLoader predictionMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(gnp.getLocation());

			final boolean fromScratch = pd.getLocation().equals("<NA>");
			
			final DoubleMatrixDataset<String, String> pathwayMatrix;
			final DoubleMatrixDataset<String, String> pathwayMatrix2;
			final ArrayList<String> overlappingGenes;
			if (fromScratch) {
				
				overlappingGenes = new ArrayList<>(predictionMatrixLoader.getOriginalRowMap());
				overlappingGenes.retainAll(genes.keySet());
				
				pathwayMatrix = new DoubleMatrixDataset<>(overlappingGenes, predictionMatrixLoader.getOriginalColMap());
				pathwayMatrix2 = pathwayMatrix;
				
			} else {
				DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pd.getLocation());
				overlappingGenes = new ArrayList<>(pathwayMatrixLoader.getOriginalRowMap());
				overlappingGenes.retainAll(genes.keySet());
				//First load all genes in gene file so that they will be in output file
				pathwayMatrix = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(overlappingGenes);

				overlappingGenes.retainAll(predictionMatrixLoader.getOriginalRowMap());

				//subset pathways to overlap with predictions
				pathwayMatrix2 = pathwayMatrix.viewRowSelection(overlappingGenes);

			}

			LOGGER.info("Term-gene annotation before expand: " + pathwayMatrix.getMatrix().cardinality());

			LinkedHashMap<String, Integer> geneCountMatrixCols = new LinkedHashMap<>();
			geneCountMatrixCols.put("Raw", 0);
			geneCountMatrixCols.put("Expanded", 1);
			final DoubleMatrixDataset<String, String> geneCountMatrix = new DoubleMatrixDataset<>(pathwayMatrix.getHashColsCopy(), geneCountMatrixCols);

			final ArrayList<String> terms = pathwayMatrix.getColObjects();

			final DoubleMatrixDataset<String, String> predictionMatrix = predictionMatrixLoader.loadSubsetOfRowsBinaryDoubleData(overlappingGenes);

			LOGGER.info("Pathway: " + pathwayMatrix.rows() + " x " + pathwayMatrix.columns());
			LOGGER.info("Predictions: " + predictionMatrix.rows() + " x " + predictionMatrix.columns());

			final double bonfThresholdP = 0.05d / (double) overlappingGenes.size();
			final double bonfThresholdZ = ZScores.pToZTwoTailed(bonfThresholdP) * -1;
			LOGGER.info("Bonf threshold: " + bonfThresholdZ);

			final int geneCount = pathwayMatrix2.rows();

			if (LOGGER.isDebugEnabled()) {
				for (int c = 0; c < 10; ++c) {
					LOGGER.debug(pathwayMatrix.getCol(c).cardinality());
				}
			}

			int termsOverlapping = 0;
			int annotationsUpdated = 0;

			final ArrayList<String> termsIncludedInOutput = new ArrayList<>();

			for (String term : terms) {

				final DoubleMatrix1D pathwayMatrixCol = pathwayMatrix2.getCol(term);

				int count = pathwayMatrixCol.cardinality();

				geneCountMatrix.setElement(term, "Raw", count);

				if (!fromScratch && count < 10) {
					geneCountMatrix.setElement(term, "Expanded", count);
					continue;
				}

				termsIncludedInOutput.add(term);

				if (predictionMatrix.containsCol(term)) {

					termsOverlapping++;

					final DoubleMatrix1D predictionMatrixCol = predictionMatrix.getCol(term);

					//due to loading order of genes must be identical here
					for (int r = 0; r < geneCount; ++r) {
						if (pathwayMatrixCol.get(r) == 0 && predictionMatrixCol.get(r) >= bonfThresholdZ) {
							pathwayMatrixCol.set(r, 1);
							annotationsUpdated++;
							count++;
						}
					}

				}

				geneCountMatrix.setElement(term, "Expanded", count);

			}

			LOGGER.info("Terms overlapping with prediction matrix: " + termsOverlapping);
			LOGGER.info("Terms included in output: " + termsIncludedInOutput.size());
			LOGGER.info("Annotation updated: " + annotationsUpdated);

			if (LOGGER.isDebugEnabled()) {
				for (int c = 0; c < 10; ++c) {
					LOGGER.debug(pathwayMatrix.getCol(c).cardinality());
				}
			}

			pathwayMatrix.viewColSelection(termsIncludedInOutput).saveBinary(options.getOutputBasePath() + "_" + pd.getName());
			geneCountMatrix.save(options.getOutputBasePath() + "_" + pd.getName() + "_counts.txt.gz");

		}

	}

}
