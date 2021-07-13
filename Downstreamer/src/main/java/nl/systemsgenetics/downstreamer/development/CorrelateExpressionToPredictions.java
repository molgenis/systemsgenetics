/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;

/**
 *
 * @author patri
 */
public class CorrelateExpressionToPredictions {

		private static final Logger LOGGER = Logger.getLogger(CorrelateExpressionToPredictions.class);

	
	public static void run(DownstreamerOptions options) throws IOException, Exception{
		
		DoubleMatrixDatasetFastSubsetLoader corePredictionZscoresLoader = new DoubleMatrixDatasetFastSubsetLoader(options.getIntermediateFolder() + "/" + options.getX() + "_Enrichment_zscoreExHla");
		DoubleMatrixDataset<String, String> samplePredictionZscores = DoubleMatrixDataset.loadDoubleBinaryData(options.getIntermediateFolder() + "/expression_Enrichment_zscoreExHla");
		
		DoubleMatrixDatasetFastSubsetLoader expressionDataLoader = new DoubleMatrixDatasetFastSubsetLoader(options.getY());
		
		LinkedHashSet<String> sharedGenes = new LinkedHashSet<>(corePredictionZscoresLoader.getOriginalRowMap().keySet());
		sharedGenes.retainAll(expressionDataLoader.getOriginalRowMap().keySet());
		
		DoubleMatrixDataset<String, String> corePredictionZscores = corePredictionZscoresLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes);
		DoubleMatrixDataset<String, String> expressionData = expressionDataLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes);
		
		LOGGER.info("corePredictionZscores: " + corePredictionZscores.getMatrix().toStringShort());
		LOGGER.info("expressionData: " + expressionData.getMatrix().toStringShort());
		LOGGER.info("samplePredictionZscores: " + samplePredictionZscores.getMatrix().toStringShort());
		
		corePredictionZscores.normalizeColumns();
		expressionData.normalizeColumns();
		
		DoubleMatrixDataset<String, String> correlations = DoubleMatrixDataset.correlateColumnsOf2ColumnNormalizedDatasets(corePredictionZscores, expressionData);
		
		LOGGER.info("correlations: " + correlations.getMatrix().toStringShort());
		
		samplePredictionZscores = samplePredictionZscores.viewRowSelection(correlations.getHashRows().keySet());
		
		LinkedHashMap<String, Integer> outputColumns = new LinkedHashMap<String, Integer>(2);
		outputColumns.put("SampleEnrichmentZscore", 0);
		outputColumns.put("CorrelationSampleExpressionTo" + options.getX(), 1);
		
		DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<String, String> (correlations.getHashRows(), outputColumns);
		
		LOGGER.info("output: " + output.getMatrix().toStringShort());
		
		output.getCol(0).assign(samplePredictionZscores.getCol(0));
		output.getCol(1).assign(correlations.getCol(0));
		
		output.save(options.getOutputBasePath() + options.getX() + "vsTissueEnrichment.txt");
		
		
	}
	
}
