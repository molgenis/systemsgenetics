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
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import org.apache.log4j.Logger;
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

	private static final Logger LOGGER = Logger.getLogger(PredictedPathwayAnnotations.class);
	
	public static void expandAnnotations(DownstreamerOptions options) throws Exception {

		List<PathwayDatabase> pathwayAnnotations = options.getPathwayDatabases();
		List<PathwayDatabase> gnPredictions = options.getPathwayDatabases2();
		HashMap<String, PathwayDatabase> gnPredictionsMap = new HashMap<>();
		for(PathwayDatabase gnp : gnPredictions){
			gnPredictionsMap.put(gnp.getName(), gnp);
		}
		
		LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
		LOGGER.info("Loaded " + genes.size() + " genes");
		
		for(PathwayDatabase pd : pathwayAnnotations){
			final PathwayDatabase gnp = gnPredictionsMap.get(pd.getName());
			if(gnp == null){
				throw new Exception("No predictions specified for: " + pd.getName());
			}
			
			LOGGER.info(pd.getName());
			
			DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pd.getLocation());
			final ArrayList<String> overlappingGenes = new ArrayList<>(pathwayMatrixLoader.getOriginalRowMap().keySet());
			overlappingGenes.retainAll(genes.keySet());
			
			
			//First load all genes in gene file so that they will be in output file
			final DoubleMatrixDataset<String, String> pathwayMatrix = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(overlappingGenes);
			
			LOGGER.info("Term-gene annotation before expand: " + pathwayMatrix.getMatrix().cardinality());
			
			DoubleMatrixDatasetFastSubsetLoader predictionMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(gnp.getLocation());
			
			overlappingGenes.retainAll(predictionMatrixLoader.getOriginalRowMap().keySet());
			
			//subset pathways to overlap with predictions
			DoubleMatrixDataset<String, String> pathwayMatrix2 = pathwayMatrix.viewRowSelection(overlappingGenes);
			final DoubleMatrixDataset<String, String> predictionMatrix = predictionMatrixLoader.loadSubsetOfRowsBinaryDoubleData(overlappingGenes);
			
			LOGGER.info("Pathway: " + pathwayMatrix.rows() + " x " + pathwayMatrix.columns());
			LOGGER.info("Predictions: " + predictionMatrix.rows() + " x " + predictionMatrix.columns());
			
			
			final double bonfThresholdP = 0.05d / (double) overlappingGenes.size();
			final double bonfThresholdZ = ZScores.pToZTwoTailed(bonfThresholdP) * -1;
			LOGGER.info("Bonf threshold: " + bonfThresholdZ);
			
			
			final int geneCount = pathwayMatrix2.rows(); 
			final ArrayList<String> terms = pathwayMatrix.getColObjects();
				
		
			if(LOGGER.isDebugEnabled()){
				for(int c = 0 ; c < 10 ; ++c){
					LOGGER.debug(pathwayMatrix.getCol(c).cardinality());
				}
			}
			
			int termsOverlapping = 0;
			int annotationsUpdated = 0;
			
			final ArrayList<String> termsIncludedInOutput = new ArrayList<>();
			
			for(String term : terms){
				
				final DoubleMatrix1D pathwayMatrixCol = pathwayMatrix2.getCol(term);
				
				if(pathwayMatrixCol.cardinality() < 10){
					continue;
				}
				
				termsIncludedInOutput.add(term);
				
				if(predictionMatrix.containsCol(term)){
					
					termsOverlapping++;
					
					final DoubleMatrix1D predictionMatrixCol = predictionMatrix.getCol(term);
					
					//due to loading order of genes must be identical here
					for(int r = 0 ; r < geneCount ; ++r){
						if(pathwayMatrixCol.get(r) == 0 && predictionMatrixCol.get(r) >= bonfThresholdZ){
							pathwayMatrixCol.set(r, 1);
							annotationsUpdated++;
						}
					}
					
				}
	
			}
			
			LOGGER.info("Terms overlapping with prediction matrix: " + termsOverlapping);
			LOGGER.info("Terms included in output: " + termsIncludedInOutput.size());
			LOGGER.info("Annotation updated: " + annotationsUpdated);
			
			if(LOGGER.isDebugEnabled()){
				for(int c = 0 ; c < 10 ; ++c){
					LOGGER.debug(pathwayMatrix.getCol(c).cardinality());
				}
			}

			pathwayMatrix.viewColSelection(termsIncludedInOutput).saveBinary(options.getOutputBasePath() + "_" + pd.getName());
			
		}

	}

}
