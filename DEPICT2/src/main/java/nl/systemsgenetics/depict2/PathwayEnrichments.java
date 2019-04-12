/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import static nl.systemsgenetics.depict2.Depict2.readMatrixAnnotations;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.WeightedCorrelations;

/**
 *
 * @author patri
 */
public class PathwayEnrichments {

	private static final Logger LOGGER = Logger.getLogger(Depict2Options.class);

	public static void doEnrichments(final DoubleMatrixDataset<String, String> genePvalues, final DoubleMatrixDataset<String, String> genePvaluesNullGwas, final DoubleMatrixDataset<String, String> geneWeights, final List<PathwayDatabase> pathwayDatabases, final String outputBasePath) {

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
					final DoubleMatrixDataset<String, String> genePvaluesNullGwasSubset = genePvaluesNullGwas.viewRowSelection(genesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> geneWeightsSubset = geneWeights.viewRowSelection(genesInPathwayMatrix);

					final DoubleMatrixDataset<String, String> pathwayMatrix = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(pathwayDatabase.getLocation(), genesInPathwayMatrix);

					//All matrices should now contain the same genes in the same order
					//Do enrichment analysis using weighted correlations
					DoubleMatrixDataset<String, String> enrichment = WeightedCorrelations.weightedCorrelationColumnsOf2Datasets(genePvaluesSubset, pathwayMatrix, geneWeightsSubset);
					DoubleMatrixDataset<String, String> enrichmentNull = WeightedCorrelations.weightedCorrelationColumnsOf2Datasets(genePvaluesNullGwasSubset, pathwayMatrix, geneWeightsSubset);

					
					PathwayAnnotations pathwayAnnotations = new PathwayAnnotations(new File(pathwayDatabase.getLocation() + ".colAnnotations.txt"));
					
					writeEnrichment(pathwayDatabase, pathwayAnnotations, enrichment, outputBasePath, "_correlations");
					writeEnrichment(pathwayDatabase, pathwayAnnotations, enrichmentNull, outputBasePath, "_correlations_null");
					
					DoubleMatrix2D enrichmentMatrix = enrichment.getMatrix();
					DoubleMatrix2D enrichmentNullMatrix = enrichmentNull.getMatrix();
					
					final int numberOfPathways = enrichmentMatrix.rows();
					final int numberOfPhenotypes = enrichmentMatrix.columns();
					final int numberOfNullGwasPhenotypes = enrichmentNullMatrix.columns();
					final double numberOfNullGwasPhenotypesPlus1Double = (double) enrichmentNullMatrix.columns() + 1;
					
					for(int r = 0 ; r < numberOfPathways ; ++r){
						
						for(int c = 0 ; c < numberOfPhenotypes ; ++c){
							
							final double corr = Math.abs(enrichmentMatrix.getQuick(r, c));
							
							int x = 0;
							for(int p = 0 ; p < numberOfNullGwasPhenotypes ; ++p){
								if(corr < Math.abs(enrichmentNullMatrix.getQuick(r, p))){
									++x;
								}
							}
							
							enrichmentMatrix.setQuick(r, c, (x/numberOfNullGwasPhenotypesPlus1Double));
							
						}					
						
					}				
					
					writeEnrichment(pathwayDatabase, pathwayAnnotations, enrichment, outputBasePath, "_pvalues");
					
					
				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}

				pb.step();

			});

		}

	}

	private static void writeEnrichment(PathwayDatabase pathwayDatabase, PathwayAnnotations pathwayAnnotations, DoubleMatrixDataset<String, String> enrichment, final String outputBasePath, final String nameSuffix) throws IOException {
		
		final CSVWriter enrichmentWriter = new CSVWriter(new FileWriter(new File(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + nameSuffix + ".txt")), '\t', '\0', '\0', "\n");
		final String[] outputLine = new String[1 + pathwayAnnotations.getMaxNumberOfAnnotations() + enrichment.columns()];
		int c = 0;
		outputLine[c++] = "Pathway";
		for (int i = 0; i < pathwayAnnotations.getMaxNumberOfAnnotations(); ++i) {
			outputLine[c++] = "PathwayAnnotation" + (i + 1);
		}
		for (String col : enrichment.getHashCols().keySet()) {
			outputLine[c++] = col;
		}
		enrichmentWriter.writeNext(outputLine);
		
		for (String pathwayKey : enrichment.getHashRows().keySet()) {
			c = 0;
			outputLine[c++] = pathwayKey;
			
			if (pathwayAnnotations.getMaxNumberOfAnnotations() > 0) {
				ArrayList<String> thisPathwayAnnotations = pathwayAnnotations.getAnnotationsForPathway(pathwayKey);
				if (thisPathwayAnnotations == null) {
					for (int i = 0; i < pathwayAnnotations.getMaxNumberOfAnnotations(); ++i) {
						outputLine[c++] = "";
					}
				} else {
					for (int i = 0; i < pathwayAnnotations.getMaxNumberOfAnnotations(); ++i) {
						outputLine[c++] = i < thisPathwayAnnotations.size() ? thisPathwayAnnotations.get(i) : "";
					}
				}
			}
			
			DoubleMatrix1D row = enrichment.viewRow(pathwayKey);
			for (int e = 0; e < row.size(); ++e) {
				outputLine[c++] = String.valueOf(row.getQuick(e));
			}
			
			enrichmentWriter.writeNext(outputLine);
		}
		enrichmentWriter.close();
		
		//enrichment.save(outputBasePath + "_" + pathwayDatabase.getName() + "Enrichment.txt");
	}

}
