/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
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

					File annotationFile = new File(pathwayDatabase.getLocation() + ".colAnnotations.txt");
					int maxNumberOfAnnotations = 0;
					final HashMap<String, ArrayList<String>> annotations;

					if (annotationFile.exists()) {

						annotations = new HashMap<>();
						final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
						final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(annotationFile))).withCSVParser(parser).withSkipLines(0).build();

						String[] nextLine;
						while ((nextLine = reader.readNext()) != null) {

							ArrayList<String> thisAnnotions = new ArrayList();

							for (int i = 1; i < nextLine.length; ++i) {
								thisAnnotions.add(nextLine[i]);
							}

							annotations.put(nextLine[0], thisAnnotions);

							if (thisAnnotions.size() > maxNumberOfAnnotations) {
								maxNumberOfAnnotations = thisAnnotions.size();
							}

						}

						LOGGER.debug("Found annotations for: " + pathwayDatabase.getName() + ". Max annotations is: " + maxNumberOfAnnotations);

					} else {
						LOGGER.debug("Did not found annotations at: " + annotationFile.getAbsolutePath());
						annotations = null;
					}

					final CSVWriter enrichmentWriter = new CSVWriter(new FileWriter(new File(outputBasePath + "_" + pathwayDatabase.getName() + "Enrichment.txt")), '\t', '\0', '\0', "\n");
					final String[] outputLine = new String[1 + maxNumberOfAnnotations + enrichment.columns()];
					int c = 0;
					outputLine[c++] = "Pathway";
					for (int i = 0; i < maxNumberOfAnnotations; ++i) {
						outputLine[c++] = "PathwayAnnotation" + (i + 1);
					}
					for (String col : enrichment.getHashCols().keySet()) {
						outputLine[c++] = col;
					}
					enrichmentWriter.writeNext(outputLine);

					for (String pathwayKey : enrichment.getHashRows().keySet()) {
						c = 0;
						outputLine[c++] = pathwayKey;

						if (annotations != null) {
							ArrayList<String> thisPathwayAnnotations = annotations.get(pathwayKey);
							if (thisPathwayAnnotations == null) {
								for (int i = 0; i < maxNumberOfAnnotations; ++i) {
									outputLine[c++] = "";
								}
							} else {
								for (int i = 0; i < maxNumberOfAnnotations; ++i) {
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
				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}

				pb.step();

			});

		}

	}

}
