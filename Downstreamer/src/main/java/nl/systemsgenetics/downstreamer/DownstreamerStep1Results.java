/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer;

import java.io.File;
import java.io.IOException;
import java.util.List;
import nl.systemsgenetics.downstreamer.gene.Gene;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class DownstreamerStep1Results {

	private final DoubleMatrixDataset<String, String> genePvalues;
	private final DoubleMatrixDataset<String, String> genePvaluesNullGwas;
	private final DoubleMatrixDataset<String, String> geneVariantCount;
	private final DoubleMatrixDataset<String, String> geneMaxSnpZscore;
	private final DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas;

	private static final Logger LOGGER = Logger.getLogger(DownstreamerStep1Results.class);

	public DownstreamerStep1Results(DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> genePvaluesNullGwas, DoubleMatrixDataset<String, String> geneVariantCount, DoubleMatrixDataset<String, String> geneMaxSnpZscore, DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas) {
		this.genePvalues = genePvalues;
		this.genePvaluesNullGwas = genePvaluesNullGwas;
		this.geneVariantCount = geneVariantCount;
		this.geneMaxSnpZscore = geneMaxSnpZscore;
		this.geneMaxSnpZscoreNullGwas = geneMaxSnpZscoreNullGwas;
	}

	public static DownstreamerStep1Results loadFromDisk(String run1BasePath) throws IOException, Exception {
		if (!new File(run1BasePath + "_genePvalues.dat").exists()) {
			throw new RuntimeException("Could not find gene pvalues at: " + run1BasePath + "_genePvalues.dat");
		}
		DoubleMatrixDataset<String, String> genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(run1BasePath + "_genePvalues");
		DoubleMatrixDataset<String, String> genePvaluesNullGwas = DoubleMatrixDataset.loadDoubleBinaryData(run1BasePath + "_genePvaluesNullGwas");

		// Always load to avoid nullpointers
		DoubleMatrixDataset<String, String> geneMaxSnpZscore = DoubleMatrixDataset.loadDoubleBinaryData(run1BasePath + "_geneMaxSnpScores");
		DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas = DoubleMatrixDataset.loadDoubleBinaryData(run1BasePath + "_geneMaxSnpZscoresNullGwas");

		DoubleMatrixDataset<String, String> geneVariantCount = DoubleMatrixDataset.loadDoubleTextData(run1BasePath + "_geneVariantCount.txt", '\t');

		return new DownstreamerStep1Results(genePvalues, genePvaluesNullGwas, geneVariantCount, geneMaxSnpZscore, geneMaxSnpZscoreNullGwas);
	}

	public DoubleMatrixDataset<String, String> getGenePvalues() {
		return genePvalues;
	}

	public DoubleMatrixDataset<String, String> getGenePvaluesNullGwas() {
		return genePvaluesNullGwas;
	}

	public DoubleMatrixDataset<String, String> getGeneVariantCount() {
		return geneVariantCount;
	}

	public DoubleMatrixDataset<String, String> getGeneMaxSnpZscore() {
		return geneMaxSnpZscore;
	}

	public DoubleMatrixDataset<String, String> getGeneMaxSnpZscoreNullGwas() {
		return geneMaxSnpZscoreNullGwas;
	}

	
	
}
