/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import nl.systemsgenetics.downstreamer.runners.options.OptionsModeEnrichment;
import nl.systemsgenetics.downstreamer.summarystatistic.LinearRegressionResult;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class DownstreamerEnrichment {

	private static final Logger LOGGER = LogManager.getLogger(DownstreamerEnrichment.class);

	public static DownstreamerStep2Results enrichmentAnalysis(OptionsModeEnrichment options) throws Exception {

		final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();

		//quick check if all pathway databases exist before loading all the other data.
		for (PathwayDatabase pd : pathwayDatabases) {
			if (!pd.exist()) {
				throw new FileNotFoundException("Could not read: " + pd.getLocation() + ".dat or .datg");
			}
		}

		// Load the genes to run the analysis on
		LinkedHashMap<String, Gene> genes = IoUtils.readGenesMap(options.getGeneInfoFile());
		LOGGER.info("Loaded " + genes.size() + " genes");

		DoubleMatrixDataset<String, String> gwasGenePvalues;
		DoubleMatrixDataset<String, String> gwasGeneVarCount;
		DoubleMatrixDataset<String, String> gwasGeneMinVarPvalue;

		if (options.getSingleGwasFile() != null) {
			// Load GWAS data Cols: genes, pvalue, nSNPs, min SNP p-value
			DoubleMatrixDataset<String, String> gwasData = DoubleMatrixDataset.loadDoubleData(options.getSingleGwasFile().getAbsolutePath());

			ArrayList<String> colNames = gwasData.getColObjects();

			gwasGenePvalues = gwasData.viewColSelection(colNames.get(0));
			gwasGeneVarCount = gwasData.viewColSelection(colNames.get(1));
			gwasGeneMinVarPvalue = gwasData.viewColSelection(colNames.get(2));

		} else {

			gwasGenePvalues = DoubleMatrixDataset.loadDoubleData(options.getGwasPvalueMatrixPath() + "_pvalues");
			gwasGeneVarCount = DoubleMatrixDataset.loadDoubleData(options.getGwasPvalueMatrixPath() + "_nvar");
			gwasGeneMinVarPvalue = DoubleMatrixDataset.loadDoubleData(options.getGwasPvalueMatrixPath() + "_minVarPvalue");

		}

		if (!options.isSkipPvalueToZscore()) {
			inplacePvalueToZscore(gwasGenePvalues);
		}

		final Set<String> hlaGenes;
		if (options.isExcludeHla()) {
			hlaGenes = new HashSet<>();
			for (Gene gene : genes.values()) {
				if (gene.overlaps(options.getHla())) {
					hlaGenes.add(gene.getGene());
				}
			}
			LOGGER.info("Excluding " + hlaGenes.size() + " genes");

		} else {
			hlaGenes = Collections.EMPTY_SET;
		}

		ArrayList<String> allGwasGenes = gwasGenePvalues.getRowObjects();

		LinkedHashSet<String> selectedGenes = new LinkedHashSet<>();

		//This loop will create a list of genes that should be selected
		genes:
		for (int g = 0; g < gwasGenePvalues.rows(); ++g) {

			gwasGenePvalues.setElementQuick(g, 0, -ZScores.pToZTwoTailed(gwasGenePvalues.getElementQuick(g, 0)));

			String gene = allGwasGenes.get(g);

			if (hlaGenes.contains(gene)) {
				//is hla gene so don't add to selected
				continue;
			}

			if (gwasGeneVarCount.getElementQuick(g, 0) < 1) {
				//if no SNPs are found around gene don't use
				//Assuses this is equal for alle GWASes
				continue;
			}

			if (!genes.containsKey(gene)) {
				//if gene is not in list of genes to use, skip
				continue;
			}

			for (int c = 0; c < gwasGenePvalues.columns(); ++g) {
				if (Double.isNaN(gwasGenePvalues.getElementQuick(g, c))) {
					//gene did not have a valid p-value
					continue genes;
				}
			}

			selectedGenes.add(gene);

		}

		// Load optinal covariates
		final DoubleMatrixDataset<String, String> covariatesToCorrectGenePvalues;
		if (options.getCovariates() != null) {
			covariatesToCorrectGenePvalues = DoubleMatrixDataset.loadDoubleData(options.getCovariates().getAbsolutePath());

			if (!covariatesToCorrectGenePvalues.getHashRows().keySet().containsAll(selectedGenes)) {
				throw new Exception("Not all genes are found in the covariate file");
			}

		} else {
			covariatesToCorrectGenePvalues = null;
		}

		final Map<String, List<Gene>> chrArmGeneMap = IoUtils.readGenesAsChrArmMap(options.getGeneInfoFile());

		final ArrayList<PathwayEnrichments> pathwayEnrichments = new ArrayList<>(pathwayDatabases.size());
		
		for (PathwayDatabase pathwayDatabase : pathwayDatabases) {

			final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase.getLocation());
			ArrayList<String> genesOverlappingWithPathwayDatabase = new ArrayList<>(selectedGenes.size());

			for (String pathwayGene : pathwayMatrixLoader.getAllRowIdentifiers()) {
				if (selectedGenes.contains(pathwayGene)) {
					genesOverlappingWithPathwayDatabase.add(pathwayGene);
				}
			}

			DoubleMatrixDataset<String, String> gwasGenePvaluesSubset;
			final DoubleMatrixDataset<String, String> covariatesToCorrectGenePvaluesSubset = covariatesToCorrectGenePvalues == null ? null : covariatesToCorrectGenePvalues.viewRowSelection(genesOverlappingWithPathwayDatabase);
					
			if(options.isForceNormalGenePvalues()){
				gwasGenePvaluesSubset = gwasGenePvalues.viewRowSelection(genesOverlappingWithPathwayDatabase).createColumnForceNormalDuplicate();
			} else {
				gwasGenePvaluesSubset = gwasGenePvalues.viewRowSelection(genesOverlappingWithPathwayDatabase);
			}
					
			final DoubleMatrixDataset<String, String> pathwayData = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(genesOverlappingWithPathwayDatabase);

			
			
			final List<int[]> blockDiagonalIndices = DownstreamerRegressionEngine.createBlockDiagonalIndexFromGenes(chrArmGeneMap, genesOverlappingWithPathwayDatabase);

			//TODO papamters
			List<LinearRegressionResult> pathwayRegeressionResults = DownstreamerRegressionEngine.performDownstreamerRegression(null, null, null, null, null, blockDiagonalIndices, true, true, true);

			final DoubleMatrixDataset<String, String> pathwayPvalues = new DoubleMatrixDataset<>(pathwayData.getColObjects(), gwasGenePvalues.getColObjects());
			
			
		}

		return null;

	}

	private static void inplacePvalueToZscore(DoubleMatrixDataset dataset) {
		final DoubleMatrix2D matrix = dataset.getMatrix();

		// Inplace convert gene p-values to z-scores
		IntStream.range(0, matrix.rows()).parallel().forEach(r -> {
			for (int c = 0; c < matrix.columns(); ++c) {
				matrix.setQuick(r, c, -ZScores.pToZTwoTailed(matrix.getQuick(r, c)));
			}
		});
	}

}
