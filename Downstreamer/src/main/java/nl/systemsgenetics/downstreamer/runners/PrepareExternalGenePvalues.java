/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.runners;

import cern.colt.GenericPermuting;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class PrepareExternalGenePvalues {

	private static final Logger LOGGER = Logger.getLogger(DownstreamerUtilities.class);

	public static void prepare(DownstreamerOptionsDeprecated options) throws Exception {

		LOGGER.info("This function will convert gene p-values from forinstance"
				+ " a rare variant burden test to files that can be used by step2"
				+ " of Downstreamer. The input file must be a matrix with ensg ids on rows"
				+ " and trait name as single column.");

		DoubleMatrixDataset<String, String> genePvalues = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
		
		List<Gene> genes = IoUtils.readGenes(options.getGeneInfoFile());
        Set<String> geneIds = new HashSet<>();
        for (Gene gene : genes) {
			if(genePvalues.containsRow(gene.getGene())){
				geneIds.add(gene.getGene());
			}
        }
		
		genePvalues = genePvalues.viewRowSelection(geneIds);

		final int numberOfPermutations = options.getPermutationGeneCorrelations() + options.getPermutationPathwayEnrichment() + options.getPermutationFDR();
		final int numberOfGenes = genePvalues.rows();
		
		LinkedHashMap<String, Integer> phenoHash = new LinkedHashMap<>();

		for (int i = 0; i < numberOfPermutations; ++i) {
			phenoHash.put("RanPheno" + (i + 1), i);
		}

		final DoubleMatrixDataset<String, String> genePvaluesNullGwas = new DoubleMatrixDataset<>(genePvalues.getHashRows(), phenoHash);
		DoubleMatrix2D genePvaluesNullGwasMatrix = genePvaluesNullGwas.getMatrix();
		
		System.out.println("Number of genes loaded: " + numberOfGenes);
		
		DoubleMatrix1D genePvaluesVector = genePvalues.getCol(0);
		for (int i = 0; i < numberOfPermutations; ++i) {
			
			int[] shuffle = GenericPermuting.permutation(i+2, numberOfGenes);
			DoubleMatrix1D copyGenePvalues = genePvaluesVector.copy();

			genePvaluesNullGwasMatrix.viewColumn(i).assign(copyGenePvalues.viewSelection(shuffle));
		}
		
		DoubleMatrixDataset<String, String> geneVariantCount = genePvalues.duplicate();
		geneVariantCount.getMatrix().assign(1);
		
		DoubleMatrixDataset<String, String> geneMaxSnpZscore = genePvalues.duplicate();
		geneMaxSnpZscore.getMatrix().assign(1);
		
		DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas = genePvaluesNullGwas.duplicate();
		geneMaxSnpZscoreNullGwas.getMatrix().assign(1);
		
		genePvalues.saveBinary(options.getOutputBasePath() + "_genePvalues");
		genePvaluesNullGwas.saveBinary(options.getOutputBasePath() + "_genePvaluesNullGwas");
		geneVariantCount.save(options.getOutputBasePath() + "_geneVariantCount.txt.gz");
		geneMaxSnpZscore.saveBinary(options.getOutputBasePath() + "_geneMaxSnpScores");
		geneMaxSnpZscoreNullGwas.saveBinary(options.getOutputBasePath() + "_geneMaxSnpZscoresNullGwas");
	}

}
