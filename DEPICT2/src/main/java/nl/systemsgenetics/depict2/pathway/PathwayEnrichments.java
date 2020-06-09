/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2.pathway;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.jet.math.tdouble.DoubleFunctions;
import com.opencsv.CSVWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.stream.IntStream;

import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.depict2.Depict2;
import nl.systemsgenetics.depict2.gene.Gene;
import nl.systemsgenetics.depict2.gene.GenePathwayAssociationStatistic;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

/**
 * @author patri
 */
public class PathwayEnrichments {

	private static final Logger LOGGER = Logger.getLogger(Depict2.class);

	private final PathwayDatabase pathwayDatabase;
	private final HashSet<String> hlaGenesToExclude;
	private final boolean ignoreGeneCorrelations;
	private final DoubleMatrixDataset<String, String> betas;
	//private final DoubleMatrixDataset<String, String> standardErrors;
	//private final DoubleMatrixDataset<String, String> tStatistics;
	//private final DoubleMatrixDataset<String, String> pValues;
	private final DoubleMatrixDataset<String, String> qValues;
	//private final DoubleMatrixDataset<String, String> pValuesNull;
	private final DoubleMatrixDataset<String, String> zscores;
	//private final DoubleMatrixDataset<String, String> betasNull;
	private final int numberOfPathways;
	private final File intermediateFolder;

	// TODO: can be locals??
	private final Set<String> excludeGenes;
	private final List<Gene> genes;
	private final String outputBasePath;

	public PathwayEnrichments(final PathwayDatabase pathwayDatabase,
			final HashSet<String> genesWithPvalue,
			final List<Gene> genes,
			final boolean forceNormalPathwayPvalues,
			final boolean forceNormalGenePvalues,
			DoubleMatrixDataset<String, String> geneZscores,
			DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelation,
			DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetas,
			final String outputBasePath,
			final HashSet<String> hlaGenesToExclude,
			final boolean ignoreGeneCorrelations,
			final double genePruningR,
			final File debugFolder,
			final File intermediateFolder,
			final boolean quantileNormalizePermutations,
			final boolean regressGeneLengths,
			DoubleMatrixDataset<String, String> geneMaxSnpZscore,
			DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwasCorrelation,
			DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwasBetas,
			int geneCorrelationWindow
	) throws Exception {

		this.pathwayDatabase = pathwayDatabase;
		this.genes = genes;
		this.outputBasePath = outputBasePath;
		this.hlaGenesToExclude = hlaGenesToExclude;
		this.ignoreGeneCorrelations = ignoreGeneCorrelations;
		this.intermediateFolder = intermediateFolder;
		final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase.getLocation());
		final LinkedHashSet<String> sharedGenes;

		if (this.hlaGenesToExclude == null) {
			excludeGenes = Collections.emptySet();
		} else {
			excludeGenes = this.hlaGenesToExclude;
		}

		// Determine final set of genes to analyze and overlap with genes in pathway matrix
		Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();
		sharedGenes = new LinkedHashSet<>();

		for (String gene : genesWithPvalue) {
			if (pathwayGenes.contains(gene) && !excludeGenes.contains(gene)) {
				sharedGenes.add(gene);
			}
		}

		// Subset the gene zscores to all overlapping, non NA genes
		geneZscores = geneZscores.viewRowSelection(sharedGenes);
		geneZscoresNullGwasCorrelation = geneZscoresNullGwasCorrelation.viewRowSelection(sharedGenes);
		geneZscoresNullGwasNullBetas = geneZscoresNullGwasNullBetas.viewRowSelection(sharedGenes);

		// Force normal gene p-values y/n
		if (forceNormalGenePvalues) {
			LOGGER.info("Force normalizing gene p-values / z-scores");
			geneZscores = createColumnForceNormalDuplicate(geneZscores, geneMaxSnpZscore);
			geneZscoresNullGwasCorrelation = createColumnForceNormalDuplicate(geneZscoresNullGwasCorrelation, geneMaxSnpZscoreNullGwasCorrelation);
			geneZscoresNullGwasNullBetas = createColumnForceNormalDuplicate(geneZscoresNullGwasNullBetas, geneMaxSnpZscoreNullGwasBetas);
		} else {
			// Do this because the regression is in place.
			geneZscores = geneZscores.duplicate();
			geneZscoresNullGwasCorrelation = geneZscoresNullGwasCorrelation.duplicate();
			geneZscoresNullGwasNullBetas = geneZscoresNullGwasNullBetas.duplicate();
		}

		// Do the regression with gene lengths and z-scores y/n
		if (regressGeneLengths) {
			LOGGER.info("Regressing gene lengths and determining residuals");

			// Determine (log10) gene lengths
			final double[] geneLengths = new double[sharedGenes.size()];
			int i = 0;
			// Ensures the geneLength vector is in the same order as the matrices
			for (String geneInZscoreMatrix : sharedGenes) {
				for (Gene geneInfo : genes) {
					if (geneInfo.getGene().equals(geneInZscoreMatrix)) {
						geneLengths[i] = Math.log(geneInfo.getLength()) / Math.log(10);
						i++;
						break;
					}
				}
			}

			// Determine residuals
			inplaceDetermineGeneLengthRegressionResiduals(geneZscores, geneLengths);
			inplaceDetermineGeneLengthRegressionResiduals(geneZscoresNullGwasCorrelation, geneLengths);
			inplaceDetermineGeneLengthRegressionResiduals(geneZscoresNullGwasNullBetas, geneLengths);
		}

		// Determine which genes will be merged to metagenes based on their genetic correlation
		final HashMap<String, ArrayList<MetaGene>> metaGenesPerArm;
		{
			DoubleMatrixDataset<String, String> tmp = geneZscoresNullGwasCorrelation.duplicate();
			tmp.normalizeColumns();
			metaGenesPerArm = groupCorrelatedGenesPerChrArm(tmp, genePruningR, genes, sharedGenes, debugFolder, pathwayDatabase.getName(), hlaGenesToExclude);
		}

		if (LOGGER.isDebugEnabled()) {
			final CSVWriter metaGeneWriter = new CSVWriter(new FileWriter(new File(debugFolder, pathwayDatabase.getName() + "_Enrichment_metaGenes.txt")), '\t', '\0', '\0', "\n");
			final String[] outputLine = new String[5];
			int c = 0;
			outputLine[c++] = "MetaGeneId";
			outputLine[c++] = "ChrArm";
			outputLine[c++] = "Number of genes";
			outputLine[c++] = "Start";
			outputLine[c++] = "Stop";
			metaGeneWriter.writeNext(outputLine);

			for (Map.Entry<String, ArrayList<MetaGene>> metaGenesEntry : metaGenesPerArm.entrySet()) {
				final String chrArm = metaGenesEntry.getKey();
				final ArrayList<MetaGene> metaGenes = metaGenesEntry.getValue();
				for (MetaGene metaGene : metaGenes) {
					c = 0;
					outputLine[c++] = metaGene.getMetaGeneId();
					outputLine[c++] = chrArm;
					outputLine[c++] = String.valueOf(metaGene.getGeneCount());
					outputLine[c++] = String.valueOf(metaGene.getStart());
					outputLine[c++] = String.valueOf(metaGene.getStop());
					metaGeneWriter.writeNext(outputLine);

				}
			}
			metaGeneWriter.close();
		}

		final int numberTraits = geneZscores.columns();
		final int numberTraitsNull = geneZscoresNullGwasNullBetas.columns();

		try (ProgressBar pb = new ProgressBar(pathwayDatabase.getName() + " enrichment analysis", metaGenesPerArm.size() + numberTraits + numberTraits + numberTraitsNull + 3, ProgressBarStyle.ASCII)) {

			// Collapse the gene p-values to these predefined metagenes
			final DoubleMatrixDataset<String, String> geneZscoresPathwayMatched;
			final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationPathwayMatched;
			final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetasPathwayMatched;

			geneZscoresPathwayMatched = collapseDatasetToMetaGenes(geneZscores, false, metaGenesPerArm.values());
			geneZscoresNullGwasCorrelationPathwayMatched = collapseDatasetToMetaGenes(geneZscoresNullGwasCorrelation, false, metaGenesPerArm.values());
			geneZscoresNullGwasNullBetasPathwayMatched = collapseDatasetToMetaGenes(geneZscoresNullGwasNullBetas, false, metaGenesPerArm.values());

			// Collapse the pathways to meta-genes
			final DoubleMatrixDataset<String, String> genePathwayZscores;
			if (forceNormalPathwayPvalues) {
				LOGGER.debug("Doing force normal pathway scores");
				genePathwayZscores = collapseDatasetToMetaGenes(pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes),
						true,
						metaGenesPerArm.values()).createColumnForceNormalDuplicate();
			} else {
				genePathwayZscores = collapseDatasetToMetaGenes(pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedGenes),
						true,
						metaGenesPerArm.values());
				LOGGER.debug("Center and scale pathway scores");
			}

			// Center and scale pathway scores and gene z-scores so the GLS implementation is valid
			genePathwayZscores.normalizeColumns();
			geneZscoresPathwayMatched.normalizeColumns();
			geneZscoresNullGwasCorrelationPathwayMatched.normalizeColumns();
			geneZscoresNullGwasNullBetasPathwayMatched.normalizeColumns();

			// Save normalized gene scores to file, to check distribution later on
			geneZscoresPathwayMatched.save(new File(intermediateFolder.getAbsolutePath() + "/", pathwayDatabase.getName() + "_Enrichment_normalizedGwasGeneScores" + (this.hlaGenesToExclude == null ? "" : "_ExHla") + ".txt").getAbsolutePath());

			if (LOGGER.isDebugEnabled()) {
				genePathwayZscores.saveBinary(new File(debugFolder, pathwayDatabase.getName() + "_Enrichment_normalizedPathwayScores" + (this.hlaGenesToExclude == null ? "" : "_ExHla")).getAbsolutePath());
				geneZscoresPathwayMatched.saveBinary(new File(debugFolder, pathwayDatabase.getName() + "_Enrichment_normalizedGwasGeneScores" + (this.hlaGenesToExclude == null ? "" : "_ExHla")).getAbsolutePath());
				geneZscoresNullGwasNullBetasPathwayMatched.save(new File(debugFolder, pathwayDatabase.getName() + "_Enrichment_normalizedNullGwasGeneScores" + (this.hlaGenesToExclude == null ? "" : "_ExHla") + ".txt"));
			}

			LinkedHashMap<String, Integer> singleColMap = new LinkedHashMap<>(1);
			singleColMap.put("B1", 0);

			// Storage for beta's per arm
			final List<DoubleMatrixDataset<String, String>> b1PerArm = Collections.synchronizedList(new ArrayList<>(metaGenesPerArm.size()));
			final List<DoubleMatrixDataset<String, String>> b2PerArm = Collections.synchronizedList(new ArrayList<>(metaGenesPerArm.size()));

			final List<DoubleMatrixDataset<String, String>> b1NullPerArm = Collections.synchronizedList(new ArrayList<>(metaGenesPerArm.size()));
			final List<DoubleMatrixDataset<String, String>> b2NullPerArm = Collections.synchronizedList(new ArrayList<>(metaGenesPerArm.size()));

			// Store the inverse gene-gene correlation matrices
			final Map<String, DoubleMatrixDataset<String, String>> inverseCorrelationMatrices = Collections.synchronizedMap(new HashMap<>(metaGenesPerArm.size()));

			// Advance progressbar
			pb.step();

			// Calculate the beta's per chromosome arm in parallel
			// for (Map.Entry<String, ArrayList<String>> chrArmMappingEntry : geneChrArmMapping.entrySet()) {
			metaGenesPerArm.entrySet().parallelStream().forEach((Map.Entry<String, ArrayList<MetaGene>> chrArmMappingEntry) -> {
				try {
					// Determine the meta genes to use
					final String chrArm = chrArmMappingEntry.getKey();
					final ArrayList<MetaGene> armGenes = chrArmMappingEntry.getValue();
					final ArrayList<String> chrArmGenesInPathwayMatrix = new ArrayList<>(armGenes.size());
					
					for (MetaGene armGene : armGenes) {
						chrArmGenesInPathwayMatrix.add(armGene.getMetaGeneId());
					}
					// Now genesInPathwayMatrix will only contain genes that are also in the gene p-value matrix
					LOGGER.debug("Number of meta genes in chr arm: " + chrArmGenesInPathwayMatrix.size());

					if (chrArmGenesInPathwayMatrix.isEmpty()) {
						throw new RuntimeException("This should not happen");
					}

					// b1 rows: traits cols: 1
					final DoubleMatrixDataset<String, String> b1Arm = new DoubleMatrixDataset<>(geneZscoresPathwayMatched.getHashCols(), singleColMap);
					// b2 rows: traits cols: genes
					final DoubleMatrixDataset<String, String> b2Arm = new DoubleMatrixDataset<>(geneZscoresPathwayMatched.getHashCols(), genePathwayZscores.getHashCols());
					// b1 rows: null traits cols: 1
					final DoubleMatrixDataset<String, String> b1NullGwasArm = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetasPathwayMatched.getHashCols(), singleColMap);
					// b2 rows: null traits cols: genes
					final DoubleMatrixDataset<String, String> b2NullGwasArm = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetasPathwayMatched.getHashCols(), genePathwayZscores.getHashCols());

					// Add to final output
					b1PerArm.add(b1Arm);
					b2PerArm.add(b2Arm);
					b1NullPerArm.add(b1NullGwasArm);
					b2NullPerArm.add(b2NullGwasArm);

					// Subset the full data to get the matched subset for the current chromosome arm
					final DoubleMatrixDataset<String, String> geneZscoresSubset = geneZscoresPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset = geneZscoresNullGwasCorrelationPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetasSubset = geneZscoresNullGwasNullBetasPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> genePathwayZscoresSubset = genePathwayZscores.viewRowSelection(chrArmGenesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> geneZscoresNullGwasSubsetGeneCorrelations;

					// Make gene-gene correlation matrix
					if (geneCorrelationWindow < 0) {
						LOGGER.debug("Creating full correlation matrix for chr arm");
						geneZscoresNullGwasSubsetGeneCorrelations = geneZscoresNullGwasCorrelationSubset.viewDice().calculateCorrelationMatrix();
					} else {
						LOGGER.debug("Creating correlation matrix in window: " + geneCorrelationWindow);
						geneZscoresNullGwasSubsetGeneCorrelations = createLocalGeneCorrelation(geneZscoresNullGwasCorrelationSubset, armGenes, geneCorrelationWindow);
					}
					// geneZscoresNullGwasSubsetGeneCorrelations = geneZscoresNullGwasCorrelationSubset.viewDice().calculateCorrelationMatrix();

					if (LOGGER.isDebugEnabled()) {
						geneZscoresNullGwasSubsetGeneCorrelations.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneCor.txt"));
						geneZscoresSubset.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneScores.txt"));
						genePathwayZscoresSubset.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_pathwayScores.txt"));
					}

					// Make the inverse of the gene-gene correlation matrix
					final DoubleMatrix2D geneInvCorMatrixSubsetMatrix;
					try {
						if (this.ignoreGeneCorrelations) {
							// Identity matrix, i.e. OLS
							LOGGER.info("Ignoring geOLSne correlations");
							geneInvCorMatrixSubsetMatrix = DoubleFactory2D.dense.identity(geneZscoresNullGwasSubsetGeneCorrelations.rows());
						} else {
							LOGGER.debug("Calculating correlation inverse");
							geneInvCorMatrixSubsetMatrix = new DenseDoubleAlgebra().inverse(geneZscoresNullGwasSubsetGeneCorrelations.getMatrix());
						}
					} catch (Exception ex) {
						LOGGER.fatal(pathwayDatabase.getName() + " " + chrArm + " number of genes: " + geneZscoresNullGwasSubsetGeneCorrelations.rows());
						throw ex;
					}

					// Convert to DoubleMatrixDataset
					DoubleMatrixDataset<String, String> geneInvCorMatrixSubset = new DoubleMatrixDataset<>(geneInvCorMatrixSubsetMatrix,
							geneZscoresNullGwasSubsetGeneCorrelations.getHashRows(),
							geneZscoresNullGwasSubsetGeneCorrelations.getHashCols());

					// Store in list
					inverseCorrelationMatrices.put(chrArm, geneInvCorMatrixSubset);

					if (LOGGER.isDebugEnabled()) {
						geneInvCorMatrixSubset.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneInvCor.txt").getAbsolutePath());
					}
					//final DoubleMatrixDataset<String, String> geneInvCorMatrixSubset = new DoubleMatrixDataset<>(geneInvCorMatrixSubsetMatrix, geneZscoresNullGwasNullBetasSubset.getHashRows(), geneZscoresNullGwasNullBetasSubset.getHashRows());
					//geneInvCorMatrixSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneInvCor.txt");

					// Determine B1 and B2
					glsStep1(geneZscoresSubset, geneInvCorMatrixSubsetMatrix, genePathwayZscoresSubset, b1Arm, b2Arm);
					glsStep1(geneZscoresNullGwasNullBetasSubset, geneInvCorMatrixSubsetMatrix, genePathwayZscoresSubset, b1NullGwasArm, b2NullGwasArm);

					if (LOGGER.isDebugEnabled()) {
						b1Arm.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_b1.txt"));
						b2Arm.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_b2.txt"));
					}

					// Advance progressbar
					pb.step();

				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}

			});
			// Closure of parallel computation

			// Combine beta's over different chromosome arms
			// b1 rows: traits cols: 1
			final DoubleMatrixDataset<String, String> b1 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), singleColMap);
			// b2 rows: traits cols: pathways
			final DoubleMatrixDataset<String, String> b2 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), genePathwayZscores.getHashCols());

			final DoubleMatrixDataset<String, String> b1NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), singleColMap);
			final DoubleMatrixDataset<String, String> b2NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), genePathwayZscores.getHashCols());

			// Combine the inverse gene-gene correlation matrix into one matrix
//			final DoubleMatrixDataset<String, String> fullInverseCorrelationMatrix = mergeCorrelationMatrices(inverseCorrelationMatrices.values());
			// Make sure the oder of metagenes is the same
			//TODO check if needed
//			fullInverseCorrelationMatrix.reorderRows(geneZscoresPathwayMatched.getHashRows());
//			fullInverseCorrelationMatrix.reorderCols(geneZscoresPathwayMatched.getHashRows());
			// Combine the results per arm can be done in parallel over the 4 different matrices
			for (DoubleMatrixDataset<String, String> b1Arm : b1PerArm) {
				b1.getMatrix().assign(b1Arm.getMatrix(), cern.jet.math.tdouble.DoubleFunctions.plus);
			}

			for (DoubleMatrixDataset<String, String> b2Arm : b2PerArm) {
				b2.getMatrix().assign(b2Arm.getMatrix(), cern.jet.math.tdouble.DoubleFunctions.plus);
			}

			for (DoubleMatrixDataset<String, String> b1Arm : b1NullPerArm) {
				b1NullGwas.getMatrix().assign(b1Arm.getMatrix(), cern.jet.math.tdouble.DoubleFunctions.plus);
			}

			for (DoubleMatrixDataset<String, String> b2Arm : b2NullPerArm) {
				b2NullGwas.getMatrix().assign(b2Arm.getMatrix(), cern.jet.math.tdouble.DoubleFunctions.plus);
			}

			if (LOGGER.isDebugEnabled()) {
				b1.save(new File(debugFolder, pathwayDatabase.getName() + "_Enrichment_b1.txt"));
				b2.save(new File(debugFolder, pathwayDatabase.getName() + "_Enrichment_b2.txt"));
			}

			LOGGER.debug("Merged b1 and b2");

			pb.step();

			// Determine final betas and p values
			betas = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());

			//standardErrors = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
			//tStatistics = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
			final DoubleMatrixDataset<String, String> pValues = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
			zscores = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
			//betasNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());
			numberOfPathways = genePathwayZscores.columns();

			final DoubleMatrixDataset<String, String> betasNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());
			final DoubleMatrixDataset<String, String> pValuesNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());

			if (LOGGER.isDebugEnabled()) {
				LOGGER.debug("betas: " + betas.rows() + " x " + betas.columns());
//				LOGGER.debug("standardErrors: " + standardErrors.rows() + " x " + standardErrors.columns());
//				LOGGER.debug("tStatistics: " + tStatistics.rows() + " x " + tStatistics.columns());
				LOGGER.debug("pValues: " + pValues.rows() + " x " + pValues.columns());
				LOGGER.debug("zValues: " + zscores.rows() + " x " + zscores.columns());
				//LOGGER.debug("betasNull: " + betasNull.rows() + " x " + betasNull.columns());
				LOGGER.debug("numberOfPathways: " + numberOfPathways);
			}

			//Now calculate beta's and residuals		
			final int df = geneZscoresPathwayMatched.rows() - 1;
			final TDistribution tdist = new TDistribution(df);

			for (int traitI = 0; traitI < numberTraits; ++traitI) {

				final int traitIb = traitI;//to allow use in parralle loop it must be final
				final double b1Trait = b1.getElementQuick(traitIb, 0);

				final DoubleMatrixDataset<String, String> residuals = genePathwayZscores.duplicate();

				IntStream.range(0, numberOfPathways).parallel().forEach(pathwayI -> {
					double beta = b2.getElementQuick(traitIb, pathwayI) / b1Trait;

					betas.setElementQuick(pathwayI, traitIb, beta);

					DoubleMatrix1D betaX = geneZscoresPathwayMatched.getCol(traitIb).copy();
					betaX.assign(DoubleFunctions.mult(beta));

					residuals.getCol(pathwayI).assign(betaX, DoubleFunctions.minus);

				});

				DoubleMatrixDataset<String, String> invCorXResiduals = new DoubleMatrixDataset<>(residuals.getHashRows(), residuals.getHashCols());

				inverseCorrelationMatrices.values().parallelStream().forEach((DoubleMatrixDataset<String, String> geneInvCorMatrixSubset) -> {

					DoubleMatrix2D residualsSubset = residuals.viewRowSelectionMatrix(geneInvCorMatrixSubset.getRowObjects());
					DoubleMatrix2D invCorXResidualsSubset = invCorXResiduals.viewRowSelectionMatrix(geneInvCorMatrixSubset.getRowObjects());

					geneInvCorMatrixSubset.getMatrix().zMult(residualsSubset, invCorXResidualsSubset);

				});

				IntStream.range(0, numberOfPathways).parallel().forEach(pathwayI -> {

					DoubleMatrix1D residualsPathway = residuals.getCol(pathwayI);
					DoubleMatrix1D invCorXResiddualsPathway = invCorXResiduals.getCol(pathwayI);

					double sigmaSqr = (residualsPathway.zDotProduct(invCorXResiddualsPathway)) / df;

					double seBeta = Math.sqrt(sigmaSqr / b1Trait);

					// Determine pvalue
					double beta = betas.getElementQuick(pathwayI, traitIb);

					double tstatistic = Math.abs(beta) / seBeta;
					double pvalue = tdist.cumulativeProbability(-tstatistic) * 2;
					double zscore = ZScores.pToZTwoTailed(pvalue);

					// The zscore returned by pToZTwoTailed is always negative, therefore, match direction on beta
					if (beta > 0) {
						zscore = -zscore;
					}

//					standardErrors.setElementQuick(pathwayI, traitI, seBeta);
//					tStatistics.setElementQuick(pathwayI, traitI, tstatistic);
					pValues.setElementQuick(pathwayI, traitIb, pvalue);
					zscores.setElementQuick(pathwayI, traitIb, zscore);

				});

				pb.step();

			}

			LOGGER.debug("Done calculating model coefficients");

			// Betas for null
			LOGGER.debug("Calculating null beta's");

			for (int traitI = 0; traitI < numberTraitsNull; ++traitI) {
				final int traitIb = traitI;//to allow use in parralle loop it must be final
				final double b1Trait = b1NullGwas.getElementQuick(traitI, 0);
				final DoubleMatrixDataset<String, String> residuals = genePathwayZscores.duplicate();

				IntStream.range(0, numberOfPathways).parallel().forEach(pathwayI -> {
					//for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

					double beta = b2NullGwas.getElementQuick(traitIb, pathwayI) / b1Trait;
					betasNull.setElementQuick(pathwayI, traitIb, beta);

					DoubleMatrix1D betaX = geneZscoresNullGwasNullBetasPathwayMatched.getCol(traitIb).copy();
					betaX.assign(DoubleFunctions.mult(beta));

					residuals.getCol(pathwayI).assign(betaX, DoubleFunctions.minus);

				});

				DoubleMatrixDataset<String, String> invCorXResiduals = new DoubleMatrixDataset<>(residuals.getHashRows(), residuals.getHashCols());

				inverseCorrelationMatrices.values().parallelStream().forEach((DoubleMatrixDataset<String, String> geneInvCorMatrixSubset) -> {

					DoubleMatrix2D residualsSubset = residuals.viewRowSelectionMatrix(geneInvCorMatrixSubset.getRowObjects());
					DoubleMatrix2D invCorXResidualsSubset = invCorXResiduals.viewRowSelectionMatrix(geneInvCorMatrixSubset.getRowObjects());

					geneInvCorMatrixSubset.getMatrix().zMult(residualsSubset, invCorXResidualsSubset);

				});

				IntStream.range(0, numberOfPathways).parallel().forEach(pathwayI -> {

					DoubleMatrix1D residualsPathway = residuals.getCol(pathwayI);
					DoubleMatrix1D invCorXResiddualsPathway = invCorXResiduals.getCol(pathwayI);

					double sigmaSqrt = (residualsPathway.zDotProduct(invCorXResiddualsPathway)) / df;

					double seBeta = Math.sqrt(sigmaSqrt / b1Trait);

					// Determine pvalue
					double beta = betasNull.getElementQuick(pathwayI, traitIb);
					double tstatistic = Math.abs(beta) / seBeta;
					double pvalue = tdist.cumulativeProbability(-tstatistic) * 2;

					pValuesNull.setElementQuick(pathwayI, traitIb, pvalue);

				});

				pb.step();

			}
			LOGGER.debug("Done calculating null pvalues");

			qValues = pValues.duplicate();
			final DoubleMatrix1D sortedNullPvalues = pValuesNull.getMatrix().vectorize().viewSorted();
			final long permutedPvalues = sortedNullPvalues.size();
			final long permutedPvaluesMin1 = permutedPvalues - 1;
			final double numberTraitsNullD = (double) numberTraitsNull;

			IntStream.range(0, numberTraits).parallel().forEach(traitI -> {
				//for (int traitI = 0; traitI < numberTraits; ++traitI) {

				DoubleMatrix1D qValuesTraitSorted = qValues.viewCol(traitI).viewSorted();

				int indexNullPvalues = -1;

				for (int i = 0; i < qValuesTraitSorted.size(); ++i) {

					//initially qvalue matrix contains the pvaluess
					final double currentP = qValuesTraitSorted.get(i);

					while (indexNullPvalues < permutedPvaluesMin1 && sortedNullPvalues.get(indexNullPvalues + 1) <= currentP) {
						indexNullPvalues++;
					}

					//at this point indexNullPvalues + 1 is number of times pvalue or something smaller is seen in the permutations
					double qvalue = ((indexNullPvalues + 1) / numberTraitsNullD) / (i + 1);

					//now overwrite this p-value with q-value
					qValuesTraitSorted.setQuick(i, qvalue);

				}

				for (int i = (int) qValuesTraitSorted.size() - 2; i >= 0; --i) {
					if (qValuesTraitSorted.get(i) > qValuesTraitSorted.get(i + 1)) {
						qValuesTraitSorted.set(i, qValuesTraitSorted.get(i + 1));
					}
				}

				pb.step();

			});

			// Write output
			betas.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_betas" : "_betasExHla"));
			//standardErrors.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_se" : "_seExHla"));
			//pValues.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_analyticalPvals" : "_analyticalPvalsExHla"));
			//betasNull.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_EnrichmentNull" + (this.hlaGenesToExclude == null ? "_betas" : "_betasExHla"));
			pValuesNull.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_EnrichmentNull" + (this.hlaGenesToExclude == null ? "_analyticalPvals" : "_analyticalPvalsExHla"));
			betas.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_betas.txt" : "_betasExHla.txt"));
			//standardErrors.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_se.txt" : "_seExHla.txt"));
			pValues.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_analyticalPvals.txt" : "_analyticalPvalsExHla.txt"));
			qValues.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_analyticalQvals.txt" : "_analyticalQvalsExHla.txt"));

			// Save as txt to avoid having to convert them later
			if (LOGGER.isDebugEnabled()) {
				//fullInverseCorrelationMatrix.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_inverseCorrelationMatrix.txt" : "_inverseCorrelationMatrixExHla.txt"));
				geneZscoresPathwayMatched.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_geneZscoresPathwayMatched.txt" : "_geneZscoresPathwayMatchedExHla.txt"));
				genePathwayZscores.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_genePathwayZscores.txt" : "_genePathwayZscoresExHla.txt"));
			}

			zscores.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla"));
			//this.clearZscoreCache();

			pb.step();
		}
	}

	private static void inplaceDetermineGeneLengthRegressionResiduals(final DoubleMatrixDataset<String, String> geneZscores, final double[] geneLengths) {

		if (LOGGER.isDebugEnabled()) {
			LOGGER.debug("Gene length vector: " + geneLengths.length);
			LOGGER.debug("Determining residuals for matrix of size: " + geneZscores.rows() + " x " + geneZscores.columns());
//
//			int maxColIdx = geneZscores.columns();
//			if (maxColIdx > 6) {
//				maxColIdx = 6;
//			}
//			geneZscores.viewSelection(
//					geneZscores.getRowObjects().subList(0, 10),
//					geneZscores.getColObjects().subList(0, maxColIdx)).printMatrix();
		}

		final DoubleMatrix2D geneZscoresMatrix = geneZscores.getMatrix();
		final int genes = geneZscores.rows();

		IntStream.range(0, geneZscores.columns()).parallel().forEach(c -> {
			// Build the regression model
			final SimpleRegression regression = new SimpleRegression();
			DoubleMatrix1D geneZscoresCol = geneZscoresMatrix.viewColumn(c);

			for (int r = 0; r < genes; ++r) {
				regression.addData(geneLengths[r], geneZscoresCol.getQuick(r));
			}

			// Apply model and inplace convert to residuals
			for (int r = 0; r < genes; ++r) {
				double predictedY = regression.predict(geneLengths[r]);
				geneZscoresCol.setQuick(r, (geneZscoresCol.getQuick(r) - predictedY));
			}
		});

//		if (LOGGER.isDebugEnabled()) {
//			LOGGER.debug("Determined residuals for matrix of size: " + geneZscores.rows() + " x " + geneZscores.columns());
//			int maxColIdx = geneZscores.columns();
//			if (maxColIdx > 6) {
//				maxColIdx = 6;
//			}
//			geneZscores.viewSelection(
//					geneZscores.getRowObjects().subList(0, 10),
//					geneZscores.getColObjects().subList(0, maxColIdx)).printMatrix();
//		}
	}

	/**
	 * Determine the Pvalue, Se and Tstat for a given gls regression. Only valid
	 * for cases where the data is centered and scaled and has no intercept. //
	 * TODO: Generalize this so an intercept term can be included
	 *
	 * @param beta
	 * @param b1
	 * @param geneInvCor
	 * @param genePvalues
	 * @param genePathwayZscores
	 * @return
	 */
	private static GenePathwayAssociationStatistic determineGlsAnalyticalPvalues(double beta,
			double b1,
			DoubleMatrix2D geneInvCor,
			DoubleMatrix1D genePvalues,
			DoubleMatrix1D genePathwayZscores) {
		// Grab the size of the array
		int n;
		if (genePvalues.size() < Integer.MAX_VALUE) {
			n = (int) genePvalues.size();
		} else {
			throw new IllegalArgumentException("Input too large for int, this shouldn't happen");
		}

		// df is n-1 as there is no intercept
		int df = n - 1;

		// Determine model residuals
		DoubleMatrix1D betaX = genePvalues;
		betaX.assign(DoubleFunctions.mult(beta));
		DoubleMatrix1D residuals = genePathwayZscores;
		residuals.assign(betaX, DoubleFunctions.minus);

		// TODO: Ugly, but dont know how I can do this better as I need a DoubleMatrix2D for matrix mult
		// Cat residuals to n * 1 DoubleMatrix2D
		DoubleMatrix2D residualMatrix = residuals.like2D(n, 1);
		for (int r = 0; r < n; ++r) {
			residualMatrix.setQuick(r, 0, residuals.get(r));
		}

		// Determine sigma squared
		// part1 =  t(residuals) * geneInvCor
		// 1 x n matrix
		DoubleMatrix2D part1 = residualMatrix.like(1, n);
		residualMatrix.zMult(geneInvCor, part1, 1, 0, true, false);
		// part2 =  t(residuals) * geneInvCor * residuals
		// 1 x 1 matrix
		DoubleMatrix2D part2 = residualMatrix.like(1, 1);
		part1.zMult(residualMatrix, part2, 1, 0, false, false);
		double sigmaSquared = part2.get(0, 0) / df;

		// Determine Se for beta
		// as b1 is not inverse divide instead of multiply. If you want to include an intercept,
		// this needs to be a matrix mult on the inverse of b1
		//double standardError = Math.sqrt(b1) * Math.sqrt(sigmaSquared);
		double standardError = Math.sqrt(sigmaSquared) / Math.sqrt(b1);

		// Determine pvalue
		double tstatistic = Math.abs(beta / standardError);
		double pvalue = new TDistribution(df).cumulativeProbability(-tstatistic) * 2;
		double zscore = ZScores.pToZTwoTailed(pvalue);

		// The zscore returned by pToZTwoTailed is always negative, therefore, match direction on beta
		if (beta > 0) {
			zscore = ((double) -1) * zscore;
		}

		return new GenePathwayAssociationStatistic(beta, standardError, tstatistic, zscore, pvalue);
	}

	public final DoubleMatrixDataset<String, String> getEnrichmentZscores() throws IOException {
		return zscores;
	}

	public DoubleMatrixDataset<String, String> getqValues() {
		return qValues;
	}

	@Deprecated
	public final DoubleMatrixDataset<String, String> getEmpericalEnrichmentZscores() throws IOException {

		throw new RuntimeException("No longer functional");

//		if (enrichmentPvalues == null) {
//
//			enrichmentPvalues = betas.duplicate();
//
//			final DoubleMatrix2D zscoreMatrix = enrichmentPvalues.getMatrix();
//			
////			final DoubleMatrix2D betasNullMatrix = betasNull.getMatrix();
//
//			final int numberOfPathways = zscoreMatrix.rows();
//			final int numberOfPhenotypes = zscoreMatrix.columns();
//			final int numberOfNullGwasPhenotypes = betasNullMatrix.columns();
//			final double numberOfNullGwasPhenotypesMin1Double = betasNullMatrix.columns() - 1;
//
//			LOGGER.debug("numberOfNullGwasPhenotypes: " + numberOfNullGwasPhenotypes);
//
////			List<String> pathwayNames;
////			if (LOGGER.isDebugEnabled()) {
////				pathwayNames = betasNull.getRowObjects();
////			} else {
////				pathwayNames = Collections.emptyList();
////			}
//			for (int r = 0; r < numberOfPathways; ++r) {
//
//				double meanNull = 0;
//				for (int p = 0; p < numberOfNullGwasPhenotypes; ++p) {
//					meanNull += betasNullMatrix.getQuick(r, p);
//				}
//				meanNull /= numberOfNullGwasPhenotypes;
//
//				double x = 0;
//				for (int p = 0; p < numberOfNullGwasPhenotypes; ++p) {
//					x += (betasNullMatrix.getQuick(r, p) - meanNull) * (betasNullMatrix.getQuick(r, p) - meanNull);
//				}
//				double sdNull = Math.sqrt(x / numberOfNullGwasPhenotypesMin1Double);
//
////				if (LOGGER.isDebugEnabled()) {750
////					LOGGER.debug(pathwayNames.get(r) + " mean: " + meanNull + " sd: " + sdNull);
////				}
//				for (int c = 0; c < numberOfPhenotypes; ++c) {
//
//					final double corr = zscoreMatrix.getQuick(r, c);
//
//					zscoreMatrix.setQuick(r, c, (corr - meanNull) / sdNull);
//
//				}
//			}
//
//			enrichmentPvalues.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla"));
//		}
//
//		return enrichmentPvalues;
	}

//	public final void clearZscoreCache() {
//		enrichmentPvalues = null;
//	}
	private static void glsStep1(DoubleMatrixDataset<String, String> geneZscoresSubset, DoubleMatrix2D geneInvCorMatrix, DoubleMatrixDataset<String, String> genePathwayZscoresSubset, DoubleMatrixDataset<String, String> b1, DoubleMatrixDataset<String, String> b2) {

		final int numberOfGenes = geneZscoresSubset.rows();
		final int numberTraits = geneZscoresSubset.columns();
		final int numberOfPathways = genePathwayZscoresSubset.columns();
		//final DoubleMatrix2D geneInvCorMatrix = geneInvCorMatrixSubset.getMatrix();
		final DoubleMatrix2D genePathwayZscoresMatrix = genePathwayZscoresSubset.getMatrix();

		// Result of transpose geneZscoresTrait times inv correlation matrix
		DoubleMatrix2D A = geneZscoresSubset.getMatrix().like(1, numberOfGenes);

		// Trait = gwas
		for (int traitI = 0; traitI < numberTraits; ++traitI) {
			try {

				DoubleMatrix2D geneZscoresTrait = geneZscoresSubset.viewColAsMmatrix(traitI);
				geneZscoresTrait.zMult(geneInvCorMatrix, A, 1, 0, true, false);
				final double x = A.viewRow(0).zDotProduct(geneZscoresTrait.viewColumn(0));

				// Col order should be the same
				b1.setElementQuick(0, traitI, x + b1.getElementQuick(0, traitI));
				DoubleMatrix2D b2Row = b2.viewRowAsMmatrix(traitI);
				A.zMult(genePathwayZscoresMatrix, b2Row, 1, 0, false, false);

			} catch (Exception e) {
				LOGGER.fatal("Number of pathways: " + numberOfPathways);
				LOGGER.fatal("Current trait index: " + traitI);
				LOGGER.fatal("Dim genePathwayZscores: " + genePathwayZscoresSubset.rows() + "x" + genePathwayZscoresSubset.columns());
				LOGGER.fatal("Dim genePathwayZscores internal: " + genePathwayZscoresSubset.getMatrix().rows() + "x" + genePathwayZscoresSubset.getMatrix().columns());
				throw (e);
			}
		}

	}

	private static Map<String, ArrayList<Gene>> createChrArmGeneMapping(List<Gene> genes, Set<String> includedGenes) {
		Map<String, ArrayList<Gene>> chrArmToGeneMapping = new HashMap<>(25);
		for (Gene gene : genes) {

			if (includedGenes.contains(gene.getGene())) {

				String chrArm = gene.getChrAndArm();

				ArrayList<Gene> armGenes = chrArmToGeneMapping.get(chrArm);
				if (armGenes == null) {
					armGenes = new ArrayList<>();
					chrArmToGeneMapping.put(chrArm, armGenes);
				}

				armGenes.add(gene);

			}

		}
		return chrArmToGeneMapping;
	}

	public PathwayDatabase getPathwayDatabase() {
		return pathwayDatabase;
	}

	public int getNumberOfPathways() {
		return numberOfPathways;
	}

	private static LinkedHashSet<String> findUncorrelatedGenes(DoubleMatrixDataset<String, String> geneZscoresNullGwas, HashSet<String> genesWithPvalue, List<Gene> genes, double maxCorrelationBetweenGenes) {

		final Map<String, ArrayList<Gene>> chrArmToGeneMapping = createChrArmGeneMapping(genes, genesWithPvalue);

		final LinkedHashSet<String> includedUncorrelatedGenes = new LinkedHashSet<>();

		chrArmToGeneMapping.keySet().parallelStream().forEach((String chrArm) -> {

			final ArrayList<Gene> armGenes = chrArmToGeneMapping.get(chrArm);
			final ArrayList<String> armGenesIds = new ArrayList<>(armGenes.size());

			final HashSet<String> includedUncorrelatedGenesArm = new HashSet<>();

			for (Gene armGene : armGenes) {
				armGenesIds.add(armGene.getGene());
			}

			final DoubleMatrixDataset<String, String> geneZscoresNullGwasArm = geneZscoresNullGwas.viewRowSelection(armGenesIds);

			//final DoubleMatrixDataset<String, String> invCorMatrixArmGenes = invCorMatrix.viewSelection(armGenes, armGenes);
			final DoubleMatrixDataset<String, String> genePvaluesNullGwasArmT = geneZscoresNullGwasArm.viewDice();

			final DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelation = genePvaluesNullGwasArmT.calculateCorrelationMatrix();

			//We need to take the inverse of the correlation matrix. To do that the correlation between genes can't be correlated
			//Simply removing highly correlated genes did not always work, therefor:
			//(1) create correlation matrix of correlations
			//(2) identifie genes that have correlated correlation
			//(3) prune gene correlation matrix
			DoubleMatrixDataset<String, String> correlationOfCorrelations = genePvaluesNullGwasGeneArmCorrelation.calculateCorrelationMatrix();

			ArrayList<String> variantNames = correlationOfCorrelations.getRowObjects();

			rows:
			for (int r = 0; r < correlationOfCorrelations.rows(); ++r) {
				cols:
				for (int c = 0; c < r; ++c) {
					if (Math.abs(correlationOfCorrelations.getElementQuick(r, c)) >= maxCorrelationBetweenGenes && includedUncorrelatedGenesArm.contains(variantNames.get(c))) {
						continue rows;
					}
				}
				includedUncorrelatedGenesArm.add(variantNames.get(r));
			}

			synchronized (includedUncorrelatedGenes) {
				includedUncorrelatedGenes.addAll(includedUncorrelatedGenesArm);
			}

		});

		return includedUncorrelatedGenes;
	}

		private static DoubleMatrixDataset<String, String> createLocalGeneCorrelation(final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset, final ArrayList<MetaGene> genes, final int correlationWindow) {

		if (genes.size() != geneZscoresNullGwasCorrelationSubset.rows()) {
			throw new RuntimeException("Genes should match geneZscoresNullGwasCorrelationSubset");
		}

		final DoubleMatrixDataset<String, String> correlations = new DoubleMatrixDataset<>(geneZscoresNullGwasCorrelationSubset.getHashRows(), geneZscoresNullGwasCorrelationSubset.getHashRows());
		final DoubleMatrix2D correlationMatrix = correlations.getMatrix();
		final int geneCount = geneZscoresNullGwasCorrelationSubset.rows();
		final int nullGwasCount = geneZscoresNullGwasCorrelationSubset.columns();
		DoubleMatrix2D geneZscoresNullGwasCorrelationSubsetMatrix = geneZscoresNullGwasCorrelationSubset.getMatrix();

		final SimpleRegression regression = new SimpleRegression();

		for (int i = geneCount; --i >= 0;) {
			for (int j = i + 1; --j >= 0;) {
				regression.clear();

				if (i == j) {
					correlationMatrix.setQuick(i, j, 1);
				} else {

					//Genes should be in the same order as the matrix
					MetaGene geneI = genes.get(i);
					MetaGene geneJ = genes.get(j);

					//Only look at position because this is done per chromosome arm
					int geneIStart = geneI.getStart();
					int geneIStop = geneI.getStop();

					int geneJStart = geneJ.getStart();
					int geneJStop = geneJ.getStop();

					if (Math.abs(geneIStart - geneJStart) <= correlationWindow
							|| Math.abs(geneIStart - geneJStop) <= correlationWindow
							|| Math.abs(geneIStop - geneJStart) <= correlationWindow
							|| Math.abs(geneIStop - geneJStop) <= correlationWindow) {
						for (int n = 0; n < nullGwasCount; ++n) {
							regression.addData(geneZscoresNullGwasCorrelationSubsetMatrix.getQuick(i, n), geneZscoresNullGwasCorrelationSubsetMatrix.getQuick(j, n));
						}

						double x = regression.getR();

						correlationMatrix.setQuick(i, j, x);
						correlationMatrix.setQuick(j, i, x); // symmetric
					}

				}
			}
		}

		return correlations;

	}

	protected static HashMap<String, ArrayList<MetaGene>> groupCorrelatedGenesPerChrArm(final DoubleMatrixDataset<String, String> geneZscoresNullGwas, final double maxCorrelationBetweenGenes, final List<Gene> genes, final Set<String> includedGenes, final File debugFolder, final String pathwayDatabaseName, final HashSet<String> hlaGenesToExclude) {

		HashMap<String, ArrayList<MetaGene>> metaGenes = new HashMap<>();

		final Map<String, ArrayList<Gene>> geneChrArmMapping = createChrArmGeneMapping(genes, includedGenes);

		geneChrArmMapping.keySet().parallelStream().forEach((String chrArm) -> {

			final ArrayList<Gene> armGenes = geneChrArmMapping.get(chrArm);
			final ArrayList<String> armGenesIds = new ArrayList<>(armGenes.size());

			HashMap<String, MetaGene> metaGenesArm = new HashMap<>(armGenes.size());

			for (int i = 0; i < armGenes.size(); i++) {
				armGenesIds.add(armGenes.get(i).getGene());
			}

			final DoubleMatrixDataset<String, String> geneZscoresNullGwasArm = geneZscoresNullGwas.viewRowSelection(armGenesIds);

			final DoubleMatrixDataset<String, String> genePvaluesNullGwasGeneArmCorrelation = geneZscoresNullGwasArm.viewDice().calculateCorrelationMatrix();

			//We need to take the inverse of the correlation matrix. To do that the correlation between genes can't be correlated
			//Simply removing highly correlated genes did not always work, therefor:
			//(1) create correlation matrix of correlations
			//(2) identifie genes that have correlated correlation
			//(3) prune gene correlation matrix
			//DoubleMatrixDataset<String, String> correlationOfCorrelations = genePvaluesNullGwasGeneArmCorrelation.calculateCorrelationMatrix();
			if (LOGGER.isDebugEnabled()) {
				try {
					genePvaluesNullGwasGeneArmCorrelation.save(new File(debugFolder, pathwayDatabaseName + "_" + chrArm + "_Enrichment_geneCorrelationsForMetaGenes" + (hlaGenesToExclude == null ? "" : "_ExHla") + ".txt").getAbsolutePath());
					//correlationOfCorrelations.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneCorrelationsOfCorrelationsForMetaGenes" + (hlaGenesToExclude == null ? "" : "_ExHla") + ".txt").getAbsolutePath());
				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}
			}

			rows:
			for (int r = 0; r < genePvaluesNullGwasGeneArmCorrelation.rows(); ++r) {

				String currentGene = armGenesIds.get(r);
				
				if(!currentGene.equals(armGenes.get(r).getGene())){
					throw new RuntimeException("Internal error");
				}

				MetaGene currentMetaGene = new MetaGene(currentGene, armGenes.get(r).getStart(), armGenes.get(r).getStop());
				metaGenesArm.put(currentGene, currentMetaGene);

				cols:
				for (int c = 0; c < r; ++c) {
					if (Math.abs(genePvaluesNullGwasGeneArmCorrelation.getElementQuick(r, c)) >= maxCorrelationBetweenGenes) {

						//Never null because c < r
						MetaGene otherMetaGene = metaGenesArm.get(armGenesIds.get(c));
						currentMetaGene.addOtherMetaGene(otherMetaGene);

						for (String otherGene : otherMetaGene.getGenes()) {
							metaGenesArm.put(otherGene, currentMetaGene);
						}

					}
				}

			}

			HashSet<MetaGene> metaGenesArmUnique = new HashSet<>(metaGenesArm.values());
			
//			Iterator<MetaGene> iterator = metaGenesArmUnique.iterator();
//			while (iterator.hasNext()) {
//				if (iterator.next().getGeneCount() > 1) {
//					iterator.remove();
//				}
//			}
			synchronized (metaGenes) {
				metaGenes.put(chrArm, new ArrayList<>(metaGenesArmUnique));
			}

		});

		return metaGenes;
	}

	protected static class MetaGene {

		final HashSet<String> genes = new HashSet<>();
		int start;
		int stop;

		public MetaGene(String gene, int start, int stop) {
			genes.add(gene);
			//make sure reverse genes have start stop relative to forward
			this.start = Math.min(start, stop);
			this.stop = Math.max(start, stop);
		}

		/**
		 * only use for unit testing
		 * @param genes 
		 */
		protected MetaGene(int start, int stop, String... genes) {
			this.genes.addAll(Arrays.asList(genes));
			this.start = 0;
			this.stop = 0;
		}

		public HashSet<String> getGenes() {
			return genes;
		}

		public void addGene(String gene) {
			genes.add(gene);
		}

		public void addOtherMetaGene(MetaGene other) {
			genes.addAll(other.getGenes());
			this.start = Math.min(this.start, other.start);
			this.stop = Math.max(this.stop, other.stop);
		}

		public int getGeneCount() {
			return genes.size();
		}

		public String getMetaGeneId() {
			return String.join("_", genes);
		}

		@Override
		public int hashCode() {
			return this.genes.hashCode();
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			final MetaGene other = (MetaGene) obj;
			if (!Objects.equals(this.genes, other.genes)) {
				return false;
			}
			return true;
		}

		public int getStart() {
			return start;
		}

		public int getStop() {
			return stop;
		}

	}

	/**
	 * @param dataset with genes on rows
	 */
	protected static DoubleMatrixDataset<String, String> collapseDatasetToMetaGenes(final DoubleMatrixDataset<String, String> dataset, final boolean zscoreSum, final Collection<ArrayList<MetaGene>> metaGenes) {

		LinkedHashMap<String, Integer> metaGenesRows = new LinkedHashMap<>();

		int r = 0;
		for (ArrayList<MetaGene> metaGenesPerArm : metaGenes) {
			for (MetaGene metaGene : metaGenesPerArm) {
				metaGenesRows.put(metaGene.getMetaGeneId(), r++);
			}
		}

		final int cols = dataset.columns();
		final DoubleMatrixDataset datasetCollapsed = new DoubleMatrixDataset<>(metaGenesRows, dataset.getHashCols());

		for (ArrayList<MetaGene> metaGenesPerArm : metaGenes) {
			for (MetaGene metaGene : metaGenesPerArm) {
				if (metaGene.getGeneCount() == 1) {
					//In this case meta gene ID is only the original gene ID.
					DoubleMatrix1D collapedRow = datasetCollapsed.getRow(metaGene.getMetaGeneId());
					DoubleMatrix1D originalRow = dataset.getRow(metaGene.getMetaGeneId());

					for (int i = 0; i < cols; ++i) {
						collapedRow.setQuick(i, originalRow.getQuick(i));
					}

				} else {

					DoubleMatrix1D collapedRow = datasetCollapsed.getRow(metaGene.getMetaGeneId());

					for (String geneInMeta : metaGene.getGenes()) {

						DoubleMatrix1D originalRow = dataset.getRow(geneInMeta);
						collapedRow.assign(originalRow, DoubleFunctions.plus);

					}

					final double denominator = zscoreSum ? Math.sqrt(metaGene.getGeneCount()) : metaGene.getGeneCount();

					for (int i = 0; i < cols; ++i) {
						collapedRow.setQuick(i, collapedRow.getQuick(i) / denominator);
					}

				}
			}
		}

		return datasetCollapsed;

	}

	/**
	 * Merges a list of correlation matrices so they are on the diagonal, and
	 * the rest is padded with 0 Matrix needs to be square
	 *
	 * @param correlationMatrices list of correlation matrices as
	 * DoubleMatrixDatasets
	 * @return One correlation matrix with the off diagonal parts not covered by
	 * a matrix as 0
	 */
	public static DoubleMatrixDataset<String, String> mergeCorrelationMatrices(Collection<DoubleMatrixDataset<String, String>> correlationMatrices) throws Exception {

		// Determine what the dimension of the output is going to be and check for matrix squareness
		int matrixDimension = 0;
		for (DoubleMatrixDataset<String, String> curMatrix : correlationMatrices) {
			if (curMatrix.columns() != curMatrix.rows()) {
				LOGGER.fatal("Input matrix not square, this should be square if it is a correlation matrix");
				throw new IllegalArgumentException("Input matrix not square, this should be square if it is a correlation matrix");
			}
			matrixDimension += curMatrix.columns();
		}

		if (LOGGER.isDebugEnabled()) {
			LOGGER.debug("Merging over " + correlationMatrices.size() + " matrices");
			LOGGER.debug("Total dimension: " + matrixDimension);
		}

		// Initialize the output storage
		DoubleMatrixDataset<String, String> fullCorrelationMatrix = new DoubleMatrixDataset<>(matrixDimension, matrixDimension);
		List<String> rowNames = new ArrayList<>(matrixDimension);
		List<String> colNames = new ArrayList<>(matrixDimension);

		// Merge the matrices
		int index = 0;
		for (DoubleMatrixDataset<String, String> curMatrix : correlationMatrices) {

			int start = index;
			// The column count is not zero based, so it should work as the end index, i.e. column count 100 is index 101
			// Since r < end, r goes to 99 not 100
			int end = index + curMatrix.columns();

			for (int r = start; r < end; r++) {
				for (int c = start; c < end; c++) {
					fullCorrelationMatrix.setElementQuick(r, c, curMatrix.getElementQuick(r - start, c - start));
				}
				rowNames.add(curMatrix.getRowObjects().get(r - start));
				colNames.add(curMatrix.getColObjects().get(r - start));
			}
			// The column count is not zero based, so it should work as the next index, i.e. column count 100 is index 101
			index += curMatrix.columns();
		}

		fullCorrelationMatrix.setRowObjects(rowNames);
		fullCorrelationMatrix.setColObjects(colNames);

		if (LOGGER.isDebugEnabled()) {
			LOGGER.debug("Index count. Should be equal to total dimension: " + index);
			LOGGER.debug("Done merging");
		}

		return fullCorrelationMatrix;
	}

	public static DoubleMatrixDataset<String, String> createColumnForceNormalDuplicate(DoubleMatrixDataset<String, String> matrix, DoubleMatrixDataset<String, String> tieBreaker) {

		DoubleMatrixDataset<String, String> newDataset = new DoubleMatrixDataset<>(matrix.getHashRows(), matrix.getHashCols());

		NaturalRankingTieFighter ranking = new NaturalRankingTieFighter(NaNStrategy.FAILED,
				TiesStrategy.AVERAGE);

		IntStream.range(0, matrix.columns()).parallel().forEach(c -> {

			double[] col = matrix.getCol(c).toArray();
			double[] colTie = tieBreaker.getCol(c).toArray();

			double mean = JSci.maths.ArrayMath.mean(col);
			double stdev = JSci.maths.ArrayMath.standardDeviation(col);

			double[] rankedValues = ranking.rank(col, colTie);

			for (int s = 0; s < matrix.rows(); s++) {
				double pValue = (0.5d + rankedValues[s] - 1d) / (double) (rankedValues.length);

				newDataset.setElementQuick(s, c, mean + cern.jet.stat.Probability.normalInverse(pValue) * stdev);
			}

		});

		return newDataset;

	}

}
