/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleSingularValueDecomposition;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.log4j.Logger;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;

/**
 *
 * @author patri
 */
public class PathwayEnrichments {

	private static final Logger LOGGER = Logger.getLogger(Depict2.class);

	private final PathwayDatabase pathwayDatabase;
	private final Set<String> excludeGenes;
	private final List<Gene> genes;
	private final String outputBasePath;
	private final HashSet<String> hlaGenesToExclude;
	private final boolean ignoreGeneCorrelations;
	private final DoubleMatrixDataset<String, String> betas;
	private final DoubleMatrixDataset<String, String> betasNull;
	private DoubleMatrixDataset<String, String> enrichmentPvalues = null;
	private final int numberOfPathways;

	public PathwayEnrichments(final PathwayDatabase pathwayDatabase, final HashSet<String> genesWithPvalue, final List<Gene> genes, final boolean forceNormalPathwayPvalues, final boolean forceNormalGenePvalues, final DoubleMatrixDataset<String, String> geneZscores, final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelation, final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetas, final String outputBasePath, final HashSet<String> hlaGenesToExclude, final boolean ignoreGeneCorrelations, final double genePruningR) throws IOException {
		this.pathwayDatabase = pathwayDatabase;
		this.genes = genes;
		this.outputBasePath = outputBasePath;
		this.hlaGenesToExclude = hlaGenesToExclude;
		this.ignoreGeneCorrelations = ignoreGeneCorrelations;

		final DoubleMatrixDatasetFastSubsetLoader pathwayMatrixLoader = new DoubleMatrixDatasetFastSubsetLoader(pathwayDatabase.getLocation());
		final LinkedHashSet<String> sharedGenes;

		if (this.hlaGenesToExclude == null) {
			excludeGenes = Collections.emptySet();
		} else {
			excludeGenes = this.hlaGenesToExclude;
		}

		Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();

		sharedGenes = new LinkedHashSet<>();

		for (String gene : genesWithPvalue) {
			if (pathwayGenes.contains(gene) && !excludeGenes.contains(gene)) {
				sharedGenes.add(gene);
			}
		}

		final Map<String, ArrayList<Gene>> geneChrArmMapping = createChrArmGeneMapping(genes, sharedGenes);

		try (ProgressBar pb = new ProgressBar(pathwayDatabase.getName() + " enrichtment analysis", geneChrArmMapping.size() + 2, ProgressBarStyle.ASCII)) {

			DoubleMatrixDataset<String, String> tmp = geneZscoresNullGwasCorrelation.viewRowSelection(sharedGenes).duplicate();
			tmp.normalizeColumns();

			LinkedHashSet<String> sharedUncorrelatedGenes = findUncorrelatedGenes(tmp, sharedGenes, genes, genePruningR);

			LOGGER.debug("Number of uncorrelated genes: " + sharedUncorrelatedGenes.size());

			final DoubleMatrixDataset<String, String> genePathwayZscores;
			if (forceNormalPathwayPvalues) {
				LOGGER.debug("Doing force normal pathway scores");
				genePathwayZscores = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
			} else {
				genePathwayZscores = pathwayMatrixLoader.loadSubsetOfRowsBinaryDoubleData(sharedUncorrelatedGenes);
				LOGGER.debug("Center and scale pathway scores");
				genePathwayZscores.normalizeColumns();
			}

			final DoubleMatrixDataset<String, String> geneZscoresPathwayMatched;
			final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationPathwayMatched;
			final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetasPathwayMatched;

			if (forceNormalGenePvalues) {
				LOGGER.debug("Doing force normal gene p-values");
				geneZscoresPathwayMatched = geneZscores.viewRowSelection(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
				geneZscoresNullGwasCorrelationPathwayMatched = geneZscoresNullGwasCorrelation.viewRowSelection(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
				geneZscoresNullGwasNullBetasPathwayMatched = geneZscoresNullGwasNullBetas.viewRowSelection(sharedUncorrelatedGenes).createColumnForceNormalDuplicate();
			} else {
				geneZscoresPathwayMatched = geneZscores.viewRowSelection(sharedUncorrelatedGenes).duplicate();
				geneZscoresNullGwasCorrelationPathwayMatched = geneZscoresNullGwasCorrelation.viewRowSelection(sharedUncorrelatedGenes).duplicate();
				geneZscoresNullGwasNullBetasPathwayMatched = geneZscoresNullGwasNullBetas.viewRowSelection(sharedUncorrelatedGenes).duplicate();
			}

			geneZscoresPathwayMatched.centerColumns();
			geneZscoresNullGwasCorrelationPathwayMatched.centerColumns();
			geneZscoresNullGwasNullBetasPathwayMatched.centerColumns();

			genePathwayZscores.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment_normalizedPathwayScores" + (this.hlaGenesToExclude == null ? "" : "_ExHla") + ".txt");
			geneZscoresPathwayMatched.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment_normalizedGwasGeneScores" + (this.hlaGenesToExclude == null ? "" : "_ExHla") + ".txt");
			geneZscoresNullGwasNullBetasPathwayMatched.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment_normalizedNullGwasGeneScores" + (this.hlaGenesToExclude == null ? "" : "_ExHla") + ".txt");

			LinkedHashMap<String, Integer> singleColMap = new LinkedHashMap<>(1);
			singleColMap.put("B1", 0);

			final List<DoubleMatrixDataset<String, String>> b1PerArm = Collections.synchronizedList(new ArrayList<>(geneChrArmMapping.size()));
			final List<DoubleMatrixDataset<String, String>> b2PerArm = Collections.synchronizedList(new ArrayList<>(geneChrArmMapping.size()));

			final List<DoubleMatrixDataset<String, String>> b1NullPerArm = Collections.synchronizedList(new ArrayList<>(geneChrArmMapping.size()));
			final List<DoubleMatrixDataset<String, String>> b2NullPerArm = Collections.synchronizedList(new ArrayList<>(geneChrArmMapping.size()));

			pb.step();

//for (Map.Entry<String, ArrayList<String>> chrArmMappingEntry : geneChrArmMapping.entrySet()) {
			geneChrArmMapping.entrySet().parallelStream().forEach((Map.Entry<String, ArrayList<Gene>> chrArmMappingEntry) -> {

				try {

					final String chrArm = chrArmMappingEntry.getKey();
					final ArrayList<Gene> armGenes = chrArmMappingEntry.getValue();

					final ArrayList<String> chrArmGenesInPathwayMatrix = new ArrayList<>(armGenes.size());
					final ArrayList<Gene> armGenesInPathwayMatrix = new ArrayList<>(armGenes.size());

					for (Gene armGene : armGenes) {
						if (genePathwayZscores.containsRow(armGene.getGene())) {
							chrArmGenesInPathwayMatrix.add(armGene.getGene());
							armGenesInPathwayMatrix.add(armGene);
						}
					}
					//Now genesInPathwayMatrix will only contain genes that are also in the gene p-value matrix

					LOGGER.debug("Number of genes in chr arm: " + chrArmGenesInPathwayMatrix.size());

					if (chrArmGenesInPathwayMatrix.isEmpty()) {
						throw new RuntimeException("This should not happen");
					}

					//b1 rows: traits cols: 1
					final DoubleMatrixDataset<String, String> b1Arm = new DoubleMatrixDataset<>(geneZscores.getHashCols(), singleColMap);
					//b2 rows: traits cols: pathways
					final DoubleMatrixDataset<String, String> b2Arm = new DoubleMatrixDataset<>(geneZscores.getHashCols(), genePathwayZscores.getHashCols());

					final DoubleMatrixDataset<String, String> b1NullGwasArm = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), singleColMap);
					final DoubleMatrixDataset<String, String> b2NullGwasArm = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), genePathwayZscores.getHashCols());

					b1PerArm.add(b1Arm);
					b2PerArm.add(b2Arm);

					b1NullPerArm.add(b1NullGwasArm);
					b2NullPerArm.add(b2NullGwasArm);

					final DoubleMatrixDataset<String, String> geneZscoresSubset = geneZscoresPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset = geneZscoresNullGwasCorrelationPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetasSubset = geneZscoresNullGwasNullBetasPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
					final DoubleMatrixDataset<String, String> genePathwayZscoresSubset = genePathwayZscores.viewRowSelection(chrArmGenesInPathwayMatrix);

					final DoubleMatrixDataset<String, String> geneZscoresNullGwasSubsetGeneCorrelations = createLocalGeneCorrelation(geneZscoresNullGwasCorrelationSubset, armGenesInPathwayMatrix, 2000000);

					//geneZscoresNullGwasSubsetGeneCorrelations = geneZscoresNullGwasCorrelationSubset.viewDice().calculateCorrelationMatrix();
					geneZscoresNullGwasSubsetGeneCorrelations.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneCor.txt");

					final DoubleMatrix2D geneInvCorMatrixSubsetMatrix;
					try {
						if (this.ignoreGeneCorrelations) {
							LOGGER.debug("Ignoring gene correlations");
							geneInvCorMatrixSubsetMatrix = DoubleFactory2D.dense.identity(geneZscoresNullGwasSubsetGeneCorrelations.rows());
						} else {
							LOGGER.debug("Calculating correlation inverse");
							geneInvCorMatrixSubsetMatrix = new DenseDoubleAlgebra().inverse(geneZscoresNullGwasSubsetGeneCorrelations.getMatrix());
						}
					} catch (Exception ex) {
						LOGGER.fatal(pathwayDatabase.getName() + " " + chrArm);
						throw ex;
					}

//final DoubleMatrixDataset<String, String> geneInvCorMatrixSubset = new DoubleMatrixDataset<>(geneInvCorMatrixSubsetMatrix, geneZscoresNullGwasNullBetasSubset.getHashRows(), geneZscoresNullGwasNullBetasSubset.getHashRows());
//geneInvCorMatrixSubset.save(outputBasePath + "_" + pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneInvCor.txt");
					glsStep1(geneZscoresSubset, geneInvCorMatrixSubsetMatrix, genePathwayZscoresSubset, b1Arm, b2Arm);
					glsStep1(geneZscoresNullGwasNullBetasSubset, geneInvCorMatrixSubsetMatrix, genePathwayZscoresSubset, b1NullGwasArm, b2NullGwasArm);

					pb.step();

				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}

			});

//b1 rows: traits cols: 1
			final DoubleMatrixDataset<String, String> b1 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), singleColMap);
//b2 rows: traits cols: pathways
			final DoubleMatrixDataset<String, String> b2 = new DoubleMatrixDataset<>(geneZscores.getHashCols(), genePathwayZscores.getHashCols());

			final DoubleMatrixDataset<String, String> b1NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), singleColMap);
			final DoubleMatrixDataset<String, String> b2NullGwas = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetas.getHashCols(), genePathwayZscores.getHashCols());

//Combine the results per arm can be done in parrale over the 4 differnt matrices
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

			betas = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
			betasNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());

			numberOfPathways = genePathwayZscores.columns();
			final int numberTraits = geneZscores.columns();

			for (int traitI = 0; traitI < numberTraits; ++traitI) {
				final double b1Trait = b1.getElementQuick(traitI, 0);
				for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

					double beta = b2.getElementQuick(traitI, pathwayI) / b1Trait;
					betas.setElementQuick(pathwayI, traitI, beta);

				}
			}

			final int numberTraitsNull = geneZscoresNullGwasNullBetas.columns();

			for (int traitI = 0; traitI < numberTraitsNull; ++traitI) {
				final double b1Trait = b1NullGwas.getElementQuick(traitI, 0);
				for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

					double beta = b2NullGwas.getElementQuick(traitI, pathwayI) / b1Trait;
					betasNull.setElementQuick(pathwayI, traitI, beta);

				}
			}

			betas.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_betas" : "_betasExHla") + ".txt");
			betasNull.save(outputBasePath + "_" + pathwayDatabase.getName() + "_EnrichmentNull" + (this.hlaGenesToExclude == null ? "_betas" : "_betasExHla") + ".txt");

			pb.step();
		}
	}

	public DoubleMatrixDataset<String, String> getEnrichmentZscores() throws IOException {

		if (enrichmentPvalues == null) {

			enrichmentPvalues = betas.duplicate();

			final DoubleMatrix2D zscoreMatrix = enrichmentPvalues.getMatrix();
			final DoubleMatrix2D betasNullMatrix = betasNull.getMatrix();

			final int numberOfPathways = zscoreMatrix.rows();
			final int numberOfPhenotypes = zscoreMatrix.columns();
			final int numberOfNullGwasPhenotypes = betasNullMatrix.columns();
			final double numberOfNullGwasPhenotypesMin1Double = betasNullMatrix.columns() - 1;

			LOGGER.debug("numberOfNullGwasPhenotypes: " + numberOfNullGwasPhenotypes);

			List<String> pathwayNames;
			if (LOGGER.isDebugEnabled()) {
				pathwayNames = betasNull.getRowObjects();
			} else {
				pathwayNames = Collections.emptyList();
			}

			for (int r = 0; r < numberOfPathways; ++r) {

				double meanNull = 0;
				for (int p = 0; p < numberOfNullGwasPhenotypes; ++p) {
					meanNull += betasNullMatrix.getQuick(r, p);
				}
				meanNull /= numberOfNullGwasPhenotypes;

				double x = 0;
				for (int p = 0; p < numberOfNullGwasPhenotypes; ++p) {
					x += (betasNullMatrix.getQuick(r, p) - meanNull) * (betasNullMatrix.getQuick(r, p) - meanNull);
				}
				double sdNull = Math.sqrt(x / numberOfNullGwasPhenotypesMin1Double);

				if (LOGGER.isDebugEnabled()) {
					LOGGER.debug(pathwayNames.get(r) + " mean: " + meanNull + " sd: " + sdNull);
				}

				for (int c = 0; c < numberOfPhenotypes; ++c) {

					final double corr = zscoreMatrix.getQuick(r, c);

					zscoreMatrix.setQuick(r, c, (corr - meanNull) / sdNull);

				}
			}

			enrichmentPvalues.save(outputBasePath + "_" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla") + ".txt");
		}

		return enrichmentPvalues;

	}

	private static DoubleMatrix2D getPseudoInverseOfSquareMatrix(DoubleMatrix2D a) {

		final DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(a, true, true);

		double[] s = svd.getSingularValues();
		final int numberOfComponents = s.length;

		int firstNonPositive;
		for (firstNonPositive = 0; firstNonPositive < numberOfComponents; ++firstNonPositive) {
			if (s[firstNonPositive] <= 0) {
				break;
			}
		}

		s = Arrays.copyOfRange(s, 0, firstNonPositive);

		final DoubleMatrix2D v = svd.getV().viewPart(0, 0, numberOfComponents, firstNonPositive);
		final DoubleMatrix2D ut = svd.getU().viewPart(0, 0, numberOfComponents, firstNonPositive).viewDice();

		for (int r = 0; r < s.length; r++) {
			final double x = 1 / s[r];
			for (int c = 0; c < firstNonPositive; c++) {
				ut.setQuick(r, c, x * ut.getQuick(r, c));
			}
		}

		return v.zMult(ut, null);

	}

	private static void glsStep1(DoubleMatrixDataset<String, String> geneZscoresSubset, DoubleMatrix2D geneInvCorMatrix, DoubleMatrixDataset<String, String> genePathwayZscoresSubset, DoubleMatrixDataset<String, String> b1, DoubleMatrixDataset<String, String> b2) {

		final int numberOfGenes = geneZscoresSubset.rows();
		final int numberTraits = geneZscoresSubset.columns();
		final int numberOfPathways = genePathwayZscoresSubset.columns();

		//final DoubleMatrix2D geneInvCorMatrix = geneInvCorMatrixSubset.getMatrix();
		final DoubleMatrix2D genePathwayZscoresMatrix = genePathwayZscoresSubset.getMatrix();

		//result of transpose geneZscoresTrait times inv correlation matrix
		DoubleMatrix2D A = geneZscoresSubset.getMatrix().like(1, numberOfGenes);

		for (int traitI = 0; traitI < numberTraits; ++traitI) {

			try {

				DoubleMatrix2D geneZscoresTrait = geneZscoresSubset.viewColAsMmatrix(traitI);
				geneZscoresTrait.zMult(geneInvCorMatrix, A, 1, 0, true, false);

				final double x = A.viewRow(0).zDotProduct(geneZscoresTrait.viewColumn(0));

				//Col order should be the same
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

	private static DoubleMatrixDataset<String, String> createLocalGeneCorrelation(final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset, final ArrayList<Gene> genes, final int correlationWindow) {

		if(genes.size() != geneZscoresNullGwasCorrelationSubset.rows()){
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
					Gene geneI = genes.get(i);
					Gene geneJ = genes.get(j);

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

}
