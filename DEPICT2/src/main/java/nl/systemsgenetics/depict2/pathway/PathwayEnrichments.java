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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.depict2.Depict2;
import nl.systemsgenetics.depict2.gene.Gene;
import nl.systemsgenetics.depict2.gene.GenePathwayAssociationStatistic;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.log4j.Logger;
import umcg.genetica.graphics.panels.HistogramPanel;
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
    private final DoubleMatrixDataset<String, String> standardErrors;
    private final DoubleMatrixDataset<String, String> tStatistics;
    private final DoubleMatrixDataset<String, String> pValues;
    private final DoubleMatrixDataset<String, String> betasNull;
    private DoubleMatrixDataset<String, String> enrichmentPvalues = null;
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
                              final DoubleMatrixDataset<String, String> geneZscores,
                              final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelation,
                              final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetas,
                              final String outputBasePath,
                              final HashSet<String> hlaGenesToExclude,
                              final boolean ignoreGeneCorrelations,
                              final double genePruningR,
                              final int geneCorrelationWindow,
                              final File debugFolder,
                              final File intermediateFolder,
                              final boolean quantileNormalizePermutations) throws Exception {

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

        // Determine final set of genes to analyze
        Set<String> pathwayGenes = pathwayMatrixLoader.getOriginalRowMap().keySet();
        sharedGenes = new LinkedHashSet<>();

        for (String gene : genesWithPvalue) {
            if (pathwayGenes.contains(gene) && !excludeGenes.contains(gene)) {
                sharedGenes.add(gene);
            }
        }


        // TODO: with a little cleanup, this can be moved to Depict2.java, saves doing this each time for each pathway enrichment
        // Determine (log10) gene lengths
        LOGGER.info("Determining gene lengths");
        final double[] geneLengths = new double[geneZscores.getRowObjects().size()];
        int i = 0;
        // TODO: very ugly
        // Ensures the geneLength vector is in the same order as the matrices
        for (String geneInZscoreMatrix : geneZscores.getRowObjects()) {
            for (Gene geneInfo : genes) {
                if (geneInfo.getGene().equals(geneInZscoreMatrix)) {
                    geneLengths[i] = Math.log(geneInfo.getLength()) / Math.log(10);
                    i++;
                    break;
                }
            }
        }

        // Do the regression with gene lengths and zscores
        LOGGER.info("Regressing gene lengths and determining residuals");
        inplaceDetermineGeneLengthRegressionResiduals(geneZscores, geneLengths);
        inplaceDetermineGeneLengthRegressionResiduals(geneZscoresNullGwasCorrelation, geneLengths);
        inplaceDetermineGeneLengthRegressionResiduals(geneZscoresNullGwasNullBetas, geneLengths);
        LOGGER.info("Done regressing gene lengths and determining residuals");

        final HashMap<String, HashSet<MetaGene>> metaGenesPerArm;
        {
            DoubleMatrixDataset<String, String> tmp = geneZscoresNullGwasCorrelation.viewRowSelection(sharedGenes).duplicate();
            tmp.normalizeColumns();
            metaGenesPerArm = groupCorrelatedGenesPerChrArm(tmp, genePruningR, genes, sharedGenes);
        }

        try (ProgressBar pb = new ProgressBar(pathwayDatabase.getName() + " enrichment analysis", metaGenesPerArm.size() + 2, ProgressBarStyle.ASCII)) {
            // Collapse genes that are highly correlated and merge them to meta-genes
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

            final DoubleMatrixDataset<String, String> geneZscoresPathwayMatched;
            final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationPathwayMatched;
            final DoubleMatrixDataset<String, String> geneZscoresNullGwasNullBetasPathwayMatched;

            if (forceNormalGenePvalues) {
                LOGGER.debug("Doing force normal gene p-values");
                geneZscoresPathwayMatched = collapseDatasetToMetaGenes(geneZscores, false, metaGenesPerArm.values()).createColumnForceNormalDuplicate();
                geneZscoresNullGwasCorrelationPathwayMatched = collapseDatasetToMetaGenes(geneZscoresNullGwasCorrelation, false, metaGenesPerArm.values()).createColumnForceNormalDuplicate();
                geneZscoresNullGwasNullBetasPathwayMatched = collapseDatasetToMetaGenes(geneZscoresNullGwasNullBetas, false, metaGenesPerArm.values()).createColumnForceNormalDuplicate();
            } else {
                geneZscoresPathwayMatched = collapseDatasetToMetaGenes(geneZscores, false, metaGenesPerArm.values());
                geneZscoresNullGwasCorrelationPathwayMatched = collapseDatasetToMetaGenes(geneZscoresNullGwasCorrelation, false, metaGenesPerArm.values());
                geneZscoresNullGwasNullBetasPathwayMatched = collapseDatasetToMetaGenes(geneZscoresNullGwasNullBetas, false, metaGenesPerArm.values());
            }

            //TODO: Is this appropriate on the this type of half normal dist?
            //TODO: also this double normalizes if forceNormalGenePvalues is true
            genePathwayZscores.normalizeColumns();
            geneZscoresPathwayMatched.normalizeColumns();
            geneZscoresNullGwasCorrelationPathwayMatched.normalizeColumns();
            geneZscoresNullGwasNullBetasPathwayMatched.normalizeColumns();

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
            final List<DoubleMatrixDataset<String, String>> inverseCorrelationMatrices = Collections.synchronizedList(new ArrayList<>(metaGenesPerArm.size()));

            // Advance progressbar
            pb.step();

            // Calculate the beta's per chromosome arm in parallel
            // for (Map.Entry<String, ArrayList<String>> chrArmMappingEntry : geneChrArmMapping.entrySet()) {
            metaGenesPerArm.entrySet().parallelStream().forEach((Map.Entry<String, HashSet<MetaGene>> chrArmMappingEntry) -> {
                try {
                    // Determine the meta genes to use
                    final String chrArm = chrArmMappingEntry.getKey();
                    final HashSet<MetaGene> armGenes = chrArmMappingEntry.getValue();
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
                    // if (geneCorrelationWindow < 0) {
                    LOGGER.debug("Creating full correlation matrix for chr arm");
                    geneZscoresNullGwasSubsetGeneCorrelations = geneZscoresNullGwasCorrelationSubset.viewDice().calculateCorrelationMatrix();
                    // } else {
                    // LOGGER.debug("Creating correlation matrix in window: " + geneCorrelationWindow);
                    // geneZscoresNullGwasSubsetGeneCorrelations = createLocalGeneCorrelation(geneZscoresNullGwasCorrelationSubset, armGenesInPathwayMatrix, geneCorrelationWindow);
                    // }
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
                            LOGGER.debug("Ignoring gene correlations");
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
                    inverseCorrelationMatrices.add(geneInvCorMatrixSubset);

                    if (LOGGER.isDebugEnabled()) {
                        geneInvCorMatrixSubset.saveBinary(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneInvCor").getAbsolutePath());
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
            final DoubleMatrixDataset<String, String> fullInverseCorrelationMatrix = mergeCorrelationMatrices(inverseCorrelationMatrices);

            // Make sure the oder of metagenes is the same
            fullInverseCorrelationMatrix.reorderRows(geneZscoresPathwayMatched.getHashRows());
            fullInverseCorrelationMatrix.reorderCols(geneZscoresPathwayMatched.getHashRows());

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

            // Determine final betas and p values
            betas = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            standardErrors = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            tStatistics = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            pValues = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            betasNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());
            numberOfPathways = genePathwayZscores.columns();

            // Should be eq to geneZscoresPathwayMatched.columns();
            final int numberTraits = geneZscores.columns();

            // Betas for actual data
            LOGGER.info("Calculating model coefficients");
            for (int traitI = 0; traitI < numberTraits; ++traitI) {
                final double b1Trait = b1.getElementQuick(traitI, 0);

                for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {
                    // Determine beta
                    double beta = b2.getElementQuick(traitI, pathwayI) / b1Trait;

                    // Determine analytical pvalues
                    GenePathwayAssociationStatistic curStat = determineGlsAnalyticalPvalues(beta,
                            b1Trait,
                            fullInverseCorrelationMatrix.getMatrix(),
                            geneZscoresPathwayMatched.getCol(traitI),
                            genePathwayZscores.getCol(pathwayI));

                    // TODO: Just make this one DoubleMatrixDataset
                    betas.setElementQuick(pathwayI, traitI, beta);
                    standardErrors.setElementQuick(pathwayI, traitI, curStat.getStandardError());
                    tStatistics.setElementQuick(pathwayI, traitI, curStat.gettStatistic());
                    pValues.setElementQuick(pathwayI, traitI, curStat.getAnalyticalPvalue());
                }
            }

            LOGGER.info("Done calculating model coefficients");

            // Betas for null
            LOGGER.info("Calculating null beta's");
            final int numberTraitsNull = geneZscoresNullGwasNullBetas.columns();

            for (int traitI = 0; traitI < numberTraitsNull; ++traitI) {
                final double b1Trait = b1NullGwas.getElementQuick(traitI, 0);
                for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {

                    double beta = b2NullGwas.getElementQuick(traitI, pathwayI) / b1Trait;
                    betasNull.setElementQuick(pathwayI, traitI, beta);

                }
            }
            LOGGER.info("Done calculating null beta's");

            // Write output
            betas.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_betas" : "_betasExHla"));
            standardErrors.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_se" : "_seExHla"));
            pValues.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_analyticalPvals" : "_analyticalPvalsExHla"));
            betasNull.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_EnrichmentNull" + (this.hlaGenesToExclude == null ? "_betas" : "_betasExHla"));

            // Save as txt to avoid having to convert them later
            if (LOGGER.isDebugEnabled()) {
                betas.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_betas.txt" : "_betasExHla.txt"));
                standardErrors.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_se.txt" : "_seExHla.txt"));
                pValues.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_analyticalPvals.txt" : "_analyticalPvalsExHla.txt"));
                fullInverseCorrelationMatrix.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_inverseCorrelationMatrix.txt" : "_inverseCorrelationMatrixExHla.txt"));
                geneZscoresPathwayMatched.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_geneZscoresPathwayMatched.txt" : "_geneZscoresPathwayMatchedExHla.txt"));
                genePathwayZscores.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_genePathwayZscores.txt" : "_genePathwayZscoresExHla.txt"));
            }


            this.getEnrichmentZscores();
            this.clearZscoreCache();

            pb.step();
        }
    }


    private static void inplaceDetermineGeneLengthRegressionResiduals(DoubleMatrixDataset<String, String> geneZscores, double[] geneLengths) {
        IntStream.range(0, geneZscores.getMatrix().columns()).parallel().forEach(c -> {
            // Build the regression model
            SimpleRegression regression = new SimpleRegression();
            for (int r = 0; r < geneZscores.rows(); ++r) {
                regression.addData(geneLengths[r], geneZscores.getElement(r, c));
            }

            // Apply model and inplace convert to residuals
            for (int r = 0; r < geneZscores.rows(); ++r) {
                double predictedY = regression.predict(geneLengths[r]);
                geneZscores.setElementQuick(r, c, geneZscores.getElement(r, c) - predictedY);
            }
        });
    }


    /**
     * Determine the Pvalue, Se and Tstat for a given gls regression. Only valid for cases where the data
     * is centered and scaled and has no intercept.
     * // TODO: Generalize this so an intercept term can be included
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

        return new GenePathwayAssociationStatistic(beta, standardError, tstatistic, pvalue);
    }

    public final DoubleMatrixDataset<String, String> getEnrichmentZscores() throws IOException {

        if (enrichmentPvalues == null) {

            enrichmentPvalues = betas.duplicate();

            final DoubleMatrix2D zscoreMatrix = enrichmentPvalues.getMatrix();
            final DoubleMatrix2D betasNullMatrix = betasNull.getMatrix();

            final int numberOfPathways = zscoreMatrix.rows();
            final int numberOfPhenotypes = zscoreMatrix.columns();
            final int numberOfNullGwasPhenotypes = betasNullMatrix.columns();
            final double numberOfNullGwasPhenotypesMin1Double = betasNullMatrix.columns() - 1;

            LOGGER.debug("numberOfNullGwasPhenotypes: " + numberOfNullGwasPhenotypes);

//			List<String> pathwayNames;
//			if (LOGGER.isDebugEnabled()) {
//				pathwayNames = betasNull.getRowObjects();
//			} else {
//				pathwayNames = Collections.emptyList();
//			}
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

//				if (LOGGER.isDebugEnabled()) {750
//					LOGGER.debug(pathwayNames.get(r) + " mean: " + meanNull + " sd: " + sdNull);
//				}
                for (int c = 0; c < numberOfPhenotypes; ++c) {

                    final double corr = zscoreMatrix.getQuick(r, c);

                    zscoreMatrix.setQuick(r, c, (corr - meanNull) / sdNull);

                }
            }

            enrichmentPvalues.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla"));
        }

        return enrichmentPvalues;

    }

    public final void clearZscoreCache() {
        enrichmentPvalues = null;
    }

    private static void glsStep1(DoubleMatrixDataset<String, String> geneZscoresSubset, DoubleMatrix2D geneInvCorMatrix, DoubleMatrixDataset<String, String> genePathwayZscoresSubset, DoubleMatrixDataset<String, String> b1, DoubleMatrixDataset<String, String> b2) {

        final int numberOfGenes = geneZscoresSubset.rows();
        final int numberTraits = geneZscoresSubset.columns();
        final int numberOfPathways = genePathwayZscoresSubset.columns();
        //final DoubleMatrix2D geneInvCorMatrix = geneInvCorMatrixSubset.getMatrix();
        final DoubleMatrix2D genePathwayZscoresMatrix = genePathwayZscoresSubset.getMatrix();

        // Result of transpose geneZscoresTrait times inv correlation matrix
        DoubleMatrix2D A = geneZscoresSubset.getMatrix().like(1, numberOfGenes);

        // Trait = gene
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

    private static DoubleMatrixDataset<String, String> createLocalGeneCorrelation(final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset, final ArrayList<Gene> genes, final int correlationWindow) {

        if (genes.size() != geneZscoresNullGwasCorrelationSubset.rows()) {
            throw new RuntimeException("Genes should match geneZscoresNullGwasCorrelationSubset");
        }

        final DoubleMatrixDataset<String, String> correlations = new DoubleMatrixDataset<>(geneZscoresNullGwasCorrelationSubset.getHashRows(), geneZscoresNullGwasCorrelationSubset.getHashRows());
        final DoubleMatrix2D correlationMatrix = correlations.getMatrix();
        final int geneCount = geneZscoresNullGwasCorrelationSubset.rows();
        final int nullGwasCount = geneZscoresNullGwasCorrelationSubset.columns();
        DoubleMatrix2D geneZscoresNullGwasCorrelationSubsetMatrix = geneZscoresNullGwasCorrelationSubset.getMatrix();

        final SimpleRegression regression = new SimpleRegression();

        for (int i = geneCount; --i >= 0; ) {
            for (int j = i + 1; --j >= 0; ) {
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

    private static HashMap<String, HashSet<MetaGene>> groupCorrelatedGenesPerChrArm(final DoubleMatrixDataset<String, String> geneZscoresNullGwas, final double maxCorrelationBetweenGenes, final List<Gene> genes, final Set<String> includedGenes) {

        HashMap<String, HashSet<MetaGene>> metaGenes = new HashMap<>();

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
            DoubleMatrixDataset<String, String> correlationOfCorrelations = genePvaluesNullGwasGeneArmCorrelation.calculateCorrelationMatrix();

            rows:
            for (int r = 0; r < correlationOfCorrelations.rows(); ++r) {

                String currentGene = armGenesIds.get(r);

                MetaGene currentMetaGene = new MetaGene(armGenesIds.get(r));
                metaGenesArm.put(currentGene, currentMetaGene);

                cols:
                for (int c = 0; c < r; ++c) {
                    if (Math.abs(correlationOfCorrelations.getElementQuick(r, c)) >= maxCorrelationBetweenGenes) {

                        //Never null because c < r
                        MetaGene otherMetaGene = metaGenesArm.get(armGenesIds.get(0));
                        currentMetaGene.addOtherMetaGene(otherMetaGene);

                        for (String otherGene : otherMetaGene.getGenes()) {
                            metaGenesArm.put(otherGene, currentMetaGene);
                        }

                    }
                }

            }

            HashSet<MetaGene> metaGenesArmUnique = new HashSet<>(metaGenesArm.values());

            synchronized (metaGenes) {
                metaGenes.put(chrArm, metaGenesArmUnique);
            }

        });

        return metaGenes;
    }

    private static class MetaGene {

        final HashSet<String> genes = new HashSet<>();

        public MetaGene(String gene) {
            genes.add(gene);
        }

        public HashSet<String> getGenes() {
            return genes;
        }

        public void addGene(String gene) {
            genes.add(gene);
        }

        public void addOtherMetaGene(MetaGene other) {
            genes.addAll(other.getGenes());
        }

        public int getGeneCount() {
            return genes.size();
        }

        public String getMetaGeneId() {
            return String.join("_", genes);
        }

    }

    /**
     * @param dataset with genes on rows
     */
    private static DoubleMatrixDataset<String, String> collapseDatasetToMetaGenes(final DoubleMatrixDataset<String, String> dataset, final boolean zscoreSum, final Collection<HashSet<MetaGene>> metaGenes) {

        LinkedHashMap<String, Integer> metaGenesRows = new LinkedHashMap<>();

        int r = 0;
        for (HashSet<MetaGene> metaGenesPerArm : metaGenes) {
            for (MetaGene metaGene : metaGenesPerArm) {
                metaGenesRows.put(metaGene.getMetaGeneId(), r++);
            }
        }

        final int cols = dataset.columns();
        final DoubleMatrixDataset datasetCollapsed = new DoubleMatrixDataset<>(metaGenesRows, dataset.getHashCols());

        for (HashSet<MetaGene> metaGenesPerArm : metaGenes) {
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
     * Merges a list of correlation matrices so they are on the diagonal, and the rest is padded with 0
     * Matrix needs to be square
     *
     * @param correlationMatrices list of correlation matrices as DoubleMatrixDatasets
     * @return One correlation matrix with the off diagonal parts not covered by a matrix as 0
     */
    public static DoubleMatrixDataset<String, String> mergeCorrelationMatrices(List<DoubleMatrixDataset<String, String>> correlationMatrices) throws Exception {

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

}
