/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleCholeskyDecomposition;
import cern.colt.matrix.tdouble.algo.decomposition.SparseDoubleCholeskyDecomposition;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.gene.GenePathwayAssociationStatistic;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.log4j.Logger;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import umcg.genetica.graphics.panels.HistogramPanel;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.ZScores;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;


@Deprecated
public class CholeskyBasedPathwayEnrichments  { //extends PathwayEnrichments

    private static final Logger LOGGER = Logger.getLogger(Downstreamer.class);

    private final PathwayDatabase pathwayDatabase;
    private final HashSet<String> hlaGenesToExclude;
    private final boolean ignoreGeneCorrelations;
    private final DoubleMatrixDataset<String, String> betas;
    private final DoubleMatrixDataset<String, String> standardErrors;
    private final DoubleMatrixDataset<String, String> pValues;
    private final DoubleMatrixDataset<String, String> zscores;
    //private final DoubleMatrixDataset<String, String> qValues;
    private final DoubleMatrixDataset<String, String> pValuesNull;
    //private final DoubleMatrixDataset<String, String> betasNull;
    private DoubleMatrixDataset<String, String> enrichmentPvalues = null;
    private final int numberOfPathways;
    private final File intermediateFolder;

    // TODO: can be locals??
    private final Set<String> excludeGenes;
    private final List<Gene> genes;
    private final String outputBasePath;

    public CholeskyBasedPathwayEnrichments(final PathwayDatabase pathwayDatabase,
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
                                           final int geneCorrelationWindow,
                                           final File debugFolder,
                                           final File intermediateFolder,
                                           final boolean quantileNormalizePermutations,
                                           final boolean regressGeneLengths,
                                           DoubleMatrixDataset<String, String> geneMaxSnpZscore,
                                           DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwasCorrelation,
                                           DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwasBetas
    ) throws Exception {
       // super();

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
            LOGGER.info("Determining gene lengths");
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
        final HashMap<String, HashSet<MetaGene>> metaGenesPerArm;
        {
            DoubleMatrixDataset<String, String> tmp = geneZscoresNullGwasCorrelation.duplicate();
            tmp.normalizeColumns();
            metaGenesPerArm = groupCorrelatedGenesPerChrArm(tmp, genePruningR, genes, sharedGenes);
        }

        try (ProgressBar pb = new ProgressBar(pathwayDatabase.getName() + " enrichment analysis", metaGenesPerArm.size() + 6, ProgressBarStyle.ASCII)) {

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


            // Store the inverse gene-gene correlation matrices
            final Map<String, DoubleMatrixDataset<String, String>> correlationMatrices = Collections.synchronizedMap(new HashMap<>(metaGenesPerArm.size()));

            // Advance progressbar
            pb.step();

            // Calculate the correlations per chromosome arm in parallel
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

                    // Subset the full data to get the matched subset for the current chromosome arm
                    final DoubleMatrixDataset<String, String> geneZscoresNullGwasCorrelationSubset = geneZscoresNullGwasCorrelationPathwayMatched.viewRowSelection(chrArmGenesInPathwayMatrix);
                    final DoubleMatrixDataset<String, String> geneZscoresNullGwasSubsetGeneCorrelations;

                    // Make gene-gene correlation matrix
                    LOGGER.debug("Creating full correlation matrix for chr arm");
                    geneZscoresNullGwasSubsetGeneCorrelations = geneZscoresNullGwasCorrelationSubset.viewDice().calculateCorrelationMatrix();

                    if (LOGGER.isDebugEnabled()) {
                        geneZscoresNullGwasSubsetGeneCorrelations.save(new File(debugFolder, pathwayDatabase.getName() + "_" + chrArm + "_Enrichment_geneCor.txt"));
                    }

                    if (this.ignoreGeneCorrelations) {
                        // Identity matrix, i.e. OLS
                        // TODO: re-implement
                        LOGGER.debug("Ignoring geOLSne correlations");
                        throw new UnsupportedOperationException("Not implemented");
                    }
                    correlationMatrices.put(chrArm, geneZscoresNullGwasSubsetGeneCorrelations);
                    // Advance progressbar
                    pb.step();

                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }

            });
            // Closure of parallel computation

            // Combine the inverse gene-gene correlation matrix into one matrix
            LOGGER.debug("Merging correlation matrices per arm");

            // TODO: Can be optimized by returning sparse matrix direcl, but does not take that long to do
            final DoubleMatrixDataset<String, String> tmpCorrelation = mergeCorrelationMatrices(correlationMatrices.values());
            // Make sure the oder of metagenes is the same
            tmpCorrelation.reorderRows(geneZscoresPathwayMatched.getHashRows());
            tmpCorrelation.reorderCols(geneZscoresPathwayMatched.getHashRows());

            //final SparseCCDoubleMatrix2D fullCorrelationMatrix = new SparseCCDoubleMatrix2D(tmpCorrelation.rows(), tmpCorrelation.columns());
            //fullCorrelationMatrix.assign(tmpCorrelation.getMatrix());
           // LOGGER.debug(fullCorrelationMatrix.getClass());

            final DoubleMatrix2D inverseCholeskyDecomp = new DenseDoubleAlgebra().inverse(new DenseDoubleCholeskyDecomposition(tmpCorrelation.getMatrix()).getL());
            LOGGER.debug("Done calculating inverse of Cholesky decomp of full correlation matrix");

            LOGGER.debug("Pre-multiplying x and y using inverse of Cholesky decomp");
            final DoubleMatrix2D x = new DoubleMatrixDataset<>(geneZscoresPathwayMatched.getHashRowsCopy(), geneZscoresPathwayMatched.getHashColsCopy()).getMatrix();
            final DoubleMatrix2D y = new DoubleMatrixDataset<>(genePathwayZscores.getHashRowsCopy(), genePathwayZscores.getHashColsCopy()).getMatrix();
            final DoubleMatrix2D xNull = new DoubleMatrixDataset<>(geneZscoresNullGwasNullBetasPathwayMatched.getHashRowsCopy(), geneZscoresNullGwasNullBetasPathwayMatched.getHashColsCopy()).getMatrix();

            LOGGER.debug("x: " + x.rows() + " x " + x.columns());
            LOGGER.debug("y: " + y.rows() + " y " + y.columns());
            LOGGER.debug("xNull: " + xNull.rows() + " xNull " + xNull.columns());

            inverseCholeskyDecomp.zMult(geneZscoresPathwayMatched.getMatrix(), x);
            LOGGER.debug("X done");
            inverseCholeskyDecomp.zMult(genePathwayZscores.getMatrix(), y);
            LOGGER.debug("Y done");
            inverseCholeskyDecomp.zMult(geneZscoresNullGwasNullBetasPathwayMatched.getMatrix(), xNull);
            LOGGER.debug("Xnull done");

            pb.step();

            // Determine final betas and p values
            betas = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            standardErrors = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            pValues = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            zscores = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscores.getHashColsCopy());
            pValuesNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());
            //betasNull = new DoubleMatrixDataset<>(genePathwayZscores.getHashColsCopy(), geneZscoresNullGwasNullBetas.getHashColsCopy());
            numberOfPathways = genePathwayZscores.columns();

            if (LOGGER.isDebugEnabled()) {
                LOGGER.debug("betas: " + betas.rows() + " x " + betas.columns());
                LOGGER.debug("standardErrors: " + standardErrors.rows() + " x " + standardErrors.columns());
                LOGGER.debug("pValues: " + pValues.rows() + " x " + pValues.columns());
                LOGGER.debug("zValues: " + zscores.rows() + " x " + zscores.columns());
                //LOGGER.debug("betasNull: " + betasNull.rows() + " x " + betasNull.columns());
                LOGGER.debug("numberOfPathways: " + numberOfPathways);
            }

            // TODO: Parallelize
            // Now calculate beta's and residuals
            final int numberTraits = geneZscoresPathwayMatched.columns();
            LOGGER.debug("Calculating regression statistics");
            for (int traitI = 0; traitI < numberTraits; ++traitI) {
                for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {
                    // New regression without intercept
                    DoubleMatrix1D curX = x.viewColumn(traitI);
                    DoubleMatrix1D curY = y.viewColumn(pathwayI);
                    SimpleRegression curRegression = new SimpleRegression(false);

                    for (int obs = 0; obs < curX.size(); obs++) {
                        curRegression.addData(curX.get(obs), curY.get(obs));
                    }

                    betas.setElementQuick(pathwayI, traitI, curRegression.getSlope());
                    pValues.setElementQuick(pathwayI, traitI, curRegression.getSignificance());
                    standardErrors.setElementQuick(pathwayI, traitI, curRegression.getSlopeStdErr());

                    // The zscore returned by pToZTwoTailed is always negative, therefore, match direction on beta
                    double zscore = ZScores.pToZTwoTailed(curRegression.getSignificance());
                    if (curRegression.getSlope() > 0) {
                        zscore = -zscore;
                    }
                    zscores.setElementQuick(pathwayI, traitI, zscore);
                }
            }

            // TODO: Parallelize
            // Betas for null
            LOGGER.debug("Calculating null beta's");
            final int numberTraitsNull = geneZscoresNullGwasNullBetas.columns();

            for (int traitI = 0; traitI < numberTraitsNull; ++traitI) {
                for (int pathwayI = 0; pathwayI < numberOfPathways; ++pathwayI) {
                    DoubleMatrix1D curX = xNull.viewColumn(traitI);
                    DoubleMatrix1D curY = y.viewColumn(pathwayI);
                    SimpleRegression curRegression = new SimpleRegression(false);

                    for (int obs = 0; obs < curX.size(); obs++) {
                        curRegression.addData(curX.get(obs), curY.get(obs));
                    }
                    pValuesNull.setElementQuick(pathwayI, traitI, curRegression.getSignificance());
                }
            }
            pb.step();


            // Write output
            betas.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_betas" : "_cholesky_betasExHla"));
            //standardErrors.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_se" : "_cholesky_seExHla"));
            //pValues.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_analyticalPvals" : "_cholesky_analyticalPvalsExHla"));
            //betasNull.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_EnrichmentNull" + (this.hlaGenesToExclude == null ? "_cholesky_betas" : "_cholesky_betasExHla"));
            pValuesNull.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_EnrichmentNull" + (this.hlaGenesToExclude == null ? "_cholesky_analyticalPvals" : "_cholesky_analyticalPvalsExHla"));
            betas.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_betas.txt" : "_cholesky_betasExHla.txt"));
            standardErrors.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_se.txt" : "_cholesky_seExHla.txt"));
            pValues.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_analyticalPvals.txt" : "_cholesky_analyticalPvalsExHla.txt"));

            // Save as txt to avoid having to convert them later
            if (LOGGER.isDebugEnabled()) {
                //fullCorrelationMatrix.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_correlationMatrix.txt" : "_cholesky_inverseCorrelationMatrixExHla.txt"));
                geneZscoresPathwayMatched.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_geneZscoresPathwayMatched.txt" : "_cholesky_geneZscoresPathwayMatchedExHla.txt"));
                genePathwayZscores.save(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_cholesky_genePathwayZscores.txt" : "_cholesky_genePathwayZscoresExHla.txt"));
            }

            this.getEnrichmentZscores();
            this.clearZscoreCache();

            pb.step();
        }
    }


    public DoubleMatrixDataset<String, String> getqValues() {
        throw new NotImplementedException();
    }

    public DoubleMatrixDataset<String, String> getEnrichmentZscores() throws IOException {
        zscores.saveBinary(intermediateFolder.getAbsolutePath() + "/" + pathwayDatabase.getName() + "_Enrichment" + (this.hlaGenesToExclude == null ? "_zscore" : "_zscoreExHla"));
        return zscores;
    }


    public void clearZscoreCache() {
        enrichmentPvalues = null;
    }

    public PathwayDatabase getPathwayDatabase() {
        return pathwayDatabase;
    }

    public int getNumberOfPathways() {
        return numberOfPathways;
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

    private static void inplaceDetermineGeneLengthRegressionResiduals(DoubleMatrixDataset<String, String> geneZscores, double[] geneLengths) {

        if (LOGGER.isDebugEnabled()) {
            LOGGER.debug("Gene length vector: " + geneLengths.length);
            LOGGER.debug("Determining residuals for matrix of size: " + geneZscores.rows() + " x " + geneZscores.columns());

            int maxColIdx = geneZscores.columns();
            if (maxColIdx > 10) {
                maxColIdx = 10;
            }
            geneZscores.viewSelection(
                    geneZscores.getRowObjects().subList(0, 10),
                    geneZscores.getColObjects().subList(0, maxColIdx)).printMatrix();
        }

        IntStream.range(0, geneZscores.columns()).parallel().forEach(c -> {
            // Build the regression model
            SimpleRegression regression = new SimpleRegression();
            for (int r = 0; r < geneZscores.rows(); ++r) {
                regression.addData(geneLengths[r], geneZscores.getElement(r, c));
            }

            // Apply model and inplace convert to residuals
            for (int r = 0; r < geneZscores.rows(); ++r) {
                double predictedY = regression.predict(geneLengths[r]);
                geneZscores.setElementQuick(r, c, (geneZscores.getElement(r, c) - predictedY));
            }
        });

        if (LOGGER.isDebugEnabled()) {
            LOGGER.debug("Determined residuals for matrix of size: " + geneZscores.rows() + " x " + geneZscores.columns());
            int maxColIdx = geneZscores.columns();
            if (maxColIdx > 10) {
                maxColIdx = 10;
            }
            geneZscores.viewSelection(
                    geneZscores.getRowObjects().subList(0, 10),
                    geneZscores.getColObjects().subList(0, maxColIdx)).printMatrix();
        }
    }

    /**
     * Merges a list of correlation matrices so they are on the diagonal, and
     * the rest is padded with 0 Matrix needs to be square
     *
     * @param correlationMatrices list of correlation matrices as
     *                            DoubleMatrixDatasets
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

        for (int c = 0; c < matrix.columns(); ++c) {

            double[] col = matrix.getCol(c).toArray();
            double[] colTie = tieBreaker.getCol(c).toArray();

            double mean = JSci.maths.ArrayMath.mean(col);
            double stdev = JSci.maths.ArrayMath.standardDeviation(col);

            double[] rankedValues = ranking.rank(col, colTie);

            for (int s = 0; s < matrix.rows(); s++) {
                double pValue = (0.5d + rankedValues[s] - 1d) / (double) (rankedValues.length);

                newDataset.setElementQuick(s, c, mean + cern.jet.stat.Probability.normalInverse(pValue) * stdev);
            }

        }

        return newDataset;

    }


}
