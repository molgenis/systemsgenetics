/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import eqtlmappingpipeline.Main;
import eqtlmappingpipeline.normalization.Normalizer;
import gnu.trove.map.hash.THashMap;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.molgenis.genotype.Allele;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.binInteraction.*;
import umcg.genetica.io.binInteraction.gene.BinaryInteractionGeneCreator;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.*;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author harm-jan Multi-threaded implementation of the OLS model
 */
public class InteractionAnalysisMultiThreaded {

    private int nrInOutput;

    public void prepareDataForCelltypeSpecificEQTLMapping(String inexpraw, String outdirectory, Double correlationThreshold, String celltypeSpecificProbeFile, String mdsComponentFile, String cellCountFile, String gte, Integer threads) throws Exception {
        String rawExpressionDataFile = inexpraw;
        Normalizer n = new Normalizer();
        if (correlationThreshold == null) {
            correlationThreshold = 0.9;
        }

        if (rawExpressionDataFile == null || rawExpressionDataFile.trim().length() == 0 || !Gpio.exists(rawExpressionDataFile)) {
            throw new IllegalArgumentException("Error: Raw gene expression file: " + rawExpressionDataFile + "  either does not exist or was not provided to the program.");
        }
        if (outdirectory == null || outdirectory.trim().length() == 0) {
            throw new IllegalArgumentException("Error: output directory not provided");
        }

        if (Math.abs(correlationThreshold) > 1 || Math.abs(correlationThreshold) < 0) {
            throw new IllegalArgumentException("Error: PC1 sample correlation threshold should be between 0 and 1");
        }

        if (celltypeSpecificProbeFile == null || celltypeSpecificProbeFile.trim().length() == 0 || !Gpio.exists(celltypeSpecificProbeFile)) {
            throw new IllegalArgumentException("Error: Cell type specific probe list has not been provided or does not exist: " + celltypeSpecificProbeFile);
        }

        if (mdsComponentFile == null || mdsComponentFile.trim().length() == 0 || !Gpio.exists(mdsComponentFile)) {
            System.err.println("Warning: will not correct for possible population stratification effects!");
            mdsComponentFile = null;
        }

        // create the output directory
        outdirectory = Gpio.formatAsDirectory(outdirectory);

        Gpio.createDir(outdirectory);
        String expressionOutputDirectory = outdirectory + "ExpressionData/";
        Gpio.createDir(expressionOutputDirectory);

        // read genotype to expression coupling for removal of gene expression samples not linked to genotypes
        HashSet<String> expressionSamplestoInclude = null;
        if (gte != null) {
            System.out.println("Loading genotype to expression coupling: " + gte);
            expressionSamplestoInclude = new HashSet<String>();
            TextFile tf = new TextFile(gte, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length > 1) {
                    expressionSamplestoInclude.add(elems[1]);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            System.out.println("Your genotype to expression coupling file has: " + expressionSamplestoInclude.size() + " individuals.");
        }

        // 7. select Cell type specific probes
        System.out.println("Loading list of cell type specific probes from: " + celltypeSpecificProbeFile);
        HashSet<String> cellTypeSpecificProbeSet = new HashSet<String>();
        TextFile cellSpecificProbeTF = new TextFile(celltypeSpecificProbeFile, TextFile.R);
        cellTypeSpecificProbeSet.addAll(cellSpecificProbeTF.readAsArrayList());
        cellSpecificProbeTF.close();

        if (cellTypeSpecificProbeSet.isEmpty()) {
            System.err.println("Error: " + celltypeSpecificProbeFile + " is empty!");
            System.exit(-1);
        } else {
            System.out.println(cellTypeSpecificProbeSet.size() + " cell type specific probes loaded.");
        }

        // 1. load gene expression data
        System.out.println("Loading gene expression data.");
        DoubleMatrixDataset<String, String> rawExpressionDataset = DoubleMatrixDataset.loadSubsetOfTextDoubleData(rawExpressionDataFile, '\t', null, expressionSamplestoInclude);
        //new DoubleMatrixDataset<String, String>(rawExpressionDataFile, null, expressionSamplestoInclude);
        // double[][] rawExpressionData = rawExpressionDataset.getMatrix().toArray();

        // determine the number of cell type specific probes in this dataset
        int probeCounter = 0;
        List<String> probes = rawExpressionDataset.getRowObjects();
        for (int i = 0; i < probes.size(); i++) {
            if (cellTypeSpecificProbeSet.contains(probes.get(i))) {
                probeCounter++;
            }
        }

        if (probeCounter == 0) {
            System.err.println("Error: none of the cell type specific probes defined in " + celltypeSpecificProbeFile + " are present in expression dataset: " + rawExpressionDataFile);
            System.exit(-1);
        } else {
            System.out.println(probeCounter + " of the cell type specific probes are in your dataset.");
        }

        // if we have filtered the gene expression data for those samples not having genotypes, save the file now.
        if (expressionSamplestoInclude != null) {
            rawExpressionDataset.save(expressionOutputDirectory + "ExpressionDataForSamplesWithGenotypes.txt.gz");
        }

        // 2. QN + Log2 transform
        QuantileNormalization.quantilenormalize(rawExpressionDataset);
        Log2Transform.log2transform(rawExpressionDataset);

        // Correct for population stratification
        if (mdsComponentFile != null) {
            String file = n.adjustCovariates(rawExpressionDataset, expressionOutputDirectory + "ExpressionDataRaw-QNormLog2Transformed", mdsComponentFile, 0);
            System.out.println("MDS component corrected file: " + file + ".txt.gz");
            rawExpressionDataset = DoubleMatrixDataset.loadDoubleData(file + ".txt.gz");
        } else {
            rawExpressionDataset.save(expressionOutputDirectory + "ExpressionDataRaw-QNormLog2Transformed.txt.gz");
        }

        rawExpressionDataset = rawExpressionDataset.viewDice(); // put the samples on the rows

        //Set gene expression mean and standard deviation to 0 and 1, respectively, to speed up sample correlation matrix calculation (covariance = correlation matrix in this case):
        for (int i = 0; i < rawExpressionDataset.rows(); i++) {
            double[] row = rawExpressionDataset.getRow(i).toArray();
            double mean = Descriptives.mean(row);
            double var = Descriptives.variance(row);
            double sd = Math.sqrt(var);

            for (int j = 0; j < rawExpressionDataset.columns(); j++) {
                double v = (rawExpressionDataset.getElementQuick(i, j) - mean) / sd;
                rawExpressionDataset.setElementQuick(i, j, v);
            }
        }

        System.out.println("Will now determine the sample correlation matrix");
        // 3. Sample correlation matrix:
        ConcurrentCorrelation correlator = null;
        if (threads != null) {
            correlator = new ConcurrentCorrelation(threads);
        } else {
            correlator = new ConcurrentCorrelation();
        }

        double[][] sampleCorrelationMatrix = correlator.pairwiseCorrelation(rawExpressionDataset.getMatrix().toArray());

        DoubleMatrixDataset<String, String> sampleCorrelationMatrixOut = new DoubleMatrixDataset<String, String>();
        sampleCorrelationMatrixOut.setMatrix(sampleCorrelationMatrix);
        sampleCorrelationMatrixOut.setColObjects(rawExpressionDataset.getRowObjects());
        sampleCorrelationMatrixOut.setRowObjects(rawExpressionDataset.getRowObjects());
        sampleCorrelationMatrixOut.save(outdirectory + "SampleCorrelationMatrix.txt");

        // 4. PCA on sample correlation matrix
        rawExpressionDataset = rawExpressionDataset.viewDice(); // put samples back on columns
        // this method returns two DoubleMatrixDatasets: left are the PC scores, right are the Eigenvalues and expects the samples to be on the columns
        Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = n.calculatePCA(rawExpressionDataset, sampleCorrelationMatrixOut, expressionOutputDirectory + "PCAResults", 2);

        // 5. Correlate samples with PC1 - scores (QC step to determine poor RNA samples)
        // This dataset needs to be transposed if rows are currently PCs, and columns contain samples.
        DoubleMatrixDataset<String, String> pcScores = PCAResults.getLeft();
        pcScores = pcScores.viewDice();
        double[] firstPC = pcScores.getRow(0).toArray();

        // iterate through the samples
        TextFile sampleToPCScoreCorrelationOut = new TextFile(outdirectory + "SamplePC1Correlations.txt", TextFile.W);
        // transfer individuals to rows..
        rawExpressionDataset = rawExpressionDataset.viewDice();
        List<String> individuals = rawExpressionDataset.getRowObjects();

        HashSet<String> individualsPassingQC = new HashSet<String>();
        sampleToPCScoreCorrelationOut.writeln("Individual\tSpearmanCorrelationWithPC1\tPearsonCorrelationWithPC1");
        for (int i = 0; i < individuals.size(); i++) {
            String individual = individuals.get(i);
            double[] x = rawExpressionDataset.getRow(i).toArray();
            double[] y = firstPC;
            double pearson = Correlation.correlate(x, y);
            double spearman = Correlation.rankCorrelate(x, y);
            if (Math.abs(pearson) > correlationThreshold) {
                individualsPassingQC.add(individual);
            } else {
                System.out.println(individual + "\tDid not pass QC. Correlation with PC1: " + Math.abs(pearson));
            }
            sampleToPCScoreCorrelationOut.writeln(individual + "\t" + spearman + "\t" + pearson);
        }
        sampleToPCScoreCorrelationOut.close();

        if (individualsPassingQC.size() < rawExpressionDataset.rows() * 0.1) {
            System.err.println("Error: QC method includes less than 10% of your samples (" + individualsPassingQC.size() + "/" + rawExpressionDataset.rows() + "). There may be something wrong with your data! Please contact us!");
            System.exit(-1);
        } else {
            System.out.println("QC method includes " + individualsPassingQC.size() + " out of " + rawExpressionDataset.rows() + " samples.");
        }

        // clear some memory
        pcScores = null;
        PCAResults = null;
        individuals = null;

        String rawExpressionFileToUseForNextStep = null;

        System.out.println("Now reloading the gene expression data for the samples that passed the QC.");
        // 6. Remove samples with r < 0.9 for PC1
        // reload expression file, include only samples that pass QC...
        rawExpressionDataset = DoubleMatrixDataset.loadSubsetOfTextDoubleData(rawExpressionDataFile, '\t', null, individualsPassingQC); // new DoubleMatrixDataset<String, String>(rawExpressionDataFile, null, individualsPassingQC);


        // quantile normalize, log2 transform again, because the number of samples might have been changed..
        QuantileNormalization.quantilenormalize(rawExpressionDataset);
        Log2Transform.log2transform(rawExpressionDataset);

        if (mdsComponentFile != null) {
            System.out.println("Correcting for MDS components..");
            String file = n.adjustCovariates(rawExpressionDataset, expressionOutputDirectory + "ExpressionDataSamplePCQC-QNormLog2Transform", mdsComponentFile, 0);
            System.out.println("MDS component corrected file: " + file + ".txt.gz");
            rawExpressionDataset = DoubleMatrixDataset.loadDoubleData(file + ".txt.gz");
            rawExpressionFileToUseForNextStep = file + ".txt.gz";
        } else {
            rawExpressionDataset.save(outdirectory + "ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.txt.gz");
            rawExpressionFileToUseForNextStep = outdirectory + "ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.txt.gz";
        }

        // collect data for cell type specific probes
        DoubleMatrixDataset<String, String> cellTypeSpecificDataset = new DoubleMatrixDataset<>(probeCounter, rawExpressionDataset.columns());
        probeCounter = 0;
        ArrayList<String> cellTypeSpecificProbeDatasetRowNames = new ArrayList<String>();
        for (int i = 0; i < probes.size(); i++) {
            if (cellTypeSpecificProbeSet.contains(probes.get(i))) {
                for (int c = 0; c < rawExpressionDataset.columns(); c++) {
                    cellTypeSpecificDataset.setElementQuick(probeCounter, c, rawExpressionDataset.getElementQuick(i, c));
                }
                cellTypeSpecificProbeDatasetRowNames.add(probes.get(i));
                probeCounter++;
            }
        }

        // initiate cell type specific probe correlation matrix
        DoubleMatrixDataset<String, String> celltypeSpecificCorrelationMatrix = new DoubleMatrixDataset<>(probeCounter, probeCounter);
        for (int i = 0; i < probeCounter; i++) {
            for (int j = i + 1; j < probeCounter; j++) {
                double r = Correlation.correlate(cellTypeSpecificDataset.getRow(i).toArray(), cellTypeSpecificDataset.getRow(j).toArray());
                celltypeSpecificCorrelationMatrix.setElementQuick(i, j, r);
                celltypeSpecificCorrelationMatrix.setElementQuick(j, i, r);
            }
            celltypeSpecificCorrelationMatrix.setElementQuick(i, i, 1);
        }

        // save the correlation matrix
        celltypeSpecificCorrelationMatrix.setColObjects(cellTypeSpecificProbeDatasetRowNames);
        celltypeSpecificCorrelationMatrix.setRowObjects(cellTypeSpecificProbeDatasetRowNames);
        celltypeSpecificCorrelationMatrix.save(outdirectory + "CelltypeSpecificProbeCorrelationMatrix.txt.gz");

        // 9. PCA over cell specific probe correlation matrix
        cellTypeSpecificDataset.setColObjects(rawExpressionDataset.getColObjects());
        cellTypeSpecificDataset.setRowObjects(cellTypeSpecificProbeDatasetRowNames);
        cellTypeSpecificDataset.save(expressionOutputDirectory + "CellTypeSpecificProbeExpression.txt.gz");
        cellTypeSpecificDataset = cellTypeSpecificDataset.viewDice();

        // calculate first Principal Component over the cell type specific probe matrix...
        PCAResults = n.calculatePCA(cellTypeSpecificDataset, celltypeSpecificCorrelationMatrix, outdirectory + "CellTypeSpecificProbePCA", cellTypeSpecificProbeDatasetRowNames.size());

        // 10. PC1 scores: cell specific proxy -- write to file for future use...
        DoubleMatrixDataset<String, String> cellSpecificPCScores = PCAResults.getLeft();

        //Ensure that the cellTypeSpecificPCScores correlate positively with the set of probes that we have used to determine this component:
        double[] pcScoresSamples = new double[cellSpecificPCScores.rows()];
        for (int i = 0; i < cellSpecificPCScores.rows(); i++) {
            pcScoresSamples[i] = cellSpecificPCScores.getElementQuick(i, 0);
        }
        cellTypeSpecificDataset = cellTypeSpecificDataset.viewDice();
        int nrProbesCorrelatingPositively = 0;
        for (int i = 0; i < cellTypeSpecificDataset.rows(); i++) {
            double corr = JSci.maths.ArrayMath.correlation(pcScoresSamples, cellTypeSpecificDataset.getRow(i).toArray());
            if (corr >= 0) {
                nrProbesCorrelatingPositively++;
            } else {
                nrProbesCorrelatingPositively--;
            }
        }
        if (nrProbesCorrelatingPositively < 0) {
            for (int i = 0; i < cellSpecificPCScores.rows(); i++) {
                cellSpecificPCScores.setElementQuick(i, 0, -cellSpecificPCScores.getElementQuick(i, 0));
            }
        }

        TextFile tfOutCellSpecific = new TextFile(outdirectory + "CellTypeProxyFile.txt", TextFile.W);
        tfOutCellSpecific.writeln("Sample\tCellCountProxyValue");
        for (int i = 0; i < cellSpecificPCScores.rows(); i++) {
            tfOutCellSpecific.writeln(cellSpecificPCScores.getRowObjects().get(i) + "\t" + cellSpecificPCScores.getElementQuick(i, 0));
        }
        tfOutCellSpecific.close();

        // 11. Correlate PC1 scores with cell counts (if any)
        if (cellCountFile != null) {
            HashMap<String, Double> cellCounts = new HashMap<String, Double>();
            TextFile cellcountfile = new TextFile(cellCountFile, TextFile.R);
            cellcountfile.readLine();
            String[] elems = cellcountfile.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length > 1) {
                    String sample = elems[0];
                    try {
                        Double d = Double.parseDouble(elems[1]);
                        cellCounts.put(sample, d);
                    } catch (NumberFormatException e) {
                        System.err.println("Error parsing number in " + cellCountFile + ": " + elems[1]);
                    }
                }
                elems = cellcountfile.readLineElems(TextFile.tab);
            }
            cellcountfile.close();

            if (cellCounts.isEmpty()) {
                System.err.println("ERROR: none of the cell counts in " + cellCountFile + " could be parsed.");

            } else {

                ArrayList<Double> x = new ArrayList<Double>();
                ArrayList<Double> y = new ArrayList<Double>();

                for (int i = 0; i < cellSpecificPCScores.rows(); i++) {
                    String sample = cellSpecificPCScores.getRowObjects().get(i);
                    if (cellCounts.containsKey(sample)) {
                        x.add(cellSpecificPCScores.getElementQuick(i, 0));
                        y.add(cellCounts.get(sample));
//                    System.out.println(sample + "\t" + cellSpecificPCScores.rawData[i][0] + "\t" + cellCounts.get(sample));
                    }
                }

                double[] xArr = toPrimitiveArr(x.toArray(new Double[0]));
                double[] yArr = toPrimitiveArr(y.toArray(new Double[0]));

                double r = Correlation.correlate(xArr, yArr);
                SpearmansCorrelation corr = new SpearmansCorrelation();
                double spearman = corr.correlation(xArr, yArr);
                for (int q = 0; q < xArr.length; q++) {
                    System.out.println(q + "\t" + xArr[q] + "\t" + yArr[q]);
                }

                ScatterPlot plot = new ScatterPlot(500, 500, xArr, yArr, ScatterPlot.OUTPUTFORMAT.PDF, outdirectory + "plot.pdf");

//                plot.draw(xArr, yArr, "Cell type specific PC Scores", "Cell counts", "Comparison between cell counts and predicted cell counts", outdirectory + "Scatterplot.png");
                TextFile tfout = new TextFile(outdirectory + "ComparisonToCellCount.txt", TextFile.W);
                System.out.println("Correlation between actual cell counts and PC1 scores: " + r + "\tr2: " + (r * r) + "\tn: " + xArr.length);
                tfout.writeln("Pearson\tSpearman\tn");
                tfout.writeln(r + "\t" + spearman + "\t" + xArr.length);
                tfout.close();
            }
        }

        System.out.println("");
        System.out.println("PLEASE NOTE:");
        System.out.println("For the next step, you can use the following file as raw expression data (--inexpraw): " + rawExpressionFileToUseForNextStep);
        System.out.println("For the cell count proxy file, please use the following file for the next step (--cellcounts): " + tfOutCellSpecific.getFileName());
        System.out.println("");

    }

    public void runInteractionAnalysis(String inExpPCCorrected, String covariateFile, String ingt,
                                       String gte, String snpprobecombinationfile, Integer nrThreads, String out,
                                       String covariateList, boolean forceNormalDistribution, boolean robustSE, boolean fullStats, boolean binaryOutput, String cohort) throws IOException, Exception {
        String probeannot = null;

        double mafthreshold = 0.05;
        double hwepthreshold = 0.001;
        double crthreshold = 0.95;

        if (snpprobecombinationfile == null || !Gpio.exists(snpprobecombinationfile)) {
            throw new IllegalArgumentException("ERROR: please provide snpprobe combination file");
        }

        if (robustSE) {
            System.out.println("Running tests for robust standard errors. Now testing R connection");
            try {
                RConnection rConnection = new RConnection();
                System.out.println("R server found: " + rConnection.getServerVersion());
                rConnection.close();
            } catch (RserveException ex) {
                System.err.println(ex.getMessage());
                System.err.println("Could not connect to RServe");
                System.exit(-1);
            }
        }

        out = Gpio.formatAsDirectory(out);
        Gpio.createDir(out);

        // read SNP-probe combinations
        HashSet<Pair<String, String>> snpprobeCombos = null;
        TextFile tf = new TextFile(snpprobecombinationfile, TextFile.R);
        snpprobeCombos = tf.readAsPairs(0, 1);
        tf.close();

        if (snpprobeCombos.isEmpty()) {
            System.err.println("Error: no SNP-probe combinations loaded from file: " + snpprobecombinationfile);
            System.exit(-1);
        } else {
            System.out.println(snpprobeCombos.size() + " SNP-Probe combinations loaded from: " + snpprobecombinationfile);
        }

        HashSet<String> includeTheseIndividuals = null;

        // load dataset
        System.out.println("Now loading eQTL dataset.");
        TriTyperGeneticalGenomicsDatasetSettings settings = new TriTyperGeneticalGenomicsDatasetSettings();
        settings.cisAnalysis = true;
        settings.transAnalysis = true;
        settings.expressionLocation = inExpPCCorrected;
        settings.expressionplatform = null;
        settings.genotypeLocation = ingt;
        settings.genotypeToExpressionCoupling = gte;
        settings.logtransform = false;
        settings.quantilenormalize = false;
        settings.name = "Dataset";
        settings.probeannotation = probeannot;

        TriTyperGeneticalGenomicsDataset ds = new TriTyperGeneticalGenomicsDataset(settings);
        TriTyperGenotypeData genotypeData = ds.getGenotypeData();
        TriTyperExpressionData pcCorrectedExpressionData = ds.getExpressionData();

        THashMap<String, String> gteHash = ds.getGenotypeToExpressionCouplings();
        HashSet<String> expressionIndividualsInPCCorrectedData = new HashSet<String>();
        for (String genotypeSample : genotypeData.getIndividuals()) {
            if (gteHash.containsKey(genotypeSample)) {
                if (includeTheseIndividuals == null || includeTheseIndividuals.contains(gteHash.get(genotypeSample))) {
                    expressionIndividualsInPCCorrectedData.add(gteHash.get(genotypeSample));
                }
            }
        }

        // load the same individuals from the raw data....
        System.out.println("Now loading covariate data file");

        Set<String> covariateHash = null;
        if (covariateList != null) {
            TextFile tfcovariatelist = new TextFile(covariateList, TextFile.R);
            covariateHash = tfcovariatelist.readAsSet(0, TextFile.tab);
            tfcovariatelist.close();
        }
        DoubleMatrixDataset<String, String> covariateData = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covariateFile, '\t', covariateHash, expressionIndividualsInPCCorrectedData);

        // check whether the samples are on the columns...
        int ctr = 0;
        for (String s : covariateData.getColObjects()) {
            if (expressionIndividualsInPCCorrectedData.contains(s)) {
                ctr++;
            }
        }
        if (ctr == 0) {
            // try to load the data again, but transpose it
            covariateData = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covariateFile, '\t', covariateHash, expressionIndividualsInPCCorrectedData);
            covariateData = covariateData.viewDice();
            ctr = 0;
            for (String s : covariateData.getColObjects()) {
                if (expressionIndividualsInPCCorrectedData.contains(s)) {
                    ctr++;
                }
            }
            if (ctr == 0) {
                System.err.println("Error: covariate sample identifiers don't match up with those in expression data");
                System.exit(-1);
            }
        }
        System.out.println(ctr + " gene expression samples have covariates.");

        // since the number of samples has changed, we might need to reperform q-norm and log2 transform...
        // it may be a good idea to remove these last steps from the normalization step..


        // investigate which SNPs to run..
        LinkedHashSet<Pair<String, String>> snpProbeCombinationsToTest = new LinkedHashSet<Pair<String, String>>();
        HashSet<String> snpsPassingQC = new HashSet<String>();
        HashSet<String> snpsVisited = new HashSet<String>();
        HashMap<String, SNP> snpStats = new HashMap<String, SNP>(); // for the binary output
        SNPLoader loader = genotypeData.createSNPLoader();

        System.out.println("Parsing SNP probe combos");
        TextFile tfOut = new TextFile(out + "eQTLsNotPassingQC.txt", TextFile.W);
        for (Pair<String, String> p : snpprobeCombos) {
            String snp = p.getLeft();
            String probe = p.getRight();
            Integer snpId = genotypeData.getSnpToSNPId().get(snp);

            Integer probeIdInPCCorrectedData = pcCorrectedExpressionData.getProbeToId().get(probe);

            if (snpId != -9 && probeIdInPCCorrectedData != -9) {
                if (snpsPassingQC.contains(snp)) {
                    snpProbeCombinationsToTest.add(p);
                } else if (!snpsVisited.contains(snp)) {
                    // snp has not been seen before.. test for QC parameters.
                    SNP snpObj = genotypeData.getSNPObject(snpId);
                    loader.loadGenotypes(snpObj);
                    if (snpObj.getMAF() >= mafthreshold && snpObj.getHWEP() >= hwepthreshold && snpObj.getCR() >= crthreshold) {
                        snpsPassingQC.add(snp);
                        snpProbeCombinationsToTest.add(p);

                        if (binaryOutput) {
                            snpStats.put(snp, snpObj);
                        }
                    } else {
                        tfOut.writeln(p.toString() + "\tSNP fails QC (MAF/HWEP/CR)\t" + snpObj.getMAF() + "\t" + snpObj.getHWEP() + "\t" + snpObj.getCR());
                    }
                    snpObj.clearGenotypes();
                }
            } else {
                tfOut.writeln(p.toString() + "\tProbe or SNP not on platform\t" + probe + " ID:(" + probeIdInPCCorrectedData + ")\t" + snp + " ID: (" + snpId + ")");
            }
            snpsVisited.add(snp);
        }
        tfOut.close();

        if (snpProbeCombinationsToTest.isEmpty()) {
            System.err.println("None of the specified SNP-probe combinations to test are present in the dataset!");
            System.exit(-1);
        } else {
            System.out.println(snpProbeCombinationsToTest.size() + " eQTLs can be tested in your dataset, using " + covariateData.rows() + " covariates.");
        }

        // make a base cellcount array
        String[] expInds = pcCorrectedExpressionData.getIndividuals();

        out = Gpio.formatAsDirectory(out);
        Gpio.createDir(out);

        ArrayList<String> rowNames = new ArrayList<String>();
        rowNames.addAll(covariateData.getRowObjects());

        Correlation.correlationToZScore(covariateData.columns());

        System.out.println("Output matrix will be " + snpProbeCombinationsToTest.size() + "(x5) x " + rowNames.size());

        double[][] expressiondata = pcCorrectedExpressionData.getMatrix();
        int[] wgaId = ds.getExpressionToGenotypeIdArray();

        if (forceNormalDistribution) {
            System.out.println("Forcing normal distribution on covariate and expression data");
            System.out.println("Warning: normal distribution is forced before covariate samples are matched to genotypes.");
            System.out.println("Make sure that the number of samples between samples and covariates are more or less equal");
            System.out.println("Currently: " + pcCorrectedExpressionData.getColNames().length + " for expression and " + covariateData.columns() + " for covariates");

            Normalizer norm = new Normalizer();

            for (int row = 0; row < expressiondata.length; row++) {
                expressiondata[row] = norm.forceNormal(expressiondata[row]);
            }

            double[][] covariates = covariateData.getMatrix().toArray();
            for (int row = 0; row < covariates.length; row++) {
                covariates[row] = norm.forceNormal(covariates[row]);
            }
            covariateData.setMatrix(covariates);
            System.out.println("Done. And you have been warned.");
        }


        TextFile snpFile = new TextFile(out + "SNPSummaryStatistics.txt", TextFile.W);
        snpFile.writeln("SNP\tChr\tChrPos\tAlleles\tMinorAllele\tMAF\tCallRate\tHWE\tGenotypesCalled");

        String[] snpsPassingQCArr = snpsPassingQC.toArray(new String[0]);
        int nrSubmitted = 0;
        if (nrThreads == null) {
            nrThreads = Runtime.getRuntime().availableProcessors();
        }


        System.out.println("Running with: " + nrThreads + " threads");

        ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
        CompletionService<InteractionAnalysisResults> pool = new ExecutorCompletionService<InteractionAnalysisResults>(threadPool);

        nrInOutput = 0;
        TextFile outputFile = null;
        BinaryInteractionFile binaryInteractionFile = null;

        // Write binary output header
        if (binaryOutput) {
            File binaryOutFile = new File(out + "InteractionResults.binary.dat");
            String description = "Genotypes: " + ingt +
                    " Expresion: " + inExpPCCorrected +
                    " GTE: " + gte +
                    " Covariates: " + covariateFile +
                    " Covariates List: " + covariateList +
                    " SNP-probes: " + snpprobecombinationfile +
                    " Software version: " + Main.VERSION;
            binaryInteractionFile = createBinaryOutputHeader(binaryOutFile, snpsPassingQCArr, snpStats,
                    snpProbeCombinationsToTest, covariateData, expressionIndividualsInPCCorrectedData, cohort, description);
        } else {
            System.out.println("Output will be written to: " + out + "InteractionResults.txt");
            outputFile = new TextFile(out + "InteractionResults.txt", TextFile.W);
            String outputheader = "SNP\tProbe\tCovariate\tZ-SNP\tZ-Cov\tZ-Interaction\tZ-Main\tZ-Interaction-Flipped\tN\tRSquared";
            if (fullStats) {
                outputheader
                        += "\tsnpBeta"
                        + "\tsnpSE"
                        + "\tcovariateBeta"
                        + "\tcovariateSE"
                        + "\tinteractionBeta"
                        + "\tinteractionSE"
                        + "\tinteractionBeta-Flipped";
            }

            outputFile.writeln(outputheader);
        }


        ProgressBar pb = new ProgressBar(snpProbeCombinationsToTest.size(), "Now testing available eQTL effects for interactions.");
        int maxbuffer = (nrThreads * 8);
        for (int i = 0; i < snpsPassingQCArr.length; i++) {
            String snp = snpsPassingQCArr[i];
            ArrayList<Pair<String, String>> eQTLsForSNP = new ArrayList<Pair<String, String>>();

            for (Pair<String, String> eQTL : snpProbeCombinationsToTest) {
                if (eQTL.getLeft().equals(snp)) {
                    eQTLsForSNP.add(eQTL);
                }
            }

            if (eQTLsForSNP.size() > 0) {
                Integer snpId = genotypeData.getSnpToSNPId().get(snp);
                SNP snpObj = genotypeData.getSNPObject(snpId);
                loader.loadGenotypes(snpObj);
                if (loader.hasDosageInformation()) {
                    loader.loadDosage(snpObj);
                }

                // push the actual work to thread..
                InteractionAnalysisTask t = new InteractionAnalysisTask(
                        snpObj,
                        eQTLsForSNP,
                        expressiondata,
                        wgaId,
                        expInds,
                        covariateData,
                        pcCorrectedExpressionData,
                        robustSE,
                        fullStats
                );
                pool.submit(t);
                nrSubmitted++;
            }

            if (nrSubmitted % maxbuffer == 0) {
//                System.out.println("Retrieving results..");
                int nrRetrieved = 0;
                while (nrRetrieved < nrSubmitted) {
                    // pick up results

                    try {
                        InteractionAnalysisResults result = pool.take().get();
                        if (result != null) {
                            if (binaryOutput)
                                binaryInteractionFile = processResultWriteBinaryOutput(result, binaryInteractionFile, snpFile, covariateData, pb, fullStats);
                            else
                                processResult(result, outputFile, snpFile, covariateData, pb, fullStats);
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                        System.exit(1);
                    }

                    nrRetrieved++;
                }

                nrSubmitted = 0;
            }
        }

        pb.close();

        // check if there's still some work to be done..
        if (nrSubmitted > 0) {
            int nrRetrieved = 0;
            while (nrRetrieved < nrSubmitted) {
                // pick up results
                try {
                    InteractionAnalysisResults result = pool.take().get();
                    if (result != null) {
                        if (binaryOutput)
                            binaryInteractionFile = processResultWriteBinaryOutput(result, binaryInteractionFile, snpFile, covariateData, pb, fullStats);
                        else
                            processResult(result, outputFile, snpFile, covariateData, pb, fullStats);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
                nrRetrieved++;
            }

            nrSubmitted = 0;
        }

        threadPool.shutdown();

        snpFile.close();

        if (binaryOutput) {
            binaryInteractionFile.finalizeWriting();
            System.out.println("Interaction results writer buffer flushed: " + binaryInteractionFile.getInteractionWriteBufferFlushed());
            System.out.println("QTL results writer buffer flushed: " + binaryInteractionFile.getQtlWriteBufferFlushed());
            System.out.println("Total number of expected interactions: " + binaryInteractionFile.getTotalNumberInteractions());
            System.out.println("Total number of writen interactions: " + binaryInteractionFile.getInteractionZscoresSet());
            System.out.println("Number of QTL z-scores: " + binaryInteractionFile.getVariantCount());
            binaryInteractionFile.close();

            if (binaryInteractionFile.getInteractionZscoresSet() != binaryInteractionFile.getTotalNumberInteractions()) {
                System.out.println("WARNING!!! written and expected interactions not the same");
                System.err.println("WARNING!!! written and expected interactions not the same");
            }

        } else {
            outputFile.close();
        }
//        datasetOut.colObjects = colNames;
//        datasetOut.recalculateHashMaps();

//        if (binaryOutput) {
//            datasetOut.save(out + "CellTypeSpecificityMatrix.binary");
//        } else {
//            datasetOut.save(out + "CellTypeSpecificityMatrix.txt");
//        }
        loader.close();

        System.out.println("Done.");
    }

    private void processResult(InteractionAnalysisResults result, TextFile outputFile, TextFile snpFile, DoubleMatrixDataset<String, String> covariateData, ProgressBar pb, boolean fullStats) throws IOException {
        double[][] interactionZScoreMatrix = result.getInteractionZScoreMatrix();
        double[][] SNPZResultMatrix = result.getSNPZResultMatrix();
        double[][] covariateZResultMatrix = result.getCovariateZResultMatrix();
        double[][] maineffectZResultMatrix = result.getMaineffectZResultMatrix();
        int[][] nMatrix = result.getnMatrix();
        ArrayList<Pair<String, String>> eqtls = result.geteQTLsTested();

        double[][] covariateBeta = result.getCovariateBeta();
        double[][] covariateSE = result.getCovariateSE();
        double[][] interactionBeta = result.getInteractionBeta();
        double[][] interactionSE = result.getInteractionSE();
        double[][] mainBeta = result.getSNPBeta();
        double[][] mainSE = result.getSNPSE();
        double[][] rsquared = result.getRsquared();

        for (int e = 0; e < eqtls.size(); e++) {
            Pair<String, String> eqtl = eqtls.get(e);
            for (int c = 0; c < SNPZResultMatrix[e].length; c++) {

                StringBuilder builder = new StringBuilder();
                builder.append(eqtl.getLeft());
                builder.append("\t");
                builder.append(eqtl.getRight());
                builder.append("\t");
                builder.append(covariateData.getRowObjects().get(c));

                builder.append("\t");
                builder.append(SNPZResultMatrix[e][c]);
                builder.append("\t");
                builder.append(covariateZResultMatrix[e][c]);
                builder.append("\t");
                double interactionZ = interactionZScoreMatrix[e][c];
                builder.append(interactionZ);
                builder.append("\t");
                double mainZ = maineffectZResultMatrix[e][c];
                builder.append(mainZ);

                if (mainZ < 0) {
                    interactionZ *= -1;
                }
                builder.append("\t");
                builder.append(interactionZ);
                builder.append("\t");
                builder.append(nMatrix[e][c]);

                builder.append("\t");
                builder.append(rsquared[e][c]);

                if (fullStats) {
                    builder.append("\t");
                    builder.append(mainBeta[e][c]);
                    builder.append("\t");
                    builder.append(mainSE[e][c]);
                    builder.append("\t");
                    builder.append(covariateBeta[e][c]);
                    builder.append("\t");
                    builder.append(covariateSE[e][c]);
                    builder.append("\t");
                    double interactionB = interactionBeta[e][c];
                    builder.append(interactionB);
                    builder.append("\t");
                    builder.append(interactionSE[e][c]);

                    if (mainZ < 0) {
                        interactionB *= -1;
                    }
                    builder.append("\t");
                    builder.append(interactionB);

                }

                outputFile.writeln(builder.toString());
            }
            nrInOutput++;
            pb.iterate();
        }
        snpFile.writeln(result.getQcString());
    }

    private BinaryInteractionFile processResultWriteBinaryOutput(InteractionAnalysisResults result, BinaryInteractionFile createdInteractions, TextFile snpFile, DoubleMatrixDataset<String, String> covariateData, ProgressBar pb, boolean fullStats) throws IOException, BinaryInteractionFileException {

        double[][] interactionZScoreMatrix = result.getInteractionZScoreMatrix();
        double[][] SNPZResultMatrix = result.getSNPZResultMatrix();
        double[][] covariateZResultMatrix = result.getCovariateZResultMatrix();
        double[][] maineffectZResultMatrix = result.getMaineffectZResultMatrix();
        int[][] nMatrix = result.getnMatrix();
        double[][] rsquared = result.getRsquared();
        ArrayList<Pair<String, String>> eqtls = result.geteQTLsTested();

        for (int e = 0; e < eqtls.size(); e++) {
            Pair<String, String> eqtl = eqtls.get(e);
            int numSamples = nMatrix[e][0];
            String snp = eqtl.getLeft();
            String gene = eqtl.getRight();

            //main effect z-score
            double mainZ = maineffectZResultMatrix[e][0];
            BinaryInteractionQtlZscores qtlZscore = new BinaryInteractionQtlZscores(new double[]{mainZ}, new int[]{numSamples});
            createdInteractions.setQtlResults(snp, gene, qtlZscore);
            for (int c = 0; c < SNPZResultMatrix[e].length; c++) {
                String covariate = covariateData.getRowObjects().get(c);

                //interaction z-scores
                double interactionZ = interactionZScoreMatrix[e][c];
                final int[] samplesInteractionCohort = {(Double.isNaN(interactionZ) ? 0 : numSamples)};
                final double[] zscoreSnpCohort = {SNPZResultMatrix[e][c]};
                final double[] zscoreCovariateCohort = {covariateZResultMatrix[e][c]};
                final double[] zscoreInteractionCohort = {interactionZ};
                final double[] rSquaredCohort = {rsquared[e][c]};
                if (mainZ < 0) {
                    interactionZ *= -1;
                }
                final double[] zscoreInteractionFlippedCohort = {interactionZ};

                BinaryInteractionZscores interactionZscores = new BinaryInteractionZscores(samplesInteractionCohort, zscoreSnpCohort, zscoreCovariateCohort, zscoreInteractionCohort, rSquaredCohort, zscoreInteractionFlippedCohort);
                createdInteractions.setInteractionResults(snp, gene, covariate, interactionZscores);
            }
            pb.iterate();
        }
        return createdInteractions;
    }


    private BinaryInteractionFile createBinaryOutputHeader(File binaryOutFile, String[] snpsPassingQCArr, HashMap<String, SNP> snpStats, LinkedHashSet<Pair<String, String>> snpProbeCombinationsToTest, DoubleMatrixDataset<String, String> covariateData, HashSet<String> expressionIndividualsInPCCorrectedData, String cohort, String description) throws BinaryInteractionFileException, IOException {
        LinkedHashSet<String> geneIds = new LinkedHashSet<String>();
        System.out.println("snpProbeCombinationsToTest size: " + snpProbeCombinationsToTest.size());
        for (Pair<String, String> snpProbePair : snpProbeCombinationsToTest) {
            String gene = snpProbePair.getRight();
            geneIds.add(gene);
        }

        int numSNPs = snpsPassingQCArr.length;
        int numGenes = geneIds.size();

        //fill variants
        BinaryInteractionVariantCreator[] variants = new BinaryInteractionVariantCreator[numSNPs];
        for (int snpIdx = 0; snpIdx < numSNPs; snpIdx++) {
            String snpId = snpsPassingQCArr[snpIdx];

            SNP snpObj = snpStats.get(snpId);
            byte[] alleles = snpObj.getAlleles();
            byte minorAllele = snpObj.getMinorAllele();
            byte majorAllele;
            if (alleles[0] == minorAllele)
                majorAllele = alleles[1];
            else
                majorAllele = alleles[0];

            variants[snpIdx] = new BinaryInteractionVariantCreator(snpId, snpObj.getChr() + "", snpObj.getChrPos(), Allele.create((char) majorAllele), Allele.create((char) minorAllele));

        }

        //fill genes
        BinaryInteractionGeneCreator[] genes = new BinaryInteractionGeneCreator[numGenes];
        int geneIdx = 0;
        for (String gene : geneIds) {
            genes[geneIdx] = new BinaryInteractionGeneCreator(gene);
            geneIdx++;
        }


        //fill covariates
        String[] covariates = covariateData.getRowObjects().toArray(new String[0]);

        //fill cohort
        int numSamples = 0;
        for (String s : expressionIndividualsInPCCorrectedData) {
            if (covariateData.getHashCols().containsKey(s))
                numSamples++;
        }
        BinaryInteractionCohort[] cohorts = new BinaryInteractionCohort[1];
        cohorts[0] = new BinaryInteractionCohort(cohort, numSamples);

        // initialize
        BinaryInteractionFileCreator creator = new BinaryInteractionFileCreator(binaryOutFile, variants, genes, cohorts, covariates, true, false, true, true);

        creator.setDescription(description);

        for (Pair<String, String> eqtl : snpProbeCombinationsToTest) {
            creator.addTestedVariantGene(eqtl.getLeft(), eqtl.getRight());
        }
        BinaryInteractionFile createdInteractions = creator.create();
        return createdInteractions;
    }

    private double[] toPrimitiveArr(Double[] toArray) {
        double[] arr = new double[toArray.length];
        for (int i = 0; i < toArray.length; i++) {
            arr[i] = toArray[i];
        }
        return arr;
    }
}
