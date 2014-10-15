/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import umcg.genetica.graphics.ScatterPlot;
import eqtlmappingpipeline.normalization.Normalizer;
import gnu.trove.map.hash.THashMap;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;

/**
 *
 * @author harm-jan Multi-threaded implementation of the OLS model
 */
public class InteractionAnalysisMultiThreaded {

    private int nrInOutput;

    public void prepareDataForCelltypeSpecificEQTLMapping(String inexpraw, String outdirectory, Double correlationThreshold, String celltypeSpecificProbeFile, String mdsComponentFile, String cellCountFile, String gte, Integer threads) throws IOException {
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
        DoubleMatrixDataset<String, String> rawExpressionDataset = new DoubleMatrixDataset<String, String>(rawExpressionDataFile, null, expressionSamplestoInclude);
        double[][] rawExpressionData = rawExpressionDataset.getRawData();

        // determine the number of cell type specific probes in this dataset
        int probeCounter = 0;
        List<String> probes = rawExpressionDataset.rowObjects;
        for (int i = 0; i < probes.size(); i++) {
            if (cellTypeSpecificProbeSet.contains(probes.get(i))) {
                probeCounter++;
            }
        }

        if (probeCounter == 0) {
            System.err.println("Error: none of the cell type specific probes defined in " + celltypeSpecificProbeFile + " are present in expression dataset: " + rawExpressionDataset.fileName);
            System.exit(-1);
        } else {
            System.out.println(probeCounter + " of the cell type specific probes are in your dataset.");
        }

        // if we have filtered the gene expression data for those samples not having genotypes, save the file now.
        if (expressionSamplestoInclude != null) {
            rawExpressionDataset.save(expressionOutputDirectory + "ExpressionDataForSamplesWithGenotypes.txt.gz");
        }

        // 2. QN + Log2 transform
        QuantileNormalization.quantilenormalize(rawExpressionData);
        Log2Transform.log2transform(rawExpressionData);

        // Correct for population stratification
        if (mdsComponentFile != null) {
            String file = n.adjustCovariates(rawExpressionDataset, expressionOutputDirectory + "ExpressionDataRaw-QNormLog2Transformed", mdsComponentFile, true, 0);
            System.out.println("MDS component corrected file: " + file + ".txt.gz");
            rawExpressionDataset = new DoubleMatrixDataset<String, String>(file + ".txt.gz");
        } else {
            rawExpressionDataset.save(expressionOutputDirectory + "ExpressionDataRaw-QNormLog2Transformed.txt.gz");
        }

        rawExpressionDataset.transposeDataset(); // put the samples on the rows
        rawExpressionData = rawExpressionDataset.getRawData();

        //Set gene expression mean and standard deviation to 0 and 1, respectively, to speed up sample correlation matrix calculation (covariance = correlation matrix in this case):
        for (int i = 0; i < rawExpressionData.length; i++) {
            double mean = Descriptives.mean(rawExpressionData[i]);
            double var = Descriptives.variance(rawExpressionData[i]);
            double sd = Math.sqrt(var);
            for (int j = 0; j < rawExpressionData[i].length; j++) {
                rawExpressionData[i][j] -= mean;
                rawExpressionData[i][j] /= sd;
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

        double[][] sampleCorrelationMatrix = correlator.pairwiseCorrelation(rawExpressionData);

        DoubleMatrixDataset<String, String> sampleCorrelationMatrixOut = new DoubleMatrixDataset<String, String>();
        sampleCorrelationMatrixOut.rawData = sampleCorrelationMatrix;
        sampleCorrelationMatrixOut.colObjects = rawExpressionDataset.rowObjects;
        sampleCorrelationMatrixOut.rowObjects = rawExpressionDataset.rowObjects;
        sampleCorrelationMatrixOut.recalculateHashMaps();
        sampleCorrelationMatrixOut.save(outdirectory + "SampleCorrelationMatrix.txt");

        // 4. PCA on sample correlation matrix
        rawExpressionDataset.transposeDataset(); // put samples back on columns
        // this method returns two DoubleMatrixDatasets: left are the PC scores, right are the Eigenvalues and expects the samples to be on the columns
        Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = n.calculatePCA(rawExpressionDataset, sampleCorrelationMatrix, expressionOutputDirectory + "PCAResults", 1);

        // 5. Correlate samples with PC1 - scores (QC step to determine poor RNA samples)
        // This dataset needs to be transposed if rows are currently PCs, and columns contain samples.
        DoubleMatrixDataset<String, String> pcScores = PCAResults.getLeft();
        pcScores.transposeDataset();
        double[] firstPC = pcScores.rawData[0];

        // iterate through the samples
        TextFile sampleToPCScoreCorrelationOut = new TextFile(outdirectory + "SamplePC1Correlations.txt", TextFile.W);
        // transfer individuals to rows..
        rawExpressionDataset.transposeDataset(); // put the samples on the rows
        rawExpressionData = rawExpressionDataset.getRawData();
        List<String> individuals = rawExpressionDataset.rowObjects;

        HashSet<String> individualsPassingQC = new HashSet<String>();
        sampleToPCScoreCorrelationOut.writeln("Individual\tSpearmanCorrelationWithPC1\tPearsonCorrelationWithPC1");
        for (int i = 0; i < individuals.size(); i++) {
            String individual = individuals.get(i);
            double[] x = rawExpressionData[i];
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

        if (individualsPassingQC.size() < rawExpressionDataset.rowObjects.size() * 0.1) {
            System.err.println("Error: QC method includes less than 10% of your samples (" + individualsPassingQC.size() + "/" + rawExpressionDataset.rowObjects.size() + "). There may be something wrong with your data! Please contact us!");
            System.exit(-1);
        } else {
            System.out.println("QC method includes " + individualsPassingQC.size() + " out of " + rawExpressionDataset.rowObjects.size() + " samples.");
        }

        // clear some memory
        pcScores = null;
        PCAResults = null;
        individuals = null;

        String rawExpressionFileToUseForNextStep = null;

        System.out.println("Now reloading the gene expression data for the samples that passed the QC.");
        // 6. Remove samples with r < 0.9 for PC1
        // reload expression file, include only samples that pass QC...      
        rawExpressionDataset = new DoubleMatrixDataset<String, String>(rawExpressionDataFile, null, individualsPassingQC);
        rawExpressionData = rawExpressionDataset.getRawData();

        // quantile normalize, log2 transform again, because the number of samples might have been changed..
        QuantileNormalization.quantilenormalize(rawExpressionData);
        Log2Transform.log2transform(rawExpressionData);

        if (mdsComponentFile != null) {
            System.out.println("Correcting for MDS components..");
            String file = n.adjustCovariates(rawExpressionDataset, expressionOutputDirectory + "ExpressionDataSamplePCQC-QNormLog2Transform", mdsComponentFile, true, 0);
            System.out.println("MDS component corrected file: " + file + ".txt.gz");
            rawExpressionDataset = new DoubleMatrixDataset<String, String>(file + ".txt.gz");
            rawExpressionFileToUseForNextStep = file + ".txt.gz";
        } else {
            rawExpressionDataset.save(outdirectory + "ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.txt.gz");
            rawExpressionFileToUseForNextStep = outdirectory + "ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.txt.gz";
        }

        rawExpressionData = rawExpressionDataset.rawData;

        // collect data for cell type specific probes
        double[][] probeData = new double[probeCounter][rawExpressionDataset.colObjects.size()];
        probeCounter = 0;
        ArrayList<String> cellTypeSpecificProbeDatasetRowNames = new ArrayList<String>();
        for (int i = 0; i < probes.size(); i++) {
            if (cellTypeSpecificProbeSet.contains(probes.get(i))) {
                probeData[probeCounter] = rawExpressionData[i];
                cellTypeSpecificProbeDatasetRowNames.add(probes.get(i));
                probeCounter++;
            }
        }

        // initiate cell type specific probe correlation matrix
        double[][] celltypeSpecificCorrelationMatrix = new double[probeCounter][probeCounter];
        for (int i = 0; i < probeCounter; i++) {
            for (int j = i + 1; j < probeCounter; j++) {
                double r = Correlation.correlate(probeData[i], probeData[j]);
                celltypeSpecificCorrelationMatrix[i][j] = r;
                celltypeSpecificCorrelationMatrix[j][i] = r;
            }
            celltypeSpecificCorrelationMatrix[i][i] = 1;
        }

        // save the correlation matrix
        DoubleMatrixDataset<String, String> probeCorrelationMatrixOut = new DoubleMatrixDataset<String, String>();
        probeCorrelationMatrixOut.colObjects = cellTypeSpecificProbeDatasetRowNames;
        probeCorrelationMatrixOut.rowObjects = cellTypeSpecificProbeDatasetRowNames;
        probeCorrelationMatrixOut.rawData = celltypeSpecificCorrelationMatrix;
        probeCorrelationMatrixOut.recalculateHashMaps();
        probeCorrelationMatrixOut.save(outdirectory + "CelltypeSpecificProbeCorrelationMatrix.txt.gz");

        // 9. PCA over cell specific probe correlation matrix
        DoubleMatrixDataset<String, String> cellTypeSpecificDataset = new DoubleMatrixDataset<String, String>(probeData);
        cellTypeSpecificDataset.colObjects = rawExpressionDataset.colObjects;
        cellTypeSpecificDataset.rowObjects = cellTypeSpecificProbeDatasetRowNames;
        cellTypeSpecificDataset.save(expressionOutputDirectory + "CellTypeSpecificProbeExpression.txt.gz");
        cellTypeSpecificDataset.transposeDataset();

        // calculate first Principal Component over the cell type specific probe matrix...
        PCAResults = n.calculatePCA(cellTypeSpecificDataset, celltypeSpecificCorrelationMatrix, outdirectory + "CellTypeSpecificProbePCA", 1);

        // 10. PC1 scores: cell specific proxy -- write to file for future use...
        DoubleMatrixDataset<String, String> cellSpecificPCScores = PCAResults.getLeft();

        //Ensure that the cellTypeSpecificPCScores correlate positively with the set of probes that we have used to determine this component:
        double[] pcScoresSamples = new double[cellSpecificPCScores.nrRows];
        for (int i = 0; i < cellSpecificPCScores.nrRows; i++) {
            pcScoresSamples[i] = cellSpecificPCScores.rawData[i][0];
        }
        cellTypeSpecificDataset.transposeDataset();
        int nrProbesCorrelatingPositively = 0;
        for (int i = 0; i < cellTypeSpecificDataset.rawData.length; i++) {
            double corr = JSci.maths.ArrayMath.correlation(pcScoresSamples, cellTypeSpecificDataset.rawData[i]);
            if (corr >= 0) {
                nrProbesCorrelatingPositively++;
            } else {
                nrProbesCorrelatingPositively--;
            }
        }
        if (nrProbesCorrelatingPositively < 0) {
            for (int i = 0; i < cellSpecificPCScores.nrRows; i++) {
                cellSpecificPCScores.rawData[i][0] = -cellSpecificPCScores.rawData[i][0];
            }
        }

        TextFile tfOutCellSpecific = new TextFile(outdirectory + "CellTypeProxyFile.txt", TextFile.W);
        tfOutCellSpecific.writeln("Sample\tCellCountProxyValue");
        for (int i = 0; i < cellSpecificPCScores.nrRows; i++) {
            tfOutCellSpecific.writeln(cellSpecificPCScores.rowObjects.get(i) + "\t" + cellSpecificPCScores.rawData[i][0]);
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

                for (int i = 0; i < cellSpecificPCScores.rowObjects.size(); i++) {
                    String sample = cellSpecificPCScores.rowObjects.get(i);
                    if (cellCounts.containsKey(sample)) {
                        x.add(cellSpecificPCScores.rawData[i][0]);
                        y.add(cellCounts.get(sample));
//                    System.out.println(sample + "\t" + cellSpecificPCScores.rawData[i][0] + "\t" + cellCounts.get(sample));
                    }
                }

                double[] xArr = toPrimitiveArr(x.toArray(new Double[0]));
                double[] yArr = toPrimitiveArr(y.toArray(new Double[0]));

                double r = Correlation.correlate(xArr, yArr);

                for (int q = 0; q < xArr.length; q++) {
                    System.out.println(q + "\t" + xArr[q] + "\t" + yArr[q]);
                }

                ScatterPlot plot = new ScatterPlot(500, 500, xArr, yArr, ScatterPlot.OUTPUTFORMAT.PDF, outdirectory + "plot.pdf");

//                plot.draw(xArr, yArr, "Cell type specific PC Scores", "Cell counts", "Comparison between cell counts and predicted cell counts", outdirectory + "Scatterplot.png");
                TextFile tfout = new TextFile(outdirectory + "ComparisonToCellCount.txt", TextFile.W);
                System.out.println("Correlation between actual cell counts and PC1 scores: " + r + "\tr2: " + (r * r) + "\tn: " + xArr.length);
                tfout.writeln("Correlation between actual cell counts and PC1 scores: " + r + "\tr2: " + (r * r) + "\tn: " + xArr.length);
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
            String covariateList) throws IOException, Exception {
        String probeannot = null;

        double mafthreshold = 0.05;
        double hwepthreshold = 0.001;
        double crthreshold = 0.95;

        if (snpprobecombinationfile == null || !Gpio.exists(snpprobecombinationfile)) {
            throw new IllegalArgumentException("ERROR: please provide snpprobe combination file");
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
        DoubleMatrixDataset<String, String> covariateData = new DoubleMatrixDataset<String, String>(covariateFile, covariateHash, expressionIndividualsInPCCorrectedData);

        // since the number of samples has changed, we might need to reperform q-norm and log2 transform...
        // it may be a good idea to remove these last steps from the normalization step..
        // investigate which SNPs to run..
        HashSet<Pair<String, String>> snpProbeCombinationsToTest = new HashSet<Pair<String, String>>();
        HashSet<String> snpsPassingQC = new HashSet<String>();
        HashSet<String> snpsVisited = new HashSet<String>();
        SNPLoader loader = genotypeData.createSNPLoader();

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
            System.out.println(snpProbeCombinationsToTest.size() + " eQTLs can be tested in your dataset, using " + covariateData.nrRows + " covariates.");
        }

        // make a base cellcount array
        String[] expInds = pcCorrectedExpressionData.getIndividuals();

        out = Gpio.formatAsDirectory(out);
        Gpio.createDir(out);

        ArrayList<String> rowNames = new ArrayList<String>();
        rowNames.addAll(covariateData.rowObjects);
//        if (cellcounts != null) {
//            rowNames.add("CellTypeSNPZScore");
//            rowNames.add("CellTypeZScore");
//            rowNames.add("CellTypeInteractionZScore");
//            rowNames.add("MainEffectZScore");
//        }
        Correlation.correlationToZScore(covariateData.nrCols);

//        DoubleMatrixDataset<String, String> datasetOut = new DoubleMatrixDataset<String, String>(rowNames.size(), snpProbeCombinationsToTest.size());
        System.out.println("Output matrix will be " + snpProbeCombinationsToTest.size() + "(x5) x " + rowNames.size());
//        datasetOut.rowObjects = rowNames;

//        ArrayList<String> colNames = new ArrayList<String>();
        double[][] expressiondata = pcCorrectedExpressionData.getMatrix();
        int[] wgaId = ds.getExpressionToGenotypeIdArray();

        TextFile snpFile = new TextFile(out + "SNPSummaryStatistics.txt", TextFile.W);
        snpFile.writeln("SNP\tChr\tChrPos\tAlleles\tMinorAllele\tMAF\tCallRate\tHWE\tGenotypesCalled");

//        TextFile proxyEffectFile = null;
//        if (cellcounts != null) {
//            proxyEffectFile = new TextFile(out + "CelltypeSpecificEQTLEffects.txt", TextFile.W);
//            proxyEffectFile.writeln("#/#\tSNP\tProbe\tnrCalled\tCorrelation\tanovaFTestP\tbetaInteraction\tseInteraction\ttInteraction\tpValueInteraction\tzScoreInteraction");
//        }
        String[] snpsPassingQCArr = snpsPassingQC.toArray(new String[0]);
        int nrSubmitted = 0;
        if (nrThreads == null) {
            nrThreads = Runtime.getRuntime().availableProcessors();
        }

        ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
        CompletionService<InteractionAnalysisResults> pool = new ExecutorCompletionService<InteractionAnalysisResults>(threadPool);

        nrInOutput = 0;

        TextFile outputFile = new TextFile(out + "InteractionResults.txt", TextFile.W);
        String outputheader = "SNP\tProbe\tCovariate\tZ-SNP\tZ-Cov\tZ-Interaction\tZ-Main\tN";

        outputFile.writeln(outputheader);
        ProgressBar pb = new ProgressBar(snpProbeCombinationsToTest.size(), "Now testing available eQTL effects for cell type specificity.");
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
                InteractionAnalysisTask t = new InteractionAnalysisTask(snpObj, eQTLsForSNP, expressiondata, wgaId, expInds, covariateData, pcCorrectedExpressionData);
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
                            processResult(result, outputFile, snpFile, covariateData, pb);
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
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
                        processResult(result, outputFile, snpFile, covariateData, pb);
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

        outputFile.close();
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

    private void processResult(InteractionAnalysisResults result, TextFile outputFile, TextFile snpFile, DoubleMatrixDataset<String, String> covariateData, ProgressBar pb) throws IOException {
        double[][] interactionZScoreMatrix = result.getInteractionZScoreMatrix();
        double[][] SNPZResultMatrix = result.getSNPZResultMatrix();
        double[][] covariateZResultMatrix = result.getCovariateZResultMatrix();
        double[][] maineffectZResultMatrix = result.getMaineffectZResultMatrix();
        int[][] nMatrix = result.getnMatrix();
        ArrayList<Pair<String, String>> eqtls = result.geteQTLsTested();

        for (int e = 0; e < eqtls.size(); e++) {
            Pair<String, String> eqtl = eqtls.get(e);
            for (int c = 0; c < SNPZResultMatrix[e].length; c++) {
                String outputForeQTL = eqtl.getLeft() + "\t" + eqtl.getRight() + "\t" + covariateData.rowObjects.get(c);

                for (int param = 0; param < 5; param++) {
                    switch (param) {
                        case 0:
                            outputForeQTL += "\t" + SNPZResultMatrix[e][c];
                            break;
                        case 1:
                            outputForeQTL += "\t" + covariateZResultMatrix[e][c];
                            break;
                        case 2:
                            outputForeQTL += "\t" + interactionZScoreMatrix[e][c];
                            break;
                        case 3:
                            outputForeQTL += "\t" + maineffectZResultMatrix[e][c];
                            break;
                        case 4:
                            outputForeQTL += "\t" + nMatrix[e][c];
                            break;
                    }
                }
                outputFile.writeln(outputForeQTL);
            }
            nrInOutput++;
            pb.iterate();
        }
        snpFile.writeln(result.getQcString());
    }

    private double[] toPrimitiveArr(Double[] toArray) {
        double[] arr = new double[toArray.length];
        for (int i = 0; i < toArray.length; i++) {
            arr[i] = toArray[i];
        }
        return arr;
    }
}
