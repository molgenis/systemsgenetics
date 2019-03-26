package eqtlmappingpipeline.normalization;

import JSci.maths.matrices.DoubleMatrix;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.PCA;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix.MatrixHandling;
import umcg.genetica.math.matrix.MatrixTools;
import umcg.genetica.math.stats.*;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;
import umcg.genetica.methylation.ConvertBetaAndMvalues;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * @author harmjan
 */
public class Normalizer {

    //nrIntermediatePCAsOverSamplesToRemoveToOutput = 5
    //nrPCAsOverSamplesToRemove = 100
    public void normalize(String expressionFile, String probeIncludeList, String sampleIncludeList, int nrPCAsOverSamplesToRemove,
                          int nrIntermediatePCAsOverSamplesToRemoveToOutput, String covariatesToRemove, boolean orthogonalizecovariates,
                          boolean useOLSforCovariates, String outdir, boolean runQQNorm, boolean runLog2Transform,
                          boolean runMTransform, boolean runCenterScale, boolean runPCA, boolean adjustCovariates,
                          boolean forceMissingValues, boolean forceReplacementOfMissingValues,
                          boolean forceReplacementOfMissingValues2, boolean treatZerosAsNulls, boolean forceNormalDistribution) throws Exception {

        System.out.println("Running normalization.");
        if (outdir != null) {
            outdir = Gpio.formatAsDirectory(outdir);
            Gpio.createDir(outdir);
        } else {
            if (Gpio.getParentDir(expressionFile) == null) {
                //This happens for relative paths in current dir
                //This happens for relative paths in current dir
                outdir = "";
            } else {
                outdir = Gpio.getParentDir(expressionFile) + Gpio.getFileSeparator();
            }

        }

        String parentDir = Gpio.getParentDir(expressionFile);
        String expressionFileName = Gpio.getFileName(expressionFile);
        if (parentDir == null) {
            parentDir = "";
        }

        if (expressionFileName.contains(".txt.gz")) {
            expressionFileName = expressionFileName.replaceAll(".txt.gz", "");
        } else {
            expressionFileName = expressionFileName.replaceAll(".txt", "");
        }

        String outputFileNamePrefix = outdir + expressionFileName;


        Set<String> s = null;
        if (sampleIncludeList != null) {
            TextFile t = new TextFile(sampleIncludeList, TextFile.R);
            s = new HashSet<String>(t.readAsArrayList());
        }
        Set<String> p = null;
        if (probeIncludeList != null) {
            TextFile t = new TextFile(probeIncludeList, TextFile.R);
            p = new HashSet<String>(t.readAsArrayList());
        }
        DoubleMatrixDataset<String, String> dataset = null;

        if (s != null || p != null) {
            dataset = DoubleMatrixDataset.loadSubsetOfTextDoubleData(expressionFile, '\t', p, s); //new DoubleMatrixDataset<String, String>(expressionFile, p, s);
            //Check if samples are correclty loaded.
            boolean breakAfterCheck = false;
            if (s != null) {
                outputFileNamePrefix = outputFileNamePrefix + ".SampleSelection";
                HashSet<String> tmpNames = new HashSet<String>();
                tmpNames.addAll(dataset.getColObjects());
                tmpNames.addAll(s);
                HashSet<String> missingNames = new HashSet<String>();
                HashSet<String> extraNames = new HashSet<String>();
                for (String colName : tmpNames) {
                    if (!s.contains(colName)) {
                        extraNames.add(colName);
                    }
                    if (dataset.getHashCols().get(colName) == null) {
                        missingNames.add(colName);
                    }
                }
                if (!missingNames.isEmpty()) {
                    System.err.println("\nMatrix does not contains desired columns, please check filtering list.");
                    System.err.println(missingNames.toString() + "\n");
                    breakAfterCheck = true;
                } else if (!extraNames.isEmpty()) {
                    System.err.println("\nMatrix contains unwanted columns, please check filtering list.");
                    System.err.println(extraNames.toString() + "\n");
                    breakAfterCheck = true;
                }
            }
            //Check if probes are correclty loaded.
            if (p != null) {
                outputFileNamePrefix = outputFileNamePrefix + ".ProbeSelection";
                HashSet<String> tmpNames = new HashSet<String>();
                tmpNames.addAll(dataset.getRowObjects());
                tmpNames.addAll(p);
                HashSet<String> missingNames = new HashSet<String>();
                HashSet<String> extraNames = new HashSet<String>();
                for (String rowName : tmpNames) {
                    if (!p.contains(rowName)) {
                        extraNames.add(rowName);
                    }
                    if (dataset.getHashRows().get(rowName) == null) {
                        missingNames.add(rowName);
                    }
                }
                if (!missingNames.isEmpty()) {
                    System.err.println("\nMatrix does not contains desired rows, please check filtering list.");
                    System.err.println(missingNames.toString() + "\n");
                    breakAfterCheck = true;
                } else if (!extraNames.isEmpty()) {
                    System.err.println("\nMatrix contains unwanted rows, please check filtering list.");
                    System.err.println(extraNames.toString() + "\n");
                    breakAfterCheck = true;
                }
            }

//            if(breakAfterCheck){
//                System.exit(-1);
//            }

            dataset.save(outputFileNamePrefix + ".txt.gz");
        } else {
            dataset = DoubleMatrixDataset.loadDoubleTextData(expressionFile, '\t');
        }


        // check for probes with zero variance, if there > 3 samples in the dataset
        if (dataset.columns() > 3) {
            outputFileNamePrefix = removeProbesWithZeroVariance(dataset, outputFileNamePrefix);
        }

        if (runQQNorm) {
            outputFileNamePrefix = quantileNormalize(dataset, outputFileNamePrefix, forceMissingValues, forceReplacementOfMissingValues, forceReplacementOfMissingValues2, treatZerosAsNulls);
        }
        if (runLog2Transform) {
            outputFileNamePrefix = log2transform(dataset, outputFileNamePrefix);
        }
        if (runMTransform) {
            outputFileNamePrefix = mValueTransform(dataset, outputFileNamePrefix);
        }
        if (runCenterScale) {
            outputFileNamePrefix = centerAndScale(dataset, outputFileNamePrefix);
        }

        if (adjustCovariates && covariatesToRemove != null) {
            outputFileNamePrefix = adjustCovariates(dataset, outputFileNamePrefix, covariatesToRemove, useOLSforCovariates, 1E-10);
        }

        if (runPCA) {
            int cores = Runtime.getRuntime().availableProcessors();
            ConcurrentCorrelation c = new ConcurrentCorrelation(cores);

            DoubleMatrixDataset<String, String> cormat = c.pairwiseCorrelation(dataset.viewDice());
            Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = calculatePCA(dataset, cormat, outputFileNamePrefix, null);
            if (nrPCAsOverSamplesToRemove != 0 || nrIntermediatePCAsOverSamplesToRemoveToOutput != 0) {
                correctDataForPCs(dataset, outputFileNamePrefix, nrPCAsOverSamplesToRemove, nrIntermediatePCAsOverSamplesToRemoveToOutput, PCAResults.getLeft(), PCAResults.getRight());
            }
        }

        if (forceNormalDistribution) {
            outputFileNamePrefix = forceNormalDistribution(dataset, outputFileNamePrefix);
        }
    }


    NaturalRanking ranking = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);


    public Triple<String, String, String> calculatePcaOnly(String expressionFile) throws Exception {
        String outdir = Gpio.getParentDir(expressionFile) + Gpio.getFileSeparator();

        String parentDir = Gpio.getParentDir(expressionFile);
        String expressionFileName = Gpio.getFileName(expressionFile);
        if (parentDir == null) {
            parentDir = "";
        }

        if (expressionFileName.contains(".txt.gz")) {
            expressionFileName = expressionFileName.replaceAll(".txt.gz", "");
        } else {
            expressionFileName = expressionFileName.replaceAll(".txt", "");
        }

        DoubleMatrixDataset<String, String> dataset = DoubleMatrixDataset.loadDoubleTextData(expressionFile, '\t');

        String outputFileNamePrefix = outdir + expressionFileName;

        ConcurrentCorrelation c = new ConcurrentCorrelation(2);
        DoubleMatrixDataset<String, String> cormat = c.pairwiseCorrelation(dataset.viewDice());

        calculatePCA(dataset, cormat, outputFileNamePrefix, null);
        return new Triple<String, String, String>(outputFileNamePrefix + ".PCAOverSamplesEigenvectorsTransposed.txt.gz", outputFileNamePrefix + ".PCAOverSamplesEigenvectors.txt.gz", outputFileNamePrefix + ".PCAOverSamplesPrincipalComponents.txt.gz");
    }


    public double[] forceNormal(double[] data) {
        double[] rankedValues = ranking.rank(data);
        for (int s = 0; s < data.length; s++) {
            //Convert the rank to a proportion, with range <0, 1>
            double pValue = (0.5d + rankedValues[s] - 1d) / (double) (rankedValues.length);
            //Convert the pValue to a Z-Score:
            data[s] = cern.jet.stat.Probability.normalInverse(pValue);
        }
        return data;
    }


    public String forceNormalDistribution(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws Exception {
        double[][] rawData = dataset.getMatrix().toArray();
        for (int p = 0; p < dataset.getRowObjects().size(); p++) {
            rawData[p] = forceNormal(rawData[p]);
        }

        DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<>();
        datasetNormalized.setRowObjects(dataset.getRowObjects());
        datasetNormalized.setColObjects(dataset.getColObjects());

        fileNamePrefix += ".ForcedNormal";
        datasetNormalized.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;


    }

    public String quantileNormalize(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix,
                                    boolean forceMissingValues, boolean forceReplacementOfMissingValues,
                                    boolean forceReplacementOfMissingValues2,
                                    boolean treatZerosAsNulls) throws IOException {
        boolean dataContainsNulls = false;
        {
            double[][] rawData = dataset.getMatrix().toArray();

            dataContainsNulls = MatrixTools.containsNaNs(rawData);

            if (treatZerosAsNulls && dataContainsNulls) {
                System.out.println("Warning: Data already contains nulls before treating zeros as nulls.\n Later on it will not be possible to distinguish between those two!");
            }
            if (treatZerosAsNulls) {
                MatrixHandling.ReplaceZerosToNull(rawData);
                dataContainsNulls = MatrixTools.containsNaNs(rawData);
            }

            dataset.setMatrix(rawData);
        }

        if (!dataContainsNulls) {
            QuantileNormalization.quantilenormalize(dataset);
        } else if (forceReplacementOfMissingValues) {
            QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, false, false, false);
        } else if (forceReplacementOfMissingValues2) {
            QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, false, true, false);
        } else if (forceMissingValues && treatZerosAsNulls) {
            QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, true, false, true);
        } else if (forceMissingValues) {
            QuantileNormalization.QuantileNormAdressingNaValuesAfterInitialQN(dataset, true, false, false);
        } else {
            System.out.println("Warning: Your data contains missing values and missing value treatment is not selected.\n"
                    + "If desired please supply additional flag: --forceMissingValues or --forceReplacementOfMissingValues");
            System.exit(0);
        }

        if (treatZerosAsNulls) {
            double[][] rawData = dataset.getMatrix().toArray();
            MatrixHandling.ReplaceNullToZero(rawData);
            dataset.setMatrix(rawData);
        }

        fileNamePrefix += ".QuantileNormalized";
        dataset.save(fileNamePrefix + ".txt.gz");

        return fileNamePrefix;
    }

    public String log2transform(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
        Log2Transform.log2transform(dataset);
        fileNamePrefix += ".Log2Transformed";
        dataset.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;
    }

    public String mValueTransform(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
        double[][] rawData = dataset.getMatrix().toArray();
        ConvertBetaAndMvalues.transformToMvalue(rawData);
        dataset.setMatrix(rawData);
        fileNamePrefix += ".MvalueTransformed";
        dataset.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;
    }

    public String centerAndScale(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
//        double[][] rawData = dataset.getMatrix().toArray();
        System.out.println("Standardizing probe mean");
        IntStream.range(0, dataset.rows()).parallel().forEach(p -> {
            double[] row = dataset.getRow(p).toArray();
            double mean = Descriptives.mean(row);
            //double stdev = Math.sqrt(Descriptives.variance(rawData[p], mean));
            for (int s = 0; s < dataset.columns(); s++) {
                double val = dataset.getElementQuick(p, s) - mean;
                dataset.setElementQuick(p, s, val);
            }
        });

        fileNamePrefix += ".ProbesCentered";
        dataset.save(fileNamePrefix + ".txt.gz");

        System.out.println("- Standardizing sample mean and standard deviation");
        IntStream.range(0, dataset.columns()).parallel().forEach(s -> {
            double[] vals = new double[dataset.rows()];
            for (int p = 0; p < dataset.rows(); p++) {
                vals[p] = dataset.getElementQuick(p, s); //getRawData()[p][s];
            }
            double mean = Descriptives.mean(vals);
            for (int p = 0; p < dataset.rows(); p++) {
                vals[p] -= mean;
            }
            double var = Descriptives.variance(vals, mean);
            double stdev = Math.sqrt(var);
            for (int p = 0; p < dataset.rows(); p++) {
                dataset.setElementQuick(p, s, (vals[p] / stdev));
            }
        });
//		for (int s = 0; s < dataset.columns(); s++) {
//			double[] vals = new double[dataset.rows()];
//			for (int p = 0; p < dataset.rows(); p++) {
//				vals[p] = dataset.getRawData()[p][s];
//			}
//			double mean = Descriptives.mean(vals);
//			for (int p = 0; p < dataset.rows(); p++) {
//				vals[p] -= mean;
//			}
//			double var = Descriptives.variance(vals, mean);
//			double stdev = Math.sqrt(var);
//			for (int p = 0; p < dataset.rows(); p++) {
//				dataset.getRawData()[p][s] = (vals[p] / stdev);
//			}
//		}

        fileNamePrefix += ".SamplesZTransformed";
        dataset.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;
    }


    public static void main(String[] args) {

        Normalizer norm = new Normalizer();
        String ds = "D:\\norm\\GD660.GeneQuantCount-EUR-CPM-TMM.txt.gz";
        DoubleMatrixDataset<String, String> dataset = null;
        try {
            dataset = DoubleMatrixDataset.loadDoubleData(ds);
            String filenameprefix = "D:\\norm\\tmp\\test";
            String covariatefile = "D:\\norm\\tmp\\covariates.txt";
            norm.adjustCovariates(dataset, filenameprefix, covariatefile, true, 0);
        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    public String adjustCovariates(DoubleMatrixDataset<String, String> traitData,
                                   String fileNamePrefix,
                                   String covariatesToRemove,
                                   double varianceExplainedCutoff) throws IOException, Exception {
        return adjustCovariates(traitData, fileNamePrefix, covariatesToRemove, false, varianceExplainedCutoff);
    }

    public String adjustCovariates(DoubleMatrixDataset<String, String> traitData,
                                   String fileNamePrefix,
                                   String covariatesToRemove,
                                   boolean useOLS,
                                   double varianceExplainedCutoff) throws IOException, Exception {

        // load covariate data, and remove samples for which there is missing covariate data.
        Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> covariateData = loadCovariateValues(covariatesToRemove, traitData);
        DoubleMatrixDataset<String, String> covariateDataset = covariateData.getLeft();
        traitData = covariateData.getRight();


        if (useOLS) {
            // use Apache multivariate OLS in stead of PCA for covariate correction

            double[][] covariateDataMatrix = new double[traitData.columns()][covariateDataset.rows()];

            int[] sampleMap = new int[traitData.columns()];
            for (int s = 0; s < covariateDataset.getColObjects().size(); s++) {
                String sample = covariateDataset.getColObjects().get(s);
                Integer id = traitData.getHashCols().get(sample);
                sampleMap[s] = id;
            }

            for (int row = 0; row < covariateDataset.rows(); row++) {
                for (int col = 0; col < covariateDataset.columns(); col++) {
                    Integer sampleid = sampleMap[col];
                    covariateDataMatrix[sampleid][row] = covariateDataset.getElementQuick(row, col);
                }
            }

            // TODO: variance inflation factor
            System.out.println("Calculating OLS residuals...");
            DoubleMatrixDataset<String, String> outputmat = new DoubleMatrixDataset<>(traitData.rows(), traitData.columns());
            DoubleMatrixDataset<String, String> finalTraitData = traitData;
            IntStream.range(0, traitData.rows()).parallel().forEach(row -> {
                double[] y = finalTraitData.getRow(row).toArray();
                OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
                ols.newSampleData(y, covariateDataMatrix);
                double[] yout = ols.estimateResiduals();

                for (int c = 0; c < yout.length; c++) {
                    outputmat.setElementQuick(row, c, yout[c]);
                }
            });

            traitData.setMatrix(outputmat.getMatrix());
            fileNamePrefix += ".CovariatesRemovedOLS";
            traitData.save(fileNamePrefix + ".txt.gz");


            return fileNamePrefix;
        } else {

            // use PCA based approach
//        double[][] covariateValues = null;
            double[] pcaExpVar = null;

            System.out.println("Covariate data has " + covariateDataset.rows() + " rows and " + covariateDataset.columns() + " columns.");

            // quick Z-transform
            for (int p = 0; p < covariateDataset.rows(); p++) {
                double[] row = covariateDataset.getRow(p).toArray();
                double mean = Descriptives.mean(row);
                double stdev = Math.sqrt(Descriptives.variance(row, mean));
                for (int s = 0; s < covariateDataset.columns(); s++) {
                    double replacement = (row[s] - mean) / stdev;
                    covariateDataset.setElementQuick(p, s, replacement);
                }
            }

            //Covariation on a centered and scaled matrix equals the correlation.
            //Covariation is faster to compute.
            DoubleMatrixDataset<String, String> correlationMatrix = new DoubleMatrixDataset(covariateDataset.rows(), covariateDataset.rows());
            DoubleMatrixDataset<String, String> finalCovariateDataset = covariateDataset;
            IntStream.range(0, covariateDataset.rows()).parallel().forEach(row -> {
                double[] rowdataA = finalCovariateDataset.getRow(row).toArray();
                for (int row2 = row + 1; row2 < finalCovariateDataset.rows(); row2++) {
                    double[] rowdataB = finalCovariateDataset.getRow(row2).toArray();
                    double r = Correlation.correlate(rowdataA, rowdataB);
                    correlationMatrix.setElementQuick(row, row2, r);
                    correlationMatrix.setElementQuick(row2, row, r);
                }
                correlationMatrix.setElementQuick(row, row, 1);
            });

            DoubleMatrixDataset<String, String> covariateDatasetTransposed = covariateDataset.viewDice();
            Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = calculatePCA(
                    covariateDatasetTransposed,
                    correlationMatrix, covariatesToRemove,
                    null);

            // replace covariateValues with orthogonal ones...
            covariateDataset = PCAResults.getLeft();
            covariateDataset = covariateDataset.viewDice();
//        covariateValues = covariateDataset.getRawData();

            System.out.println(covariateDataset.rows() + " covariates finally loaded.");

            // load the eigenvalues
            pcaExpVar = new double[covariateDataset.rows()];
            System.out.println("Loading eigenvalues from: " + covariatesToRemove + ".PCAOverSamplesEigenvalues.txt.gz");
            TextFile tf = new TextFile(covariatesToRemove + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.R); //
            // skip header
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length > 2) {
                    int pcanr = Integer.parseInt(elems[0]);
                    double expvar = Double.parseDouble(elems[1]);
                    pcaExpVar[pcanr - 1] = expvar;
                    System.out.println(pcanr + "\t" + expvar);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
//        } else {
//            // PCA has been performed a-priori. Just check whether the user has supplied proper covariates.
//            if (covariateValues.length > 1) {
//                // check whether the covariates are orthogonal, by calculating the sum of products (inner product)
//                System.out.println("Determining whether covariates are orthogonal, since you defined > 1 covariate:");
//                System.out.println("Covariate1\tCovariate2\tInnerProduct\tCorrelation");
//                double dotproductthreshold = 1E-5;
//                if (covariateValues.length < 100) {
//                    dotproductthreshold = 0.05;
//                }
//                for (int i = 0; i < covariateValues.length; i++) {
//
//                    for (int j = i + 1; j < covariateValues.length; j++) {
//                        double dotproduct = 0;
//                        for (int v = 0; v < covariateValues[i].length; v++) {
//                            dotproduct += covariateValues[i][v] * covariateValues[j][v];
//                        }
//                        double corr = JSci.maths.ArrayMath.correlation(covariateValues[i], covariateValues[j]);
//
//                        if (Math.abs(dotproduct) > dotproductthreshold) {
//                            System.out.println("Innerproduct > 1E-5 for covariates " + covariateDataset.rowObjects.get(i) + " and " + covariateDataset.rowObjects.get(j) + ", InnerProduct: " + Math.abs(dotproduct) + "\tCorrelation: " + corr);
//                            System.out.println("If you want, we can orthogonalize the covariates for you: use --covpca in your command line.");
//                            System.exit(0);
//                        }
//
//                        System.out.println(covariateDataset.rowObjects.get(i) + "\t" + covariateDataset.rowObjects.get(j) + "\t" + dotproduct + "\t" + corr);
//                    }
//
//                }
//
//                System.out.println("Covariates are orthogonal. Now adjusting for covariates.");
//            }
//        }


//        double[][] rawdata = traitData.getRawData();
            for (int i = 0; i < covariateDataset.rows(); i++) {
                if (pcaExpVar == null || pcaExpVar[i] > varianceExplainedCutoff) {
                    correctForCovariate(traitData, covariateDataset, i);
                } else {
                    System.out.println("Not regressing covariate: " + i + " because explained variance < " + varianceExplainedCutoff + ": " + pcaExpVar[i]);
                }
            }

            fileNamePrefix += ".CovariatesRemoved";
            traitData.save(fileNamePrefix + ".txt.gz");


            return fileNamePrefix;
        }


    }

//    /**
//     * Calculate correlation over columns in DoubleMatrixDataset. WARNING: this
//     * method assumes that SD == 1 and mean == 0 (which makes the covariance
//     * equal to the correlation).
//     *
//     * @param dataset
//     * @return
//     */
//    private double[][] correlateSamples(DoubleMatrixDataset<String, String> dataset) {
//        double[][] correlationMatrix = new double[dataset.columns()][dataset.columns()];
//        double probeCountMinusOne = dataset.rows() - 1;
//
//        ProgressBar pb = new ProgressBar(dataset.columns(), "- Calculating correlations: " + dataset.columns() + " x " + dataset.columns());
//
//        for (int f = 0; f < dataset.columns(); f++) {
//
//
//            for (int g = f; g < dataset.columns(); g++) {
//                double covarianceInterim = 0;
//
//                for (int p = 0; p < dataset.rows(); p++) {
//                    covarianceInterim += dataset.getRawData()[p][f] * dataset.getRawData()[p][g];
//                }
//
//                double covariance = covarianceInterim / probeCountMinusOne;
//                correlationMatrix[f][g] = covariance;
//                correlationMatrix[g][f] = covariance;
////                System.out.println(f + "\t" + g + "\t" + covariance);
//            }
//            pb.iterate();
//        }
//        pb.close();
//        return correlationMatrix;
//    }

//    public DoubleMatrixDataset<String, String> correlateProbes(DoubleMatrixDataset<String, String> dataset) {
//
//        ConcurrentCorrelation c = new ConcurrentCorrelation();
//        return c.pairwiseCorrelationDoubleM(dataset);
//        DoubleMatrixDataset<String, String> correlationMatrix = new DoubleMatrixDataset<>(dataset.rows(), dataset.rows());
//        double probeCountMinusOne = dataset.rows() - 1;
//
//        ProgressBar pb = new ProgressBar(dataset.rows(), "- Calculating correlations: " + dataset.rows() + " x " + dataset.rows());
//        for (int f = 0; f < dataset.rows(); f++) {
//            for (int g = f; g < dataset.row; g++) {
//                double covarianceInterim = 0;
//                for (int p = 0; p < dataset.rows(); p++) {
//                    covarianceInterim += dataset.getRawData()[p][f] * dataset.getRawData()[p][g];
//                }
//                double covariance = covarianceInterim / probeCountMinusOne;
//                correlationMatrix[f][g] = covariance;
//                correlationMatrix[g][f] = covariance;
//                System.out.println(f + "\t" + g + "\t" + covariance);
//            }
//            pb.iterate();
//        }
//        pb.close();
//        return correlationMatrix;
//    }


    public Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> calculatePCA(DoubleMatrixDataset<String, String> dataset,
                                                                                                       DoubleMatrixDataset<String, String> correlationMatrix,
                                                                                                       String fileNamePrefix,
                                                                                                       Integer nrOfPCsToCalculate) throws Exception {
        String expressionFile = fileNamePrefix;
        System.out.println("Calculating PCA over file: " + fileNamePrefix);
        System.out.println("- Performing PCA over correlation matrix of size: " + correlationMatrix.rows() + "x" + correlationMatrix.rows());
        Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationMatrix.getMatrix().toArray());


        if (nrOfPCsToCalculate == null || nrOfPCsToCalculate > dataset.columns()) {
            nrOfPCsToCalculate = dataset.columns();
        } else if (nrOfPCsToCalculate < 1) {
            throw new IllegalArgumentException("Number of PCs to calculate should be at least 1");
        }

        DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<String, String>(dataset.columns(), nrOfPCsToCalculate);
        datasetEV.setRowObjects(dataset.getColObjects());
        datasetEV.setColObjects(new ArrayList<>());

        double[] eigenValues = eig.getRealEigenvalues();
        System.out.println("Eigenvalue results:");

        System.out.println("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");

        TextFile out = new TextFile(expressionFile + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.W);
        double cumExpVarPCA = 0;

        out.writeln("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");

        ArrayList<String> evcolnames = new ArrayList<>();
        for (int pca = 0; pca < nrOfPCsToCalculate; pca++) {
            double expVarPCA = PCA.getEigenValueVar(eigenValues, pca);
            double[] pca1ExpEigenVector = PCA.getEigenVector(eig, eigenValues, pca);
            for (int s = 0; s < dataset.columns(); s++) {
                datasetEV.setElementQuick(s, pca, pca1ExpEigenVector[s]);
            }
            int pcaNr = pca + 1;
            cumExpVarPCA += expVarPCA;
            out.write(pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
            evcolnames.add(pca, "Comp" + pcaNr);
            System.out.println("PCA:\t" + pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA);
        }
        out.close();
        datasetEV.setColObjects(evcolnames);
        datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectors.txt.gz");

        datasetEV = datasetEV.viewDice();

        datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectorsTransposed.txt.gz");

        datasetEV = datasetEV.viewDice();
        System.out.println("Calculating PCs");
        System.out.println("Initializing PCA matrix");
        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>(dataset.rows(), nrOfPCsToCalculate);
        datasetPCAOverSamplesPCAs.setRowObjects(dataset.getRowObjects());
        ArrayList<String> pcacolnames = new ArrayList<>();
        for (int s = 0; s < nrOfPCsToCalculate; s++) {
            pcacolnames.add("Comp" + (s + 1));
        }
        datasetPCAOverSamplesPCAs.setColObjects(pcacolnames);

        for (int p = 0; p < dataset.rows(); p++) {
            for (int t = 0; t < nrOfPCsToCalculate; t++) {
                datasetPCAOverSamplesPCAs.setElementQuick(p, t, 0);
            }
        }


        int nrprobes = dataset.rows();
        int nrsamples = dataset.columns();

//        double[][] evrawdata = datasetEV.getRawData();
//        double[][] datasetPCAOverSamplesPCAsrawdata = datasetPCAOverSamplesPCAs.getRawData();
//        double[][] datasetrawdata = dataset.getRawData();

        // multithread
        ProgressBar pb = new ProgressBar(dataset.rows(), "Calculating the PCA scores per probe: ");
        Integer finalNrOfPCsToCalculate = nrOfPCsToCalculate;
        DoubleMatrixDataset<String, String> finalDatasetEV = datasetEV;
        IntStream.range(0, nrprobes).parallel().forEach(probe -> {
//			double[] probePCAs = datasetPCAOverSamplesPCAsrawdata[probe];
            double[] probedata = dataset.getRow(probe).toArray(); // datasetrawdata[probe];
            for (int pc = 0; pc < finalNrOfPCsToCalculate; pc++) {
                for (int sample = 0; sample < nrsamples; sample++) {
                    double probeCoefficient = finalDatasetEV.getElementQuick(sample, pc); //evrawdata[sample][pc];
                    double v = datasetPCAOverSamplesPCAs.getElementQuick(probe, pc) + (probeCoefficient * probedata[sample]);
                    datasetPCAOverSamplesPCAs.setElementQuick(probe, pc, v);
                }
            }
            pb.iterateSynched();
        });
        pb.close();

//		for (int probe = 0; probe < dataset.rows(); probe++) {
//			for (int sample1 = 0; sample1 < nrOfPCsToCalculate; sample1++) {
//				for (int sample2 = 0; sample2 < dataset.columns(); sample2++) {
//					double probeCoefficient = datasetEV.getRawData()[sample2][sample1];
//					datasetPCAOverSamplesPCAs.getRawData()[probe][sample1] += probeCoefficient * dataset.getRawData()[probe][sample2];
//				}
//			}
//			pb.iterate();
//		}

        String outfilename = expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz";
        System.out.println("Saving PCA scores: " + outfilename);
        datasetPCAOverSamplesPCAs.save(outfilename);

        return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(datasetPCAOverSamplesPCAs, datasetEV);
    }

    public void correctDataForPCs(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput,
                                  DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs, DoubleMatrixDataset<String, String> datasetEV) throws IOException {
        String expressionFile = fileNamePrefix;

        System.out.println("\nInitializing residual gene expression matrix");

        if (dataset.columns() < nrPCAsOverSamplesToRemove) {
            int remainder = dataset.columns() % nrIntermediatePCAsOverSamplesToRemoveToOutput;
            nrPCAsOverSamplesToRemove = dataset.columns() - remainder;
        }

        for (int t = 0; t < nrPCAsOverSamplesToRemove; t++) {
            for (int p = 0; p < dataset.rows(); p++) {
                for (int s = 0; s < dataset.columns(); s++) {
                    double v = dataset.getElementQuick(p, s) - datasetPCAOverSamplesPCAs.getElementQuick(p, t) * datasetEV.getElementQuick(s, t);
                    dataset.setElementQuick(p, s, v);
                }
            }
            int nrPCAs = t + 1;
            if (nrIntermediatePCAsOverSamplesToRemoveToOutput > 0 && nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
                dataset.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
                System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
            }

        }
        dataset.save(expressionFile + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved.txt.gz");
    }

    public void repeatPCAOmitCertainPCAs(HashSet<Integer> pcasNotToRemove, String parentDir, String expressionFile,
                                         String eigenVectorFile, String pcaFile, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput) throws Exception {

        System.out.println("Will write output to: " + parentDir);
        String[] files = Gpio.getListOfFiles(parentDir);
        String startExpressionFileName = expressionFile;
        File st = new File(startExpressionFileName);

        // strip the parent dir name
        parentDir += Gpio.getFileSeparator();
        String minimalFilename = st.getName();
        String[] expressionFileNameElems = minimalFilename.split("\\.");
        String eigenvectorFile = null;
        String principalComponentsFile = null;

        if (minimalFilename.contains("PCAsOverSamplesRemoved")) {
            System.out.println("Warning, it seems like this data is already normalized for PCA's over samples.");
        }

        eigenvectorFile = eigenVectorFile;
        principalComponentsFile = pcaFile;

        System.out.println("Detected core file name to be: " + minimalFilename);

//        DoubleMatrixDataset<String, String> expressionDataset = new DoubleMatrixDataset<String, String>(parentDir + minimalFilename);
//        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>(principalComponentsFile);
//        DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<String, String>(eigenvectorFile);

        DoubleMatrixDataset<String, String> expressionDataset = DoubleMatrixDataset.loadDoubleData(parentDir + minimalFilename); // new DoubleMatrixDataset<String, String>(parentDir + minimalFilename);
        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = DoubleMatrixDataset.loadDoubleBinaryData(principalComponentsFile); //new DoubleMatrixDataset<String, String>(principalComponentsFile);
        DoubleMatrixDataset<String, String> datasetEV = DoubleMatrixDataset.loadDoubleData(eigenvectorFile); // new DoubleMatrixDataset<String, String>(eigenvectorFile);


        if (expressionDataset.columns() < nrPCAsOverSamplesToRemove) {
            int remainder = expressionDataset.columns() % nrIntermediatePCAsOverSamplesToRemoveToOutput;
            nrPCAsOverSamplesToRemove = expressionDataset.columns() - remainder;
        }

//        DoubleMatrixDataset<String, String> datasetResidualExpressionBasedOnPCAOverSamples = new DoubleMatrixDataset<String, String>(expressionDataset.rows(), expressionDataset.columns());
//        datasetResidualExpressionBasedOnPCAOverSamples.rowObjects = expressionDataset.rowObjects;
//        datasetResidualExpressionBasedOnPCAOverSamples.colObjects = expressionDataset.colObjects;
//
//        for (int p = 0; p < expressionDataset.rows(); p++) {
//            System.arraycopy(expressionDataset.getRawData()[p], 0, datasetResidualExpressionBasedOnPCAOverSamples.getRawData()[p], 0, expressionDataset.columns());
//        }

        if (minimalFilename.endsWith(".txt")) {
            minimalFilename = minimalFilename.substring(0, minimalFilename.length() - 4);
        } else if (minimalFilename.endsWith(".txt.gz")) {
            minimalFilename = minimalFilename.substring(0, minimalFilename.length() - 7);
        }

        System.out.println("Will not remove " + pcasNotToRemove.size() + " PCs");

        for (int t = 0; t < nrPCAsOverSamplesToRemove; t++) {
            if (!pcasNotToRemove.contains(t + 1)) {
//				System.out.println("Removing pc " + t + " or " + (t + 1));

                for (int p = 0; p < expressionDataset.rows(); p++) {
                    for (int s = 0; s < expressionDataset.columns(); s++) {
                        //datasetResidualExpressionBasedOnPCAOverSamples.rawData[p][s]-= datasetPCAOverSamplesPCAs.rawData[p][t] * datasetEV.rawData[s][t];
                        double v = expressionDataset.getElementQuick(p, s) - datasetPCAOverSamplesPCAs.getElementQuick(p, t) * datasetEV.getElementQuick(s, t);
                        expressionDataset.setElementQuick(p, s, v);
                    }
                }
            } else {
                System.out.println("Omitting PCA: " + (t + 1) + " since this component is under genetic control");
            }

            int nrPCAs = t + 1;

            if (nrIntermediatePCAsOverSamplesToRemoveToOutput > 0 && nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
                //datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt");
                expressionDataset.save(parentDir + minimalFilename + "." + nrPCAs + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");
                System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + minimalFilename + "." + nrPCAs + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");
            }

        }
        //datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved.txt");
        expressionDataset.save(parentDir + minimalFilename + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");

        System.out.println("Done\n");
    }

    private void correctForCovariate(DoubleMatrixDataset<String, String> rawdata,
                                     DoubleMatrixDataset<String, String> covariateValues, int covariateToCorrect) {

        IntStream.range(0, rawdata.rows()).parallel().forEach(probe -> {
            double[] y = rawdata.getRow(probe).toArray();
            double meanY = JSci.maths.ArrayMath.mean(y);
            double varianceY = JSci.maths.ArrayMath.variance(y);
            double[] x = covariateValues.getRow(covariateToCorrect).toArray();

            double[] rc = Regression.getLinearRegressionCoefficients(x, y);
            double correlation = JSci.maths.ArrayMath.correlation(x, y);
            double propExplainedVarianceTrait = correlation * correlation - 1.0d / (double) y.length;

            if (propExplainedVarianceTrait < 0) {
                propExplainedVarianceTrait = 0;
            }

//            explainedVariancePerEQTLProbe[d][(int) Math.round(propExplainedVarianceTrait * 100d)]++;
            double[] rawDataUpdated = new double[x.length];
            for (int s = 0; s < x.length; s++) {
                double residual = y[s] - x[s] * rc[0];
                rawDataUpdated[s] = residual;
            }

            double meanUpdated = JSci.maths.ArrayMath.mean(rawDataUpdated);
            double stdDevRatio = JSci.maths.ArrayMath.standardDeviation(rawDataUpdated) / Math.sqrt(varianceY);
            for (int s = 0; s < x.length; s++) {
                double val = ((rawDataUpdated[s] - meanUpdated) / stdDevRatio) + meanY;
                rawdata.setElementQuick(probe, s, val);
            }
        });
    }

    // NOTE: this new code switches around columns and rows for the covariate matrix
    private Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> loadCovariateValues(String covariatesToRemove, DoubleMatrixDataset<String, String> dataset) throws Exception {
        System.out.println("- Removing covariates as defined in: " + covariatesToRemove);
        TextFile covariates = new TextFile(covariatesToRemove, TextFile.R);
        int numRows = covariates.countLines() - 1; // minus the header :)
        int numCols = covariates.countCols(TextFile.tab) - 1; // minus the header's row identifier (if any)

        if (numRows == 0 || numCols == 0) {
            System.err.println("Covariate file is empty, but no covariates found in file! Is your file format correct?");
            System.err.println("The program is expecting the following: tab separated, one covariate per row, one sample per column, with sample identifiers identical to your --in file.");
            System.exit(0);
        } else {
            System.out.println("Covariate file has " + numRows + " rows and " + numCols + " columns");
        }


        // first hash up which samples are in the dataset
        HashMap<String, Integer> samplesInDatasetIndex = new HashMap<String, Integer>();
        String[] allSamplesInDataset = dataset.getColObjects().toArray(new String[0]);
        for (int i = 0; i < allSamplesInDataset.length; i++) {
            samplesInDatasetIndex.put(allSamplesInDataset[i], i);
        }

        // read the column names from the covariate file
        // expect the samples on the columns
        String[] elems = covariates.readLineElemsReturnReference(TextFile.tab); // header

        int ctr = 0;
        boolean[] sampleInDatasetIncludedInCovariates = new boolean[dataset.columns()];
        ArrayList<String> columnNames = new ArrayList<String>();
        for (int i = 1; i < elems.length; i++) {
            Integer index = samplesInDatasetIndex.get(elems[i]);
            columnNames.add(elems[i]);
            if (index != null) {
                sampleInDatasetIncludedInCovariates[index] = true;
                ctr++;
            }
        }

        // read the covariate names, expect them to be on the rows
        ArrayList<String> rowNames = new ArrayList<String>();
        elems = covariates.readLineElemsReturnReference(TextFile.tab); // first line
        while (elems != null) {
            rowNames.add(elems[0]);
            elems = covariates.readLineElemsReturnReference(TextFile.tab);
        }
        covariates.close();

        boolean isTransposed = false;
        if (ctr == 0) {
            System.err.println("No matching samples detected between covariate file and dataset. Maybe your covariate file needs to be transposed? Will test that for you now:");
            for (String rowName : rowNames) {
                Integer index = samplesInDatasetIndex.get(rowName);
                if (index != null) {
                    sampleInDatasetIncludedInCovariates[index] = true;
                    ctr++;
                }
            }

            if (ctr == 0) {
                System.err.println("Transposing the data does not seem to resolve the issue. Please check your sample identifiers.");
                System.exit(0);
            } else {
                System.out.println("Transposing the covariate file reveals: " + ctr + " samples present.");
                isTransposed = true;

            }


        }

//        if (dataset.columns() != numSamples) {
//            System.out.println("Covariates loaded from: " + covariatesToRemove + ", but the number of samples does not correspond! " + numSamples + " in covariates file, " + dataset.columns() + " in dataset...");
//            System.out.println("Please note that missing samples will be removed from your eventual corrected --in file.");
//        }
        if (!isTransposed && ctr <= numRows + 2 || isTransposed && ctr <= numCols + 2) {
            System.err.println("Less samples present than minimaly required for the normalization, (minimaly covariats+3 samples needed).");
            System.exit(0);
        }
        if (ctr < dataset.columns()) {
            System.err.println("Covariates loaded from: " + covariatesToRemove + ", but not all samples present in covariates file! " + ctr + " present in covariates file, out of " + dataset.columns() + " in dataset...");
            System.out.println("Your dataset will be adjusted accordingly.");
        }
        int nrCovariates = numRows;
        if (isTransposed) {
            nrCovariates = numCols;
        }

        // make matrix with equal sample size
//        double[][] covariateValues = new double[nrCovariates][dataset.columns()];
        DoubleMatrixDataset<String, String> covariateValues = new DoubleMatrixDataset<>(nrCovariates, dataset.columns());
        covariateValues.getMatrix().assign(Double.NaN);

        int lineCtr = 0;
        covariates.open();
        String[] headerElems = covariates.readLineElemsReturnReference(TextFile.tab); // header
        elems = covariates.readLineElemsReturnReference(TextFile.tab);
        while (elems != null) {
            if (isTransposed) {
                String sampleName = elems[0];
                Integer sampleIdInDataset = samplesInDatasetIndex.get(sampleName);
                if (sampleIdInDataset != null) {
                    for (int i = 1; i < elems.length; i++) {
                        try {
                            covariateValues.setElementQuick(i - 1, sampleIdInDataset, Double.parseDouble(elems[i]));
                        } catch (NumberFormatException e) {
//                            System.out.println("WARNING: " + elems[i] + " is not a numeric value! in " + covariatesToRemove + " at line: " + (lineCtr + 1) + ".");
//                            covariateValues[i - 1][sampleIdInDataset] = Double.NaN;
//                            sampleInDatasetIncludedInCovariates[sampleIdInDataset] = false;
                        }
                    }
                }
            } else {
                for (int i = 1; i < elems.length; i++) {
                    String sampleName = headerElems[i];
                    Integer sampleIdInDataset = samplesInDatasetIndex.get(sampleName);
                    if (sampleIdInDataset != null) {
                        try {
                            covariateValues.setElementQuick(lineCtr, sampleIdInDataset, Double.parseDouble(elems[i]));
                        } catch (NumberFormatException e) {
//                            System.out.println("WARNING: " + elems[i] + " is not a numeric value at line: " + (lineCtr + 1) + "\tcolumn: " + i);
                        }
                    }
                }
            }
            elems = covariates.readLineElemsReturnReference(TextFile.tab);
            lineCtr++;
        }
        covariates.close();

        // investigate how many covariates there actually is data for.
        int covariateCtr = 0;

        boolean[] includeCovariate = new boolean[covariateValues.rows()];
        for (int row = 0; row < covariateValues.rows(); row++) {
            int nrColsFilled = 0;
            for (int col = 0; col < covariateValues.columns(); col++) {
                if (!Double.isNaN(covariateValues.getElementQuick(row, col))) {
                    nrColsFilled++;
                }
            }

            if (nrColsFilled == 0) {
                // there's no data for this covariate....
                includeCovariate[row] = false;
            } else {
                includeCovariate[row] = true;
                covariateCtr++;
            }
        }

        if (covariateCtr == 0) {
            System.err.println("ERROR: none of your covariates seem to have valid numerical values.. Please check your covariate file.");
            System.exit(0);
        } else {
            System.out.println("After removing covariates without data, your dataset will have " + covariateCtr + " covariates (out of: " + covariateValues.rows() + ") .");
        }

        ArrayList<String> covariateNames = null;
        if (isTransposed) {
            covariateNames = columnNames;
        } else {
            covariateNames = rowNames;
        }

        if (covariateCtr != covariateValues.rows()) {
            // remove covariates with missing values
            System.out.println("Removing covariates that have no data at all.");
//            double[][] newCovariateData = new double[covariateCtr][dataset.columns()];
            DoubleMatrixDataset<String, String> newCovariateData = new DoubleMatrixDataset<>(covariateCtr, dataset.columns());
            ArrayList<String> newCovariateNames = new ArrayList<String>();
            int newCovariateCTR = 0;
            for (int row = 0; row < covariateValues.rows(); row++) {
                if (includeCovariate[row]) {
                    newCovariateNames.add(covariateNames.get(row));

                    for (int col = 0; col < covariateValues.columns(); col++) {
                        double val = covariateValues.getElementQuick(row, col);
                        newCovariateData.setElementQuick(newCovariateCTR, col, val);

                        // check whether we should include all samples, but don't remove yet: sync this with the expression/whatever dastaset
                        if (Double.isNaN(val)) {
                            sampleInDatasetIncludedInCovariates[col] = false;
                        }
                    }
                    newCovariateCTR++;
                } else {
                    System.out.println(covariateNames.get(row) + " removed.");
                }
            }


            nrCovariates = newCovariateCTR;
            covariateValues = newCovariateData;
            covariateNames = newCovariateNames;
        }
        System.out.println("");
        System.out.println("Remaining covariates: ");
        for (String s : covariateNames) {
            System.out.println(s);
        }
        System.out.println("");
        // investigate how many samples there actually is data for.
        for (int row = 0; row < covariateValues.rows(); row++) {
            for (int col = 0; col < covariateValues.columns(); col++) {
                if (Double.isNaN(covariateValues.getElementQuick(row, col))) {
                    sampleInDatasetIncludedInCovariates[col] = false;
                }
            }
        }

        int sampleCtr = 0;
        for (int q = 0; q < sampleInDatasetIncludedInCovariates.length; q++) {
            if (sampleInDatasetIncludedInCovariates[q]) {
                sampleCtr++;
            }
        }

        // remove samples that have a missing value for at least one covariate
//        if (sampleCtr == sampleInDatasetIncludedInCovariates.length) {
//            System.out.println("There were no missing values or samples in your covariate file. Sample size will remain unchanged.");
//            DoubleMatrixDataset<String, String> covariateDataset = new DoubleMatrixDataset<String, String>(covariateValues, dataset.rowObjects, covariateNames);
//            return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(covariateDataset, dataset);
//        } else {

        System.out.println("Your covariate corrected dataset will have " + sampleCtr + " samples, after removing samples with missing covariate values.");


        DoubleMatrixDataset<String, String> finalData = new DoubleMatrixDataset<>(dataset.rows(), sampleCtr);
        DoubleMatrixDataset<String, String> finalCovariates = new DoubleMatrixDataset<>(nrCovariates, sampleCtr);
        ArrayList<String> newColObjects = new ArrayList<String>();

        for (int col = 0; col < dataset.columns(); col++) {
            if (sampleInDatasetIncludedInCovariates[col]) {
                newColObjects.add(dataset.getColObjects().get(col));
            }
        }

        for (int row = 0; row < dataset.rows(); row++) {
            int includedSampleCtr = 0;
            for (int col = 0; col < dataset.columns(); col++) {
                if (sampleInDatasetIncludedInCovariates[col]) {
                    // include sample
                    finalData.setElementQuick(row, includedSampleCtr, dataset.getElementQuick(row, col));
                    includedSampleCtr++;
                }
            }
        }

        for (int row = 0; row < covariateValues.rows(); row++) {
            int includedCovariateSampleCtr = 0;
            for (int col = 0; col < dataset.columns(); col++) {
                // replace covariate data...
                if (sampleInDatasetIncludedInCovariates[col]) {
                    finalCovariates.setElementQuick(row, includedCovariateSampleCtr, covariateValues.getElementQuick(row, col));
                    includedCovariateSampleCtr++;
                }
            }
        }

        finalCovariates.setRowObjects(covariateNames);
        finalCovariates.setColObjects(newColObjects);

//        finalCovariates.save(covariatesToRemove + "-asLoadedByNormalizer.txt");
        finalData.setRowObjects(dataset.getRowObjects());
        finalData.setColObjects(newColObjects);
//        finalData.save(dataset + "-SampleSizeCorrectedForCovariates.txt");
        return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(finalCovariates, finalData);
//        }
    }

    private String removeProbesWithZeroVariance(DoubleMatrixDataset<String, String> dataset, String outputFileNamePrefix) throws Exception {
        boolean[] dataHasZeroVariance = new boolean[dataset.rows()];
        AtomicInteger nrRowsWithZeroVariance = new AtomicInteger();
        IntStream.range(0, dataset.rows()).parallel().forEach(row -> {
                    DoubleMatrix1D data = dataset.getRow(row);
                    double var = JSci.maths.ArrayMath.variance(data.toArray());
                    if (var == 0d) {
                        System.out.println("Removing probe with zero variance: " + dataset.getRowObjects().get(row) + " on line " + (row + 1));
                        nrRowsWithZeroVariance.getAndIncrement();
                        dataHasZeroVariance[row] = true;
                    }
                }
        );

        if (nrRowsWithZeroVariance.get() > 0) {
            int newNrRows = dataset.rows() - nrRowsWithZeroVariance.get();
            if (newNrRows == 0) {
                System.err.println("ERROR: all probes have zero variance!");
                System.exit(-1);
            }


            DoubleMatrixDataset<String, String> newData = new DoubleMatrixDataset<String, String>(newNrRows, dataset.columns());
//            double[][] newData = new double[newNrRows][dataset.nrCols];
            int ctr = 0;
            ArrayList<String> newRowHeader = new ArrayList<String>();
            for (int row = 0; row < dataset.rows(); row++) {
                if (!dataHasZeroVariance[row]) {
                    for (int c = 0; c < dataset.columns(); c++) {
                        newData.setElementQuick(ctr, c, dataset.getElementQuick(row, c));
                    }

                    newRowHeader.add(dataset.getRowObjects().get(row));
                    ctr++;
                }
            }
            newData.setRowObjects(newRowHeader);
            String outputFileName = outputFileNamePrefix + ".ProbesWithZeroVarianceRemoved";
            newData.save(outputFileName + ".txt.gz");
            return outputFileName;
        }
//
        return outputFileNamePrefix;
//        return null;
    }


}
