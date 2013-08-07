package eqtlmappingpipeline.normalization;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.PCA;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;
import umcg.genetica.math.stats.concurrent.ConcurrentCovariation;
import umcg.genetica.methylation.ConvertBetaToMvalue;

/**
 *
 * @author harmjan
 */
public class Normalizer {

    //nrIntermediatePCAsOverSamplesToRemoveToOutput = 5
    //nrPCAsOverSamplesToRemove = 100
    public void normalize(String expressionFile, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput, String covariatesToRemove, boolean orthogonalizecovariates, String outdir,
            boolean runQQNorm, boolean runLog2Transform, boolean runMTransform, boolean runCenterScale, boolean runPCA, boolean adjustCovariates) throws IOException {
        System.out.println("Running normalization.");
        if (outdir != null) {
            outdir = Gpio.formatAsDirectory(outdir);
            Gpio.createDir(outdir);
        } else {
            outdir = Gpio.getParentDir(expressionFile) + Gpio.getFileSeparator();
        }

        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(expressionFile);

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

        if (runQQNorm) {
            outputFileNamePrefix = quantileNormalize(dataset, outputFileNamePrefix);
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
            adjustCovariates(dataset, outputFileNamePrefix, covariatesToRemove, orthogonalizecovariates, 1E-10);
        }

        if (runPCA) {
            ConcurrentCorrelation c = new ConcurrentCorrelation(2);
            double[][] correlationMatrix = c.pairwiseCorrelation(dataset.getRawDataTransposed());
            expressionFileName = expressionFile.replace(parentDir, "");
            //outputFileNamePrefix = outdir + expressionFileName;
            Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = calculatePCA(dataset, correlationMatrix, outputFileNamePrefix, null);
            correctDataForPCs(dataset, outputFileNamePrefix, nrPCAsOverSamplesToRemove, nrIntermediatePCAsOverSamplesToRemoveToOutput, PCAResults.getLeft(), PCAResults.getRight());
        }
    }

    public String quantileNormalize(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
        double[][] rawData = dataset.getRawData();
        QuantileNormalization.quantilenormalize(rawData);
        DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
        fileNamePrefix += ".QuantileNormalized";
        datasetNormalized.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;
    }

    public String log2transform(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
        double[][] rawData = dataset.getRawData();
        Log2Transform.log2transform(rawData);
        DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
        fileNamePrefix += ".Log2Transformed";
        datasetNormalized.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;
    }

    public String mValueTransform(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
        double[][] rawData = dataset.getRawData();
        ConvertBetaToMvalue.transToMvalue(rawData);
        DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
        fileNamePrefix += ".MvalueTransformed";
        datasetNormalized.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;
    }

    public String centerAndScale(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix) throws IOException {
        double[][] rawData = dataset.getRawData();
        System.out.println("Standardizing probe mean and standard deviation");
        for (int p = 0; p < dataset.rowObjects.size(); p++) {
            double mean = Descriptives.mean(rawData[p]);
            double stdev = Math.sqrt(Descriptives.variance(rawData[p], mean));
            for (int s = 0; s < dataset.colObjects.size(); s++) {
                rawData[p][s] -= mean;
            }
        }

        dataset.setRawData(rawData);
        fileNamePrefix += ".ProbesCentered";
        dataset.save(fileNamePrefix + ".txt.gz");

        System.out.println("- Standardizing sample mean and standard deviation");
        for (int s = 0; s < dataset.colObjects.size(); s++) {
            double[] vals = new double[dataset.rowObjects.size()];
            for (int p = 0; p < dataset.rowObjects.size(); p++) {
                vals[p] = dataset.getRawData()[p][s];
            }
            double mean = Descriptives.mean(vals);
            for (int p = 0; p < dataset.rowObjects.size(); p++) {
                vals[p] -= mean;
            }
            double var = Descriptives.variance(vals, mean);
            double stdev = Math.sqrt(var);
            for (int p = 0; p < dataset.rowObjects.size(); p++) {
                dataset.getRawData()[p][s] = (vals[p] / stdev);
            }
        }

        DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawData, dataset.rowObjects, dataset.colObjects);
        fileNamePrefix += ".SamplesZTransformed";
        datasetNormalized.save(fileNamePrefix + ".txt.gz");
        return fileNamePrefix;
    }

    public String adjustCovariates(DoubleMatrixDataset<String, String> traitData, String fileNamePrefix, String covariatesToRemove, boolean orthogonalizecovariates, double varianceExplainedCutoff) throws IOException {
        // load covariate data, and remove samples for which there is missing covariate data.
        Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> covariateData = loadCovariateValues(covariatesToRemove, traitData);
        DoubleMatrixDataset<String, String> covariateDataset = covariateData.getLeft();
        traitData = covariateData.getRight();
        double[][] covariateValues = covariateDataset.getRawData();
        double[] pcaExpVar = null;
        if (orthogonalizecovariates) {
            // run PCA over the covariates.
            // we need to transpose the dataset...

            for (int p = 0; p < covariateDataset.rowObjects.size(); p++) {
                double mean = Descriptives.mean(covariateDataset.getRawData()[p]);
                double stdev = Math.sqrt(Descriptives.variance(covariateDataset.getRawData()[p], mean));
                for (int s = 0; s < covariateDataset.colObjects.size(); s++) {
                    covariateDataset.getRawData()[p][s] -= mean;
                    covariateDataset.getRawData()[p][s] /= stdev;
                }
            }

            //Covariation on a centered and scaled matrix equels the correlation.
            //Covariation is faster to compute.
            ConcurrentCovariation c = new ConcurrentCovariation(2);
            double[][] correlationMatrix = c.pairwiseCovariation(covariateDataset.getRawDataTransposed());

//            DoubleMatrixDataset<String, String> correlationMatrixDs = new DoubleMatrixDataset<String, String>(correlationMatrix, covariateDataset.colObjects, covariateDataset.colObjects);
//            correlationMatrixDs.save(covariatesToRemove+"-CorrelationMatrix.txt");

            Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = calculatePCA(covariateDataset, correlationMatrix, covariatesToRemove, null);
            // replace covariateValues with orthogonal ones...
            covariateDataset = PCAResults.getLeft();
            covariateDataset.transposeDataset();

            covariateValues = covariateDataset.getRawData();
            System.out.println(covariateValues.length + " covariates finally loaded.");
            // load the eigenvalues
            pcaExpVar = new double[covariateValues.length];
            TextFile tf = new TextFile(covariatesToRemove + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.R); // 
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
        } else {
            // PCA has been performed a-priori. Just check whether the user has supplied proper covariates.
            if (covariateValues.length > 1) {
                // check whether the covariates are orthogonal, by calculating the sum of products (inner product)
                System.out.println("Determining whether covariates are orthogonal, since you defined > 1 covariate:");
                System.out.println("Covariate1\tCovariate2\tInnerProduct\tCorrelation");
                double dotproductthreshold = 1E-5;
                if (covariateValues.length < 100) {
                    dotproductthreshold = 0.05;
                }
                for (int i = 0; i < covariateValues.length; i++) {

                    for (int j = i + 1; j < covariateValues.length; j++) {
                        double dotproduct = 0;
                        for (int v = 0; v < covariateValues[i].length; v++) {
                            dotproduct += covariateValues[i][v] * covariateValues[j][v];
                        }
                        double corr = JSci.maths.ArrayMath.correlation(covariateValues[i], covariateValues[j]);

                        if (Math.abs(dotproduct) > dotproductthreshold) {
                            System.out.println("Innerproduct > 1E-5 for covariates " + covariateDataset.rowObjects.get(i) + " and " + covariateDataset.rowObjects.get(j) + ", InnerProduct: " + Math.abs(dotproduct) + "\tCorrelation: " + corr);
                            System.out.println("If you want, we can orthogonalize the covariates for you: use --covpca in your command line.");
                            System.exit(0);
                        }

                        System.out.println(covariateDataset.rowObjects.get(i) + "\t" + covariateDataset.rowObjects.get(j) + "\t" + dotproduct + "\t" + corr);
                    }

                }

                System.out.println("Covariates are orthogonal. Now adjusting for covariates.");
            }
        }


        double[][] rawdata = traitData.getRawData();
        for (int i = 0; i < covariateValues.length; i++) {
            if (pcaExpVar == null || pcaExpVar[i] > varianceExplainedCutoff) {
                correctForCovariate(rawdata, covariateValues, i);
            } else {
                System.out.println("Not regressing covariate: " + i + " because explained variance < " + varianceExplainedCutoff + ": " + pcaExpVar[i]);
            }
        }

        DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String>(rawdata, traitData.rowObjects, traitData.colObjects);
        fileNamePrefix += ".CovariatesRemoved";
        datasetNormalized.save(fileNamePrefix + ".txt.gz");



        return fileNamePrefix;
    }

    /**
     * Calculate correlation over columns in DoubleMatrixDataset. WARNING: this
     * method assumes that SD == 1 and mean == 0 (which makes the covariance
     * equal to the correlation).
     *
     * @param dataset
     * @return
     */
    private double[][] correlateSamples(DoubleMatrixDataset<String, String> dataset) {
        double[][] correlationMatrix = new double[dataset.colObjects.size()][dataset.colObjects.size()];
        double probeCountMinusOne = dataset.rowObjects.size() - 1;

        ProgressBar pb = new ProgressBar(dataset.colObjects.size(), "- Calculating correlations: " + dataset.colObjects.size() + " x " + dataset.colObjects.size());

        for (int f = 0; f < dataset.colObjects.size(); f++) {



            for (int g = f; g < dataset.colObjects.size(); g++) {
                double covarianceInterim = 0;
                for (int p = 0; p < dataset.rowObjects.size(); p++) {
                    covarianceInterim += dataset.getRawData()[p][f] * dataset.getRawData()[p][g];
                }
                double covariance = covarianceInterim / probeCountMinusOne;
                correlationMatrix[f][g] = covariance;
                correlationMatrix[g][f] = covariance;
//                System.out.println(f + "\t" + g + "\t" + covariance);
            }
            pb.iterate();
        }
        pb.close();
        return correlationMatrix;
    }

    public double[][] correlateProbes(DoubleMatrixDataset<String, String> dataset) {

        double[][] correlationMatrix = new double[dataset.rowObjects.size()][dataset.rowObjects.size()];
        double probeCountMinusOne = dataset.rowObjects.size() - 1;

        ProgressBar pb = new ProgressBar(dataset.rowObjects.size(), "- Calculating correlations: " + dataset.rowObjects.size() + " x " + dataset.rowObjects.size());
        for (int f = 0; f < dataset.rowObjects.size(); f++) {
            for (int g = f; g < dataset.rowObjects.size(); g++) {
                double covarianceInterim = 0;
                for (int p = 0; p < dataset.rowObjects.size(); p++) {
                    covarianceInterim += dataset.getRawData()[p][f] * dataset.getRawData()[p][g];
                }
                double covariance = covarianceInterim / probeCountMinusOne;
                correlationMatrix[f][g] = covariance;
                correlationMatrix[g][f] = covariance;
                System.out.println(f + "\t" + g + "\t" + covariance);
            }
            pb.iterate();
        }
        pb.close();
        return correlationMatrix;
    }

    public Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> calculatePCA(DoubleMatrixDataset<String, String> dataset, double[][] correlationMatrix, String fileNamePrefix, Integer nrOfPCsToCalculate) throws IOException {
        String expressionFile = fileNamePrefix;
        System.out.println("- Performing PCA over correlation matrix of size: " + correlationMatrix.length + "x" + correlationMatrix.length);
        Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationMatrix);

        if (nrOfPCsToCalculate == null || nrOfPCsToCalculate > dataset.colObjects.size()) {
            nrOfPCsToCalculate = dataset.colObjects.size();
        } else if (nrOfPCsToCalculate < 1) {
            throw new IllegalArgumentException("Number of PCs to calculate should be at least 1");
        }

        DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<String, String>(dataset.colObjects.size(), nrOfPCsToCalculate);
        datasetEV.rowObjects = dataset.colObjects;
        double[] eigenValues = eig.getRealEigenvalues();
        System.out.println("Eigenvalue results:");

        System.out.println("PCA\tPCANr\tEigenValue\t\tExplainedVariance");

        TextFile out = new TextFile(expressionFile + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.W);
        double cumExpVarPCA = 0;




        for (int pca = 0; pca < nrOfPCsToCalculate; pca++) {
            double expVarPCA = PCA.getEigenValueVar(eigenValues, pca);
            double[] pca1ExpEigenVector = PCA.getEigenVector(eig, eigenValues, pca);
            for (int s = 0; s < dataset.colObjects.size(); s++) {
                datasetEV.getRawData()[s][pca] = pca1ExpEigenVector[s];
            }
            int pcaNr = pca + 1;
            cumExpVarPCA += expVarPCA;
            out.write(pcaNr + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
            datasetEV.colObjects.set(pca, "Comp" + String.valueOf(pcaNr));
            System.out.println("PCA:\t" + pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA);
        }
        out.close();

        datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectors.txt.gz");

        datasetEV.transposeDataset();

        datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectorsTransposed.txt.gz");

        datasetEV.transposeDataset();
        System.out.println("Calculating PCs");
        System.out.println("Initializing PCA matrix");
        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>(dataset.rowObjects.size(), nrOfPCsToCalculate);
        datasetPCAOverSamplesPCAs.rowObjects = dataset.rowObjects;
        for (int s = 0; s < nrOfPCsToCalculate; s++) {
            datasetPCAOverSamplesPCAs.colObjects.set(s, "Comp" + String.valueOf(s + 1));
        }
        for (int p = 0; p < dataset.rowObjects.size(); p++) {
            for (int t = 0; t < nrOfPCsToCalculate; t++) {
                datasetPCAOverSamplesPCAs.getRawData()[p][t] = 0;
            }
        }


        ProgressBar pb = new ProgressBar(dataset.rowObjects.size(), "Calculating the PCA scores per probe: ");
        for (int probe = 0; probe < dataset.rowObjects.size(); probe++) {
            for (int sample1 = 0; sample1 < nrOfPCsToCalculate; sample1++) {
                for (int sample2 = 0; sample2 < dataset.colObjects.size(); sample2++) {
                    double probeCoefficient = datasetEV.getRawData()[sample2][sample1];
                    datasetPCAOverSamplesPCAs.getRawData()[probe][sample1] += probeCoefficient * dataset.getRawData()[probe][sample2];
                }
            }
            pb.iterate();
        }
        pb.close();

        String outfilename = expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz";
        System.out.println("Saving PCA scores: " + outfilename);
        datasetPCAOverSamplesPCAs.save(outfilename);

        return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(datasetPCAOverSamplesPCAs, datasetEV);
    }

    public void correctDataForPCs(DoubleMatrixDataset<String, String> dataset, String fileNamePrefix, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput,
            DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs, DoubleMatrixDataset<String, String> datasetEV) throws IOException {
        String expressionFile = fileNamePrefix;
        System.out.println("\nInitializing residual gene expression matrix");
        DoubleMatrixDataset<String, String> datasetResidualExpressionBasedOnPCAOverSamples = new DoubleMatrixDataset<String, String>(dataset.rowObjects.size(), dataset.colObjects.size());
        datasetResidualExpressionBasedOnPCAOverSamples.rowObjects = dataset.rowObjects;
        datasetResidualExpressionBasedOnPCAOverSamples.colObjects = dataset.colObjects;
        for (int p = 0; p < dataset.rowObjects.size(); p++) {
            System.arraycopy(dataset.getRawData()[p], 0, datasetResidualExpressionBasedOnPCAOverSamples.getRawData()[p], 0, dataset.colObjects.size());
        }

        if (dataset.colObjects.size() < nrPCAsOverSamplesToRemove) {
            int remainder = dataset.colObjects.size() % nrIntermediatePCAsOverSamplesToRemoveToOutput;
            nrPCAsOverSamplesToRemove = dataset.colObjects.size() - remainder;
        }

        for (int t = 0; t < nrPCAsOverSamplesToRemove; t++) {
            for (int p = 0; p < dataset.rowObjects.size(); p++) {
                for (int s = 0; s < dataset.colObjects.size(); s++) {
                    //datasetResidualExpressionBasedOnPCAOverSamples.rawData[p][s]-= datasetPCAOverSamplesPCAs.rawData[p][t] * datasetEV.rawData[s][t];
                    dataset.getRawData()[p][s] -= datasetPCAOverSamplesPCAs.getRawData()[p][t] * datasetEV.getRawData()[s][t];
                }
            }
            int nrPCAs = t + 1;
            if (nrIntermediatePCAsOverSamplesToRemoveToOutput > 0 && nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
                //datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt");
                dataset.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
                System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt.gz");
            }

        }
        //datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved.txt");
        dataset.save(expressionFile + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved.txt.gz");

    }

    public void repeatPCAOmitCertainPCAs(HashSet<Integer> pcasNotToRemove, String expressionFile, int nrPCAsOverSamplesToRemove, int nrIntermediatePCAsOverSamplesToRemoveToOutput) throws IOException {
        String parentDir = Gpio.getParentDir(expressionFile);
        String[] files = Gpio.getListOfFiles(parentDir);
        String startExpressionFileName = expressionFile;

        // strip the parent dir name
        parentDir += "/";
        String minimalFilename = startExpressionFileName.replaceAll(parentDir, "");
        String[] expressionFileNameElems = minimalFilename.split("\\.");
        String eigenvectorFile = null;
        String principalComponentsFile = null;

        for (String file : files) {
            if (file.length() < minimalFilename.length() && file.contains(expressionFileNameElems[0])) {
                minimalFilename = file;
            } else if (file.toLowerCase().contains("pcaoversampleseigenvectors")) {
                eigenvectorFile = parentDir + "" + file;
            } else if (file.toLowerCase().contains("pcaoversamplesprincipalcomponents")) {
                principalComponentsFile = parentDir + "" + file;
            }
        }

        boolean fileFound = true;
        if (eigenvectorFile == null) {
            System.err.println("Could not find file containing 'PCAOverSamplesEigenvectors' in directory: " + parentDir);
            fileFound = false;
        }

        if (eigenvectorFile == null) {
            System.err.println("Could not find file containing 'PCAOverSamplesPrincipalComponents' in directory: " + parentDir);
            fileFound = false;
        }

        if (!fileFound) {
            System.exit(0);
        }

        System.out.println("Detected core file name to be: " + minimalFilename);

        DoubleMatrixDataset<String, String> expressionDataset = new DoubleMatrixDataset<String, String>(startExpressionFileName);
        DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new DoubleMatrixDataset<String, String>(principalComponentsFile);
        DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<String, String>(eigenvectorFile);

        if (expressionDataset.colObjects.size() < nrPCAsOverSamplesToRemove) {
            int remainder = expressionDataset.colObjects.size() % nrIntermediatePCAsOverSamplesToRemoveToOutput;
            nrPCAsOverSamplesToRemove = expressionDataset.colObjects.size() - remainder;
        }


        DoubleMatrixDataset<String, String> datasetResidualExpressionBasedOnPCAOverSamples = new DoubleMatrixDataset<String, String>(expressionDataset.rowObjects.size(), expressionDataset.colObjects.size());
        datasetResidualExpressionBasedOnPCAOverSamples.rowObjects = expressionDataset.rowObjects;
        datasetResidualExpressionBasedOnPCAOverSamples.colObjects = expressionDataset.colObjects;

        for (int p = 0; p < expressionDataset.rowObjects.size(); p++) {
            System.arraycopy(expressionDataset.getRawData()[p], 0, datasetResidualExpressionBasedOnPCAOverSamples.getRawData()[p], 0, expressionDataset.colObjects.size());
        }


        for (int t = 0; t < nrPCAsOverSamplesToRemove; t++) {
            if (!pcasNotToRemove.contains(t + 1)) {

                for (int p = 0; p < expressionDataset.rowObjects.size(); p++) {
                    for (int s = 0; s < expressionDataset.colObjects.size(); s++) {
                        //datasetResidualExpressionBasedOnPCAOverSamples.rawData[p][s]-= datasetPCAOverSamplesPCAs.rawData[p][t] * datasetEV.rawData[s][t];
                        expressionDataset.getRawData()[p][s] -= datasetPCAOverSamplesPCAs.getRawData()[p][t] * datasetEV.getRawData()[s][t];
                    }
                }
            } else {
                System.out.println("Omitting PCA: " + (t + 1) + " since this component is under genetic control");
            }

            int nrPCAs = t + 1;
            if (nrIntermediatePCAsOverSamplesToRemoveToOutput > 0 && nrPCAs % nrIntermediatePCAsOverSamplesToRemoveToOutput == 0) {
                //datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAs + "PCAsOverSamplesRemoved.txt");
                expressionDataset.save(minimalFilename + "." + nrPCAs + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");
                System.out.println("Removed\t" + nrPCAs + "\tPCs. File:\t" + minimalFilename + "." + nrPCAs + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");
            }

        }
        //datasetResidualExpressionBasedOnPCAOverSamples.save(expressionFile + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved.txt");
        expressionDataset.save(minimalFilename + "." + nrPCAsOverSamplesToRemove + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz");

        System.out.println("Done\n");
    }

    private void correctForCovariate(double[][] rawdata, double[][] covariateValues, int covariateToCorrect) {
        for (int probe = 0; probe < rawdata.length; probe++) {
            double[] y = rawdata[probe];
            double meanY = JSci.maths.ArrayMath.mean(y);
            double varianceY = JSci.maths.ArrayMath.variance(y);
            double[] x = covariateValues[covariateToCorrect];



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
                rawDataUpdated[s] -= meanUpdated;
                rawDataUpdated[s] /= stdDevRatio;
                rawDataUpdated[s] += meanY;
            }
            System.arraycopy(rawDataUpdated, 0, rawdata[probe], 0, x.length);
        }
    }

    // NOTE: this new code switches around columns and rows for the covariate matrix
    private Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> loadCovariateValues(String covariatesToRemove, DoubleMatrixDataset<String, String> dataset) throws IOException {
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
        String[] allSamplesInDataset = dataset.colObjects.toArray(new String[0]);
        for (int i = 0; i < allSamplesInDataset.length; i++) {
            samplesInDatasetIndex.put(allSamplesInDataset[i], i);
        }

        // read the column names from the covariate file
        // expect the samples on the columns
        String[] elems = covariates.readLineElemsReturnReference(TextFile.tab); // header

        int ctr = 0;
        boolean[] sampleInDatasetIncludedInCovariates = new boolean[dataset.colObjects.size()];
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

//        if (dataset.colObjects.size() != numSamples) {
//            System.out.println("Covariates loaded from: " + covariatesToRemove + ", but the number of samples does not correspond! " + numSamples + " in covariates file, " + dataset.colObjects.size() + " in dataset...");
//            System.out.println("Please note that missing samples will be removed from your eventual corrected --in file.");
//        }

        if (ctr < dataset.colObjects.size()) {
            System.err.println("Covariates loaded from: " + covariatesToRemove + ", but not all samples present in covariates file! " + ctr + " present in covariates file, out of " + dataset.colObjects.size() + " in dataset...");
            System.out.println("Your dataset will be adjusted accordingly.");
        }
        int nrCovariates = numRows;
        if (isTransposed) {
            nrCovariates = numCols;
        }

        // make matrix with equal sample size
        double[][] covariateValues = new double[nrCovariates][dataset.colObjects.size()];
        for (int row = 0; row < covariateValues.length; row++) {
            for (int col = 0; col < covariateValues[row].length; col++) {
                covariateValues[row][col] = Double.NaN;
            }
        }

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
                            covariateValues[i - 1][sampleIdInDataset] = Double.parseDouble(elems[i]);
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
                            covariateValues[lineCtr][sampleIdInDataset] = Double.parseDouble(elems[i]);
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
        boolean[] includeCovariate = new boolean[covariateValues.length];
        for (int row = 0; row < covariateValues.length; row++) {
            int nrColsFilled = 0;
            for (int col = 0; col < covariateValues[row].length; col++) {
                if (!Double.isNaN(covariateValues[row][col])) {
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
            System.out.println("After removing covariates without data, your dataset will have " + covariateCtr + " covariates (out of: " + covariateValues.length + ") .");
        }

        ArrayList<String> covariateNames = null;
        if (isTransposed) {
            covariateNames = columnNames;
        } else {
            covariateNames = rowNames;
        }

        if (covariateCtr != covariateValues.length) {
            // remove covariates with missing values
            System.out.println("Removing covariates that have no data at all.");
            double[][] newCovariateData = new double[covariateCtr][dataset.colObjects.size()];
            ArrayList<String> newCovariateNames = new ArrayList<String>();
            int newCovariateCTR = 0;
            for (int row = 0; row < covariateValues.length; row++) {
                if (includeCovariate[row]) {
                    newCovariateNames.add(covariateNames.get(row));

                    for (int col = 0; col < covariateValues[row].length; col++) {
                        newCovariateData[newCovariateCTR][col] = covariateValues[row][col];

                        // check whether we should include all samples, but don't remove yet: sync this with the expression/whatever dastaset
                        if (Double.isNaN(covariateValues[row][col])) {
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
        for (int row = 0; row < covariateValues.length; row++) {
            for (int col = 0; col < covariateValues[row].length; col++) {
                if (Double.isNaN(covariateValues[row][col])) {
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
        double[][] rawData = dataset.getRawData();
        double[][] newRawData = new double[rawData.length][sampleCtr];
        double[][] finalCovariateData = new double[nrCovariates][sampleCtr];
        ArrayList<String> newColObjects = new ArrayList<String>();

        for (int col = 0; col < dataset.colObjects.size(); col++) {
            if (sampleInDatasetIncludedInCovariates[col]) {
                newColObjects.add(dataset.colObjects.get(col));
            }
        }

        for (int row = 0; row < rawData.length; row++) {
            int includedSampleCtr = 0;
            for (int col = 0; col < dataset.colObjects.size(); col++) {
                if (sampleInDatasetIncludedInCovariates[col]) {
                    // include sample
                    newRawData[row][includedSampleCtr] = rawData[row][col];
                    includedSampleCtr++;
                }
            }
        }

        for (int row = 0; row < covariateValues.length; row++) {
            int includedCovariateSampleCtr = 0;
            for (int col = 0; col < dataset.colObjects.size(); col++) {
                // replace covariate data...
                if (sampleInDatasetIncludedInCovariates[col]) {
                    finalCovariateData[row][includedCovariateSampleCtr] = covariateValues[row][col];
                    includedCovariateSampleCtr++;
                }
            }
        }

        DoubleMatrixDataset<String, String> covariateDataset = new DoubleMatrixDataset<String, String>(finalCovariateData, covariateNames, newColObjects);
        covariateDataset.save(covariatesToRemove + "-asLoadedByNormalizer.txt");
        DoubleMatrixDataset<String, String> newDataset = new DoubleMatrixDataset<String, String>(newRawData, dataset.rowObjects, newColObjects);
        newDataset.save(dataset.fileName + "-SampleSizeCorrectedForCovariates.txt");
        return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(covariateDataset, newDataset);
//        }
    }
}
