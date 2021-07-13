/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.console.MultiThreadProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.*;
import umcg.genetica.math.PCA;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.math.stats.VIF;
import umcg.genetica.text.Strings;
import umcg.genetica.util.RankDoubleArray;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

/**
 * @author harmjan
 */
public class EQTLRegression {

    public void regressOutEQTLEffects(ArrayList<Pair<String, String>> eqtls, TriTyperGeneticalGenomicsDataset[] gg) throws Exception {
        regressOutEQTLEffects(eqtls, gg, true);
    }

    public void regressOutEQTLEffects(ArrayList<Pair<String, String>> eqtls, TriTyperGeneticalGenomicsDataset[] gg, boolean useOLS) throws IOException {


        EQTL[] eqtlsToRegressOut = new EQTL[eqtls.size()];
        for (int q = 0; q < eqtls.size(); q++) {
            eqtlsToRegressOut[q] = new EQTL();
            eqtlsToRegressOut[q].setRsName(eqtls.get(q).getLeft());
            eqtlsToRegressOut[q].setProbe(eqtls.get(q).getRight());
        }
        System.out.println("About to regress out: " + eqtls.size() + " QTLs from data.");
        if (useOLS) {
            regressOLS(eqtlsToRegressOut, gg);
        } else {
            regress(eqtlsToRegressOut, gg);
        }

    }

    public void regressOutEQTLEffects(EQTL[] eqtls, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {
        regressOutEQTLEffects(eqtls, gg, true);
    }

    public void regressOutEQTLEffects(EQTL[] eqtlsToRegressOut, TriTyperGeneticalGenomicsDataset[] gg, boolean useOLS) throws IOException {

        System.out.println("About to regress out: " + eqtlsToRegressOut.length + " QTLs from data.");
        if (useOLS) {
            regressOLS(eqtlsToRegressOut, gg);
        } else {
            regress(eqtlsToRegressOut, gg);
        }
    }

    public void regressOutEQTLEffects(String regressOutEQTLEffectFileName, boolean outputfiles, TriTyperGeneticalGenomicsDataset[] gg) throws Exception {
        regressOutEQTLEffects(regressOutEQTLEffectFileName, outputfiles, gg, true);
    }

    public void regressOutEQTLEffects(String regressOutEQTLEffectFileName, boolean outputfiles, TriTyperGeneticalGenomicsDataset[] gg, boolean useOLS) throws Exception {

        System.out.println("\n\n\nRemoving eQTL effects from the following eQTL file: '" + regressOutEQTLEffectFileName);

        // determine whether we're looking at an eQTL file or at a SNP/gene combo file
        TextFile tf = new TextFile(regressOutEQTLEffectFileName, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        EQTL[] eqtlsToRegressOut = null;
        if (elems != null && elems.length == 2) {
            ArrayList<EQTL> toRegress = new ArrayList<>();
            while (elems != null) {
                EQTL e = new EQTL();
                e.setProbe(elems[1]);
                e.setRsName(elems[0]);
                toRegress.add(e);
                elems = tf.readLineElems(TextFile.tab);
            }
            eqtlsToRegressOut = toRegress.toArray(new EQTL[0]);
        } else if (elems != null && elems.length > 2) {
            // assume eqtl file
            elems = tf.readLineElems(TextFile.tab); // skip header
            ArrayList<EQTL> toRegress = new ArrayList<>();
            while (elems != null) {
                EQTL e = new EQTL();
                e.setProbe(elems[4]);
                e.setRsName(elems[1]);
                toRegress.add(e);
                elems = tf.readLineElems(TextFile.tab);
            }
            eqtlsToRegressOut = toRegress.toArray(new EQTL[0]);
        }
        tf.close();


        System.out.println("Number of eQTLs to regress out found in file:\t" + eqtlsToRegressOut.length);
        if (useOLS) {
            regressOLS(eqtlsToRegressOut, gg);
        } else {
            regress(eqtlsToRegressOut, gg);
        }
        if (outputfiles) {
            for (int d = 0; d < gg.length; d++) {
                TriTyperGeneticalGenomicsDataset ds = gg[d];
                TriTyperExpressionData dsexp = ds.getExpressionData();
                double[][] matrix = dsexp.getMatrix();
                String[] probes = dsexp.getProbes();
                String[] individuals = dsexp.getIndividuals();
                String filename = ds.getSettings().expressionLocation;

                DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>(matrix, Arrays.asList(probes), Arrays.asList(individuals));

                System.out.println("Saving expression file after removal of eQTL effects: " + filename + "-EQTLEffectsRemoved.txt.gz");
                dsout.save(filename + "-EQTLEffectsRemoved.txt.gz");
            }
        }

    }


    // regress out eQTLs, replace missing values using genotype probabilities given gene expression values
    private void regressWithmissingData(EQTL[] eqtlsToRegressOut, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {
        //Inventorize whether for a particular probe there are multiple SNPs that we want to regress out:
        HashMap<String, ArrayList<EQTL>> qtlsPerGene = new HashMap<String, ArrayList<EQTL>>();
        HashMap<EQTL, Integer> hashEQTLIds = new HashMap<EQTL, Integer>();
        int nrProbesWithMultipleCovariates = 0;

        for (int v = 0; v < eqtlsToRegressOut.length; v++) {
            EQTL current = eqtlsToRegressOut[v];
            hashEQTLIds.put(current, v);
            String probe = current.getProbe();

            if (!qtlsPerGene.containsKey(probe)) {
                ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
                eqtls.add(current);
                qtlsPerGene.put(probe, eqtls);
            } else {
                qtlsPerGene.get(probe).add(current);
                nrProbesWithMultipleCovariates++;
            }
        }

        if (nrProbesWithMultipleCovariates > 0) {
            System.out.println("There are:\t" + nrProbesWithMultipleCovariates + "\tprobes for which we want to regress out multiple SNPs. This will be conducted through multiple regression employing PCA.");
        }

        // remove the eqtl effects
        System.out.println("Removing eQTLs:");
        int[] nrEQTLsRegressedOut = new int[gg.length];
        int[][] explainedVariancePerEQTLProbe = new int[gg.length][101];

        // set up IO
        SNPLoader[] ggSNPLoaders = new SNPLoader[gg.length];
        for (int d = 0; d < gg.length; d++) {
            ggSNPLoaders[d] = gg[d].getGenotypeData().createSNPLoader(1);
        }

        // iterate eQTL genes
        for (String gene : qtlsPerGene.keySet()) {
            ArrayList<EQTL> qtls = qtlsPerGene.get(gene);

            EQTL maxeqtl = null;
            Double maxZ = null;

            // select the strongest effect
            for (EQTL e : qtls) {
                double absz = e.getZscoreAbs();
                if (maxZ == null) {
                    maxZ = absz;
                    maxeqtl = e;
                }
            }


            // 1. fit model for which there is data
            // 2. impute missing genotypes using
            ArrayList<Double> allY = new ArrayList<>();
            ArrayList<Double> allX = new ArrayList<>();
            ArrayList<Integer> datasetNr = new ArrayList<>();


            // load the data
            ArrayList<HashMap<Integer, Boolean>> snpPassesQCPerDs = new ArrayList<HashMap<Integer, Boolean>>();
            for (int d = 0; d < gg.length; d++) {
                TriTyperGeneticalGenomicsDataset currentDataset = gg[d];
                Integer snpId = currentDataset.getGenotypeData().getSnpToSNPId().get(maxeqtl.getRsName());
                HashMap<Integer, Boolean> snpPassesQC = snpPassesQCPerDs.get(d);
                if (snpPassesQC == null) {
                    snpPassesQCPerDs.add(new HashMap<>());
                }

                if (snpId == -9 || (snpPassesQC.get(snpId) == null && !snpPassesQC.get(snpId))) {
                    // fill x and y with -1?
                } else {
                    SNP currentSNP = currentDataset.getGenotypeData().getSNPObject(snpId);
                    try {
                        ggSNPLoaders[d].loadGenotypes(currentSNP);
                    } catch (IOException ex) {
                        ex.printStackTrace();
                        System.exit(-1);
                    }

                    if (ggSNPLoaders[d].hasDosageInformation()) {
                        try {
                            ggSNPLoaders[d].loadDosage(currentSNP);
                        } catch (IOException ex) {
                            ex.printStackTrace();
                            System.exit(-1);
                        }
                    }

                    if (!currentSNP.passesQC()) {
                        // fill x and y with -1?
                        snpPassesQC.put(snpId, false);
                        currentSNP.clearGenotypes();
                    } else {
                        int[] indWGA = currentDataset.getExpressionToGenotypeIdArray();
                        double[] x = currentSNP.selectGenotypes(indWGA, true, true); // we want missing genotypes to be included in this case
                        double meanX = 0;
                        int ctr = 0;

                        // calculate mean and variance for genotype, taking into account missing values, if any
                        for (int i = 0; i < x.length; i++) {
                            System.out.println(i + "\t" + x[i]);
                            if (x[i] != -1) {
                                meanX += x[i];
                                ctr++;
                            }
                        }

                        System.out.println();
                        meanX /= ctr;
                        double varianceX = 0.0;
                        for (int i = 0; i < x.length; i++) {
                            if (x[i] != -1) {
                                varianceX += (x[i] - meanX) * (x[i] - meanX);
                            }
                        }
                        varianceX /= (ctr - 1);

                        if (varianceX != 0 && currentDataset.getTotalGGSamples() == x.length) {
                            for (int i = 0; i < x.length; i++) {
                                if (x[i] != -1) {
                                    // add to x and y
                                } else {

                                }
                            }
                            snpPassesQC.put(snpId, true);
                        } else {
                            // fill x and y with -1?
                            snpPassesQC.put(snpId, false);
                        }
                    }
                }


            }

        }


    }

    /**
     * Removes the effect of a supplied list of eQTL from the datasets by use of
     * regression
     */
    private void regress(EQTL[] eqtlsToRegressOut, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {

        //Inventorize whether for a particular probe there are multiple SNPs that we want to regress out:
        HashMap<String, ArrayList<EQTL>> hashProbesCovariates = new HashMap<String, ArrayList<EQTL>>();
        HashMap<EQTL, Integer> hashEQTLIds = new HashMap<EQTL, Integer>();
        int nrProbesWithMultipleCovariates = 0;

        for (int v = 0; v < eqtlsToRegressOut.length; v++) {
            EQTL current = eqtlsToRegressOut[v];
            hashEQTLIds.put(current, v);
            String probe = current.getProbe();

            if (!hashProbesCovariates.containsKey(probe)) {
                ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
                eqtls.add(current);
                hashProbesCovariates.put(probe, eqtls);
            } else {
                hashProbesCovariates.get(probe).add(current);
                nrProbesWithMultipleCovariates++;
            }
        }

        if (nrProbesWithMultipleCovariates > 0) {
            System.out.println("There are:\t" + nrProbesWithMultipleCovariates + "\tprobes for which we want to regress out multiple SNPs. This will be conducted through multiple regression employing PCA.");
        }

        // remove the eqtl effects
        System.out.println("Removing eQTLs:");
        int[] nrEQTLsRegressedOut = new int[gg.length];
        int[][] explainedVariancePerEQTLProbe = new int[gg.length][101];

        SNPLoader[] ggSNPLoaders = new SNPLoader[gg.length];
        boolean dosageInformationPresentForAllDatasets = true;
        for (int d = 0; d < gg.length; d++) {
            ggSNPLoaders[d] = gg[d].getGenotypeData().createSNPLoader(1);
            if (!ggSNPLoaders[d].hasDosageInformation()) {
                dosageInformationPresentForAllDatasets = false;
            }
        }

        //Remove multiple SNPs acting on one single probe:
        MultiThreadProgressBar pb = new MultiThreadProgressBar(gg.length);
        IntStream.range(0, gg.length).parallel().forEach(d -> {
            HashSet<EQTL> hashEQTLsMultipleRegressionRegressedOut = new HashSet<EQTL>();
            HashMap<Integer, Boolean> snpPassesQC = new HashMap<Integer, Boolean>();

            TriTyperGeneticalGenomicsDataset currentDataset = gg[d];
            String[] probes = gg[d].getExpressionData().getProbes();
            System.out.print("Dataset:\t" + gg[d].getSettings().name);
            pb.setSubtasks(d, probes.length);
            for (int p = 0; p < probes.length; p++) {


                ArrayList<EQTL> covariatesForThisProbe = hashProbesCovariates.get(probes[p]);

                if (covariatesForThisProbe != null) {
                    ArrayList<EQTL> eventualListOfEQTLs = new ArrayList<EQTL>();
                    ArrayList<SNP> snpsForProbe = new ArrayList<SNP>();
                    ArrayList<double[]> xs = new ArrayList<double[]>();


                    // preload SNPs
                    for (EQTL e : covariatesForThisProbe) {
                        if (!hashEQTLsMultipleRegressionRegressedOut.contains(e)) {
                            Integer snpId = gg[d].getGenotypeData().getSnpToSNPId().get(e.getRsName());

                            if (snpId != -9 && (snpPassesQC.get(snpId) == null || snpPassesQC.get(snpId))) {
                                // load SNP

                                SNP currentSNP = currentDataset.getGenotypeData().getSNPObject(snpId);
                                try {
                                    ggSNPLoaders[d].loadGenotypes(currentSNP);
                                } catch (IOException ex) {
                                    ex.printStackTrace();
                                    System.exit(-1);
                                }

                                if (ggSNPLoaders[d].hasDosageInformation()) {
                                    try {
                                        ggSNPLoaders[d].loadDosage(currentSNP);
                                    } catch (IOException ex) {
                                        ex.printStackTrace();
                                        System.exit(-1);
                                    }
                                }

                                if (currentSNP.passesQC()) {

                                    int[] indWGA = currentDataset.getExpressionToGenotypeIdArray();
                                    double[] x = currentSNP.selectGenotypes(indWGA, true, true); // we want missing genotypes to be included in this case
                                    double meanX = 0;
                                    int ctr = 0;

                                    // calculate mean and variance for genotype, taking into account missing values, if any
                                    for (int i = 0; i < x.length; i++) {
                                        System.out.println(i + "\t" + x[i]);
                                        if (x[i] != -1) {
                                            meanX += x[i];
                                            ctr++;
                                        }
                                    }

                                    System.out.println();
                                    meanX /= ctr;
                                    double varianceX = 0.0;
                                    for (int i = 0; i < x.length; i++) {
                                        if (x[i] != -1) {
                                            varianceX += (x[i] - meanX) * (x[i] - meanX);
                                        }
                                    }
                                    varianceX /= (ctr - 1);

                                    if (varianceX != 0 && currentDataset.getTotalGGSamples() == x.length) {
                                        for (int i = 0; i < x.length; i++) {
                                            if (x[i] != -1) {
                                                x[i] -= meanX;
                                            } else {
                                                x[i] = 0; // replace missing values with mean (which is now 0)
                                            }
                                        }
                                        eventualListOfEQTLs.add(e);
                                        snpsForProbe.add(currentSNP);
                                        xs.add(x);
                                        snpPassesQC.put(snpId, true);
                                    } else {
                                        snpPassesQC.put(snpId, false);
                                    }
                                } else {
                                    snpPassesQC.put(snpId, false);
                                    currentSNP.clearGenotypes();
                                }
                            }
                        }
                    }

                    // regress out single effects
                    if (eventualListOfEQTLs.size() == 1) {
                        SNP currentSNP = snpsForProbe.get(0);

                        int[] expressionToGenotypeId = currentDataset.getExpressionToGenotypeIdArray();
                        double[] x = xs.get(0);

                        //Get the expression data:
                        double[][] rawData = currentDataset.getExpressionData().getMatrix();
                        double meanY;
                        double varianceY;

                        //Check what the number of samples is with genotype data available:
                        int nrSamplesWGenotypeData = x.length;

                        double[] y = new double[nrSamplesWGenotypeData];
                        int totalGGSamples = currentDataset.getTotalGGSamples();

                        // Copy expression data.
                        int itr = 0;
                        for (int s = 0; s < rawData[p].length; s++) {
                            int genotypeId = expressionToGenotypeId[s];
                            if (genotypeId != -1) {
                                byte genotype = currentSNP.getGenotypes()[genotypeId];
                                if (genotype != -1 && currentDataset.getGenotypeData().getIsIncluded()[genotypeId]) {
                                    double dVal = rawData[p][s];
                                    y[itr] = dVal;
                                    itr++;
                                }
                            }
                        }
                        //Normalize subset of data:
                        meanY = JSci.maths.ArrayMath.mean(y);
                        varianceY = JSci.maths.ArrayMath.variance(y);
                        for (int i = 0; i < y.length; i++) {
                            y[i] -= meanY;
                        }


                        //Get regression coefficient:

                        double[] rc = Regression.getLinearRegressionCoefficients(x, y);

                        double correlation = JSci.maths.ArrayMath.correlation(x, y);
                        double propExplainedVarianceTrait = correlation * correlation - 1.0d / (double) y.length;
                        if (propExplainedVarianceTrait < 0) {
                            propExplainedVarianceTrait = 0;
                        }
                        explainedVariancePerEQTLProbe[d][(int) Math.round(propExplainedVarianceTrait * 100d)]++;

                        //Make copy of this particular eQTL:
                        double[] rawDataUpdated = new double[totalGGSamples];

                        if (nrSamplesWGenotypeData == totalGGSamples) {

                            //Regress out eQTL affect in linear regression way:
                            for (int s = 0; s < totalGGSamples; s++) {
                                double residual = y[s] - x[s] * rc[0];
                                rawDataUpdated[s] = residual;
                            }

                        } else {

                            //Correct for missing genotypes, first determine average genotype of called samples:
//                           int[] indWGA = currentDataset.getWGAGenotypeIDs();
                            double[] genotypes = xs.get(0);
                            for (int s = 0; s < totalGGSamples; s++) {
                                int ind = expressionToGenotypeId[s];
                                if (ind != -1) {
                                    double valX = genotypes[ind]; // currentSNP.getGenotypes()[ind];
//									if (valX == -1) {
//										valX = 0;
//									} else {
//										valX -= meanX;
//									}
                                    rawDataUpdated[s] = (double) rawData[p][s] - valX * rc[0];
                                }
                            }

                        }

                        //Make mean and standard deviation of residual gene expression identical to what it was before:
                        double meanUpdated = JSci.maths.ArrayMath.mean(rawDataUpdated);
                        double stdDevRatio = JSci.maths.ArrayMath.standardDeviation(rawDataUpdated) / Math.sqrt(varianceY);
                        for (int s = 0; s < totalGGSamples; s++) {
                            rawDataUpdated[s] -= meanUpdated;
                            rawDataUpdated[s] /= stdDevRatio;
                            rawDataUpdated[s] += meanY;
                        }
                        System.arraycopy(rawDataUpdated, 0, rawData[p], 0, totalGGSamples);
                        nrEQTLsRegressedOut[d]++;
//                    } else if (eventualListOfEQTLs.size() > 1 && !dosageInformationPresentForAllDatasets) {
//                        System.err.println("Multiple linear regression is not supported for datasets that do not have dosage information.");
//                        System.exit(-1);
//                    } else if (eventualListOfEQTLs.size() > 1 && dosageInformationPresentForAllDatasets) {

                    } else if (eventualListOfEQTLs.size() > 1) {

                        // use multiple linear regression via PCA
                        hashEQTLsMultipleRegressionRegressedOut.addAll(eventualListOfEQTLs);

                        int nrSNPs = snpsForProbe.size();
                        int totalGGSamples = currentDataset.getTotalGGSamples();

                        // use PCA
                        // Multiple SNPs need to be regressed out. Get SNP genotype values, standardize mean and std dev for each of these:
                        double[][] dataMatrix = new double[nrSNPs][0];
                        for (int i = 0; i < dataMatrix.length; i++) {
                            dataMatrix[i] = xs.get(i);
                        }

                        //Calculate covariance matrix:
                        double[][] correlationMatrix = new double[nrSNPs][nrSNPs];
                        int[][] sampleSizeMatrix = new int[nrSNPs][nrSNPs];
                        double sampleCountMinusOne = totalGGSamples - 1;

                        for (int f = 0; f < nrSNPs; f++) {
                            for (int g = f; g < nrSNPs; g++) {
                                // correct for missing values
                                int nrcalledforbothsnps = 0;

                                double covarianceInterim = 0;
                                for (int h = 0; h < totalGGSamples; h++) {

                                    // assume missing values are corrected for here.
//									if (dataMatrix[f][h] > -1 && dataMatrix[g][h] > -1) {
                                    covarianceInterim += dataMatrix[f][h] * dataMatrix[g][h];
                                    nrcalledforbothsnps++;
//									}
                                }

                                double covariance = covarianceInterim / (nrcalledforbothsnps - 1);
                                correlationMatrix[f][g] = covariance;
                                correlationMatrix[g][f] = covariance;
                                sampleSizeMatrix[f][g] = nrcalledforbothsnps;
                                sampleSizeMatrix[g][f] = nrcalledforbothsnps;
//								System.out.println(f + "\t" + g + "\t" + covariance + "\t" + nrcalledforbothsnps);
                            }
                        }

//						System.exit(-1);
                        //Perform eigenvalue decomposition:
                        Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationMatrix);
                        double[][] eigenArrayLists = new double[correlationMatrix.length][correlationMatrix.length];
                        for (int pca = 0; pca < nrSNPs; pca++) {
                            eigenArrayLists[pca] = PCA.getEigenVector(eig, pca);
                        }

                        //Calculate principal component scores:
                        double[][] dataMatrixPCScores = new double[nrSNPs][totalGGSamples];
                        for (int sample = 0; sample < totalGGSamples; sample++) {
                            for (int pca = 0; pca < nrSNPs; pca++) {
                                for (int snp = 0; snp < nrSNPs; snp++) {
                                    double probeCoefficient = eigenArrayLists[pca][snp];
                                    dataMatrixPCScores[pca][sample] += dataMatrix[snp][sample] * probeCoefficient;
                                }
                            }
                        }

                        //Orthogonal PCAs have been determined, use these to regress out the effects on gene expression levels:
                        TriTyperExpressionData expresionData = currentDataset.getExpressionData();
                        double[][] rawData = currentDataset.getExpressionData().getMatrix();

                        //Check what the number of samples is with genotype data available:
                        double[] y = new double[totalGGSamples];

                        //All genotypes have been succesfully called:
                        double meanYOriginal = expresionData.getProbeMean()[p];
                        double varianceYOriginal = expresionData.getProbeVariance()[p];

                        System.arraycopy(rawData[p], 0, y, 0, totalGGSamples);

                        boolean[] regressOutPCA = new boolean[nrSNPs];
                        double[] eigenValues = eig.getRealEigenvalues();
                        boolean atLeastOnePCANotRegressedOut = false;
                        for (int pca = 0; pca < nrSNPs; pca++) {
                            regressOutPCA[pca] = true;
                            if (PCA.getEigenValueVar(eigenValues, pca) < 1e-10) {
                                regressOutPCA[pca] = false;
                                atLeastOnePCANotRegressedOut = true;
                            }
                        }

                        if (atLeastOnePCANotRegressedOut) {
                            //Provide information on the PCAs:
                            System.out.println("There is at least one PCA that has not been regressed out as it does not explain a lot of genetic variation!:");
                            for (int pca = 0; pca < nrSNPs; pca++) {
                                double[] x = dataMatrixPCScores[pca];
                                double correlation = JSci.maths.ArrayMath.correlation(x, y);
                                double r2 = correlation * correlation;
                                int pcaNr = pca + 1;
                                String snpsStronglyCorrelatedWithPCA = "";
                                for (int snp = 0; snp < nrSNPs; snp++) {
                                    double correlationPCASNP = Math.abs(JSci.maths.ArrayMath.correlation(x, dataMatrix[snp]));
                                    double r2PCASNP = correlationPCASNP * correlationPCASNP;
                                    if (r2PCASNP > 0.1) {
                                        snpsStronglyCorrelatedWithPCA += "\t" + snpsForProbe.get(snp).getName() + ", " + r2PCASNP;
                                    }
                                }
                                System.out.println(probes[p] + "\tPCA" + pcaNr + "\tExplainedVariance:\t" + PCA.getEigenValueVar(eigenValues, pca) + "\tEigenvalue:\t" + eigenValues[eigenValues.length - 1 - pca] + "\tPCATraitR2:\t" + r2 + "\tSNPsStronglyCorrelatedWithPCA:\t" + snpsStronglyCorrelatedWithPCA);
                                System.out.println("PC" + pcaNr + " correlation with phenotype: " + JSci.maths.ArrayMath.correlation(y, dataMatrix[pca]));
                            }
                            System.out.println("");
                        }

                        //Process each PC, determine total amount of variation explained by the combination of PCs:
                        double propExplainedVarianceTrait = 0;
                        for (int pca = 0; pca < nrSNPs; pca++) {
                            if (regressOutPCA[pca]) {
                                //Get PCA scores:
                                double[] x = dataMatrixPCScores[pca];
                                //Get correlation coefficient with trait:
                                double correlation = JSci.maths.ArrayMath.correlation(x, y);
                                propExplainedVarianceTrait += correlation * correlation - 1.0d / (double) y.length;
                            }
                        }
                        if (propExplainedVarianceTrait < 0) {
                            propExplainedVarianceTrait = 0;
                        } else if (propExplainedVarianceTrait > 1) {

                            System.out.println("Warning: proportion of explained variance > 100%\n" +
                                    "Explained variance: " + propExplainedVarianceTrait + "\tDataset: " + gg[d].getSettings().name + "\tTrait: " + probes[p]);
                            System.out.println("Affected QTL: ");
                            for (EQTL s : covariatesForThisProbe) {
                                System.out.println(s.getRsName() + "-" + s.getProbe());
                            }
                            System.exit(-1);
                            propExplainedVarianceTrait = 1;
                        }
                        try {
                            explainedVariancePerEQTLProbe[d][(int) Math.round(propExplainedVarianceTrait * 100d)]++;
                        } catch (ArrayIndexOutOfBoundsException e) {
                            e.printStackTrace();
                            System.out.println("Error: propExplainedVarianceTrait == " + propExplainedVarianceTrait);
                            System.exit(-1);
                        }

                        //Regress out PC effects on trait:
                        for (int pca = 0; pca < nrSNPs; pca++) {
                            if (regressOutPCA[pca]) {
                                //Get PC scores:
                                double[] x = dataMatrixPCScores[pca];

                                //Get regression coefficient:
                                double[] rc = Regression.getLinearRegressionCoefficients(x, y);

                                //Regress out eQTL affect in linear regression way:
                                for (int s = 0; s < totalGGSamples; s++) {
                                    y[s] = y[s] - x[s] * rc[0];
                                }
                            }
                        }

                        double meanYUpdated = JSci.maths.ArrayMath.mean(y);
                        double varianceYUpdated = JSci.maths.ArrayMath.variance(y);

                        //Make mean and standard deviation of residual gene expression identical to what it was before:
                        double stdDevRatio = Math.sqrt(varianceYUpdated) / Math.sqrt(varianceYOriginal);
                        for (int s = 0; s < totalGGSamples; s++) {
                            y[s] -= meanYUpdated;
                            y[s] /= stdDevRatio;
                            y[s] += meanYOriginal;
                        }

                        //Replace original expression data with updated residual gene expression data:
                        for (int s = 0; s < totalGGSamples; s++) {
                            if (Double.isNaN(y[s])) {
                                System.out.println("Error!:\t" + probes[p] + "\t" + gg[d].getSettings().name + "\t" + s + "\t" + meanYUpdated + "\t" + stdDevRatio + "\t" + meanYOriginal);
                            }
                            rawData[p][s] = y[s];
                        }


                        if (false == true) {
                            for (int pca = 0; pca < nrSNPs; pca++) {
                                double[] x = dataMatrixPCScores[pca];
                                double correlation = JSci.maths.ArrayMath.correlation(x, y);
                                double r2 = correlation * correlation;
                                int pcaNr = pca + 1;
                                String snpsStronglyCorrelatedWithPCA = "";
                                for (int snp = 0; snp < nrSNPs; snp++) {
                                    double correlationPCASNP = Math.abs(JSci.maths.ArrayMath.correlation(x, dataMatrix[snp]));
                                    double r2PCASNP = correlationPCASNP * correlationPCASNP;
                                    if (r2PCASNP > 0.1) {
                                        snpsStronglyCorrelatedWithPCA += "\t" + snpsForProbe.get(snp).getName() + ", " + r2PCASNP;
                                    }
                                }
                                System.out.println(probes[p] + "\tPCA" + pcaNr
                                        + "\tExplainedVariance:\t" + PCA.getEigenValueVar(eigenValues, pca)
                                        + "\tEigenvalue:\t" + eigenValues[eigenValues.length - 1 - pca]
                                        + "\tPCATraitR2:\t" + r2
                                        + "\tSNPsStronglyCorrelatedWithPCA:\t" + snpsStronglyCorrelatedWithPCA);

                                System.out.println("PC" + pcaNr + " correlation with phenotype: " + JSci.maths.ArrayMath.correlation(y, dataMatrix[pca]));
                                System.out.println("SNP" + pcaNr + " correlation with phenotype: " + JSci.maths.ArrayMath.correlation(y, dataMatrix[pca]));
                            }
                        }
                        nrEQTLsRegressedOut[d]++;

                    }

                    for (SNP s : snpsForProbe) {
                        s.clearGenotypes();
                    }
                }
                pb.iterate(d);
                if (p % 10000 == 0) {
                    pb.display();
                }
            }
            pb.complete(d);
            System.out.println("");
        });
        pb.allCompleted();

        for (int ds = 0; ds < gg.length; ds++) {
//            gg[ds].getExpressionData().calcMeanAndVariance();
            ggSNPLoaders[ds].close();
            ggSNPLoaders[ds] = null;
        }

        System.out.println("\n");
        System.out.println("eQTLs regressed per dataset:");
        for (int d = 0; d < gg.length; d++) {
            System.out.println(gg[d].getSettings().name + "\t" + nrEQTLsRegressedOut[d]);
        }

        String output;
        System.out.println("\n");
        System.out.println("Proportion explained variance of genotypic variation on eQTLs per dataset:");


        output = "r2";
        for (TriTyperGeneticalGenomicsDataset gg1 : gg) {
            output += "\t" + gg1.getSettings().name;
        }

        System.out.println(output);
        for (int e = 0; e <= 100; e++) {
            double r2 = (double) e / 100;
            output = String.valueOf(r2);
            for (int d = 0; d < gg.length; d++) {
                output += "\t" + explainedVariancePerEQTLProbe[d][e];
            }
            System.out.println(output);
        }

    }


    private int logiter = 1;
    private String logdir = null;

    public void setLog(String logdir, int iteration) {
        logiter = iteration;
        this.logdir = logdir;
    }

    private void regressOLS(EQTL[] eqtlsToRegressOut, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {
        //Inventorize whether for a particular probe there are multiple SNPs that we want to regress out:
        HashMap<String, ArrayList<EQTL>> hashProbesCovariates = new HashMap<String, ArrayList<EQTL>>();
        HashMap<EQTL, Integer> hashEQTLIds = new HashMap<EQTL, Integer>();
        int nrProbesWithMultipleCovariates = 0;

        for (int v = 0; v < eqtlsToRegressOut.length; v++) {
            EQTL current = eqtlsToRegressOut[v];
            hashEQTLIds.put(current, v);
            String probe = current.getProbe();

            if (!hashProbesCovariates.containsKey(probe)) {
                ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
                eqtls.add(current);
                hashProbesCovariates.put(probe, eqtls);
            } else {
                hashProbesCovariates.get(probe).add(current);

            }
        }

        if (nrProbesWithMultipleCovariates > 0) {
            System.out.println("There are:\t" + hashProbesCovariates.size() + "\tprobes with > 1 covariate. Regression will be conducted through multiple regression employing OLS.");
        }

        // remove the eqtl effects
        System.out.println("Removing eQTLs:");
        int[] nrEQTLGenesRegressedOut = new int[gg.length];
        int[] nrEQTLsRegressedOut = new int[gg.length];
        int[][] explainedVariancePerEQTLProbe = new int[gg.length][101];

        SNPLoader[] ggSNPLoaders = new SNPLoader[gg.length];
        for (int d = 0; d < gg.length; d++) {
            ggSNPLoaders[d] = gg[d].getGenotypeData().createSNPLoader(1);
        }

        //Remove multiple SNPs acting on one single probe:
        MultiThreadProgressBar pb = new MultiThreadProgressBar(gg.length);


        // parallelize when multiple datasets present
        IntStream.range(0, gg.length).parallel().forEach(d -> {
//			HashSet<EQTL> hashEQTLsMultipleRegressionRegressedOut = new HashSet<EQTL>();
            HashMap<Integer, Boolean> snpPassesQC = new HashMap<Integer, Boolean>();

            SNPLoader currentloader = ggSNPLoaders[d];
            TriTyperGeneticalGenomicsDataset currentDataset = gg[d];
            String[] probes = currentDataset.getExpressionData().getProbes();
            System.out.print("Dataset:\t" + currentDataset.getSettings().name);
            pb.setSubtasks(d, hashProbesCovariates.size());

            TextFile logout = null;
            if (logdir != null) {
                try {
                    logout = new TextFile(logdir + currentDataset.getSettings().name + "-RegressionLog-Iteration" + logiter + ".txt.gz", TextFile.W);
                    System.out.println("Logging dataset " + currentDataset.getSettings().name + " to " + logout.getFileName());
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

//			for (int p = 0; p < probes.length; p++) {
            int qtlctr = 0;
            for (Map.Entry<String, ArrayList<EQTL>> eqtl : hashProbesCovariates.entrySet()) {
                String gene = eqtl.getKey();

                // check whether gene is in dataset
                Integer geneId = currentDataset.getExpressionData().getProbeToId().get(gene);
                if (geneId == -9) {
                    if (logout != null) {
                        try {
                            ArrayList<EQTL> covariatesForThisProbe = eqtl.getValue();
                            for (EQTL e : covariatesForThisProbe) {
                                logout.writeln(gene + "\t" + e.getRsName() + "\tGene not present");
                            }
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                } else {
                    ArrayList<EQTL> covariatesForThisProbe = eqtl.getValue();
                    ArrayList<EQTL> eventualListOfEQTLs = new ArrayList<EQTL>();
                    ArrayList<SNP> snpsForProbe = new ArrayList<SNP>();
                    ArrayList<double[]> xs = new ArrayList<double[]>();

                    // preload SNPs
                    for (EQTL e : covariatesForThisProbe) {
                        Integer snpId = currentDataset.getGenotypeData().getSnpToSNPId().get(e.getRsName());
                        if (snpId <= -1) {
                            if (logout != null) {
                                try {
                                    logout.writeln(gene + "\t" + e.getRsName() + "\tSNP not present");
                                } catch (IOException ex) {
                                    ex.printStackTrace();
                                }
                            }
                        } else if (snpId >= 0 && (snpPassesQC.get(snpId) == null || snpPassesQC.get(snpId))) {
                            // load SNP, if we have not seen it before
                            SNP currentSNP = currentDataset.getGenotypeData().getSNPObject(snpId);
                            try {
                                currentloader.loadGenotypes(currentSNP);
                            } catch (IOException ex) {
                                ex.printStackTrace();
                                System.exit(-1);
                            }

                            if (currentloader.hasDosageInformation()) {
                                try {
                                    currentloader.loadDosage(currentSNP);
                                } catch (IOException ex) {
                                    ex.printStackTrace();
                                    System.exit(-1);
                                }
                            }

                            if (!currentSNP.passesQC()) {
                                if (logout != null) {
                                    try {
                                        logout.writeln(gene + "\t" + e.getRsName() + "\tSNP failed QC.\tMAF: " + currentSNP.getMAF() + "\tHWEP: " + currentSNP.getHWEP() + "\tCR: " + currentSNP.getCR());
                                    } catch (IOException ex) {
                                        ex.printStackTrace();
                                    }
                                }
                                snpPassesQC.put(snpId, false);
                                currentSNP.clearGenotypes();

                            } else {
                                int[] indWGA = currentDataset.getExpressionToGenotypeIdArray();

                                double[] x = currentSNP.selectGenotypes(indWGA, true, true);
                                double meanX = 0;
                                int ctr = 0;

                                // calculate mean and variance for genotype, taking into account missing values, if any

//								double[] xorig = new double[x.length];
                                for (int i = 0; i < x.length; i++) {
//									xorig[i] = x[i];
                                    if (x[i] != -1) {
                                        meanX += x[i];
                                        ctr++;
                                    }
                                }

                                meanX /= ctr;

                                double varianceX = 0.0;
                                for (int i = 0; i < x.length; i++) {
                                    if (x[i] != -1) {
                                        double xmeanx = x[i] - meanX;
                                        varianceX += (xmeanx) * (xmeanx);
                                    }
                                }
                                varianceX /= (ctr - 1);


                                if (varianceX == 0 || currentDataset.getTotalGGSamples() != x.length) {
                                    if (logout != null) {
                                        try {
                                            logout.writeln(gene + "\t" + e.getRsName() + "\tSNP has zero variance or wrong nr of inds.\tVariance: " + varianceX + "\tinds: " + x.length + ", expected: " + currentDataset.getTotalGGSamples() + "\tMAF: " + currentSNP.getMAF() + "\tHWEP: " + currentSNP.getHWEP() + "\tCR: " + currentSNP.getCR());
                                        } catch (IOException ex) {
                                            ex.printStackTrace();
                                        }
                                    }
                                    snpPassesQC.put(snpId, false);
                                    currentSNP.clearGenotypes();
                                } else {
                                    // replace missing values with mean
//									System.out.println("Meanx: " + meanX + "\tVarX: " + varianceX);
//									System.out.println("Genotypes:");
                                    for (int i = 0; i < x.length; i++) {
                                        if (x[i] == -1) {
                                            x[i] = meanX;
                                        }
//										System.out.println(i + "\t" + xorig[i] + "\t" + x[i]);
//											if (x[i] != -1) {
//												x[i] -= meanX;
//											}
//											else {
//												x[i] = 0; // replace missing values with mean (which is now 0)
//											}
                                    }
//									System.out.println();
                                    eventualListOfEQTLs.add(e);
                                    snpsForProbe.add(currentSNP);
                                    xs.add(x);
                                    snpPassesQC.put(snpId, true);
                                }
                            }
                        }

                    }

                    // use OLS for both single variants as well as multiple variants.
                    if (xs.isEmpty()) {
                        if (logout != null) {
                            try {
                                logout.writeln(gene + "\thas 0 SNPs");
                            } catch (IOException ex) {
                                ex.printStackTrace();
                            }
                        }
                    } else {
                        // setup design matrix
                        int nrIndividuals = xs.get(0).length;
                        DoubleMatrixDataset<String, String> xcovars = new DoubleMatrixDataset<>(xs.size(), xs.get(0).length);
                        for (int i = 0; i < nrIndividuals; i++) {
                            xcovars.getHashCols().put("Ind" + i, i);
                        }

                        for (int i = 0; i < xs.size(); i++) {
                            double[] vals = xs.get(i);
                            xcovars.getRow(i).assign(vals);
                            String snpid = snpsForProbe.get(i).getName();
                            xcovars.getHashRows().put(snpid, i);
                        }

                        // transpose (samples should be on rows)
                        xcovars = xcovars.viewDice();

                        // check whether there are more predictors than data rows
                        if (xcovars.rows() < xcovars.columns()) {
                            // remove the rows with lowest variance
                            int toRemove = (xcovars.columns() - xcovars.rows()) + 1;
                            System.out.println("\nWarning: " + gg[d].getSettings().name + " has more predictors than datapoints for gene " + gene + ": " + xcovars.rows() + "x" + xcovars.columns() + " removing " + toRemove + " lowest variance covars");
                            xcovars = removeCovarWithLowestVariance(xcovars, toRemove);
                            System.out.println("\nWarning: " + gg[d].getSettings().name + " gene had few covars for " + gene + ". Remaining covars: " + xcovars.rows() + "x" + xcovars.columns());
                        }

                        try {
                            // prevent aliasing; correct for variance inflation.
                            if (xcovars.columns() > 1) {
                                VIF vif = new VIF();
                                int prevCovars = xcovars.columns();
                                xcovars = vif.vifCorrect(xcovars, (1 - 1E-4));
                                int currentCovars = xcovars.columns();
                                if (logout != null) {
                                    logout.writeln(gene + "\t had " + prevCovars + " before VIF, and " + currentCovars + " after.");
                                }
                            }


                            // prepare expression
                            int[] expressionToGenotypeId = currentDataset.getExpressionToGenotypeIdArray();
                            double[][] rawData = currentDataset.getExpressionData().getMatrix();
                            double meanY;
                            double varianceY;
                            int nrSamplesWGenotypeData = nrIndividuals;

                            double[] y = new double[nrSamplesWGenotypeData];
                            int totalGGSamples = currentDataset.getTotalGGSamples();

                            // Copy expression data.
                            int itr = 0;
                            for (int s = 0; s < rawData[geneId].length; s++) {
                                int genotypeId = expressionToGenotypeId[s];

                                // there should not be any missing values at this point..
                                if (currentDataset.getGenotypeData().getIsIncluded()[genotypeId]) {
                                    double dVal = rawData[geneId][s];
                                    y[itr] = dVal;
                                    itr++;
                                }
                            }

                            //Normalize/center subset of data:
                            meanY = JSci.maths.ArrayMath.mean(y);
                            varianceY = JSci.maths.ArrayMath.variance(y);
                            if (Double.isNaN(meanY) || Double.isNaN(varianceY)) {

                                System.err.println("ERROR: variance " + varianceY + " mean " + meanY + " for gene " + gene);
                                if (logout != null) {
                                    logout.writeln("ERROR: variance " + varianceY + " mean " + meanY + " for gene " + gene);
                                }
                                System.exit(-1);
                            }

                            for (int i = 0; i < y.length; i++) {
                                y[i] -= meanY;
                            }

                            OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();


                            double[] rawDataUpdated = null;
                            boolean singular = true;
                            double rsq = 0;
                            double[][] covars = xcovars.getMatrixAs2dDoubleArray();
                            while (singular) {
                                if (covars[0].length > 0) {
                                    ols.newSampleData(y, covars);
                                    try {
                                        // use OLS to determine regression coefficients
                                        rawDataUpdated = ols.estimateResiduals();
                                        singular = false;
                                        rsq = ols.calculateRSquared(); // I'm assuming this is an appropriate approximation of the explained variance.
                                        // debug: check whether residuals are correlated to genotype?
                                    } catch (SingularMatrixException e) {
                                        // remove lowest variance covariate
                                        // covars has samples on rows, covars on cols

                                        System.err.println("WARNING: singular matrix exception when regressing eQTLs for: " + gene + " with " + covars[0].length + " covariates (variants). Removing lowest variance covariate.");

                                        if (covars[0].length > 1) {
                                            covars = removeCovarWithLowestVariance(covars, 1);

                                        } else {
                                            System.err.println("WARNING: could not resolve covariate issue for: " + gene + " keeping original data.");
                                            if (logout != null) {
                                                logout.writeln("WARNING: could not resolve covariate issue for: " + gene + " keeping original data.");
                                            }
                                            singular = false;
                                            rsq = 0;
                                        }
                                    }
                                } else {
                                    // nothing more to do, all covariates have some issue or another
                                    singular = false;
                                }
                            }

                            if (covars[0].length > 0) {
                                if (rsq < 0) {
                                    if (rsq < -1E-9) {
                                        System.out.println("Warning: large negative r-squared: " + rsq + ". MeanY: " + meanY + ", varY: " + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
                                        if (logout != null) {
                                            logout.writeln("Warning: large negative r-squared: " + rsq + ". MeanY: " + meanY + ", varY: " + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
                                        }
                                    }
                                    rsq = 0d;
                                } else if (rsq > 1) {
                                    System.out.println("Warning: r-squared > 1.0: " + rsq + ". MeanY: " + meanY + ", varY" + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
                                    if (logout != null) {
                                        logout.writeln("Warning: r-squared > 1.0: " + rsq + ". MeanY: " + meanY + ", varY" + varianceY + ", SumSqTotal: " + ols.calculateTotalSumOfSquares() + ", SumSqResid: " + ols.calculateResidualSumOfSquares());
                                    }
                                    rsq = 1d;
                                }

                                explainedVariancePerEQTLProbe[d][(int) Math.round(rsq * 100d)]++;

                                if (logout != null) {
                                    SpearmansCorrelation sp = new SpearmansCorrelation();

                                    RankDoubleArray rda = new RankDoubleArray();
                                    double[] ry = rda.rank(y);
                                    double[] correlcoeff = new double[xcovars.columns()];

//								String[] indsY = currentDataset.getExpressionData().getIndividuals();
//								String[] indsX = currentDataset.getGenotypeData().getIndividuals();

//								for (int q = 0; q < y.length; q++) {
//									String outln = indsY[q] + "\t" + indsX[q] + "\t" + indsX[currentDataset.getExpressionToGenotypeIdArray()[q]] + "\t" + y[q] + "\t" + ry[q];
//									for (int xc = 0; xc < xcovars.columns(); xc++) {
//										outln += "\t" + xcovars.getElementQuick(q, xc);
//									}
//									System.out.println(q + "\t" + outln);
//								}
                                    for (int c = 0; c < xcovars.columns(); c++) {
                                        correlcoeff[c] = sp.correlation(xcovars.getCol(c).toArray(), ry);
                                    }
                                    String logln = gene + "\tNr SNPs: " + xcovars.columns() + "\tMeanY: " + meanY + "\tVarY: " + varianceY + "\trsq: " + rsq + "\tcorrel: " + Strings.concat(correlcoeff, Strings.tab);
                                    logout.writeln(logln);
//								System.out.println(logln);
//								System.exit(-1);
                                }


                                //Make mean and standard deviation of residual gene expression identical to what it was before:
                                double meanUpdated = JSci.maths.ArrayMath.mean(rawDataUpdated);
                                double stdDevRatio = JSci.maths.ArrayMath.standardDeviation(rawDataUpdated) / Math.sqrt(varianceY);

                                if (!Double.isNaN(meanUpdated) && !Double.isNaN(stdDevRatio) && stdDevRatio > 0) {
                                    for (int s = 0; s < totalGGSamples; s++) {
                                        rawDataUpdated[s] -= meanUpdated;
                                        rawDataUpdated[s] /= stdDevRatio;
                                        rawDataUpdated[s] += meanY;
                                    }
                                    System.arraycopy(rawDataUpdated, 0, rawData[geneId], 0, totalGGSamples);
                                    nrEQTLGenesRegressedOut[d]++;
                                    nrEQTLsRegressedOut[d] += xcovars.columns();
                                } else {
                                    String logln = "Error: " + gene + "\tNr SNPs: " + xcovars.columns() + "\tMeanY: " + meanY + "\tVarY: " + varianceY + "\trsq: " + rsq + "\tmeanUpdated: " + meanUpdated + "\tstdevRatio: " + stdDevRatio;
                                    if (logout != null) {
                                        logout.writeln(logln);
                                    }
                                }
                            }


                        } catch (Exception e) {
                            e.printStackTrace();
                        }

                        for (SNP s : snpsForProbe) {
                            s.clearGenotypes();
                        }
                    }

                }
                pb.iterate(d);
                if (qtlctr % 10000 == 0) {
                    pb.display();
                }
                qtlctr++;
            }
            pb.complete(d);
            System.out.println("");
            if (logout != null) {
                try {
                    logout.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        });
        pb.allCompleted();


        for (int ds = 0; ds < gg.length; ds++) {
//            gg[ds].getExpressionData().calcMeanAndVariance();
            ggSNPLoaders[ds].close();
            ggSNPLoaders[ds] = null;
        }

        System.out.println("\n");
        System.out.println("eQTLs regressed per dataset:");
        for (int d = 0; d < gg.length; d++) {
            System.out.println(gg[d].getSettings().name + "\tGenes: " + nrEQTLGenesRegressedOut[d] + "\tTotal eQTLs: " + nrEQTLsRegressedOut[d]);
        }

        String output;
        System.out.println("\n");
        System.out.println("Proportion explained variance of genotypic variation on eQTLs per dataset:");


        output = "r2";
        for (TriTyperGeneticalGenomicsDataset gg1 : gg) {
            output += "\t" + gg1.getSettings().name;
        }

        System.out.println(output);
        for (int e = 0; e <= 100; e++) {
            double r2 = (double) e / 100;
            output = String.valueOf(r2);
            for (int d = 0; d < gg.length; d++) {
                output += "\t" + explainedVariancePerEQTLProbe[d][e];
            }
            System.out.println(output);
        }
    }

    private DoubleMatrixDataset<String, String> removeCovarWithLowestVariance(DoubleMatrixDataset<String, String> covars, int nrToRemove) {
        ArrayList<Pair<Integer, Double>> pairs = new ArrayList<>();
        for (int col = 0; col < covars.columns(); col++) {
            DoubleMatrix1D colobj = covars.getCol(col);
            double var = Descriptives.variance(colobj.toArray());
            pairs.add(new Pair<>(col, var, Pair.SORTBY.RIGHT));
        }

        // sort ascending
        Collections.sort(pairs);
        double[][] output = new double[covars.rows()][covars.columns() - nrToRemove];
        for (int r = 0; r < covars.rows(); r++) {
            int n = 0;
            // start iterating from nrToRemove
            for (int c = nrToRemove; c < pairs.size(); c++) {
                output[r][n] = covars.getElement(r, pairs.get(c).getLeft());
                n++;
            }
        }

        ArrayList<String> newCols = new ArrayList<>();
        for (int c = nrToRemove; c < pairs.size(); c++) {
            newCols.add(covars.getColObjects().get(pairs.get(c).getLeft()));
        }

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<>();
        ds.setMatrix(output);
        try {
            ds.setRowObjects(covars.getRowObjects());
            ds.setColObjects(newCols);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ds;
    }

    private double[][] removeCovarWithLowestVariance(double[][] covars, int nrToRemove) {
        int minCovar = -1;
        double minVar = Double.MAX_VALUE;

        int startlen = covars[0].length;
        int endlen = startlen - nrToRemove;
        if (endlen <= 1) {
            return covars;
        }
        while (covars[0].length > endlen) {
            for (int c = 0; c < covars[0].length; c++) {
                double[] col = new double[covars.length];
                for (int r = 0; r < covars[0].length; r++) {
                    col[r] = covars[r][c];
                }
                double var = Descriptives.variance(col);
                if (var < minVar) {
                    minVar = var;
                    minCovar = c;
                }
            }

            double[][] tmpcovars = new double[covars.length][covars[0].length - 1];
            int cctr = 0;
            for (int c = 0; c < covars[0].length; c++) {
                if (c != minCovar) {
                    for (int r = 0; r < covars[0].length; r++) {
                        tmpcovars[r][cctr] = covars[r][c];
                    }
                    cctr++;
                }
            }
            covars = tmpcovars;
        }
        return covars;
    }
}
