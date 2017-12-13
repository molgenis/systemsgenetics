/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.PCA;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Regression;

/**
 *
 * @author harmjan
 */
public class EQTLRegression {

    TriTyperGeneticalGenomicsDataset[] gg;
    EQTL[] eqtlsToRegressOut;

    public void regressOutEQTLEffects(ArrayList<Pair<String, String>> eqtls, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {
        this.gg = gg;

        this.eqtlsToRegressOut = new EQTL[eqtls.size()];
        for (int q = 0; q < eqtls.size(); q++) {
            eqtlsToRegressOut[q] = new EQTL();
            eqtlsToRegressOut[q].setRsName(eqtls.get(q).getLeft());
            eqtlsToRegressOut[q].setProbe(eqtls.get(q).getRight());
        }
        System.out.println("About to regress out: " + eqtls.size() + " QTLs from data.");
        regressOutEQTLEffects();
    }

    public void regressOutEQTLEffects(EQTL[] eqtls, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {
        this.gg = gg;
        this.eqtlsToRegressOut = eqtls;
        System.out.println("About to regress out: " + eqtls.length + " QTLs from data.");
        regressOutEQTLEffects();
    }

    public void regressOutEQTLEffects(String regressOutEQTLEffectFileName, boolean outputfiles, TriTyperGeneticalGenomicsDataset[] gg) throws IOException {


        this.gg = gg;
        System.out.println("\n\n\nRemoving eQTL effects from the following eQTL file: '" + regressOutEQTLEffectFileName);
        QTLTextFile in = new QTLTextFile(regressOutEQTLEffectFileName, QTLTextFile.R);
        eqtlsToRegressOut = in.read();
        in.close();
        System.out.println("Number of eQTLs to regress out found in file:\t" + eqtlsToRegressOut.length);
        regressOutEQTLEffects();
        if (outputfiles) {
            for (int d = 0; d < gg.length; d++) {
                TriTyperGeneticalGenomicsDataset ds = gg[d];
                TriTyperExpressionData dsexp = ds.getExpressionData();
                double[][] matrix = dsexp.getMatrix();
                String[] probes = dsexp.getProbes();
                String[] individuals = dsexp.getIndividuals();
                String filename = ds.getSettings().expressionLocation;
                DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>(matrix, Arrays.asList(probes), Arrays.asList(individuals));
                dsout.recalculateHashMaps();
                System.out.println("Saving expression file after removal of eQTL effects: " + filename + "-EQTLEffectsRemoved.txt.gz");
                dsout.save(filename + "-EQTLEffectsRemoved.txt.gz");
            }
        }

    }

    /**
     * Removes the effect of a supplied list of eQTL from the datasets by use of
     * regression
     *
     */
    private void regressOutEQTLEffects() throws IOException {

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
        for (int d = 0; d < gg.length; d++) {
            HashSet<EQTL> hashEQTLsMultipleRegressionRegressedOut = new HashSet<EQTL>();
            HashMap<Integer, Boolean> snpPassesQC = new HashMap<Integer, Boolean>();

            TriTyperGeneticalGenomicsDataset currentDataset = gg[d];
            String[] probes = gg[d].getExpressionData().getProbes();
            System.out.print("Dataset:\t" + gg[d].getSettings().name);
            ProgressBar pgb = new ProgressBar(probes.length);

            for (int p = 0; p < probes.length; p++) {

                ArrayList<EQTL> covariatesForThisProbe = hashProbesCovariates.get(probes[p]);
    
                if (covariatesForThisProbe != null) {
                    ArrayList<EQTL> eventualListOfEQTLs = new ArrayList<EQTL>();
                    ArrayList<SNP> snpsForProbe = new ArrayList<SNP>();
                    ArrayList<double[]> xs = new ArrayList<double[]>();
                    ArrayList<Double> meanxs = new ArrayList<Double>();


                    for (EQTL e : covariatesForThisProbe) {
                        if (!hashEQTLsMultipleRegressionRegressedOut.contains(e)) {
                            Integer snpId = gg[d].getGenotypeData().getSnpToSNPId().get(e.getRsName());

                            if (snpId != -9 && (snpPassesQC.get(snpId) == null || snpPassesQC.get(snpId))) {
                                // load SNP

                                SNP currentSNP = currentDataset.getGenotypeData().getSNPObject(snpId);
                                ggSNPLoaders[d].loadGenotypes(currentSNP);

                                if (ggSNPLoaders[d].hasDosageInformation()) {
                                    ggSNPLoaders[d].loadDosage(currentSNP);
                                }
                                
                                if (currentSNP.passesQC()) {
                                    int[] indWGA = currentDataset.getExpressionToGenotypeIdArray();
                                    double[] x = currentSNP.selectGenotypes(indWGA);
                                    double meanX = JSci.maths.ArrayMath.mean(x);
                                    double varianceX = JSci.maths.ArrayMath.variance(x);
                                    
                                    if (varianceX != 0 && currentDataset.getTotalGGSamples()==x.length) {
                                        for (int i = 0; i < x.length; i++) {
                                            x[i] -= meanX;
                                        }
                                        eventualListOfEQTLs.add(e);
                                        snpsForProbe.add(currentSNP);
                                        xs.add(x);
                                        meanxs.add(meanX);
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
                        double meanX = meanxs.get(0);

                        //Get the expression data:
                        double[][] rawData = currentDataset.getExpressionData().getMatrix();
                        double meanY;
                        double varianceY;

                        //Check what the number of samples is with genotype data available:
                        int nrSamplesWGenotypeData = x.length;

                        double[] y = new double[nrSamplesWGenotypeData];
                        int totalGGSamples = currentDataset.getTotalGGSamples();
                        if (nrSamplesWGenotypeData == totalGGSamples) {
                            //All genotypes have been succesfully called, use quick approach:
                            meanY = currentDataset.getExpressionData().getProbeMean()[p];
                            varianceY = currentDataset.getExpressionData().getProbeVariance()[p];
                            for (int s = 0; s < totalGGSamples; s++) {
                                y[s] = rawData[p][s] - meanY;
                            }
                        } else {

                            //Not all genotypes have been succesfully called, use slow approach:
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
                            for (int s = 0; s < totalGGSamples; s++) {
                                int ind = expressionToGenotypeId[s];
                                if (ind != -1) {
                                    double valX = currentSNP.getGenotypes()[ind];
                                    if (valX == -1) {
                                        valX = 0;
                                    } else {
                                        valX -= meanX;
                                    }
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
                    } else if (eventualListOfEQTLs.size() > 1 && !dosageInformationPresentForAllDatasets) {
                        System.err.println("Multiple linear regression is not supported for datasets that do not have dosage information.");
                        System.exit(-1);
                    } else if (eventualListOfEQTLs.size() > 1 && dosageInformationPresentForAllDatasets) {
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
                        double sampleCountMinusOne = totalGGSamples - 1;

                        for (int f = 0; f < nrSNPs; f++) {
                            for (int g = f; g < nrSNPs; g++) {
                                double covarianceInterim = 0;
                                for (int h = 0; h < totalGGSamples; h++) {
                                    covarianceInterim += dataMatrix[f][h] * dataMatrix[g][h];
                                }
                                double covariance = covarianceInterim / (double) (sampleCountMinusOne);
                                correlationMatrix[f][g] = covariance;
                                correlationMatrix[g][f] = covariance;
                            }
                        }
                        
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
                        }
                        explainedVariancePerEQTLProbe[d][(int) Math.round(propExplainedVarianceTrait * 100d)]++;

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

                        nrEQTLsRegressedOut[d]++;

                    }

                    for (SNP s : snpsForProbe) {
                        s.clearGenotypes();
                    }
                }
                pgb.iterate();
            }
            pgb.print();
            pgb.close();
            System.out.println("");
        }

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
}
