/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.qcpca;

import eqtlmappingpipeline.graphics.ScatterPlot;
import umcg.genetica.math.PCA;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.QuantileNormalization;

/**
 * @author harmjan
 */
public class QCPCA {

    private boolean useCorrelationMatrix = false;
    private boolean LDpruning = false;

    // test
    public void run(String expressionLoc, String expressionPlatform, String genotypeLoc, String gte, String probeannotation, String outputdirectory, String prunedSNPListFile) {
        try {
            if (!outputdirectory.endsWith("/")) {
                outputdirectory += "/";
            }
            Gpio.createDir(outputdirectory);

            TriTyperGeneticalGenomicsDatasetSettings settings = new TriTyperGeneticalGenomicsDatasetSettings();
            settings.expressionLocation = expressionLoc;
            settings.expressionplatform = expressionPlatform;
            settings.genotypeLocation = genotypeLoc;
            settings.genotypeToExpressionCoupling = gte;
            settings.probeannotation = probeannotation;


            TriTyperGeneticalGenomicsDataset ds = new TriTyperGeneticalGenomicsDataset(settings);
            SNPLoader loader = ds.getGenotypeData().createSNPLoader();


            ArrayList<Integer> ldSNPs;
            int numsamples = ds.getTotalGGSamples();
            int[] indWGA = ds.getExpressionToGenotypeIdArray();

            if (prunedSNPListFile == null && LDpruning) {
                ldSNPs = pruneSNPsByLDThreshold(ds, loader);
            } else if (prunedSNPListFile != null) {
                ldSNPs = loadPrunedSNPListFromFile(prunedSNPListFile, ds);
            } else {
                // PCA based SNP window multiple regression pruning
                ldSNPs = pruneSNPsByMLRegressionPCA(ds, loader);
            }


            System.out.println("Copying data to array");
            ProgressBar pb = new ProgressBar(ldSNPs.size());
            double[][] datatmp = new double[numsamples][ldSNPs.size()];
            HashSet<Integer> snpsWoData = new HashSet<Integer>();
            for (int i = 0; i < ldSNPs.size(); i++) {
                Integer snpID = ldSNPs.get(i);
//                System.out.println(snpID);
                SNP snpObj = ds.getGenotypeData().getSNPObject(snpID);
                loader.loadGenotypes(snpObj);
                if (loader.hasDosageInformation()) {
                    loader.loadDosage(snpObj);
                }
                double[] snpdata = getSNPData(loader, numsamples, indWGA, snpObj, false);

                if (snpdata != null) {
                    for (int j = 0; j < snpdata.length; j++) {
                        datatmp[j][i] = snpdata[j];
                    }
                } else {
                    snpsWoData.add(i);
                }
                pb.iterate();

            }
            pb.close();

            double[][] datafinal = null;
            if (snpsWoData.size() > 0) {
                System.out.println("Detected " + snpsWoData.size() + " SNPs not passing QC, out of " + ldSNPs.size());
                int numsnpswdata = ldSNPs.size() - snpsWoData.size();
                datafinal = new double[numsamples][numsnpswdata];
                int snpcounter = 0;
                for (int i = 0; i < ldSNPs.size(); i++) {
                    if (!snpsWoData.contains(i)) {
                        for (int j = 0; j < datatmp.length; j++) {
                            datatmp[j][snpcounter] = datatmp[j][i];
                        }
                        snpcounter++;
                    }
                }
            } else {
                datafinal = datatmp;
            }

            double[][] correlationmatrix = new double[numsamples][numsamples];


            if (useCorrelationMatrix) {
                correlationmatrix = calculatecorrelationmatrix(datafinal, true);
            } else {

                for (int i = 0; i < numsamples; i++) {
                    double[] snpsi = datafinal[i];
                    for (int j = i + 1; j < numsamples; j++) {
                        double[] snpsj = datafinal[j];
                        int nrSNPsWithGenotypeDataAvailableForBothSamples = 0;
                        double ibsCount = 0;
                        for (int s = 0; s < snpsi.length; s++) {
                            if (snpsi[s] != -1 && snpsj[s] != -1) {
                                double ibsVal = 0;
                                if (snpsi[s] == snpsj[s]) {
                                    ibsVal = 1.0d;
                                } else {
                                    if (Math.abs(snpsi[s] - snpsj[s]) == 1) {
                                        ibsVal = 0.5d;
                                    }
                                }
                                ibsCount += ibsVal;
                                nrSNPsWithGenotypeDataAvailableForBothSamples++;
                            }
                        }


                        double corr = (double) ibsCount / (double) nrSNPsWithGenotypeDataAvailableForBothSamples;

                        correlationmatrix[i][j] = corr;
                        correlationmatrix[j][i] = corr;

                        pb.iterate();

                    }

                    correlationmatrix[i][i] = 1.0;
                }
            }

            TextFile corMat = new TextFile(outputdirectory + "snpcorrmat.txt", TextFile.W);
            for (int i = 0; i < correlationmatrix.length; i++) {
                String output = "";
                for (int j = 0; j < correlationmatrix.length; j++) {
                    output += "\t" + correlationmatrix[i][j];
                }
                corMat.write(output + "\n");
            }
            corMat.close();
            pb.close();

            for (int i = 0; i < 10; i++) {
                String output = "";
                for (int j = 0; j < 10; j++) {
                    output += "\t" + correlationmatrix[i][j];
                }
                System.out.println(output);
            }

            System.out.println("Performing eigenvalue decomposition");
            Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationmatrix);
            System.out.println("Getting eigenvalues");
            double[] eigenValues = PCA.getRealEigenvalues(eig);
            System.out.println("Getting eigenvariance");
            double genVarPC1 = PCA.getEigenValueVar(eigenValues, 1);
            System.out.println("Getting eigenvector");
            double[] PC1GenEigenVector = PCA.getEigenVector(eig, eigenValues, 1);
            System.out.println("Getting eigenvector");
            double[] PC2GenEigenVector = PCA.getEigenVector(eig, eigenValues, 2);

            TextFile eigenvectorsout = new TextFile(outputdirectory + "PCAOverSamplesEigenvalues.txt.gz", TextFile.W);

            double cumVarPCA = 0;
            for (int pca = 0; pca < numsamples; pca++) {
                double varPCA = PCA.getEigenValueVar(eigenValues, pca);
                int pcaNr = pca + 1;
                cumVarPCA += varPCA;
                eigenvectorsout.write(pcaNr + "\t" + varPCA + "\t" + cumVarPCA + "\n");

                System.out.println("PCA:\t" + pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + cumVarPCA);
            }
            eigenvectorsout.close();

            System.out.println("Done");
            System.out.println(genVarPC1);

            for (int i = 1; i < 11; i++) {
                ScatterPlot scat = new ScatterPlot();
                scat.draw(PCA.getEigenVector(eig, eigenValues, i), PCA.getEigenVector(eig, eigenValues, i + 1), "PC" + i, "PC" + (i + 1), "Genetic Eigenvalues", outputdirectory + "SNP-");

            }

            TextFile out = new TextFile(outputdirectory + "EigenVectors-SNPs.txt", TextFile.W);
            for (int i = 0; i < numsamples; i++) {
                String probeCoefficients = "";
                for (int pc = 1; pc <= numsamples - 1; pc++) {
                    probeCoefficients += "\t" + PCA.getEigenVector(eig, eigenValues, pc)[i];
                }
//                System.out.println(ds.getExpressionData().getIndividuals()[i]+"\t"+ds.getGenotypeData().getIndividuals()[indWGA[i]]+probeCoefficients);
                out.write(ds.getExpressionData().getIndividuals()[i] + "\t" + ds.getGenotypeData().getIndividuals()[indWGA[i]] + probeCoefficients + "\n");
            }
            out.close();

// EXPRESSION DATA!    

            DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>();
            dataset.setMatrix(ds.getExpressionData().getMatrix());
            dataset.setColObjects(Arrays.asList(ds.getExpressionData().getIndividuals()));
            dataset.setRowObjects(Arrays.asList(ds.getExpressionData().getProbes()));
            QuantileNormalization.quantilenormalize(dataset);
            Log2Transform.log2transform(dataset);

            int nrProbes = dataset.rows();
            int nrSamples = dataset.columns();
            System.out.println("Standardizing probe mean and standard deviation");
            for (int p = 0; p < dataset.rows(); p++) {
                double[] row = dataset.getRow(p).toArray();
                double mean = Descriptives.mean(row);
                double stdev = Math.sqrt(Descriptives.variance(row, mean));
                for (int s = 0; s < dataset.columns(); s++) {
                    double v = row[s];
                    dataset.setElementQuick(p,s,v-mean);
                    //                rawData[p][s] /= stdev;   // do not scale each probe for stdev: this will remove the variation captured by the
                }
            }

            System.out.println("- Standardizing sample mean and standard deviation");
            for (int s = 0; s < nrSamples; s++) {
                double[] vals = new double[nrProbes];
                for (int p = 0; p < nrProbes; p++) {
                    vals[p] = dataset.getElementQuick(p,s);
                }
                double mean = Descriptives.mean(vals);
                for (int p = 0; p < nrProbes; p++) {
                    vals[p] -= mean;
                }
                double var = Descriptives.variance(vals, mean);
                double stdev = Math.sqrt(var);
                for (int p = 0; p < nrProbes; p++) {
                    dataset.setElementQuick(p,s,(vals[p] / stdev));
                }
            }

            System.out.print("- Calculating correlations between all " + nrSamples + " samples: ");
            double[][] correlationMatrix = new double[nrSamples][nrSamples];
            double probeCountMinusOne = nrProbes - 1;

            ProgressBar pv2 = new ProgressBar(nrSamples * nrSamples);
            for (int f = 0; f < nrSamples; f++) {
                for (int g = f; g < nrSamples; g++) {

                    double covarianceInterim = 0;
                    for (int p = 0; p < nrProbes; p++) {
                        covarianceInterim += dataset.getElementQuick(p,f)* dataset.getElementQuick(p,g);
                    }
                    double covariance = covarianceInterim / probeCountMinusOne;
                    correlationMatrix[f][g] = covariance;
                    correlationMatrix[g][f] = covariance;
                    pv2.iterate();
                    pv2.iterate();
                }
            }
            pv2.close();
            System.out.println("100%");

            System.out.println("Performing eigenvalue decomposition");
            Jama.EigenvalueDecomposition eigExp = PCA.eigenValueDecomposition(correlationmatrix);
            System.out.println("Getting eigenvalues");
            double[] eigenValuesExp = PCA.getRealEigenvalues(eigExp);
            System.out.println("Getting eigenvariance");
            double expVarPC1 = PCA.getEigenValueVar(eigenValuesExp, 1);
            System.out.println("Getting eigenvector");
            double[] PC1ExpEigenVector = PCA.getEigenVector(eigExp, eigenValuesExp, 1);
            System.out.println("Getting eigenvector");
            double[] PC2ExpEigenVector = PCA.getEigenVector(eigExp, eigenValuesExp, 2);

            double[][] correlationmatrix2 = new double[numsamples][numsamples];
            pb = new ProgressBar(numsamples * numsamples);
            pb.print();

            for (int i = 1; i < 11; i++) {
                ScatterPlot scat = new ScatterPlot();
                scat.draw(PCA.getEigenVector(eigExp, eigenValuesExp, i), PCA.getEigenVector(eigExp, eigenValuesExp, i + 1), "PC" + i, "PC" + (i + 1), "Expression Eigenvalues", outputdirectory + "Exp-");
            }

            if (numsamples > 100) {
                numsamples = 100;
            }

            double bonferroni = 0.05 / (numsamples * numsamples);

            System.out.println("Determining significant correlations between genetic PCs and expression PCs");
            System.out.println("Threshold: " + bonferroni);
            for (int pc = 1; pc <= numsamples - 1; pc++) {
                double[] genEig = PCA.getEigenVector(eig, eigenValues, pc);
                for (int pc2 = pc; pc2 <= numsamples - 1; pc2++) {
                    double[] expEig = PCA.getEigenVector(eigExp, eigenValuesExp, pc2);
                    double corr = JSci.maths.ArrayMath.correlation(genEig, expEig);
                    correlationmatrix2[pc][pc2] = corr;
                    correlationmatrix2[pc2][pc] = corr;
                    ScatterPlot scat = new ScatterPlot();
                    int df = numsamples - 2;


                    double t = corr / Math.sqrt(((1 - (corr * corr)) / df));

                    // sqrt[(1—r2)/(N—2)]

                    cern.jet.random.tdouble.StudentT tDistColt = new cern.jet.random.tdouble.StudentT(df, (new cern.jet.random.tdouble.engine.DRand()));

                    double tTestPValue1 = tDistColt.cdf(t);
                    if (tTestPValue1 < bonferroni) {
                        System.out.println(corr + "\t" + t + "\t" + tTestPValue1);
                        scat.draw(genEig, expEig, "SNP" + pc, "EXP" + pc2, "SNP vs Gene expression PC" + pc + ", corr: " + corr + ", pval: " + tTestPValue1, outputdirectory + "SNPvsEXP");
                    }
                    pb.iterate();
                }
            }
            pb.close();


            out = new TextFile(outputdirectory + "SNP-PCvsExp-PC.txt", TextFile.W);
            for (int i = 0; i < numsamples; i++) {
                String probeCoefficients = "";
                for (int pc = 1; pc <= numsamples - 1; pc++) {
                    probeCoefficients += "\t" + correlationmatrix2[i][pc];
                }
//                System.out.println(ds.getExpressionData().getIndividuals()[i]+"\t"+ds.getGenotypeData().getIndividuals()[indWGA[i]]+probeCoefficients);
                out.write(i + probeCoefficients + "\n");
            }

            out.close();

        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    private double[] getSampleData(double[][] expressiondata, int sample) {
        double[] data = new double[expressiondata.length];
        for (int i = 0; i < expressiondata.length; i++) {
            data[i] = expressiondata[i][sample];
        }
        return data;
    }

    public double[] getSNPData(SNPLoader loader, int numsamples, int[] indWGA, SNP snpObj, boolean normalizegenotypedata) throws IOException {
//        System.out.println("Loading snp\t" +snpObj.getId()+"\t"+snpObj.getName());
        loader.loadGenotypes(snpObj);
        if (loader.hasDosageInformation()) {
            loader.loadDosage(snpObj);
        }
//        if (snpObj.passesQC() && snpObj.getMAF() > 0.05 && snpObj.getCR() > 0.95 && snpObj.getHWEP() > 0.0001) {
        double[] tmpData = new double[numsamples];
        if (normalizegenotypedata) {
            if (loader.hasDosageInformation()) {
                double[] genotypes = snpObj.getDosageValues();
                int numinds = 0;
                for (int i = 0; i < indWGA.length; i++) {
                    if (indWGA[i] != -1) {
                        tmpData[numinds] = genotypes[indWGA[i]];
                        numinds++;
                    }
                }


                double mean = JSci.maths.ArrayMath.mean(tmpData);
                double stdev = JSci.maths.ArrayMath.standardDeviation(tmpData);

                for (int i = 0; i < tmpData.length; i++) {
                    tmpData[i] -= mean;
                    tmpData[i] /= stdev;
                }

            } else {
                byte[] genotypes = snpObj.getGenotypes();
                int numinds = 0;

                int numIndsWithoutGenotypes = 0;

                for (int i = 0; i < indWGA.length; i++) {
                    if (indWGA[i] != -1) {
                        if (genotypes[indWGA[i]] == -1) {
                            numIndsWithoutGenotypes++;
                        }
                        tmpData[numinds] = genotypes[indWGA[i]];
                        numinds++;
                    }
                }

                double[] tmpData2 = new double[numsamples - numIndsWithoutGenotypes];
                int j = 0;
                for (int i = 0; i < tmpData.length; i++) {
                    if (tmpData[i] > 0) {
                        tmpData2[j] = tmpData[i];
                        j++;
                    }
                }

                double mean = JSci.maths.ArrayMath.mean(tmpData2);
                double stdev = JSci.maths.ArrayMath.standardDeviation(tmpData2);

                for (int i = 0; i < tmpData.length; i++) {
                    if (tmpData[i] < 0) {
                        tmpData[i] = mean;
                    }
                    tmpData[i] -= mean;
                    tmpData[i] /= stdev;
                }

            }
        } else {
            byte[] genotypes = snpObj.getGenotypes();
            int numinds = 0;

            int numIndsWithoutGenotypes = 0;

            for (int i = 0; i < indWGA.length; i++) {
                if (indWGA[i] != -1) {
                    if (genotypes[indWGA[i]] == -1) {
                        numIndsWithoutGenotypes++;
                    }
                    tmpData[numinds] = genotypes[indWGA[i]];
                    numinds++;
                }
            }

        }
        snpObj.clearGenotypes();
        return tmpData;
//        } else {
//            return null;
//        }
    }

    private ArrayList<Integer> loadPrunedSNPListFromFile(String prunedSNPListFile, TriTyperGeneticalGenomicsDataset ds) throws IOException {
        System.out.println("Loading list of pruned SNPs from text file: " + prunedSNPListFile);
        TextFile tf = new TextFile(prunedSNPListFile, TextFile.R);
        String[] list = tf.readAsArray();
        tf.close();

        ArrayList<Integer> ldSNPs = new ArrayList<Integer>();
        for (String s : list) {
            Integer snpId = ds.getGenotypeData().getSnpToSNPId().get(s);
            if (snpId != -9) {
                ldSNPs.add(snpId);
            }
        }

        System.out.println(ldSNPs.size() + " out of " + list.length + " SNPs in the pruned SNP list detected.");
        return ldSNPs;
    }

    private ArrayList<Integer> pruneSNPsByLDThreshold(TriTyperGeneticalGenomicsDataset ds, SNPLoader loader) {
        System.out.println("Pruning SNPs for LD");
        ArrayList<Integer> ldSNPs = new ArrayList<Integer>();
        DetermineLD ldcalc = new DetermineLD();
        try {
            for (int chr = 1; chr < 23; chr++) {
                HashSet<Integer> snpsVisited = new HashSet<Integer>();

                ArrayList<Integer> snpsForChr = getSortedListOfSNPsForChr(chr, ds);

                int numSNPsAfterPruning = 0;
                if (snpsForChr == null) {
                    System.out.println("No SNPs for Chr: " + chr);
                } else {
                    int startsnpnum = 0;
                    // sort SNPs..

                    ProgressBar pb = new ProgressBar(snpsForChr.size());
                    for (startsnpnum = 0; startsnpnum < snpsForChr.size(); startsnpnum++) {

                        int snpID = snpsForChr.get(startsnpnum);
                        //                        int snpID  = snpsForChr.get(s);
                        if (snpsVisited.contains(snpID)) {
                            // skip this SNP: it belongs to a different LD block
                        } else {
                            SNP snpObj = ds.getGenotypeData().getSNPObject(snpID);
                            loader.loadGenotypes(snpObj);
                            if (loader.hasDosageInformation()) {
                                loader.loadDosage(snpObj);
                            }
                            if (snpObj.passesQC() && snpObj.getMAF() > 0.05 && snpObj.getCR() > 0.95 && snpObj.getHWEP() > 0.0001) {

                                // add this SNP to the dataset for calculating correlations
                                //                                data.add();
                                ldSNPs.add(snpID);

                                for (int querysnpnum = startsnpnum + 1; querysnpnum < snpsForChr.size(); querysnpnum++) {
                                    int snpID2 = snpsForChr.get(querysnpnum);
                                    SNP snpObj2 = ds.getGenotypeData().getSNPObject(snpID2);
                                    loader.loadGenotypes(snpObj2);
                                    if (loader.hasDosageInformation()) {
                                        loader.loadDosage(snpObj2);
                                    }
                                    if (snpObj2.passesQC() && snpObj2.getMAF() > 0.05 && snpObj2.getCR() > 0.95 && snpObj2.getHWEP() > 0.0001) {
                                        double r2 = ldcalc.getRSquared(snpObj, snpObj2, ds.getGenotypeData(), ldcalc.RETURN_R_SQUARED, ldcalc.INCLUDE_CASES_AND_CONTROLS, false);
                                        double dp = ldcalc.getRSquared(snpObj, snpObj2, ds.getGenotypeData(), ldcalc.RETURN_D_PRIME, ldcalc.INCLUDE_CASES_AND_CONTROLS, false);
                                        if (r2 < 0.1 && dp < 0.5) {
                                            // 
                                            snpObj2.clearGenotypes();
                                            startsnpnum = querysnpnum;
                                            break;
                                        } else {
                                            snpObj2.clearGenotypes();
                                            snpsVisited.add(snpID2);
                                        }
                                    } else {
                                        snpObj2.clearGenotypes();
                                        snpsVisited.add(snpID2);
                                    }

                                    pb.iterate();
                                } // pairwise comparison
                            }
                            snpObj.clearGenotypes();
                        } // end HWE, MAF and CR filter
                        snpsVisited.add(snpID);
                        numSNPsAfterPruning++;
                        pb.iterate();
                    }
                    pb.close();

                    System.out.println(numSNPsAfterPruning + " SNPs left after pruning, out of " + snpsForChr.size() + "\t" + ldSNPs.size() + " total.");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(ldSNPs.size() + " pruned SNPs");
        return ldSNPs;
    }

    private ArrayList<Integer> getSortedListOfSNPsForChr(int chr, TriTyperGeneticalGenomicsDataset ds) {
        String[] snps = ds.getGenotypeData().getSNPs();
        ArrayList<Integer> snpsOnChr = new ArrayList<Integer>();
        int numsnps = 0;
        for (int i = 0; i < snps.length; i++) {
            byte snpchr = ds.getGenotypeData().getChr(i);
            if (snpchr == chr) {
                snpsOnChr.add(i);
            }
        }

        ArrayList<SortableSNP> snpsSorted = new ArrayList<SortableSNP>();
        for (Integer i : snpsOnChr) {
            int chrPos = ds.getGenotypeData().getChrPos(i);
            if (chrPos > -1) {
                snpsSorted.add(new SortableSNP(null, i, (byte) chr, chrPos, SortableSNP.SORTBY.CHRPOS));
            }
        }
        int numTotalOnChr = snpsOnChr.size();

        Collections.sort(snpsSorted);
        //check sorting
        snpsOnChr = new ArrayList<Integer>();
        for (SortableSNP s : snpsSorted) {
            snpsOnChr.add(s.id);
        }

        // check sorting
        int prevpos = -1;
        for (Integer i : snpsOnChr) {
            Integer chrPos = ds.getGenotypeData().getChrPos(i);
            if (prevpos == -1) {
                prevpos = chrPos;
            } else {
                if (chrPos >= prevpos) {
                    // its ok
                } else {
                    System.out.println("SNPs are not sorted!!");
                    for (int j = 0; j < snpsOnChr.size(); j++) {
                        Integer snpid = snpsOnChr.get(j);
                        System.out.println(j + "\t" + snpsOnChr.get(j) + "\t" + ds.getGenotypeData().getChrPos(snpid));
                    }
                    System.exit(0);
                }
            }

        }
        System.out.println("Chr " + chr + " has " + snpsOnChr.size() + " SNPs with annotation, out of " + numTotalOnChr);
        return snpsOnChr;
    }

    private ArrayList<Integer> pruneSNPsByMLRegressionPCA(TriTyperGeneticalGenomicsDataset ds, SNPLoader loader) throws IOException {
        ArrayList<Integer> ldSNPs = new ArrayList<Integer>();

        int numsamples = ds.getTotalGGSamples();
        int[] indWGA = ds.getExpressionToGenotypeIdArray();

        int windowsize = 50;
        int windowshift = 5;
        int vifthreshold = 2;

        double[][] snpdata = new double[windowsize][numsamples];
        int totalafterpruning = 0;
        for (int chr = 1; chr < 23; chr++) {
            HashSet<Integer> visitedSNPs = new HashSet<Integer>();
            ArrayList<Integer> sortedSNPs = getSortedListOfSNPsForChr(chr, ds);
            if (sortedSNPs.size() < windowsize) {
                System.out.println("Chromosome " + chr + " has less than " + windowsize + " SNPs for pruning");
            } else {

                int numwindowsremainder = sortedSNPs.size() % windowsize;
                int numwindows = (sortedSNPs.size() - numwindowsremainder) / windowsize;

                System.out.println("Pruning SNPs for chromosome: " + chr);
                ProgressBar pb = new ProgressBar(sortedSNPs.size());
//                double[][] correlationmatrix = new double[windowsize][windowsize];
                int window = 0;
                for (int startsnp = 0; startsnp + windowsize < sortedSNPs.size(); startsnp += windowshift) {
//                    System.out.println(startsnp+"\t"+sortedSNPs.size());
                    ArrayList<Integer> snpsInThisWindow = new ArrayList<Integer>();
                    int s = 0;
                    int currentsnp = startsnp;
                    boolean fullwindow = true;
                    while (s < windowsize) {    //for(int s = startsnp; s<stopsnp; s++){
                        if (currentsnp == sortedSNPs.size()) {
                            fullwindow = false;
                            break;
                        } else {
                            Integer snpid = sortedSNPs.get(currentsnp);
                            SNP snpObj = ds.getGenotypeData().getSNPObject(snpid);
                            double[] snpData = getSNPData(loader, numsamples, indWGA, snpObj, true);
                            if (snpData == null) {
                                startsnp++; // also move the next' window start one forward.
                            } else {
                                snpdata[s] = snpData;
                                snpsInThisWindow.add(snpid);
                                s++;
                            }
                            currentsnp++;
                        }
                    }
                    if (fullwindow) {
//                            double[][] pcscores = new double[windowsize][numsamples];
//                        
////                        for(int query=0; query<windowsize; query++) {
//                            double[][] correlationmatrix     = calculatecorrelationmatrix(snpdata, false);
//                            Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationmatrix);
//                            double[] eigenValues             = eig.getRealEigenvalues();
//                            double[][] eigenvectors          = new double[windowsize][windowsize];
//
//                            for(int i=0; i<windowsize; i++){
//                                eigenvectors[i] = PCA.getEigenVector(eig, i);
//                            }
//
//                            for(int pc=0; pc<windowsize; pc++){
//                                for(int snp=0; snp<windowsize; snp++){
//                                    for(int sample=0; sample<numsamples; sample++){      
//                                        double probecoefficient = eigenvectors[pc][snp];
////                                        if(snp >= query){
////                                            pcscores[pc][sample] += snpdata[snp+1][sample] * probecoefficient;
////                                        } else {
//                                            pcscores[pc][sample] += snpdata[snp][sample] * probecoefficient;
////                                        }
//                                    }
//                                }
//                            }
//                            
//                            
//                            double sum = 0;
//    //                        for(int sample=0; sample<10; sample++){
//    //                            System.out.println(snp1+"\t"+sample+"\t"+snpdata[snp1][sample]+"\t"+pcscores[snp1][sample]);
//    //                        }
//                            for(int snp1=0; snp1<windowsize; snp1++){
//                                sum = 0;
//                                for(int snp2=0; snp2<windowsize; snp2++){
//                                    if(snp2!=snp1){
//                                        double corr = JSci.maths.ArrayMath.correlation(snpdata[snp1], pcscores[snp2]);
//                                        sum+= (corr*corr);
//                                    }
//                                }
//                                System.out.println(snp1+"\t"+sum);
//                                double vif = ((double)1/Math.abs(1-sum));  // use abs to prevent rounding mistake influence
////                                if(vif < vifthreshold){
////                                    visitedSNPs.add(snpsInThisWindow.get(query));
////                                } else {
////                                    visitedSNPs.remove(snpsInThisWindow.get(query));
////                                }
//                            }


                        int tmpwindow = windowsize - 1;
                        for (int snp1 = 0; snp1 < 2; snp1++) {

                            double[][] pcscores = new double[tmpwindow][numsamples];
                            double[][] correlationmatrix = calculatecorrelationmatrix(snpdata, false, snp1);

                            Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationmatrix);
                            double[] eigenValues = eig.getRealEigenvalues();
                            double[][] eigenvectors = new double[tmpwindow][tmpwindow];
                            for (int i = 0; i < tmpwindow; i++) {
                                eigenvectors[i] = PCA.getEigenVector(eig, i);
                            }


                            for (int pc = 0; pc < tmpwindow; pc++) {
                                for (int snp = 0; snp < tmpwindow; snp++) {
                                    for (int sample = 0; sample < tmpwindow; sample++) {
                                        double probecoefficient = eigenvectors[pc][snp];
                                        if (snp >= snp1) {
                                            pcscores[pc][sample] += snpdata[snp + 1][sample] * probecoefficient;
                                        } else {
                                            pcscores[pc][sample] += snpdata[snp][sample] * probecoefficient;
                                        }
                                    }
                                }
                            }

                            double sum = 0;
//                            for(int snp1=0; snp1<windowsize; snp1++){
//                                sum = 0;
                            for (int snp2 = 0; snp2 < tmpwindow; snp2++) {

                                double corr = JSci.maths.ArrayMath.correlation(snpdata[snp1], pcscores[snp2]);
                                sum += (corr * corr);
                                System.out.println(snp1 + "\t" + snp2 + "\t" + corr + "\t" + (corr * corr));
                            }

                            System.out.println(snp1 + "\t" + sum);
                            double vif = ((double) 1 / Math.abs(1 - sum));  // use abs to prevent rounding mistake influence
//                                if(vif < vifthreshold){
//                                    visitedSNPs.add(snpsInThisWindow.get(query));
//                                } else {
//                                    visitedSNPs.remove(snpsInThisWindow.get(query));
//                                }
//                            }


                        }


//                        for(int query=0; query<windowsize; query++) {


                        //                        for(int sample=0; sample<10; sample++){
                        //                            System.out.println(snp1+"\t"+sample+"\t"+snpdata[snp1][sample]+"\t"+pcscores[snp1][sample]);
                        //                        }


//                        }
                        System.exit(0);
//                        double[][] pcscores = new double[windowsize][numsamples];
//                        
//                        for(int query=0; query<windowsize; query++){
//                            
//                            
//                            double[][] correlationmatrix     = calculatecorrelationmatrix(snpdata, false, query);
//                            Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationmatrix);
//                            double[] eigenValues             = eig.getRealEigenvalues();
//
//                            double[][] eigenvectors = new double[windowsize][windowsize];
//
//                            for(int i=0; i<windowsize-1; i++){
//                                eigenvectors[i] = PCA.getEigenVector(eig, i);
//                            }
//
//                            for(int pc=0; pc<windowsize-1; pc++){
//                                for(int snp=0; snp<windowsize-1; snp++){
//                                    for(int sample=0; sample<numsamples; sample++){
//                                        
//                                        double probecoefficient = eigenvectors[pc][snp];
//                                        if(snp >= query){
//                                            pcscores[pc][sample] += snpdata[snp+1][sample] * probecoefficient;
//                                        } else {
//                                            pcscores[pc][sample] += snpdata[snp][sample] * probecoefficient;
//                                        }
//                                        
//                                    }
//                                }
//                            }
//                            
//                            
//                            double sum = 0;
//    //                        for(int sample=0; sample<10; sample++){
//    //                            System.out.println(snp1+"\t"+sample+"\t"+snpdata[snp1][sample]+"\t"+pcscores[snp1][sample]);
//    //                        }
//                            for(int snp2=0; snp2<windowsize-1; snp2++){
////                                if(snp2!=query){
//                                    double corr = JSci.maths.ArrayMath.correlation(snpdata[query], pcscores[snp2]);
//                                    sum+= (corr*corr);
////                                }
//                            }
//                            System.out.println(query+"\t"+sum);
//                            double vif = ((double)1/Math.abs(1-sum));  // use abs to prevent rounding mistake influence
//                            if(vif < vifthreshold){
//                                visitedSNPs.add(snpsInThisWindow.get(query));
//                            } else {
//                                visitedSNPs.remove(snpsInThisWindow.get(query));
//                            }
//
//                        }
//                        System.exit(0);


                        //                    // check for orthogonality
                        //                    for(int i=0;i<pcscores.length; i++){
                        //                    
                        //                        for(int j=i+1; j<pcscores.length; j++){
                        //                            double sum = 0;
                        //                            for(int v=0; v<pcscores[i].length; v++){
                        //                                sum+= pcscores[i][v] * pcscores[j][v];
                        //                            }
                        //                            double corr = JSci.maths.ArrayMath.correlation(pcscores[i], pcscores[j]);
                        //                            if(sum > 1E-5){
                        //                                System.out.println("Innerproduct > 1E-5 for covariates "+i+" and "+j);
                        //                                
                        //                            }
                        //
                        //                            System.out.println((i+1)+"\t"+(j+1)+"\t"+sum+"\t"+corr);
                        //                        }
                        //
                        //                    }


                    }

                    pb.set(startsnp);
                }
                pb.close();
            }
            System.out.println("SNPs after pruning: " + visitedSNPs.size());
            totalafterpruning += visitedSNPs.size();
        }
        System.out.println(totalafterpruning);
        System.exit(0);
        return ldSNPs;
    }

    private double[][] calculatecorrelationmatrix(double[][] data, boolean verbose) {

        return calculatecorrelationmatrix(data, verbose, null);

    }

    private double[][] calculatecorrelationmatrix(double[][] data, boolean verbose, Integer skip) {


        ProgressBar pb = null;
        if (verbose) {
            System.out.println("Calculating correlation matrix");
            pb = new ProgressBar(data.length);
            pb.print();
        }
        int numrows = data.length;
        double[][] correlationmatrix = null;
        if (skip != null) {
            correlationmatrix = new double[numrows - 1][numrows - 1];
            int rows = 0;

            for (int i = 0; i < numrows; i++) {
                if (i != skip) {
                    int cols = 0;
                    for (int j = i + 1; j < numrows; j++) {
                        if (j != skip) {
                            double corr = JSci.maths.ArrayMath.correlation(data[i], data[j]);
                            correlationmatrix[rows][cols] = corr;
                            correlationmatrix[cols][rows] = corr;
                            cols++;
                        }
                    }
                    correlationmatrix[rows][rows] = 1.0;
                    rows++;
                }
            }
        } else {
            correlationmatrix = new double[numrows][numrows];

            for (int i = 0; i < numrows; i++) {
                for (int j = i + 1; j < numrows; j++) {
                    double corr = JSci.maths.ArrayMath.correlation(data[i], data[j]);
                    correlationmatrix[i][j] = corr;
                    correlationmatrix[j][i] = corr;
                }
                if (verbose) {
                    pb.iterate();
                }
                correlationmatrix[i][i] = 1.0;
            }
            if (verbose) {
                pb.close();
            }
        }

        return correlationmatrix;
    }
}
// case vs control analysis
//            if (1==2) {
//                
//                boolean[] isCase = new boolean[ds.getGenotypeData().getIndividuals().length];
//                int nrCases = 0;
//                int nrControls = 0;
//                for (int s=0; s<ds.getGenotypeData().getIndividuals().length; s++) {
//                    String sampleName = ds.getGenotypeData().getIndividuals()[s];
//                    System.out.println(s + "\t" + sampleName);
//                    if (sampleName.startsWith("Celiac")) {
//                        isCase[s] = true;
//                        nrCases++;
//                    } else {
//                        nrControls++;
//                    }
//                }
//                int nrSamples = nrCases + nrControls;
//                System.out.println(nrCases + "\t"+ nrControls + "\t"+ nrSamples);
//                
//                String[] allsnps = ds.getGenotypeData().getSNPs();
//                for(int a=0; a<allsnps.length; a++){
//                //for(int a=1457327; a<allsnps.length; a++){
//                    SNP snpObj = ds.getGenotypeData().getSNPObject(a);
//                    loader.loadGenotypes(snpObj);
//                    if(snpObj.passesQC() && snpObj.getMAF() > 0.05 && snpObj.getCR() > 0.95 && snpObj.getHWEP() > 0.0001){
//                        int[][] twoByTwo = new int[2][2];
//                        for (int s=0; s<nrSamples; s++) {
//                            int x = 0;
//                            if (isCase[s]) x = 1;
//                            twoByTwo[x][0]+=snpObj.getGenotypes()[s];
//                            twoByTwo[x][1]+=2 - snpObj.getGenotypes()[s];    
//                        }
//                        
//                        boolean twoByTwoPassesQC = true;
//                        for(int l=0; l<2;l++){
//                            for(int m=0; m<2;m++){
//                                if(twoByTwo[l][m] == 0){
//                                    twoByTwoPassesQC = false;
//                                }
//                            }
//                        }
//                        if(twoByTwoPassesQC){
//                            FisherExactTest fisher = new FisherExactTest();
//                            double pValue = fisher.getFisherPValue(twoByTwo[0][0], twoByTwo[0][1], twoByTwo[1][0], twoByTwo[1][1]);
//                            if (pValue < 1E-10) {
//
//
//                                System.out.println(a + "\t"+ snpObj.getName() + "\t"+ pValue + "\t" + twoByTwo[0][0] + "\t"+ twoByTwo[0][1] + "\t"+ twoByTwo[1][0] + "\t"+ twoByTwo[1][1]);
//
//                                if (a==-1) { //1457327) {
//
//                                    for (int s=0; s<nrSamples; s++) {
//
//                                        System.out.println(a + "\t" + s + "\t"+ ds.getGenotypeData().getIndividuals()[s] + "\t" + snpObj.getGenotypes()[s]);
//                                    }
//
//                                    System.exit(0);
//                                }
//
//                                short[] genotypes  = snpObj.getGenotypes();
//                                double[] tmpData = new double[genotypes.length];
//                                for (int t=0; t<numsamples; t++) {
//                                    tmpData[t] = genotypes[t];
//                                }
//                                double mean  = JSci.maths.ArrayMath.mean(tmpData);
//                                double stdev = JSci.maths.ArrayMath.standardDeviation(tmpData);
//                                for (int t=0; t<numsamples; t++) {
//                                    tmpData[t]-=mean;
//                                    tmpData[t]/=stdev;
//                                }
//                                data.add(tmpData);
//                                /*
//                                double[] tmpData = new double[numsamples];
//                                 short[] genotypes  = snpObj.getGenotypes();
//                                        int numinds = 0;
//
//                                        int numIndsWithoutGenotypes = 0;
//
//                                        for(int i=0; i<indWGA.length;i++){
//                                            if(indWGA[i] != -1){
//                                                if(genotypes[indWGA[i]] == -1){
//                                                    numIndsWithoutGenotypes++;
//                                                }
//                                                tmpData[numinds] = genotypes[indWGA[i]];
//                                                numinds++;
//                                            }
//                                        }
//
//                                        double[] tmpData2 = new double[numsamples - numIndsWithoutGenotypes];
//                                        int j=0;
//                                        for(int i=0; i<tmpData.length; i++){
//                                            if(tmpData[i] >= 0){
//                                                tmpData2[j] = tmpData[i];
//                                                j++;
//                                            }
//                                        }
//
//                                        double mean  = JSci.maths.ArrayMath.mean(tmpData2);
//                                        double stdev = JSci.maths.ArrayMath.standardDeviation(tmpData2);
//
//                                        for(int i=0; i<tmpData.length; i++){
//                                            if(tmpData[i] < 0){
//                                                tmpData[i] = mean;
//                                            }
//                                            tmpData[i] -= mean;
//                                            tmpData[i] /= stdev;
//                                        }
//                                        data.add(tmpData);
//
//                                 * 
//                                 */
//
//
//
//
//
//
//
//
//
//
//
//                            }
//                        }
//                    }
//                }
//                
//            }