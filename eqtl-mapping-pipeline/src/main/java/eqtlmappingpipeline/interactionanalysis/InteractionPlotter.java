/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import java.awt.Color;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.apache.commons.math3.linear.SingularMatrixException;
import umcg.genetica.containers.Triple;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class InteractionPlotter {

    public InteractionPlotter(String interactionFile, String genotypeDir, String expressionDataFile, String covariateDataFile, String gteFile, String outdir) throws IOException {
        outdir = Gpio.formatAsDirectory(outdir);
        Gpio.createDir(outdir);
        Map<String, String> gte = null;
        if (gteFile != null) {
            TextFile tf = new TextFile(gteFile, TextFile.R);
            gte = tf.readAsHashMap(0, 1);
            tf.close();
        }

        HashSet<String> expressionProbes = new HashSet<String>();
        TextFile tf = new TextFile(interactionFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<Triple<String, String, String>> triples = new ArrayList<Triple<String, String, String>>();
        while (elems != null) {
            if (elems.length == 2) {
                String snp = elems[0];
                String probe = elems[1];
                expressionProbes.add(probe);
                triples.add(new Triple<String, String, String>(snp, null, probe));
            } else if (elems.length == 3) {
                String snp = elems[0];
                String covariate = elems[1];
                String probe = elems[2];
                expressionProbes.add(probe);
                triples.add(new Triple<String, String, String>(snp, covariate, probe));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(triples.size() + " SNP - covariate - probe combinations read from: " + interactionFile);

        DoubleMatrixDataset<String, String> expressionData = new DoubleMatrixDataset<String, String>(expressionDataFile, expressionProbes);
        DoubleMatrixDataset<String, String> covariateData = new DoubleMatrixDataset<String, String>(covariateDataFile);

        int samplesHaveCovariatesOnCols = 0;
        int samplesHaveCovariatesOnRows = 0;
        for (int i = 0; i < expressionData.colObjects.size(); i++) {
            String expSample = expressionData.colObjects.get(i);
            Integer id1 = covariateData.hashCols.get(expSample);
            Integer id2 = covariateData.hashRows.get(expSample);
            if (id1 != null) {
                samplesHaveCovariatesOnCols++;
            }
            if (id2 != null) {
                samplesHaveCovariatesOnRows++;
            }
        }
        if (samplesHaveCovariatesOnRows > samplesHaveCovariatesOnCols) {
            System.out.println("Rows contain covariate samples in covariate file. Transposing covariates.");
            covariateData.transposeDataset();
        }

        TriTyperGenotypeData geno = new TriTyperGenotypeData(genotypeDir);
        SNPLoader loader = geno.createSNPLoader();

        int[] genotypeToCovariate = new int[geno.getIndividuals().length];
        int[] genotypeToExpression = new int[geno.getIndividuals().length];
        String[] genoIndividuals = geno.getIndividuals();
        for (int i = 0; i < genotypeToCovariate.length; i++) {
            String genoSample = genoIndividuals[i];
            if (geno.getIsIncluded()[i] != null && geno.getIsIncluded()[i]) {
                if (gte != null) {
                    genoSample = gte.get(genoSample);
                }

                Integer covariateSample = covariateData.hashCols.get(genoSample);
                Integer expressionSample = expressionData.hashCols.get(genoSample);
                if (genoSample != null && covariateSample != null && expressionSample != null) {
                    genotypeToCovariate[i] = covariateSample;
                    genotypeToExpression[i] = expressionSample;
                } else {
                    genotypeToCovariate[i] = -9;
                    genotypeToExpression[i] = -9;

                }
            } else {
                genotypeToCovariate[i] = -9;
                genotypeToExpression[i] = -9;
            }
        }

        OLSMultipleLinearRegression regressionFullWithInteraction = new OLSMultipleLinearRegression();
        cern.jet.random.tdouble.StudentT tDistColt = null;
        org.apache.commons.math3.distribution.FDistribution fDist = null;
        cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = null;

        Color[] colorarray = new Color[3];
        colorarray[0] = new Color(171, 178, 114);
        colorarray[1] = new Color(98, 175, 255);
        colorarray[2] = new Color(204, 86, 78);

        DecimalFormat decFormat = new DecimalFormat("#.###");
        DecimalFormat decFormatSmall = new DecimalFormat("0.#E0");
        for (Triple<String, String, String> triple : triples) {

            String snp = triple.getLeft();
            String covariate = triple.getMiddle();
            String probe = triple.getRight();

            Integer snpId = geno.getSnpToSNPId().get(snp);

            Integer probeId = expressionData.hashRows.get(probe);

            int startCovariate = -1;
            int endCovariate = -1;

            if (covariate == null) {
                startCovariate = 0;
                endCovariate = covariateData.nrRows;
            } else {
                Integer covariateId = covariateData.hashRows.get(covariate);
                if (covariateId != null) {
                    startCovariate = covariateId;
                    endCovariate = covariateId + 1;
                }
            }

            if (snpId >= 0 && probeId != null && startCovariate >= 0) {

                SNP snpObj = geno.getSNPObject(snpId);
                loader.loadGenotypes(snpObj);
                if (loader.hasDosageInformation()) {
                    loader.loadDosage(snpObj);
                }

                double signInteractionEffectDirection = 1;
                String[] genotypeDescriptions = new String[3];
                if (snpObj.getAlleles()[1] == snpObj.getMinorAllele()) {
                    signInteractionEffectDirection = -1;
                    genotypeDescriptions[2] = BaseAnnot.toString(snpObj.getAlleles()[0]) + "" + BaseAnnot.toString(snpObj.getAlleles()[0]);
                    genotypeDescriptions[1] = BaseAnnot.toString(snpObj.getAlleles()[0]) + "" + BaseAnnot.toString(snpObj.getAlleles()[1]);
                    genotypeDescriptions[0] = BaseAnnot.toString(snpObj.getAlleles()[1]) + "" + BaseAnnot.toString(snpObj.getAlleles()[1]);
                } else {
                    genotypeDescriptions[0] = BaseAnnot.toString(snpObj.getAlleles()[0]) + "" + BaseAnnot.toString(snpObj.getAlleles()[0]);
                    genotypeDescriptions[1] = BaseAnnot.toString(snpObj.getAlleles()[0]) + "" + BaseAnnot.toString(snpObj.getAlleles()[1]);
                    genotypeDescriptions[2] = BaseAnnot.toString(snpObj.getAlleles()[1]) + "" + BaseAnnot.toString(snpObj.getAlleles()[1]);

                }

                for (int q = startCovariate; q < endCovariate; q++) {
                    System.out.println("Plotting: " + snp + "\t" + covariateData.rowObjects.get(q) + "\t" + probe);
                    
                    
                    TextFile interactionOut = new TextFile(outdir + snp + "-" + probe + "-" + covariateData.rowObjects.get(q) + ".txt", TextFile.W);
                    interactionOut.writeln("Individual\tAllele1\tAllele2\tGenotype\tGenotypeFlipped\tCovariate\tExpression");
                    
                    byte[] alleles1 = snpObj.getAllele1();
                    byte[] alleles2 = snpObj.getAllele2();
                    byte[] genotypes = snpObj.getGenotypes();
                    ArrayList<Byte> genotypeArr = new ArrayList<Byte>();
                    ArrayList<Double> covariateArr = new ArrayList<Double>();
                    ArrayList<Double> expressionArr = new ArrayList<Double>();

                    int nrCalled = 0;

                    for (int i = 0; i < genoIndividuals.length; i++) {
                        if (genotypes[i] != -1 && genotypeToCovariate[i] != -9 && genotypeToExpression[i] != -9) {

                            if (!Double.isNaN(covariateData.rawData[q][genotypeToCovariate[i]])) {

                                int genotypeflipped = genotypes[i];
                                if (signInteractionEffectDirection == -1) {
                                    genotypeflipped = 2 - genotypeflipped;
                                }

                                String output = genoIndividuals[i]
                                        + "\t" + BaseAnnot.toString(alleles1[i])
                                        + "\t" + BaseAnnot.toString(alleles2[i])
                                        + "\t" + genotypes[i]
                                        + "\t" + genotypeflipped
                                        + "\t" + covariateData.rawData[q][genotypeToCovariate[i]]
                                        + "\t" + expressionData.rawData[probeId][genotypeToExpression[i]];
                                interactionOut.writeln(output);

                                genotypeArr.add(genotypes[i]);

                                covariateArr.add(covariateData.rawData[q][genotypeToCovariate[i]]);
                                expressionArr.add(expressionData.rawData[probeId][genotypeToExpression[i]]);
                                nrCalled++;
                            }

                        }
                        
                    }
                    interactionOut.close();
                    
                    System.out.println("");
                    //Fill arrays with data in order to be able to perform the ordinary least squares analysis:
                    double[] olsY = new double[nrCalled]; //Ordinary least squares: Our gene expression
                    double[][] olsXFullWithInteraction = new double[nrCalled][3];       //With interaction term, linear model: y ~ a * SNP + b * CellCount + c + d * SNP * CellCount
                    int itr = 0;

                    double[] dataExp = new double[nrCalled];
                    double[] dataCov = new double[nrCalled];
                    int[] dataGen = new int[nrCalled];

                    for (int s = 0; s < nrCalled; s++) {
                        byte originalGenotype = genotypeArr.get(s);
                        int genotype = originalGenotype;
                        if (signInteractionEffectDirection == -1) {
                            genotype = 2 - genotype;
                        }

                        olsY[s] = expressionArr.get(s);

                        olsXFullWithInteraction[s][0] = genotype;
                        olsXFullWithInteraction[s][1] = covariateArr.get(s);
                        olsXFullWithInteraction[s][2] = olsXFullWithInteraction[s][0] * olsXFullWithInteraction[s][1];

                        dataExp[s] = olsY[s];
                        dataGen[s] = genotype;
                        dataCov[s] = covariateArr.get(s);

                        itr++;
                    }

                    regressionFullWithInteraction.newSampleData(olsY, olsXFullWithInteraction);

                    try {
                        double rss2 = regressionFullWithInteraction.calculateResidualSumOfSquares();
                        double[] regressionParameters = regressionFullWithInteraction.estimateRegressionParameters();

                        double[] regressionStandardErrors = regressionFullWithInteraction.estimateRegressionParametersStandardErrors();

                        double betaInteraction = regressionParameters[3];
                        double seInteraction = regressionStandardErrors[3];
                        double tInteraction = betaInteraction / seInteraction;
                        double pValueInteraction = 1;
                        double zScoreInteraction = 0;

                        if (fDist == null) {
                            fDist = new org.apache.commons.math3.distribution.FDistribution((int) (3 - 2), (int) (olsY.length - 3));
                            randomEngine = new cern.jet.random.tdouble.engine.DRand();
                            tDistColt = new cern.jet.random.tdouble.StudentT(olsY.length - 4, randomEngine);
                        }

                        if (tInteraction < 0) {
                            pValueInteraction = tDistColt.cdf(tInteraction);
                            if (pValueInteraction < 2.0E-323) {
                                pValueInteraction = 2.0E-323;
                            }
                            zScoreInteraction = cern.jet.stat.tdouble.Probability.normalInverse(pValueInteraction);
                        } else {
                            pValueInteraction = tDistColt.cdf(-tInteraction);
                            if (pValueInteraction < 2.0E-323) {
                                pValueInteraction = 2.0E-323;
                            }

                            zScoreInteraction = -cern.jet.stat.tdouble.Probability.normalInverse(pValueInteraction);
                        }
                        pValueInteraction *= 2;
                        String pvalFormatted = "";
                        if (pValueInteraction >= 0.001) {
                            pvalFormatted = decFormat.format(pValueInteraction);
                        } else {
                            pvalFormatted = decFormatSmall.format(pValueInteraction);
                        }
                        ScatterPlot scatterPlot = new ScatterPlot(500, 500, dataCov, dataExp, dataGen, genotypeDescriptions, colorarray, ScatterPlot.OUTPUTFORMAT.PDF,
                                "Interaction between SNP " + snp + ", probe " + probe + " and covariate " + covariateData.rowObjects.get(q),
                                "Z: " + decFormat.format(zScoreInteraction) + " Pvalue: " + pvalFormatted + " n: " + nrCalled,
                                outdir + snp + "-" + probe + "-" + covariateData.rowObjects.get(q) + ".pdf", false);

                    } catch (SingularMatrixException ex) {
                        ex.printStackTrace();
                        System.out.println("\tMatrix is singular, skipping\n");
                    }
                }

                snpObj.clearGenotypes();
            }
        }

        loader.close();
    }
}
