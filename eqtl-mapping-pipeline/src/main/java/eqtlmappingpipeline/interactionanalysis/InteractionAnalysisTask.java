/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;

/**
 *
 * @author harmjan
 */
public class InteractionAnalysisTask implements Callable<InteractionAnalysisResults> {

    private SNP eQTLSNPObj;
    private double[][] pcCorrectedData;
//    private double[] cellcounts;
//    private String[] covariatesToUse;
    private int[] wgaId;
    private String[] expInds;
    private DoubleMatrixDataset<String, String> covariateData;
    private TriTyperExpressionData expressionDataPCCorrected;
    private ArrayList<Pair<String, String>> eQTLsForSNP;

    public InteractionAnalysisTask(SNP snpObj, ArrayList<Pair<String, String>> eQTLsForSNP, double[][] pcCorrectedData,
            int[] wgaId,
            String[] expInds, DoubleMatrixDataset<String, String> expressionDataRaw, TriTyperExpressionData expressionDataPCCorrected) {
        this.eQTLSNPObj = snpObj;
        this.eQTLsForSNP = eQTLsForSNP;
        this.pcCorrectedData = pcCorrectedData;
//        this.cellcounts = cellcounts;
//        this.covariatesToUse = probesToUseAsCovariateArr;
        this.wgaId = wgaId;
        this.expInds = expInds;
        this.expressionDataPCCorrected = expressionDataPCCorrected;
        this.covariateData = expressionDataRaw;
    }

    @Override
    public InteractionAnalysisResults call() throws Exception {

        ArrayList<Pair<String, String>> eQTLsTested = new ArrayList<Pair<String, String>>();

        int nrTotalCovariates = covariateData.nrRows;
//        if (cellcounts != null) {
//            nrRows += 4;
//        }

        double[][] interactionZScoreMatrix = new double[eQTLsForSNP.size()][nrTotalCovariates];
        double[][] SNPZResultMatrix = new double[eQTLsForSNP.size()][nrTotalCovariates];
        double[][] covariateZResultMatrix = new double[eQTLsForSNP.size()][nrTotalCovariates];
        double[][] maineffectZResultMatrix = new double[eQTLsForSNP.size()][nrTotalCovariates];
        int[][] nMatrix = new int[eQTLsForSNP.size()][nrTotalCovariates];

        //We are using a coding system that uses the minor allele. If allele2 is not the minor allele, change the sign of the results we will output.
        double signInteractionEffectDirection = 1;
        if (eQTLSNPObj.getAlleles()[1] == eQTLSNPObj.getMinorAllele()) {
            signInteractionEffectDirection = -1;
        }

        String qcString = null;
        Integer nrGenotypesCalled = null;

        org.apache.commons.math3.distribution.FDistribution fDist = null;

        cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = null;
        cern.jet.random.tdouble.StudentT tDistColt = null;

        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        OLSMultipleLinearRegression regressionFullWithInteraction = new OLSMultipleLinearRegression();

        int startWithCovariate = 0;
//        if(cellcounts != null){
//            startWithCovariate = -1;
//        }

        for (int e = 0; e < eQTLsForSNP.size(); e++) {
            Pair<String, String> eqtl = eQTLsForSNP.get(e);
            String eQTLProbeName = eqtl.getRight();

            eQTLsTested.add(eqtl);

            Integer eQTLProbeId = expressionDataPCCorrected.getProbeToId().get(eQTLProbeName);

            double[] valsX = eQTLSNPObj.selectGenotypes(wgaId, true, true);
            double[] valsY = pcCorrectedData[eQTLProbeId]; //Expression level

            for (int covariate = startWithCovariate; covariate < nrTotalCovariates; covariate++) {
                double[] covariateValues;
//              

                double[] tmpVarCelCount = null;
//                    if (covariateIdInRawData != null) {
                tmpVarCelCount = new double[valsY.length];
                for (int i = 0; i < tmpVarCelCount.length; i++) {
                    String sampleName = expInds[i];
                    Integer individualIdInCovariateData = covariateData.hashCols.get(sampleName);
                    if (individualIdInCovariateData != null) {
                        // presorting greatly speeds this stuff up
                        tmpVarCelCount[i] = covariateData.rawData[covariate][individualIdInCovariateData];
                    } else {
                        tmpVarCelCount[i] = Double.NaN;
                    }
                }
//                    } else {
//                        System.err.println("Covariate: " + covariateName + " not present in RAW data!");
//                    }
                covariateValues = tmpVarCelCount;
//                }

//                if (valsCellCount != null) {
                //Check whether all the expression samples have a genotype and a cell count...
                int nrCalled = 0;
                for (int i = 0; i < wgaId.length; i++) {
                    if (wgaId[i] != -1 && !Double.isNaN(covariateValues[i]) && valsX[i] != -1) {
                        nrCalled++;
                    }
                }

                //Fill arrays with data in order to be able to perform the ordinary least squares analysis:
                double[] olsY = new double[nrCalled]; //Ordinary least squares: Our gene expression
                double[] genotypesCalled = new double[nrCalled];
                double[][] olsX = new double[nrCalled][2];                          //No interaction term, linear model: y ~ a * SNP + b * CellCount + c
                double[][] olsXFullWithInteraction = new double[nrCalled][3];       //With interaction term, linear model: y ~ a * SNP + b * CellCount + c + d * SNP * CellCount
                int itr = 0;
                for (int s = 0; s < valsX.length; s++) {
                    double genotype = valsX[s];
                    if (genotype != -1 && !Double.isNaN(covariateValues[s])) {
                        if (signInteractionEffectDirection == -1) {
                            genotype = 2 - genotype;
                        }
                        genotypesCalled[itr] = genotype;
                        olsY[itr] = valsY[s];
                        olsX[itr][0] = genotype;
                        olsXFullWithInteraction[itr][0] = genotype;
                        olsX[itr][1] = covariateValues[s];
                        olsXFullWithInteraction[itr][1] = covariateValues[s];
                        olsXFullWithInteraction[itr][2] = olsXFullWithInteraction[itr][0] * olsXFullWithInteraction[itr][1];
                        itr++;
                    }
                }

                regression.newSampleData(olsY, olsX);
                regressionFullWithInteraction.newSampleData(olsY, olsXFullWithInteraction);
                double rss1 = regression.calculateResidualSumOfSquares();
                double rss2 = regressionFullWithInteraction.calculateResidualSumOfSquares();
                double anovaF = ((rss1 - rss2) / (3 - 2)) / (rss2 / (olsY.length - 3));
                // Changed this to apache maths 3, was apache maths 1.0
                if (fDist == null) {
                    fDist = new org.apache.commons.math3.distribution.FDistribution((int) (3 - 2), (int) (olsY.length - 3));
                    randomEngine = new cern.jet.random.tdouble.engine.DRand();
                    tDistColt = new cern.jet.random.tdouble.StudentT(olsY.length - 4, randomEngine);
                }

                double anovaFTestP = -1;
                try {
                    anovaFTestP = 1 - fDist.cumulativeProbability(anovaF);
                    if (anovaFTestP < 1E-16) {
                        anovaFTestP = 1E-16;
                    }
                } catch (Exception err) {
                }

                // Get the regression parameters and R-square value and print it.
                double[] regressionParameters = regressionFullWithInteraction.estimateRegressionParameters();
                double[] regressionStandardErrors = regressionFullWithInteraction.estimateRegressionParametersStandardErrors();

                double betaInteraction = regressionParameters[3];
                double seInteraction = regressionStandardErrors[3];
                double tInteraction = betaInteraction / seInteraction;
                double pValueInteraction = 1;
                double zScoreInteraction = 0;
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

                // THIS WILL GIVE ERRONEOUS VALUES WHEN THERE ARE MISSING 
                // VALUES IN VALSY THE NEXT TIME THIS SNP IS TESTED!!
                // this value is required for subsequent meta-analysis.. fix for altering sample sizes (take smallest size / omit missing values)
                if (qcString == null) {
                    qcString = eQTLSNPObj.getName() + "\t" + ChrAnnotation.parseByte(eQTLSNPObj.getChr()) + "\t" + eQTLSNPObj.getChrPos() + "\t" + BaseAnnot.toString(eQTLSNPObj.getAlleles()[0]) + "/" + BaseAnnot.toString(eQTLSNPObj.getAlleles()[1]) + "\t" + BaseAnnot.toString(eQTLSNPObj.getMinorAllele()) + "\t" + eQTLSNPObj.getMAF() + "\t" + eQTLSNPObj.getCR() + "\t" + eQTLSNPObj.getHWEP() + "\t" + genotypesCalled.length;
                    nrGenotypesCalled = genotypesCalled.length;
                } else if (genotypesCalled.length != nrGenotypesCalled) {

                    System.err.println("ERROR: the number of available values has changed. Does your gene expression data or cell count file contain missing values?");
                    System.exit(0);
                }

                double corr = JSci.maths.ArrayMath.correlation(genotypesCalled, olsY);
                double mainZ = Correlation.convertCorrelationToZScore(genotypesCalled.length, corr);

                // Get the regression parameters and R-square value and print it.
                double betaSNP = regressionParameters[1];
                double seSNP = regressionStandardErrors[1];
                double tSNP = betaSNP / seSNP;
                double pValueSNP = 1;
                double zScoreSNP = 0;
                if (tSNP < 0) {
                    pValueSNP = tDistColt.cdf(tSNP);
                    if (pValueSNP < 2.0E-323) {
                        pValueSNP = 2.0E-323;
                    }
                    zScoreSNP = cern.jet.stat.tdouble.Probability.normalInverse(pValueSNP);
                } else {
                    pValueSNP = tDistColt.cdf(-tSNP);
                    if (pValueSNP < 2.0E-323) {
                        pValueSNP = 2.0E-323;
                    }
                    zScoreSNP = -cern.jet.stat.tdouble.Probability.normalInverse(pValueSNP);
                }
                pValueSNP *= 2;

                // Get the regression parameters and R-square value and print it.
                double betaCellType = regressionParameters[2];
                double seCellType = regressionStandardErrors[2];
                double tCellType = betaCellType / seCellType;
                double pValueCellType = 1;
                double zScoreCovariate = 0;
                if (tCellType < 0) {
                    pValueCellType = tDistColt.cdf(tCellType);
                    if (pValueCellType < 2.0E-323) {
                        pValueCellType = 2.0E-323;
                    }
                    zScoreCovariate = cern.jet.stat.tdouble.Probability.normalInverse(pValueCellType);
                } else {
                    pValueCellType = tDistColt.cdf(-tCellType);
                    if (pValueCellType < 2.0E-323) {
                        pValueCellType = 2.0E-323;
                    }
                    zScoreCovariate = -cern.jet.stat.tdouble.Probability.normalInverse(pValueCellType);
                }
                pValueCellType *= 2;

                interactionZScoreMatrix[e][covariate] = zScoreInteraction;
                SNPZResultMatrix[e][covariate] = zScoreSNP;
                covariateZResultMatrix[e][covariate] = zScoreCovariate;
                maineffectZResultMatrix[e][covariate] = mainZ;
                nMatrix[e][covariate] = nrCalled;

//                      cellcountInterActionOutput = eQTLSNPObj.getName() + "\t" + eQTLProbeName + "\t" + nrCalled + "\t" + corr + "\t" + anovaFTestP + "\t" + betaInteraction + "\t" + seInteraction + "\t" + tInteraction + "\t" + pValueInteraction + "\t" + zScoreInteraction;
//                }
            }
        }

        // get genotypes, include missing ones
        InteractionAnalysisResults result = new InteractionAnalysisResults(
                
                qcString,
                eQTLsTested,
                interactionZScoreMatrix,
                SNPZResultMatrix,
                covariateZResultMatrix,
                maineffectZResultMatrix,
                nMatrix);

        eQTLSNPObj.clearGenotypes();
        eQTLSNPObj = null;

        return result;
    }
}
