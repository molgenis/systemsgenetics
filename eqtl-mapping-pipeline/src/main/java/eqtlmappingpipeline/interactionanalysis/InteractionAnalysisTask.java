/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import cern.jet.random.tdouble.StudentT;
import java.util.ArrayList;
import java.util.concurrent.Callable;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
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
    private double[][] pcCorrectedExpressionData;
//    private double[] cellcounts;
//    private String[] covariatesToUse;
    private int[] wgaId;
    private String[] expInds;
    private DoubleMatrixDataset<String, String> covariateData;
    private TriTyperExpressionData expressionData;
    private ArrayList<Pair<String, String>> eQTLsForSNP;
    private RConnection rConnection; // this requires the Rserve package in R..
    private final boolean sandwich;
    private final boolean provideFullStats;

    public InteractionAnalysisTask(SNP snpObj, ArrayList<Pair<String, String>> eQTLsForSNP, double[][] pcCorrectedData,
            int[] wgaId,
            String[] expInds, DoubleMatrixDataset<String, String> covariateData,
            TriTyperExpressionData expressionData, boolean robustSE, boolean provideFullStats) {
        this.eQTLSNPObj = snpObj;
        this.eQTLsForSNP = eQTLsForSNP;
        this.pcCorrectedExpressionData = pcCorrectedData;
//        this.cellcounts = cellcounts;
//        this.covariatesToUse = probesToUseAsCovariateArr;
        this.wgaId = wgaId;
        this.expInds = expInds;
        this.expressionData = expressionData;
        this.covariateData = covariateData;

        if (robustSE) {
            // this code is very suboptimal and is here for validation purposes only
            try {
                rConnection = new RConnection();
//                rConnection.voidEval("install.packages('sandwich')");
                rConnection.voidEval("library(sandwich)");
            } catch (RserveException ex) {
                System.err.println(ex.getMessage());
                rConnection = null;
            }
        }
        this.sandwich = robustSE;
        this.provideFullStats = provideFullStats;

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
        double[][] interactionBeta = null;
        double[][] interactionSE = null;
        double[][] mainBeta = null;
        double[][] mainSE = null;
        double[][] covariateBeta = null;
        double[][] covariateSE = null;

        if (provideFullStats) {
            interactionBeta = new double[eQTLsForSNP.size()][nrTotalCovariates];
            interactionSE = new double[eQTLsForSNP.size()][nrTotalCovariates];
            mainBeta = new double[eQTLsForSNP.size()][nrTotalCovariates];
            mainSE = new double[eQTLsForSNP.size()][nrTotalCovariates];
            covariateBeta = new double[eQTLsForSNP.size()][nrTotalCovariates];
            covariateSE = new double[eQTLsForSNP.size()][nrTotalCovariates];
        }
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

//        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        OLSMultipleLinearRegression regressionFullWithInteraction = new OLSMultipleLinearRegression();

        int startWithCovariate = 0;
//        if(cellcounts != null){
//            startWithCovariate = -1;
//        }

        for (int e = 0; e < eQTLsForSNP.size(); e++) {
            Pair<String, String> eqtl = eQTLsForSNP.get(e);
            String eQTLProbeName = eqtl.getRight();

            eQTLsTested.add(eqtl);

            Integer eQTLProbeId = expressionData.getProbeToId().get(eQTLProbeName);

            double[] valsX = eQTLSNPObj.selectGenotypes(wgaId, true, true);
            double[] valsY = pcCorrectedExpressionData[eQTLProbeId]; //Expression level

            for (int covariate = startWithCovariate; covariate < nrTotalCovariates; covariate++) {
                double[] covariateValues;
//              

                double[] tmpVarCelCount = new double[valsY.length];
//                    if (covariateIdInRawData != null) {

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

                // THIS WILL GIVE ERRONEOUS VALUES WHEN THERE ARE MISSING 
                // VALUES IN VALSY THE NEXT TIME THIS SNP IS TESTED!!
                // this value is required for subsequent meta-analysis.. fix for altering sample sizes (take smallest size / omit missing values)
                // in stead use the value for N that is now in the standard output.
                double[] genotypesCalled = new double[nrCalled];
                if (qcString == null) {
                    qcString = eQTLSNPObj.getName() + "\t" + ChrAnnotation.parseByte(eQTLSNPObj.getChr()) + "\t" + eQTLSNPObj.getChrPos() + "\t" + BaseAnnot.toString(eQTLSNPObj.getAlleles()[0]) + "/" + BaseAnnot.toString(eQTLSNPObj.getAlleles()[1]) + "\t" + BaseAnnot.toString(eQTLSNPObj.getMinorAllele()) + "\t" + eQTLSNPObj.getMAF() + "\t" + eQTLSNPObj.getCR() + "\t" + eQTLSNPObj.getHWEP() + "\t" + genotypesCalled.length;
                    nrGenotypesCalled = genotypesCalled.length;
                } else if (genotypesCalled.length != nrGenotypesCalled) {

                    System.err.println("ERROR: the number of available values has changed. Does your gene expression data or cell count file contain missing values?");
                    System.exit(0);
                }

                double zScoreInteraction = 0;
                double zScoreSNP = 0;
                double zScoreCovariate = 0;
                double mainZ = 0;

                double betaInteraction = 0;
                double seInteraction = 0;
                double betaSNP = 0;
                double seSNP = 0;
                double betaCovariate = 0;
                double seCovariate = 0;

                if (sandwich) {
                    if (rConnection == null) {
                        System.err.println("Error: using R connection but none found");
                    }
                    if (rConnection != null) {
                        try {
                            if (rConnection.isConnected()) {
                                double[] olsY = new double[nrCalled]; //Ordinary least squares: Our gene expression

                                double[] olsX = new double[nrCalled];                          //No interaction term, linear model: y ~ a * SNP + b * CellCount + c
//                                double[][] olsXFullWithInteraction = new double[nrCalled][3];       //With interaction term, linear model: y ~ a * SNP + b * CellCount + c + d * SNP * CellCount
                                int itr = 0;
                                for (int s = 0; s < valsX.length; s++) {
                                    double genotype = valsX[s];
                                    if (genotype != -1 && !Double.isNaN(covariateValues[s])) {
                                        if (signInteractionEffectDirection == -1) {
                                            genotype = 2 - genotype;
                                        }

                                        olsY[itr] = valsY[s];
                                        olsX[itr] = genotype;
                                        itr++;
                                    }
                                }

                                rConnection.assign("y", olsY);
                                rConnection.assign("x", olsX);
                                rConnection.assign("z", covariateValues);
                                rConnection.voidEval("interaction <- x*z");
                                rConnection.voidEval("m <- lm(y ~ x + z + interaction)");
                                rConnection.voidEval("modelsummary <- summary(m)");

                                rConnection.voidEval("m2 <- sqrt(diag(vcovHC(m, type = 'HC0')))"); // robust covariance model

                                if (tDistColt == null) {
                                    randomEngine = new cern.jet.random.tdouble.engine.DRand();
                                    tDistColt = new cern.jet.random.tdouble.StudentT(olsY.length - 4, randomEngine);
                                }

                                betaInteraction = rConnection.eval("modelsummary$coefficients[4,1]").asDouble();
                                seInteraction = rConnection.eval("as.numeric(m2[4])").asDouble();
                                betaSNP = rConnection.eval("modelsummary$coefficients[2,1]").asDouble();
                                seSNP = rConnection.eval("modelsummary$coefficients[2,2]").asDouble();
                                betaCovariate = rConnection.eval("modelsummary$coefficients[3,1]").asDouble();
                                seCovariate = rConnection.eval("modelsummary$coefficients[3,2]").asDouble();
                            } else {
                                System.err.println("ERROR: R is not connected.");
                            }

                        } catch (REngineException ex) {
                            System.err.println(ex.getMessage());
                        } catch (REXPMismatchException ex) {
                            System.err.println(ex.getMessage());
                        }

                    }

                } else {

                    //Fill arrays with data in order to be able to perform the ordinary least squares analysis:
                    double[] olsY = new double[nrCalled]; //Ordinary least squares: Our gene expression

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

//                    regression.newSampleData(olsY, olsX);
                    regressionFullWithInteraction.newSampleData(olsY, olsXFullWithInteraction);

                    // not sure if this is needed right now, but I will keep it in for later use.
//                    double rss1 = regression.calculateResidualSumOfSquares();
//                    double rss2 = regressionFullWithInteraction.calculateResidualSumOfSquares();
//                    double anovaF = ((rss1 - rss2) / (3 - 2)) / (rss2 / (olsY.length - 3));
//                    // Changed this to apache maths 3, was apache maths 1.0
//                    if (fDist == null) {
//                        fDist = new org.apache.commons.math3.distribution.FDistribution((int) (3 - 2), (int) (olsY.length - 3));
//                  randomEngine = new cern.jet.random.tdouble.engine.DRand();
//                    tDistColt = new cern.jet.random.tdouble.StudentT(olsY.length - 4, randomEngine);
//                }
//
//                    double anovaFTestP = -1;
//                    try {
//                        anovaFTestP = 1 - fDist.cumulativeProbability(anovaF);
//                        if (anovaFTestP < 1E-16) {
//                            anovaFTestP = 1E-16;
//                        }
//                    } catch (Exception err) {
//                    }
                    if (tDistColt == null) {
                        randomEngine = new cern.jet.random.tdouble.engine.DRand();
                        tDistColt = new cern.jet.random.tdouble.StudentT(olsY.length - 4, randomEngine);
                    }

                    // double intersect = regressionParameters[0]; 
                    double corr = JSci.maths.ArrayMath.correlation(genotypesCalled, olsY);
                    mainZ = Correlation.convertCorrelationToZScore(genotypesCalled.length, corr);
                    // Get the regression parameters and R-square value and print it.
                    double[] regressionParameters = regressionFullWithInteraction.estimateRegressionParameters();
                    double[] regressionStandardErrors = regressionFullWithInteraction.estimateRegressionParametersStandardErrors();

                    betaInteraction = regressionParameters[3];
                    seInteraction = regressionStandardErrors[3];

                    // Get the regression parameters and R-square value and print it.
                    betaSNP = regressionParameters[1];
                    seSNP = regressionStandardErrors[1];

                    betaCovariate = regressionParameters[2];
                    seCovariate = regressionStandardErrors[2];

                }

                double pValueInteraction = convertBetaToP(betaInteraction, seInteraction, tDistColt);
                zScoreInteraction = convertPToZ(pValueInteraction);

                double pValueSNP = convertBetaToP(betaSNP, seSNP, tDistColt);
                zScoreSNP = convertPToZ(pValueSNP);

                // Get the regression parameters and R-square value and print it.
                double pValueCovariate = convertBetaToP(betaCovariate, seCovariate, tDistColt);
                zScoreCovariate = convertPToZ(pValueCovariate);

                interactionZScoreMatrix[e][covariate] = zScoreInteraction;
                SNPZResultMatrix[e][covariate] = zScoreSNP;
                covariateZResultMatrix[e][covariate] = zScoreCovariate;
                maineffectZResultMatrix[e][covariate] = mainZ;
                nMatrix[e][covariate] = nrCalled;

                if (provideFullStats) {
                    interactionBeta[e][covariate] = betaInteraction;
                    interactionSE[e][covariate] = seInteraction;
                    mainBeta[e][covariate] = betaSNP;
                    mainSE[e][covariate] = seSNP;
                    covariateBeta[e][covariate] = betaCovariate;
                    covariateSE[e][covariate] = seCovariate;
                }
            }
        }

        eQTLSNPObj.clearGenotypes();
        eQTLSNPObj = null;
        // get genotypes, include missing ones

        if (provideFullStats) {

            return new InteractionAnalysisResults(
                    qcString,
                    eQTLsTested,
                    interactionZScoreMatrix,
                    SNPZResultMatrix,
                    covariateZResultMatrix,
                    maineffectZResultMatrix,
                    interactionBeta,
                    interactionSE,
                    mainBeta,
                    mainSE,
                    covariateBeta,
                    covariateSE,
                    nMatrix);
        } else {
            return new InteractionAnalysisResults(
                    qcString,
                    eQTLsTested,
                    interactionZScoreMatrix,
                    SNPZResultMatrix,
                    covariateZResultMatrix,
                    maineffectZResultMatrix,
                    nMatrix);

        }

    }

    private double convertBetaToP(double beta, double se, StudentT tDistColt) {
        double t = beta / se;
        double p = 1;
        if (t < 0) {
            p = tDistColt.cdf(t);

        } else {
            p = tDistColt.cdf(-t);
        }
        if (p < 2.0E-323) {
            p = 2.0E-323;
        }
        return p;
    }

    private double convertPToZ(double p) {
        return cern.jet.stat.tdouble.Probability.normalInverse(p);
    }
}
