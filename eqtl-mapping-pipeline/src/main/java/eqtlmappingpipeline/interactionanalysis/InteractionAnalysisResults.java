/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import java.util.ArrayList;
import umcg.genetica.containers.Pair;

/**
 *
 * @author harmjan
 */
public class InteractionAnalysisResults {

    private final String qcString;

    private final ArrayList<Pair<String, String>> eQTLsTested;
    private final double[][] interactionZScoreMatrix;
    private final double[][] SNPZResultMatrix;
    private final double[][] covariateZResultMatrix;
    private final double[][] maineffectZResultMatrix;
    private final int[][] nMatrix;

    private final double[][] interactionBeta;
    private final double[][] interactionSE;
    private final double[][] mainBeta;
    private final double[][] mainSE;
    private final double[][] covariateBeta;
    private final double[][] covariateSE;

    

    InteractionAnalysisResults(String qcString,
            ArrayList<Pair<String, String>> eQTLsTested,
            double[][] interactionZScoreMatrix,
            double[][] SNPZResultMatrix,
            double[][] covariateZResultMatrix,
            double[][] maineffectZResultMatrix,
            int[][] nMatrix) {
        this.qcString = qcString;
        this.eQTLsTested = eQTLsTested;
        this.interactionZScoreMatrix = interactionZScoreMatrix;
        this.SNPZResultMatrix = SNPZResultMatrix;
        this.covariateZResultMatrix = covariateZResultMatrix;
        this.maineffectZResultMatrix = maineffectZResultMatrix;
        this.nMatrix = nMatrix;
        interactionBeta = null;
        interactionSE = null;
        mainBeta = null;
        mainSE = null;
        covariateBeta = null;
        covariateSE = null;
    }

    InteractionAnalysisResults(String qcString,
            ArrayList<Pair<String, String>> eQTLsTested,
            double[][] interactionZScoreMatrix,
            double[][] SNPZResultMatrix,
            double[][] covariateZResultMatrix,
            double[][] maineffectZResultMatrix,
            double[][] interactionBeta,
            double[][] interactionSE,
            double[][] mainBeta,
            double[][] mainSE,
            double[][] covariateBeta,
            double[][] covariateSE,
            int[][] nMatrix) {
        this.qcString = qcString;
        this.eQTLsTested = eQTLsTested;
        this.interactionZScoreMatrix = interactionZScoreMatrix;
        this.SNPZResultMatrix = SNPZResultMatrix;
        this.covariateZResultMatrix = covariateZResultMatrix;
        this.maineffectZResultMatrix = maineffectZResultMatrix;
        this.nMatrix = nMatrix;

        this.interactionBeta = interactionBeta;
        this.interactionSE = interactionSE;
        this.mainBeta = mainBeta;
        this.mainSE = mainSE;
        this.covariateBeta = covariateBeta;
        this.covariateSE = covariateSE;
    }

    public String getQcString() {
        return qcString;
    }

    public ArrayList<Pair<String, String>> geteQTLsTested() {
        return eQTLsTested;
    }

    public double[][] getInteractionZScoreMatrix() {
        return interactionZScoreMatrix;
    }

    public double[][] getSNPZResultMatrix() {
        return SNPZResultMatrix;
    }

    public double[][] getCovariateZResultMatrix() {
        return covariateZResultMatrix;
    }

    public double[][] getMaineffectZResultMatrix() {
        return maineffectZResultMatrix;
    }

    public int[][] getnMatrix() {
        return nMatrix;
    }
    
    public double[][] getInteractionBeta() {
        return interactionBeta;
    }

    public double[][] getInteractionSE() {
        return interactionSE;
    }

    public double[][] getMainBeta() {
        return mainBeta;
    }

    public double[][] getMainSE() {
        return mainSE;
    }

    public double[][] getCovariateBeta() {
        return covariateBeta;
    }

    public double[][] getCovariateSE() {
        return covariateSE;
    }

}
