/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.SNP;

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
    private final double[][] snpBeta;
    private final double[][] snpSE;
    private final double[][] covariateBeta;
    private final double[][] covariateSE;
    private final double[][] rsquared;

    InteractionAnalysisResults(String qcString,
            ArrayList<Pair<String, String>> eQTLsTested,
            double[][] interactionZScoreMatrix,
            double[][] SNPZResultMatrix,
            double[][] covariateZResultMatrix,
            double[][] maineffectZResultMatrix,
            int[][] nMatrix,
            double[][] rsquaredMatrix) {
        this.qcString = qcString;
        this.eQTLsTested = eQTLsTested;
        this.interactionZScoreMatrix = interactionZScoreMatrix;
        this.SNPZResultMatrix = SNPZResultMatrix;
        this.covariateZResultMatrix = covariateZResultMatrix;
        this.maineffectZResultMatrix = maineffectZResultMatrix;
        this.nMatrix = nMatrix;
        this.rsquared = rsquaredMatrix;
        interactionBeta = null;
        interactionSE = null;
        snpBeta = null;
        snpSE = null;
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
            int[][] nMatrix,
            double[][] rsquaredMatrix) {
        this.qcString = qcString;
        this.eQTLsTested = eQTLsTested;
        this.interactionZScoreMatrix = interactionZScoreMatrix;
        this.SNPZResultMatrix = SNPZResultMatrix;
        this.covariateZResultMatrix = covariateZResultMatrix;
        this.maineffectZResultMatrix = maineffectZResultMatrix;
        this.nMatrix = nMatrix;

        this.interactionBeta = interactionBeta;
        this.interactionSE = interactionSE;
        this.snpBeta = mainBeta;
        this.snpSE = mainSE;
        this.covariateBeta = covariateBeta;
        this.covariateSE = covariateSE;
        this.rsquared = rsquaredMatrix;
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

    public double[][] getSNPBeta() {
        return snpBeta;
    }

    public double[][] getSNPSE() {
        return snpSE;
    }

    public double[][] getRsquared() {
        return rsquared;
    }

    public double[][] getCovariateBeta() {
        return covariateBeta;
    }

    public double[][] getCovariateSE() {
        return covariateSE;
    }

	public ArrayList<String> getProbeIds() {
		ArrayList<String> probeIds = new ArrayList<String>();

		for (Pair<String, String> eqtl : eQTLsTested){
			String gene = eqtl.getRight();
			if (! probeIds.contains(gene))
				probeIds.add(gene);
		}
		return probeIds;
	}
	public ArrayList<String> getSNPIds() {
		ArrayList<String> snpIds = new ArrayList<String>();

		for (Pair<String, String> eqtl : eQTLsTested){
			String snp = eqtl.getLeft();
			if (! snpIds.contains(snp))
				snpIds.add(snp);
		}
		return snpIds;
	}
}
