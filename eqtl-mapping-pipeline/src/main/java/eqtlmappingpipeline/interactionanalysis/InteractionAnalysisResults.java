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
    
    /*
    //            rowNames.add("CellTypeSNPZScore");
//            rowNames.add("CellTypeZScore");
//            rowNames.add("CellTypeInteractionZScore");
//            rowNames.add("MainEffectZScore");
    */

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
  
}
