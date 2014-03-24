/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import java.util.ArrayList;

/**
 *
 * @author harmjan
 */
public class InteractionAnalysisResults {
    private final String qcString;

    
    private final String cellcountInterActionOutput;
    private final ArrayList<String> eQTLsTested;
    private final double[][] interactionVector;

    InteractionAnalysisResults(String qcString, String cellcountInterActionOutput, ArrayList<String> eQTLsTested, double[][] interactionVector) {
        this.qcString = qcString;
        this.cellcountInterActionOutput = cellcountInterActionOutput;
        this.eQTLsTested = eQTLsTested;
        this.interactionVector = interactionVector;
    }
    
    public String getQcString() {
        return qcString;
    }

    public String getCellcountInterActionOutput() {
        return cellcountInterActionOutput;
    }

    public ArrayList<String> geteQTLsTested() {
        return eQTLsTested;
    }

    public double[][] getInteractionVector() {
        return interactionVector;
    }
    
}
