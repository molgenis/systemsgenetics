/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

/**
 *
 * @author harmjan
 */
public class Result {

    int snp;
    Double[] finalzscores;
    Double[] finalpvalues;
    Integer[][] numSamples; // [ds][probe]
    boolean[][] dspassingqc;
    boolean passesQC;
    Double[][] datasetZScores;
    String[] datasets;

    public void clearData() {
	for (int i = 0; i < finalzscores.length; i++) {
	    finalzscores[i] = null;
	}
	finalzscores = null;

	for (int i = 0; i < finalpvalues.length; i++) {
	    finalpvalues[i] = null;
	}
	finalpvalues = null;

	for (int i = 0; i < numSamples.length; i++) {
	    for (int j = 0; j < numSamples[i].length; j++) {
		numSamples[i][j] = null;
	    }
	    numSamples[i] = null;
	}
	numSamples = null;
	
//	for (int i = 0; i < datasetZScores.length; i++) {
//	    for (int j = 0; j < datasetZScores[i].length; j++) {
//		datasetZScores[i][j] = null;
//	    }
//	    datasetZScores[i] = null;
//	}
//	datasetZScores = null;

	for(int i=0;i<datasets.length;i++){
	    datasets[i] = null;
	}
	datasets = null;
	
	for(int i=0;i<dspassingqc.length;i++){
	    dspassingqc[i] = null;
	}
	dspassingqc = null;
    }
}
