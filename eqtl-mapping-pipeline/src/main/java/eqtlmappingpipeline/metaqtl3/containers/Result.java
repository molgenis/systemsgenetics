/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.containers;

/**
 *
 * @author harm-jan
 */
public class Result {
    public double[][]  correlations;
    public double[][]  zscores;
    public int[] numSamples;
    public int[] nrAllelesA;
    public double[] pvalues;
//    public WorkPackage workpackage;
    public double[] finalZScore;
    public boolean poison = false;
   
    public byte alleles;
    public byte[] assessedAllele;
    
    public boolean processed = false;
    public int wpid;
    public double[][] beta;
    public double[][] se;
    public double[][] fc;
    public double[] finalBetaSe;
    public double[] finalBeta;
    
    
    
    public Result(boolean poison){
        this.poison = poison;
    }
    
    public Result(int numDs, int numProbes, int workpackageid) {
        correlations = new double[numDs][numProbes];
        zscores      = new double[numDs][numProbes];
        numSamples   = new int[numDs];
        pvalues      = new double[numProbes];
        finalZScore  = new double[numProbes];
        finalBeta    = new double[numProbes];
        finalBetaSe    = new double[numProbes];
        beta         = new double[numDs][numProbes];
        se           = new double[numDs][numProbes];
        fc           = new double[numDs][numProbes];
        wpid         = workpackageid;
    }

//    public void clearValues(int p) {
//        pvalues[p]     = Double.NaN;
//        finalZScore[p] = null;
//        for(int d=0; d<correlations.length; d++){
//            correlations[d][p] = null;
//            zscores[d][p] = null;
//        }
//    }

   
    
}
