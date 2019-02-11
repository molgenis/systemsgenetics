/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.systemsgenetics.depict2.originalLude;

import java.util.concurrent.Callable;
/**
 *
 * @author lude
 */
public class EstimateChi2SumDistUsingCorrelatedVariablesThread implements Callable<double[]> {

    public double[] eigenValues = null;
    public int nrPerms = 0;
    public int threadID = 0;
    
    public EstimateChi2SumDistUsingCorrelatedVariablesThread(double[] eigenValues, int nrPerms, int threadID) {
        this.eigenValues = eigenValues;
        this.nrPerms = nrPerms;
        this.threadID = threadID;
    }
    
    @Override
    public double[] call() throws Exception {
        java.util.Random rnd = new java.util.Random();
        double[] geneChi2SumNull = new double[nrPerms];
        for (int perm=0; perm<nrPerms; perm++) {
            double weightedChi2Perm = 0;
            for (int g=0; g<eigenValues.length; g++) {
                double randomZ = rnd.nextGaussian();
                weightedChi2Perm+=eigenValues[g] * randomZ * randomZ;
            }
            geneChi2SumNull[perm]=weightedChi2Perm;
        } 
        return geneChi2SumNull;
    }    
 
   
}
