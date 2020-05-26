/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.containers;

import umcg.genetica.io.trityper.SNP;

/**
 *
 * @author harmjan
 */
public class WorkPackage implements Comparable<WorkPackage> {

    private SNP[] snps;         // 1 snp / dataset
    private int[] probes; // multi probes / dataset or null, for all probes
//    private Boolean[] passesQC;
    private Boolean[] flipSNPAlleles;
    private boolean poison;
    private short datasetsPassingQC;
    private int sortSNPsByDataset = 0;
    public Result results;
    private int id;
    private int numtested = 0;
    private int metaSNPId = -1;
    private boolean hasResults;

    /**
     * @return the snps
     */
    public SNP[] getSnps() {
        return snps;
    }

    public void setNumTested(int num) {
        numtested = num;
    }

    public int getNumTested() {
        return numtested;
    }

    /**
     * @param snps the snps to set
     */
    public void setSnps(SNP[] snps) {
        this.snps = snps;
    }

    /**
     * @return the probes
     */
    public int[] getProbes() {
        return probes;
    }

    /**
     * @param probes the probes to set
     */
    public void setProbes(int[] probes) {
        this.probes = probes;
    }

    public void setIsKillPackage(boolean b) {
        this.poison = b;
    }

    public boolean getPoison() {
        return poison;
    }

    /**
     * @return the datasetsPassingQC
     */
    public short getDatasetsPassingQC() {
        return datasetsPassingQC;
    }

    /**
     * @param datasetsPassingQC the datasetsPassingQC to set
     */
    public void setDatasetsPassingQC(short datasetsPassingQC) {
        this.datasetsPassingQC = datasetsPassingQC;
    }

    public void setDatasetToSortSNPs(int d) {
        this.sortSNPsByDataset = d;
    }

    @Override
    public int compareTo(WorkPackage o) {
        SNP otherSNP = o.getSnps()[sortSNPsByDataset];
        SNP currentSNP = snps[sortSNPsByDataset];

        if (this.equals(o)) {
            return 0;
        } else {
            if (currentSNP != null && otherSNP != null) {
                if (currentSNP.getId() > otherSNP.getId()) {
                    return 1;
                } else {
                    return -1;
                }
            } else {
                return 2;
            }
        }

    }

    public boolean equals(WorkPackage o) {
        SNP otherSNP = o.getSnps()[sortSNPsByDataset];
        SNP currentSNP = snps[sortSNPsByDataset];
        if (currentSNP != null && otherSNP != null) {
            return currentSNP.getId() == otherSNP.getId();
        } else {
            return false;
        }
    }

    public Boolean[] getFlipSNPAlleles() {
        return flipSNPAlleles;
    }

    public void setFlipSNPAlleles(Boolean[] b) {
        this.flipSNPAlleles = b;
    }

    public void setResult(Result dsResults) {
        results = dsResults;
    }

    public void clearResults() {
//        if(results.pvalues != null){
//            int nrPvalues = results.pvalues.length;
//            for(int p=0; p<nrPvalues; p++){
//                results.clearValues(p);
//            }
//        } 
//        results.pvalues = null;
        results = null;
        hasResults = false;
        numtested = 0;
        
    }

    public int getId() {
        return id;
    }

    public void setId(int s) {
        id = s;
    }

    public int getMetaSNPId() {
        return metaSNPId;
    }

    public void setMetaSNPId(int metaSNPId) {
        this.metaSNPId = metaSNPId;
    }

    public void setHasResults(boolean hasResults) {
        this.hasResults = hasResults;
    }
    
    public boolean getHasResults(){
        return hasResults;
    }
    
    synchronized public void incrementDatasetsPassingQC() {
        datasetsPassingQC++;
    }
}
