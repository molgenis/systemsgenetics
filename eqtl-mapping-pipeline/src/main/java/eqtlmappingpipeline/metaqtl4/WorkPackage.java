/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import umcg.genetica.io.trityper.SNP;

/**
 *
 * @author harmjan
 */
public class WorkPackage implements Comparable<WorkPackage> {

    int snp = -1;
    int[] probes = null;
    SNP[] snps = null;
    int workPackageId = -1;
    private Boolean[] flipSNPAlleles;
    private short datasetsPassingQC;
    private int sortSNPsByDataset = 0;
    private int metaSNPId = -1;
    private Integer[][] pvalueIndexes;

    public WorkPackage(int snp, int[] probes, SNP[] snps) {
        this.metaSNPId = snp;
        this.probes = probes;
        this.snps = snps;
    }

    public void setId(int id) {
        this.workPackageId = id;
    }

    public int[] getProbes() {
        return probes;
    }

    public SNP[] getSNPs() {
        return snps;
    }

    public Boolean[] getFlipSNPAlleles() {
        return flipSNPAlleles;
    }

    public void setFlipSNPAlleles(Boolean[] b) {
        this.flipSNPAlleles = b;
    }

    /**
     * @return the datasetsPassingQC
     */
    public short getDatasetsPassingQC() {
        return datasetsPassingQC;
    }

    public void setDatasetToSortSNPs(int d) {
        this.sortSNPsByDataset = d;
    }

    public int getMetaSNPId() {
        return metaSNPId;
    }
    
    public void setMetaSNPId(int metaSNPId) {
        this.metaSNPId = metaSNPId;
    }

    /**
     * @param datasetsPassingQC the datasetsPassingQC to set
     */
    public void setDatasetsPassingQC(short datasetsPassingQC) {
        this.datasetsPassingQC = datasetsPassingQC;
    }

    public int compareTo(WorkPackage o) {
        SNP otherSNP = o.getSNPs()[sortSNPsByDataset];
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
        SNP otherSNP = o.getSNPs()[sortSNPsByDataset];
        SNP currentSNP = snps[sortSNPsByDataset];
        if (currentSNP != null && otherSNP != null) {
            if (currentSNP.getId() == otherSNP.getId()) {
                return true;
            } else {
                return false;
            }

        } else {
            return false;
        }
    }

    public void setResults(Integer[][] pvalueIndexes) {
        this.pvalueIndexes = pvalueIndexes;
    }
    
    public Integer[][] getResults() {
        return pvalueIndexes;
    }

    void clearResults() {
        this.pvalueIndexes = null;
    }
}
