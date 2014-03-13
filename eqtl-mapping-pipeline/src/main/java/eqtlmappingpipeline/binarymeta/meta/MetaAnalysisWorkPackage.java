/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta;

import java.util.HashSet;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;

/**
 *
 * @author harmjan
 */
public class MetaAnalysisWorkPackage implements Comparable<MetaAnalysisWorkPackage> {

    private boolean poison = false;
    private EQTL[] result;
    private byte[][] datasetdata;
    private int sortForDataset = 0;
    private Integer[] snpids;
    private BinaryResultSNP[] snpObjects;
    private int snpNumber;
    private String zscoretabeline;
    private int numberOfTestedProbes;
    private Integer[] testedProbesList;
    private boolean passedQC;

    public MetaAnalysisWorkPackage(int snpnumber, int numDatasets) {
        snpNumber = snpnumber;
        snpids = new Integer[numDatasets];
        snpObjects = new BinaryResultSNP[numDatasets];
        datasetdata = new byte[numDatasets][0];
    }

    boolean getPoison() {
        return poison;
    }

    public void setSNPObject(int d, BinaryResultSNP snpObj) {
        snpObjects[d] = snpObj;
    }

    public BinaryResultSNP getSNPObject(int d) {
        return snpObjects[d];
    }

    public void setResult(EQTL[] result) {
        this.result = result;
    }

    public void poisonTheWell() {
        poison = true;
    }

    public EQTL[] getResult() {
        return result;
    }

    public Integer getSNPId(int d) {
        return snpids[d];
    }

    @Override
    public int compareTo(MetaAnalysisWorkPackage o) {
        Integer otherSNPId = o.getSNPId(sortForDataset);
        Integer currentSNPId = snpids[sortForDataset];
        if (currentSNPId != null && otherSNPId != null) {
            return currentSNPId - otherSNPId;
        } else if (currentSNPId == null) {
            return Integer.MAX_VALUE;
        } else {
            return Integer.MIN_VALUE;
        }
    }

    public boolean equals(MetaAnalysisWorkPackage o) {
        Integer otherSNPId = o.getSNPId(sortForDataset);
        Integer currentSNPId = snpids[sortForDataset];
        if (currentSNPId != null && otherSNPId != null) {
            if (currentSNPId.equals(otherSNPId)) {
                return true;
            } else {
                return false;
            }

        } else {
            return false;
        }
    }

    public void setSNPId(int d, Integer snpid) {
        snpids[d] = snpid;
    }

    public void setSortByDataset(int d) {
        sortForDataset = d;
    }

    public byte[] getData(int d) {
        return datasetdata[d];
    }

    public void setData(int d, byte[] datasetdata) {
        this.datasetdata[d] = datasetdata;
    }

    public int getSNPNum() {
        return snpNumber;
    }

    public void setZScoreOut(String toString) {
        zscoretabeline = toString;
    }

    public String getZScoreOut() {
        return zscoretabeline;
    }

    public void clearData() {
        for (int i = 0; i < result.length; i++) {
            result[i] = null;
        }
        result = null;
        for (int i = 0; i < snpObjects.length; i++) {

            snpids[i] = null;
            snpObjects[i] = null;
        }
        snpids = null;

    }

    public void clearByteData() {
        for (int i = 0; i < snpObjects.length; i++) {
            datasetdata[i] = null;
        }
        datasetdata = null;
    }

    void setNumOfTestedProbes(int probesTested) {
        this.numberOfTestedProbes = probesTested;
    }

    public int getNumberOfTestedProbes() {
        return this.numberOfTestedProbes;
    }

    public void setProbesTestedHash(HashSet<Integer> probesTestedHash) {
        this.testedProbesList = probesTestedHash.toArray(new Integer[0]);
    }

    public Integer[] getListOfTestedProbes() {
        return testedProbesList;
    }

    public boolean getPassedQC() {
        return passedQC;
    }

    public void setPassedQC(boolean b) {
        passedQC = b;
    }
}
