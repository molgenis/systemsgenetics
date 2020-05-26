/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class QTL implements Comparable<QTL> {

    private double pvalue = Double.MAX_VALUE;
    private MetaQTL4MetaTrait trait;
    private int sid = -1;
    private byte alleleAssessed;
    private double zscore = 0;
    private byte[] alleles;
    private float[] datasetZScores;
    private int[] datasetsSamples;

    public QTL(int datasets) {
        alleles = null;
        datasetZScores = null;
        datasetsSamples = null;
    }

    public QTL() {
    }

    public QTL(double pval, MetaQTL4MetaTrait t, int sid, byte assessedAllele, double zscore, byte[] alleles, float[] zscores, int[] numSamples) {
        this.pvalue = pval;
        this.trait = t;
        this.sid = sid;
        this.alleleAssessed = assessedAllele;
        this.zscore = zscore;
        this.alleles = alleles;
        this.datasetZScores = zscores;
        this.datasetsSamples = numSamples;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + (int) (Double.doubleToLongBits(this.pvalue) ^ (Double.doubleToLongBits(this.pvalue) >>> 32));
        hash = 97 * hash + (int) (Double.doubleToLongBits(Math.abs(this.zscore)) ^ (Double.doubleToLongBits(Math.abs(this.zscore)) >>> 32));
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final QTL other = (QTL) obj;

        if (Double.doubleToLongBits(this.pvalue) != Double.doubleToLongBits(other.pvalue)) {
            return false;
        }
        return Double.doubleToLongBits(Math.abs(this.zscore)) == Double.doubleToLongBits(Math.abs(other.zscore));
    }

    @Override
    public String toString() {
        return "QTL{" + "pvalue=" + pvalue + ", zscore=" + zscore + '}';
    }

    @Override
    public int compareTo(QTL o) {
        int returnval = -1;
        if (pvalue == o.pvalue) {
            if (Math.abs(zscore) == Math.abs(o.zscore)) {
                returnval = 0;
            } else if (Math.abs(zscore) > Math.abs(o.zscore)) {
                returnval = -1;
            } else {
                returnval = 1;
            }
        } else if (pvalue > o.pvalue) {
            returnval = 1;
        } else {
            returnval = -1;
        }
        return returnval;
    }

//    public boolean equals(QTL o) {
//        if (pvalue == o.pvalue) {
//            return Math.abs(zscore) == Math.abs(o.zscore);
//        } else {
//            return false;
//        }
//    }
    public void cleanUp() {

        alleles = null;
        if (datasetZScores != null) {
            for (int i = 0; i < datasetZScores.length; i++) {
                datasetZScores[i] = Float.NaN;
            }
            datasetZScores = null;
        }
        if (datasetsSamples != null) {
            for (int i = 0; i < datasetsSamples.length; i++) {
                datasetsSamples[i] = -9;
            }
            datasetsSamples = null;
        }

    }

    public double getPvalue() {
        return pvalue;
    }

    public double getZscore() {
        return zscore;
    }

    public int getSNPId() {
        return sid;
    }

    public String getAlleles() {
        return BaseAnnot.getAllelesDescription(alleles);
    }

    public String getAlleleAssessed() {
        return BaseAnnot.toString(alleleAssessed);
    }

    public float[] getDatasetZScores() {
        return datasetZScores;
    }

    public int[] getDatasetSampleSizes() {
        return datasetsSamples;
    }

    public MetaQTL4MetaTrait getMetaTrait() {
        return trait;
    }

}
