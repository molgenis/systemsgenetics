/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import cern.colt.matrix.tint.IntMatrix2D;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class QTL implements Comparable<QTL> {

    private double pvalue = Double.MAX_VALUE;
    private int pid = -1;
    private int sid = -1;
    private byte alleleAssessed;
    private double zscore = 0;
    private byte[] alleles;
    private double[] datasetZScores;
    private int[] datasetsSamples;
    private double[] correlations;
    

    public QTL(int datasets) {
        alleles = null;
        datasetZScores = null;
        datasetsSamples = null;
        correlations = null;
    }

    public QTL() {
    }

    public QTL(double pval, int pid, int sid, byte assessedAllele, double zscore, byte[] alleles, double[] zscores, int[] numSamples) {
        this.pvalue = pval;
        this.pid = pid;
        this.sid = sid;
        this.alleleAssessed = assessedAllele;
        this.zscore = zscore;
        this.alleles = alleles;
        this.datasetZScores = zscores;
        this.datasetsSamples = numSamples;
    }

    @Override
    public int compareTo(QTL o) {
        if (pvalue == o.pvalue) {
            if (Math.abs(zscore) == Math.abs(o.zscore)) {
                return 0;
            } else if (Math.abs(zscore) < Math.abs(o.zscore)) {
                return 1;
            } else {
                return -1;
            }
        } else if (pvalue > o.pvalue) {
            return 1;
        } else {
            return -1;
        }

    }

    public boolean equals(QTL o) {
        if (pvalue == o.pvalue) {
            if (Math.abs(zscore) == Math.abs(o.zscore)) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    public void cleanUp() {

        alleles = null;
        if (datasetZScores != null) {
            for (int i = 0; i < datasetZScores.length; i++) {
                datasetZScores[i] = Double.NaN;
            }
            datasetZScores = null;
        }
        if (datasetsSamples != null) {
            for (int i = 0; i < datasetsSamples.length; i++) {
                datasetsSamples[i] = -9;
            }
            datasetsSamples = null;
        }

        if (correlations != null) {
            for (int i = 0; i < correlations.length; i++) {
                correlations[i] = Double.NaN;
            }
            correlations = null;
        }
    }

    public double getPvalue() {
        return pvalue;
    }

    public double getZscore() {
        return zscore;
    }

    public double[] getCorrelations() {
        return correlations;
    }

    

    
}
