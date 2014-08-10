/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

/**
 *
 * @author harmjan
 */
public class EQTL implements Comparable<EQTL> {

    private Double pvalue = 1d;
    private Double pvalueAbs = 1d;
    private String rsName;
    private Byte rsChr;
    private Integer rsChrPos;
    private String probe;
    private Byte probeChr;
    private Integer probeChrPos;
    private String eQTLType;
    private String alleles;
    private String alleleAssessed;
    private String[] datasets;
    private Double zscore;
    private Double zscoreAbs;
    private Double[] datasetZScores;
    private Integer[] datasetsSamples;
    private Double[] probeMeans;
    private Double[] probeVariance;
    private String probeHUGO;
    private Double[] correlations;
    private Double FDR;
    private int pid;
    private String metabeta;
    private String beta;
    private String fc;
    private boolean useAbsoluteZScore = false;

    public EQTL() {
    }

//    public EQTL(Double pval, Result r, int p, boolean cisOnly, Integer[][] probeTranslation, TriTyperGeneticalGenomicsDataset[] gg, int midpointprobedist) {
//        Double[][] tmpcorrelations = r.correlations;
//        Integer[] numSamples = r.numSamples;
//        Double[][] zscores = r.zscores;
//        Double[] finalZscores = r.finalZScore;
//
//        int[] probes = r.workpackage.getProbes();
//        //SNP[] snps = r.workpackage.getSnps();
//        int numDatasets = zscores.length;
//
//        this.pvalue = pval;
//        this.zscore = finalZscores[p];
//
//        pid = p;
//
//        if (cisOnly && probes != null) {
//            pid = probes[p];
//        }
//
//        correlations    = new Double[numDatasets];
////        datasets        = new String[numDatasets];
//        datasetZScores  = new Double[numDatasets];
//        datasetsSamples = new Integer[numDatasets];
////        probeMeans      = new Double[numDatasets];
////        probeVariance   = new Double[numDatasets];
////        for(int d=0; d<numDatasets; d++){
////            if (r.workpackage.getSnps()[d] != null) {
//////                this.rsName     = r.workpackage.getSnps()[d].getName();
//////                this.rsChr      = r.workpackage.getSnps()[d].getChr();
//////                this.rsChrPos   = r.workpackage.getSnps()[d].getChrPos();
//////                byte[] tmpalleles = snps[d].getAlleles();
//////                String allele0 = BaseAnnot.toString(tmpalleles[0]);
//////                String allele1 = BaseAnnot.toString(tmpalleles[1]);
//////                alleleAssessed = BaseAnnot.toString(snps[d].getMinorAllele());
//////                alleles = allele0 + "/" + allele1;
////                break;
////            }
////        }
////
////        for(int d=0; d<numDatasets; d++){
////            Integer probeId = new Integer (probeTranslation[d][pid].intValue());
////            if (probeId != null && snps[d] != null) {
////                this.probe      = gg[d].getExpressionData().getProbes()[probeId];
////                this.probeChr   = gg[d].getExpressionData().getChr()[probeId];
////                this.probeChrPos= (int) Math.floor((double) (gg[d].getExpressionData().getChrStart()[probeId] + gg[d].getExpressionData().getChrStop()[probeId]) / 2);
////                this.probeHUGO  = gg[d].getExpressionData().getAnnotation()[probeId];
////                if (Math.abs(this.probeChrPos - this.rsChrPos) < midpointprobedist) {
////                    eQTLType = "cis";
////                } else {
////                    eQTLType = "trans";
////                }
////                break;
////            }
////        }
////
////
////        for(int d=0; d<numDatasets; d++){
////            correlations[d]     = tmpcorrelations[d][p];
////            if(correlations[d] != null){
//////                datasets[d]         = gg[d].getSettings().name;
////                datasetZScores[d]   = zscores[d][p];
////                datasetsSamples[d]  = numSamples[d];
////            }
////        }
//
//
////        for(int d=0; d<numDatasets; d++){
////            Integer probeId = probeTranslation[d][pid];
////            if (probeId != null && correlations[d] != null) {
////                probeMeans[d]       = gg[d].getExpressionData().getOriginalProbeMean()[probeId];
////                probeVariance[d]    = gg[d].getExpressionData().getOriginalProbeVariance()[probeId];
////            }
////        }
//    }
    /**
     * @return the pvalue
     */
    public Double getPvalue() {
        return pvalue;
    }

    /**
     * @param pvalue the pvalue to set
     */
    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    /**
     * @return the rsName
     */
    public String getRsName() {
        return rsName;
    }

    /**
     * @param rsName the rsName to set
     */
    public void setRsName(String rsName) {
        this.rsName = rsName;
    }

    /**
     * @return the rsChr
     */
    public Byte getRsChr() {
        return rsChr;
    }

    /**
     * @param rsChr the rsChr to set
     */
    public void setRsChr(byte rsChr) {
        this.rsChr = rsChr;
    }

    /**
     * @return the rsChrPos
     */
    public Integer getRsChrPos() {
        return rsChrPos;
    }

    /**
     * @param rsChrPos the rsChrPos to set
     */
    public void setRsChrPos(int rsChrPos) {
        this.rsChrPos = rsChrPos;
    }

    /**
     * @return the probe
     */
    public String getProbe() {
        return probe;
    }

    /**
     * @param probe the probe to set
     */
    public void setProbe(String probe) {
        this.probe = probe;
    }

    /**
     * @return the probeChr
     */
    public Byte getProbeChr() {
        return probeChr;
    }

    /**
     * @param probeChr the probeChr to set
     */
    public void setProbeChr(byte probeChr) {
        this.probeChr = probeChr;
    }

    /**
     * @return the probeChrPos
     */
    public Integer getProbeChrPos() {
        return probeChrPos;
    }

    /**
     * @param probeChrPos the probeChrPos to set
     */
    public void setProbeChrPos(int probeChrPos) {
        this.probeChrPos = probeChrPos;
    }

    /**
     * @return the type
     */
    public String getType() {
        return eQTLType;
    }

    /**
     * @param type the type to set
     */
    public void setType(String type) {
        this.eQTLType = type;
    }

    /**
     * @return the alleles
     */
    public String getAlleles() {
        return alleles;
    }

    /**
     * @param alleles the alleles to set
     */
    public void setAlleles(String alleles) {
        this.alleles = alleles;
    }

    /**
     * @return the alleleAssessed
     */
    public String getAlleleAssessed() {
        return alleleAssessed;
    }

    /**
     * @param alleleAssessed the alleleAssessed to set
     */
    public void setAlleleAssessed(String alleleAssessed) {
        this.alleleAssessed = alleleAssessed;
    }

    /**
     * @return the datasets
     */
    public String[] getDatasets() {
        return datasets;
    }

    /**
     * @param datasets the datasets to set
     */
    public void setDatasets(String[] datasets) {
        this.datasets = datasets;
    }

    /**
     * @return the zscore
     */
    public double getZscore() {
        return zscore;
    }

    /**
     * @param zscore the zscore to set
     */
    public void setZscore(double zscore) {
        this.zscore = zscore;
    }

    /**
     * @return the datasetZScores
     */
    public Double[] getDatasetZScores() {
        return datasetZScores;
    }

    /**
     * @param datasetZScores the datasetZScores to set
     */
    public void setDatasetZScores(Double[] datasetZScores) {
        this.datasetZScores = datasetZScores;
    }

    /**
     * @return the datasetsSamples
     */
    public Integer[] getDatasetsSamples() {
        return datasetsSamples;
    }

    /**
     * @param datasetsSamples the datasetsSamples to set
     */
    public void setDatasetsSamples(Integer[] datasetsSamples) {
        this.datasetsSamples = datasetsSamples;
    }

    /**
     * @return the probeMeans
     */
    public Double[] getProbeMeans() {
        return probeMeans;
    }

    /**
     * @param probeMeans the probeMeans to set
     */
    public void setProbeMeans(Double[] probeMeans) {
        this.probeMeans = probeMeans;
    }

    /**
     * @return the probeVariance
     */
    public Double[] getProbeVariance() {
        return probeVariance;
    }

    /**
     * @param probeVariance the probeVariance to set
     */
    public void setProbeVariance(Double[] probeVariance) {
        this.probeVariance = probeVariance;
    }

    /**
     * @return the probeHUGO
     */
    public String getProbeHUGO() {
        return probeHUGO;
    }

    /**
     * @param probeHUGO the probeHUGO to set
     */
    public void setProbeHUGO(String probeHUGO) {
        this.probeHUGO = probeHUGO;
    }

    /**
     * @return the correlations
     */
    public Double[] getCorrelations() {
        return correlations;
    }

    /**
     * @param correlations the correlations to set
     */
    public void setCorrelations(Double[] correlations) {
        this.correlations = correlations;
    }

    public void setFDR(double d) {
        this.FDR = d;
    }

    public Double getFDR() {
        return FDR;
    }

    public String compare(EQTL test) {
        boolean identical = true;
        String reason = "";

        if (!test.getProbe().equals(this.probe)) {
            reason += "Diff probes:\t" + test.getProbe() + "\t" + getProbe() + "\t";
            identical = false;
        } else {
            if (probeChr != null && probeChrPos != null && test.getProbeChr() != null && test.getProbeChrPos() != null) {
//                if(!test.getProbeChr().equals(this.probeChr)){
//                    reason += "Diff probe chr:\t"+test.getProbeChr()+"\t"+getProbeChr()+"\t";
//                    identical = false;
//                }
//                if(!test.getProbeChrPos().equals(this.probeChrPos)){
//                    reason += "Diff probe chrPos:\t"+test.getProbeChrPos()+"\t"+getProbeChrPos()+"\t";
//                    identical = false;
//                }
            }
        }

        if (!test.getRsName().equals(this.rsName)) {
            reason += "Diff rsName:\t" + test.getRsName() + "\t" + getRsName() + "\t";
            identical = false;
        } else {
            if (rsChr != null && rsChrPos != null && test.getRsChrPos() != null && test.getRsChr() != null) {
                // test position etc
//                if(!test.getRsChr().equals(rsChr)){
//                    reason += "Diff rsChr:\t"+test.getRsChr()+"\t"+getRsChr()+"\t";
//                    identical = false;
//                }
//                if(!test.getRsChrPos().equals(rsChrPos)){
//                    reason += "Diff rsChrPos:\t"+test.getRsChrPos()+"\t"+getRsChrPos()+"\t";
//                    identical = false;
//                }
            }
        }

        if (test.getPvalue() != pvalue) {
            reason += "Diff pval:\t" + test.getPvalue() + "\t" + getPvalue() + "\t";
            identical = false;
        } else {
            reason += "";
        }

        if (test.getZscore() != zscore.doubleValue()) {
            if (Math.abs(test.getZscore()) - Math.abs(zscore.doubleValue()) > 0.0001) {
                if (test.getAlleleAssessed().equals(alleleAssessed)) {
                    reason += "Diff zscore:\t" + test.getZscore() + "\t" + getZscore() + "\t";
                    identical = false;
                }
            }
        }

        if (!identical) {
            return reason;
        } else {
            return null;
        }
    }

    @Override
    public int hashCode() {
        int hash = 1;
        hash += probe.hashCode();
        hash += rsName.hashCode();
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
        final EQTL other = (EQTL) obj;
        if ((this.rsName == null) ? (other.rsName != null) : !this.rsName.equals(other.rsName)) {
            return false;
        }
        if ((this.probe == null) ? (other.probe != null) : !this.probe.equals(other.probe)) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        String sepStr = ";";
        String nullstr = "-";
        char tabStr = '\t';
        EQTL e = this;
        StringBuilder out = new StringBuilder();

        if (useAbsoluteZScore) {
            if (pvalueAbs == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(e.getPvalueAbs());
                out.append(tabStr);
            }

        } else {
            if (pvalue == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(e.getPvalue());
                out.append(tabStr);
            }

        }
        if (rsName == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getRsName());
            out.append(tabStr);
        }

        if (rsChr == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getRsChr());
            out.append(tabStr);
        }

        if (rsChrPos == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getRsChrPos());
            out.append(tabStr);
        }

        if (probe == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getProbe());
            out.append(tabStr);
        }

        if (probeChr == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getProbeChr());
            out.append(tabStr);
        }

        if (probeChrPos == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getProbeChrPos());
            out.append(tabStr);
        }

        if (eQTLType == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getType());
            out.append(tabStr);
        }

        if (alleles == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getAlleles());
            out.append(tabStr);
        }

        if (alleleAssessed == null) {
            out.append(nullstr);
            out.append(tabStr);
        } else {
            out.append(e.getAlleleAssessed());
            out.append(tabStr);
        }

        if (useAbsoluteZScore) {
            if (zscoreAbs == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(e.getZscoreAbs());
                out.append(tabStr);
            }
        } else {
            if (zscore == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(e.getZscore());
                out.append(tabStr);
            }
        }

        String[] ds = e.getDatasets();
//        Double[] corrs = e.getCorrelations();
//        Double[] zscores = e.getDatasetZScores();
//        Double[] probevars = e.getProbeVariance();
        Double[] probemeans = e.getProbeMeans();
//        Integer[] numsamples = e.getDatasetsSamples();

        StringBuilder outcorrs = new StringBuilder();
        StringBuilder outzscores = new StringBuilder();
        StringBuilder outsamples = new StringBuilder();
        StringBuilder outmeans = new StringBuilder();
        StringBuilder outvars = new StringBuilder();

        if (ds == null) {
            out.append(nullstr);
            out.append(tabStr);
            out.append(nullstr);
            out.append(tabStr);
            out.append(nullstr);
            out.append(tabStr);
            out.append(nullstr);
            out.append(tabStr);
            out.append(nullstr);
            out.append(tabStr);
            out.append(nullstr);
            out.append(tabStr);
            out.append(nullstr);
        } else {
            for (int d = 0; d < ds.length; d++) {
                if (d == 0) {
                    if (datasets[d] == null) {
                        out.append(nullstr);
                    } else {
                        out.append(ds[d]);
                    }

                    if (correlations == null || correlations[d] == null) {
                        outcorrs.append(nullstr);
                    } else {
                        outcorrs.append(correlations[d]);
                    }

                    if (datasetZScores == null || datasetZScores[d] == null) {
                        outzscores.append(nullstr);
                    } else {
                        outzscores.append(datasetZScores[d]);
                    }

                    if (datasetsSamples == null || datasetsSamples[d] == null) {
                        outsamples.append(nullstr);
                    } else {
                        outsamples.append(datasetsSamples[d]);
                    }

                    if (probeMeans == null || probeMeans[d] == null) {
                        outmeans.append(nullstr);
                    } else {
                        outmeans.append(probemeans[d]);
                    }

                    if (probeVariance == null || probeVariance[d] == null) {
                        outvars.append(nullstr);
                    } else {
                        outvars.append(probeVariance[d]);
                    }

                } else {
                    if (datasets[d] == null) {
                        out.append(sepStr).append(nullstr);
                    } else {
                        out.append(sepStr).append(ds[d]);
                    }

                    if (correlations == null || correlations[d] == null) {
                        outcorrs.append(sepStr).append(nullstr);
                    } else {
                        outcorrs.append(sepStr).append(correlations[d]);
                    }

                    if (datasetZScores == null || datasetZScores[d] == null) {
                        outzscores.append(sepStr).append(nullstr);
                    } else {
                        outzscores.append(sepStr).append(datasetZScores[d]);
                    }

                    if (datasetsSamples == null || datasetsSamples[d] == null) {
                        outsamples.append(sepStr).append(nullstr);
                    } else {
                        outsamples.append(sepStr).append(datasetsSamples[d]);
                    }

                    if (probeMeans == null || probeMeans[d] == null) {
                        outmeans.append(sepStr).append(nullstr);
                    } else {
                        outmeans.append(sepStr).append(probemeans[d]);
                    }

                    if (probeVariance == null || probeVariance[d] == null) {
                        outvars.append(sepStr).append(nullstr);
                    } else {
                        outvars.append(sepStr).append(probeVariance[d]);
                    }
                }
            }

            out.append(tabStr);
            out.append(outzscores.toString());
            out.append(tabStr);
            out.append(outsamples.toString());
            out.append(tabStr);
            out.append(outmeans.toString());
            out.append(tabStr);
            out.append(outvars.toString());
            out.append(tabStr);
            out.append(e.getProbeHUGO());
            out.append(tabStr);
            out.append(outcorrs.toString());
            out.append(tabStr);

            if (metabeta == null) {
                out.append(nullstr);
            } else {
                out.append(e.getMetaBeta());
            }
            out.append(tabStr);
            if (beta == null) {
                out.append(nullstr);
            } else {
                out.append(e.getBeta());
            }
            out.append(tabStr);
            if (fc == null) {
                out.append(nullstr);

            } else {
                out.append(e.getFC());
            }
            out.append(tabStr);
            if (FDR == null) {
                out.append(nullstr);
            } else {
                out.append(e.getFDR());
            }

        }

        return out.toString();
    }

    @Override
    public int compareTo(EQTL o) {
        if (useAbsoluteZScore) {
            if (pvalueAbs.doubleValue() == o.pvalueAbs.doubleValue()) {
                if (zscoreAbs.doubleValue() == o.zscoreAbs.doubleValue()) {
                    if (rsName.compareTo(o.rsName) == 0) {
                        return probe.compareTo(o.probe);
                    } else {
                        return rsName.compareTo(o.rsName);
                    }
                } else if (zscoreAbs.doubleValue() < o.zscoreAbs.doubleValue()) {
                    return 1;
                } else {
                    return -1;
                }
            } else if (pvalueAbs.doubleValue() > o.pvalueAbs.doubleValue()) {
                return 1;
            } else {
                return -1;
            }
        } else {
            if (pvalue.doubleValue() == o.pvalue.doubleValue()) {
                if (Math.abs(zscore.doubleValue()) == Math.abs(o.zscore.doubleValue())) {
                    if (rsName.compareTo(o.rsName) == 0) {
                        return probe.compareTo(o.probe);
                    } else {
                        return rsName.compareTo(o.rsName);
                    }
                } else if (Math.abs(zscore.doubleValue()) < Math.abs(o.zscore.doubleValue())) {
                    return 1;
                } else {
                    return -1;
                }
            } else if (pvalue.doubleValue() > o.pvalue.doubleValue()) {
                return 1;
            } else {
                return -1;
            }
        }
    }

    public int compareToVerbose(EQTL o) {
        if (useAbsoluteZScore) {
            if (pvalueAbs.doubleValue() == o.pvalueAbs.doubleValue()) {

                if (zscoreAbs.doubleValue() == o.zscoreAbs.doubleValue()) {
                    if (rsName.compareTo(o.rsName) == 0) {
                        return probe.compareTo(o.probe);
                    } else {
                        return rsName.compareTo(o.rsName);
                    }
                } else if (zscoreAbs.doubleValue() < o.zscoreAbs.doubleValue()) {
                    return 1;
                } else {
                    return -1;
                }
            } else if (pvalueAbs.doubleValue() > o.pvalueAbs.doubleValue()) {
                return 1;
            } else {
                return -1;
            }
        } else {
            if (pvalue.doubleValue() == o.pvalue.doubleValue()) {
                System.out.println("p value identical: " + pvalue + "\t" + o.pvalue);
                if (Math.abs(zscore.doubleValue()) == Math.abs(o.zscore.doubleValue())) {
                    System.out.println("Z value identical: " + zscore + "\t" + o.zscore);
                    if (rsName.compareTo(o.rsName) == 0) {
                        System.out.println("rs value identical: " + rsName + "\t" + o.rsName);

                        System.out.println("probe: " + probe + "\t" + o.probe + "\t" + probe.compareTo(o.probe));
                        return probe.compareTo(o.probe);
                    } else {
                        System.out.println("probe: " + rsName + "\t" + o.rsName + "\t" + rsName.compareTo(o.rsName));
                        return rsName.compareTo(o.rsName);
                    }
                } else if (Math.abs(zscore.doubleValue()) < Math.abs(o.zscore.doubleValue())) {
                    return 1;
                } else {
                    return -1;
                }
            } else if (pvalue.doubleValue() > o.pvalue.doubleValue()) {
                return 1;
            } else {
                return -1;
            }
        }
    }

    public boolean equals(EQTL o) {
        if (useAbsoluteZScore) {
            if (pvalueAbs.doubleValue() == o.pvalueAbs.doubleValue()) {
                if (zscoreAbs.doubleValue() == o.zscoreAbs.doubleValue()) {
                    if (probe.equals(o.probe)) {
                        return rsName.equals(o.rsName);
                    } else {
                        return false;
                    }
                } else {
                    return false;
                }
            } else {
                return false;
            }
        } else {
            if (pvalue.doubleValue() == o.pvalue.doubleValue()) {
                if (Math.abs(zscore.doubleValue()) == Math.abs(o.zscore.doubleValue())) {
                    if (probe.equals(o.probe)) {
                        return rsName.equals(o.rsName);
                    } else {
                        return false;
                    }
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
    }

    public boolean sameQTL(EQTL o) {
        return (zscore.equals(o.zscore) && probe.equals(o.probe) && rsName.equals(o.rsName) && FDR.equals(o.FDR) && beta.equals(o.beta));
    }

//    @Override
//    public String toString(){
//        return pvalue+"\t"+description;
//    }
    public void setMetaBeta(String string) {
        metabeta = string;
    }

    public void setBeta(String string) {
        beta = string;
    }

    public void setFC(String string) {
        fc = string;
    }

    public String getMetaBeta() {
        return metabeta;
    }

    public String getBeta() {
        return beta;
    }

    public String getFC() {
        return fc;
    }

    public void clearData() {
//	pvalue = null;
        rsName = null;
        rsChr = null;
        rsChrPos = null;
        probe = null;
        probeChr = null;
        probeChrPos = null;
        eQTLType = null;
        alleles = null;
        alleleAssessed = null;
        datasets = null;
        zscore = null;
        datasetZScores = null;
        datasetsSamples = null;
        probeMeans = null;
        probeVariance = null;
        probeHUGO = null;
        correlations = null;
        FDR = null;
        metabeta = null;
        beta = null;
        fc = null;
    }

    public void setZscoreAbs(double zScoreAbs) {
        this.zscoreAbs = zScoreAbs;
    }

    public void setPvalueAbs(double pValueOverallAbs) {
        this.pvalueAbs = pValueOverallAbs;
    }

    public Double getZscoreAbs() {
        return zscoreAbs;
    }

    public Double getPvalueAbs() {
        return pvalueAbs;
    }

    public void setUseAbsoluteZScore() {
        this.useAbsoluteZScore = true;
    }

    public boolean getUseAbsoluteZScore() {
        return useAbsoluteZScore;
    }

    public String getDiff(EQTL o) {

        return "z: " + zscore + " - " + o.zscore + "\t"
                + "probe: " + probe + " - " + o.probe + "\t"
                + "rsName: " + rsName + " - " + o.rsName + "\t"
                + "FDR: " + FDR + " - " + o.FDR + "\t"
                + "Beta: " + beta + " - " + o.beta + "\t";

    }

}
