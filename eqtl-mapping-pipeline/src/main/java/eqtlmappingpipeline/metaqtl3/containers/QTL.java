/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.containers;

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
    private double[] datasetfc;
    private double[] datasetbeta;
    private double[] datasetbetase;
    private double finalbeta;
    private double finalbetase;

    public QTL(int datasets) {
        alleles = null;
        datasetZScores = null;
        datasetsSamples = null;
        correlations = null;
    }

    public QTL() {
    }

    public QTL(double pval, int pid, int sid, byte assessedAllele, double zscore, byte[] alleles, double[] zscores, int[] numSamples, double[] correlations, double[] fc, double[] beta, double[] betase, double finalbeta, double finalbetase) {
        this.pvalue = pval;
        this.pid = pid;
        this.sid = sid;
        this.alleleAssessed = assessedAllele;
        this.zscore = zscore;
        this.alleles = alleles;
        this.datasetZScores = zscores;
        this.datasetsSamples = numSamples;
        this.correlations = correlations;
        this.datasetfc = fc;
        this.datasetbeta = beta;
        this.datasetbetase = betase;
        this.finalbeta = finalbeta;
        this.finalbetase = finalbetase;
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

//    void copy(QTL eQTL) {
//        this.pvalue = eQTL.pvalue;
//        this.pid = eQTL.pid;
//        this.sid = eQTL.sid;
//        this.alleles = eQTL.alleles;
//        this.datasetZScores = eQTL.datasetZScores;
//        this.datasetsSamples = eQTL.datasetsSamples;
//        this.correlations = eQTL.correlations;
//        this.alleleAssessed = eQTL.alleleAssessed;
//
//        this.datasetfc = eQTL.datasetfc;
//        this.datasetbeta = eQTL.datasetbeta;
//        this.datasetbetase = eQTL.datasetbetase;
//
//        this.finalbeta = eQTL.finalbeta;
//        this.finalbetase = eQTL.finalbetase;
//
//        this.zscore = eQTL.zscore;
//    }
    public String getDescription(WorkPackage[] workPackages, IntMatrix2D probeTranslation, TriTyperGeneticalGenomicsDataset[] gg,
            int maxCisDistance) {
        String sepStr = ";";
        String nullstr = "-";
        char tabStr = '\t';

        StringBuilder out = new StringBuilder();

        out.append(pvalue);
        out.append(tabStr);

        String rsName = null;
        Byte rsChr = null;
        Integer rsChrPos = null;

        if (sid == -1 && pid == -1) {
            return null;
        } else {
            WorkPackage currentWP = workPackages[sid];
            SNP[] snps = workPackages[sid].getSnps();

            for (int d = 0; d < snps.length; d++) {
                if (snps[d] != null) {
                    rsName = snps[d].getName();
                    rsChr = snps[d].getChr();
                    rsChrPos = snps[d].getChrPos();
                    break;
                }
            }

            String probe = null;
            Byte probeChr = null;
            Integer probeChrPos = null;

            for (int d = 0; d < snps.length; d++) {
                if (probeTranslation.get(d, pid) != -9) {
                    int probeId = probeTranslation.get(d, pid);
                    probe = gg[d].getExpressionData().getProbes()[probeId];
                    probeChr = gg[d].getExpressionData().getChr()[probeId];
                    probeChrPos = (gg[d].getExpressionData().getChrStart()[probeId] + gg[d].getExpressionData().getChrStop()[probeId]) / 2;
                    break;
                }
            }

            out.append(rsName);
            out.append(tabStr);

            if (rsChr == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(rsChr);
                out.append(tabStr);
            }

            if (rsChrPos == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(rsChrPos);
                out.append(tabStr);
            }

            if (probe == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(probe);
                out.append(tabStr);
            }

            if (probeChr == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(probeChr);
                out.append(tabStr);
            }

            if (probeChrPos == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(probeChrPos);
                out.append(tabStr);
            }

            String eQTLType = "trans";
            if (rsChr != null && probeChr != null && probeChrPos != null && rsChrPos != null) {
                if (rsChr.equals(probeChr) && Math.abs(probeChrPos - rsChrPos) < maxCisDistance) {
                    eQTLType = "cis";
                }
            }

//            if (eQTLType == null) {
//                out.append(nullstr);
//                out.append(tabStr);
//            } else {
            out.append(eQTLType);
            out.append(tabStr);
//            }

            if (alleles == null) {
                System.err.println(rsName + " has null alleles..?'\n" + out.toString());
                return null;
            }

            out.append(BaseAnnot.toString(alleles[0])).append("/").append(BaseAnnot.toString(alleles[1]));
            out.append(tabStr);
            out.append(BaseAnnot.toString(alleleAssessed));
            out.append(tabStr);
            out.append(zscore);
            out.append(tabStr);

            String[] ds = new String[gg.length];
            Double[] probevars = new Double[gg.length];
            Double[] probemeans = new Double[gg.length];

            String hugo = nullstr;
            for (int d = 0; d < gg.length; d++) {
                if (!Double.isNaN(correlations[d])) {
                    ds[d] = gg[d].getSettings().name;
                    if (probeTranslation.get(d, pid) != -9) {
                        int probeId = probeTranslation.get(d, pid);
                        probevars[d] = gg[d].getExpressionData().getOriginalProbeVariance()[probeId];
                        probemeans[d] = gg[d].getExpressionData().getOriginalProbeMean()[probeId];
                        hugo = gg[d].getExpressionData().getAnnotation()[probeId];
                    } else {
                        System.out.println("ERROR!!!");
                    }
                } else {
                    ds[d] = null;
                    probevars[d] = null;
                    probemeans[d] = null;
                }
            }

            StringBuilder outcorrs = new StringBuilder();
            StringBuilder outzscores = new StringBuilder();
            StringBuilder outsamples = new StringBuilder();
            StringBuilder outmeans = new StringBuilder();
            StringBuilder outvars = new StringBuilder();
            StringBuilder outfc = new StringBuilder();
            StringBuilder outbeta = new StringBuilder();

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
                        sepStr = "";
//                    if(ds[d] == null){
//                        out.append(nullstr);
//                    } else {
//                        out.append(ds[d]);
//                    }
//
//                    if(correlations == null || correlations[d] == null){
//                        outcorrs.append(nullstr);
//                    } else {
//                        outcorrs.append(correlations[d]);
//                        if(currentWP.getFlipSNPAlleles()[d]){
//                            outcorrs.append(-correlations[d]);
//                        } else {
//                            outcorrs.append(correlations[d]);
//                        }
//                    }
//
//                    if(datasetZScores == null || datasetZScores[d] == null){
//                        outzscores.append(nullstr);
//                    } else {
////                        outzscores.append(datasetZScores[d]);
//                        if(currentWP.getFlipSNPAlleles()[d]){
//                            outzscores.append(-datasetZScores[d]);
//                        } else {
//                            outzscores.append(datasetZScores[d]);
//                        }
//                    }
//
//                    if(datasetsSamples == null || datasetsSamples[d] == null){
//                        outsamples.append(nullstr);
//                    } else {
//                        outsamples.append(datasetsSamples[d]);
//                    }
//
//                    if(probemeans == null || probemeans[d] == null){
//                        outmeans.append(nullstr);
//                    } else {
//                        outmeans.append(probemeans[d]);
//                    }
//
//                    if(probevars == null || probevars[d] == null){
//                        outvars.append(nullstr);
//                    } else {
//                        outvars.append(probevars[d]);
//                    }

                    } else {
                        sepStr = ";";
                    }

                    if (ds[d] == null) {
                        out.append(sepStr).append(nullstr);
                    } else {
                        out.append(sepStr).append(ds[d]);
                    }

                    if (correlations == null || Double.isNaN(correlations[d])) {
                        outcorrs.append(sepStr).append(nullstr);
                    } else {
                        if (currentWP.getFlipSNPAlleles()[d]) {
                            outcorrs.append(sepStr).append(-correlations[d]);
                        } else {
                            outcorrs.append(sepStr).append(correlations[d]);
                        }

                    }

                    if (datasetZScores == null || Double.isNaN(datasetZScores[d])) {
                        outzscores.append(sepStr).append(nullstr);
                    } else {
                        if (currentWP.getFlipSNPAlleles()[d]) {
                            outzscores.append(sepStr).append(-datasetZScores[d]);
                        } else {
                            outzscores.append(sepStr).append(datasetZScores[d]);
                        }

                    }

                    if (datasetsSamples == null || datasetsSamples[d] == -9) {
                        outsamples.append(sepStr).append(nullstr);
                    } else {
                        outsamples.append(sepStr).append(datasetsSamples[d]);
                    }

                    if (probemeans == null || probemeans[d] == null) {
                        outmeans.append(sepStr).append(nullstr);
                    } else {
                        outmeans.append(sepStr).append(probemeans[d]);
                    }

                    if (probevars == null || probevars[d] == null) {
                        outvars.append(sepStr).append(nullstr);
                    } else {
                        outvars.append(sepStr).append(probevars[d]);
                    }

                    if (datasetfc == null || Double.isNaN(datasetfc[d])) {
                        outfc.append(sepStr).append(nullstr);
                    } else {
                        outfc.append(sepStr).append(datasetfc[d]);
                    }

                    if (datasetbeta == null || Double.isNaN(datasetbeta[d])) {
                        outbeta.append(sepStr).append(nullstr);
                    } else {
                        if (currentWP.getFlipSNPAlleles()[d]) {
                            outbeta.append(sepStr).append((-datasetbeta[d])).append(" (").append(datasetbetase[d]).append(")");
                        } else {
                            outbeta.append(sepStr).append((datasetbeta[d])).append(" (").append(datasetbetase[d]).append(")");
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
                out.append(hugo);
                out.append(tabStr);
                out.append(outcorrs.toString());
                out.append(tabStr);
                out.append(finalbeta).append(" (").append(finalbetase).append(")");
                out.append(tabStr);
                out.append(outbeta.toString());
                out.append(tabStr);
                out.append(outfc.toString());

            }


            return out.toString();
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

    public String getPermutationDescription(WorkPackage[] workPackages, IntMatrix2D probeTranslation, TriTyperGeneticalGenomicsDataset[] gg, int maxCisDistance) {
        String nullstr = "-";
        char tabStr = '\t';

        StringBuilder out = new StringBuilder();

        out.append(pvalue);
        out.append(tabStr);

        String rsName = null;

        if (sid == -1 && pid == -1) {
            return null;
        } else {
            SNP[] snps = workPackages[sid].getSnps();

            for (int d = 0; d < snps.length; d++) {
                if (snps[d] != null) {
                    rsName = snps[d].getName();
                    break;
                }
            }

            String probe = null;
            for (int d = 0; d < snps.length; d++) {
                if (probeTranslation.get(d, pid) != -9) {
                    int probeId = probeTranslation.get(d, pid);
                    probe = gg[d].getExpressionData().getProbes()[probeId];
                    break;
                }
            }

            out.append(rsName);
            out.append(tabStr);


            if (probe == null) {
                out.append(nullstr);
                out.append(tabStr);
            } else {
                out.append(probe);
                out.append(tabStr);
            }
            
            String hugo = nullstr;
            for (int d = 0; d < gg.length; d++) {
                if (!Double.isNaN(correlations[d])) {
                    if (probeTranslation.get(d, pid) != -9) {
                        int probeId = probeTranslation.get(d, pid);
                        hugo = gg[0].getExpressionData().getAnnotation()[probeId];
                        if(hugo==null){
                            hugo = nullstr;
                        }
                    }
                }
            }

            out.append(hugo);
            out.append(tabStr);
            
            if (alleles == null) {
                System.err.println(rsName + " has null alleles..?'\n" + out.toString());
                return null;
            }

            out.append(BaseAnnot.toString(alleles[0])).append("/").append(BaseAnnot.toString(alleles[1]));
            out.append(tabStr);

            out.append(BaseAnnot.toString(alleleAssessed));
            out.append(tabStr);

            out.append(zscore);
            
            return out.toString();
        }
    }
}
