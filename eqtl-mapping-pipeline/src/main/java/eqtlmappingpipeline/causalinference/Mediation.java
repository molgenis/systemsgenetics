/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.causalinference;

import java.io.IOException;
import java.util.Iterator;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Regression;

/**
 *
 * @author harm-jan
 */
public class Mediation extends IVAnalysis {

    public Mediation(String xmlSettingsFile,
            String ingt, String inexp, String inexpplatform, String inexpannot,
            String gte, String out, int perm, String snpProbeCombinationList) throws IOException, Exception {
        super(xmlSettingsFile, ingt, inexp, inexpplatform, inexpannot, gte, out, perm, snpProbeCombinationList);
    }

    @Override
    public void run() throws IOException {
        for (int d = 0; d < m_gg.length; d++) {
            // now test all triples
            SNPLoader snpLoader = m_gg[d].getGenotypeData().createSNPLoader();
            int[] indWGA = m_gg[d].getExpressionToGenotypeIdArray();

            for (int perm = 0; perm < m_settings.nrPermutationsFDR + 1; perm++) {
                String outfile = null;
                if (perm == 0) {
                    outfile = outDir + m_gg[d].getSettings().name + "_IVAnalysis-RealData.txt";
                } else {
                    outfile = outDir + m_gg[d].getSettings().name + "_IVAnalysis-PermutationRound-" + perm + ".txt";
                    m_gg[d].permuteSampleLables(m_settings.r);
                }
                TextFile out = new TextFile(outfile, TextFile.W);
                Iterator<Triple<String, String, String>> it = snpProbeCombos.iterator();
                Triple<String, String, String> next = it.next();
                ProgressBar pb = new ProgressBar(snpProbeCombos.size(), "Running Mediation Analysis - Permutation " + perm);

                out.writeln("SNP\tSNP Chr\tSNP ChrPos\t"
                        + "Alleles\tDirectionAllele\t"
                        + "N\t"
                        + "CisArrayAddress\tCisProbe Chr\tCisProbe ChrPos\t"
                        + "CisGeneName\t"
                        + "TransArrayAddress\tTransProbe Chr\tTransProbe ChrPos\t"
                        + "TransGeneName\t"
                        + "CisTrans-Correlation\t"
                        + "Cis-eQTL-Beta\t"
                        + "Cis-eQTL-SE\t"
                        + "CisTrans-Beta\t"
                        + "CisTrans-SE\t"
                        + "Trans-eQTL-Beta\t"
                        + "Trans-eQTL-SE\t"
                        + "CisTrans-Residual-Correlation\t"
                        + "CisTrans-Residual-Beta\t"
                        + "CisTrans-Residual-SE\t"
                        + "Trans-eQTL-Residual-Beta\t"
                        + "Trans-eQTL-Residual-SE\t"
                        + "Beta-Ratio");

                while (next != null) {
                    String snp = next.getLeft();
                    String cisprobe = next.getMiddle();
                    String transprobe = next.getRight();

                    Integer snpId = m_gg[d].getGenotypeData().getSnpToSNPId().get(snp);
                    Integer cisProbeId = m_gg[d].getExpressionData().getProbeToId().get(cisprobe);
                    Integer transProbeId = m_gg[d].getExpressionData().getProbeToId().get(transprobe);

                    if (snpId == -9 || cisProbeId == null || transProbeId == null) {
//                        out.writeln(snp + "\t" + snpId + "\t" + cisprobe + "\t" + cisProbeId + "\t" + null + "\t" + transprobe + "\t" + transProbeId + "\t" + null + "\t" + null + "\t" + null + "\t" + null + "\t" + null + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                    } else {

                        SNP snpObj = m_gg[d].getGenotypeData().getSNPObject(snpId);
                        snpLoader.loadGenotypes(snpObj);
                        if (snpLoader.hasDosageInformation()) {
                            snpLoader.loadDosage(snpObj);
                        }
                        double[] origCisVals = m_gg[d].getExpressionData().getMatrix()[cisProbeId];
                        double[] origTransVals = m_gg[d].getExpressionData().getMatrix()[transProbeId];

                        int calledGenotypes = 0;
                        for (int i = 0; i < m_gg[d].getExpressionData().getIndividuals().length; i++) {
                            int genotypeId = indWGA[i];
                            short gt = snpObj.getGenotypes()[genotypeId];
                            if (genotypeId > -1 && gt > -1) {
                                calledGenotypes++;
                            }
                        }

                        double[] genotypes = new double[calledGenotypes];
                        double[] cisvals = new double[calledGenotypes];
                        double[] transvals = new double[calledGenotypes];

                        calledGenotypes = 0;
                        for (int i = 0; i < m_gg[d].getExpressionData().getIndividuals().length; i++) {
                            int genotypeId = indWGA[i];
                            short gt = snpObj.getGenotypes()[genotypeId];
                            if (genotypeId > -1 && gt > -1) {
                                genotypes[calledGenotypes] = snpObj.getDosageValues()[genotypeId];
                                cisvals[calledGenotypes] = origCisVals[i];
                                transvals[calledGenotypes] = origTransVals[i];
                                calledGenotypes++;
                            }
                        }

                        // normalize genotype and cis + trans to get beta's equal to the correlation coefficient
                        genotypes = normalize(genotypes);
                        cisvals = normalize(cisvals);
                        transvals = normalize(transvals);

                        double corrCisTrans = JSci.maths.ArrayMath.correlation(cisvals, transvals); // for code validation
                        double[] cisTransRCs = Regression.getLinearRegressionCoefficients(cisvals, transvals); // returns beta, alpha, se, t
                        double[] snpCisRCs = Regression.getLinearRegressionCoefficients(genotypes, cisvals); // returns beta, alpha, se, t
                        double[] snpTransRCs = Regression.getLinearRegressionCoefficients(genotypes, transvals);

                        // remove correlation between cis and trans probe
//                        double[] resCis = new double[cisvals.length];
                        double[] resTransVals = new double[cisvals.length];
                        for (int i = 0; i < resTransVals.length; i++) {
//                            resCis[i] = cisvals[i] - snpCisRCs[0] * genotypes[i];
                            resTransVals[i] = transvals[i] - cisTransRCs[0] * cisvals[i];
                        }

                        resTransVals = normalize(resTransVals);

                        double[] cisResTransRCs = Regression.getLinearRegressionCoefficients(cisvals, resTransVals); // returns beta, alpha, se, t
                        double[] snpResTransRCs = Regression.getLinearRegressionCoefficients(genotypes, resTransVals);

                        double rescorr = JSci.maths.ArrayMath.correlation(cisvals, resTransVals); // for code validation

                        out.writeln(snp
                                + "\t" + snpObj.getChr()
                                + "\t" + snpObj.getChrPos()
                                + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0]) + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                                + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0])
                                + "\t" + transvals.length
                                + "\t" + cisprobe
                                + "\t" + m_gg[d].getExpressionData().getChr()[cisProbeId]
                                + "\t" + m_gg[d].getExpressionData().getChrStart()[cisProbeId]
                                + ":" + m_gg[d].getExpressionData().getChrStop()[cisProbeId]
                                + "\t" + m_gg[d].getExpressionData().getAnnotation()[cisProbeId]
                                + "\t" + transprobe
                                + "\t" + m_gg[d].getExpressionData().getChr()[transProbeId]
                                + "\t" + m_gg[d].getExpressionData().getChrStart()[transProbeId]
                                + ":" + m_gg[d].getExpressionData().getChrStop()[transProbeId]
                                + "\t" + m_gg[d].getExpressionData().getAnnotation()[transProbeId]
                                + "\t" + corrCisTrans
                                + "\t" + snpCisRCs[0]
                                + "\t" + snpCisRCs[2]
                                + "\t" + cisTransRCs[0]
                                + "\t" + cisTransRCs[2]
                                + "\t" + snpTransRCs[0]
                                + "\t" + snpTransRCs[2]
                                + "\t" + rescorr
                                + "\t" + cisResTransRCs[0]
                                + "\t" + cisResTransRCs[2]
                                + "\t" + snpResTransRCs[0]
                                + "\t" + snpResTransRCs[2]
                                + "\t" + (snpResTransRCs[0] / snpTransRCs[0]));
                        snpObj.clearGenotypes();
                    }

                    if (it.hasNext()) {
                        next = it.next();
                    } else {
                        next = null;
                    }
                    pb.iterate();
                }
                pb.close();
                out.close();
            }
            snpLoader.close();
        }

    }

    private double[] normalize(double[] x) {
        double mean = Descriptives.mean(x);
        double sd = Math.sqrt(Descriptives.variance(x, mean));
        double[] xcorr = new double[x.length];
        for (int i = 0; i < xcorr.length; i++) {
            xcorr[i] = (x[i] - mean) / sd;
        }
        return xcorr;
    }
}
