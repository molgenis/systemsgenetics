/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.textmeta;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class FixedEffectMetaAnalysisTask implements Callable<String> {

    Pair<String, String> eqtl;
    String[] filesInDir;
    EQTL[][] allEQTLs;
    private final HashMap<String, Integer> index;
    private final int minimalNrDatasets;
    private final int minimalNrSamples;

    public FixedEffectMetaAnalysisTask(HashMap<String, Integer> index, Pair<String, String> eqtl, String[] files, EQTL[][] allEQTLs, int minimalNrDatasets, int minimalNrSamples) {
        this.eqtl = eqtl;
        this.index = index;
        this.filesInDir = files;
        this.allEQTLs = allEQTLs;
        this.minimalNrDatasets = minimalNrDatasets;
        this.minimalNrSamples = minimalNrSamples;
    }

    @Override
    public String call() throws Exception {
        ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
        String snp = eqtl.getLeft();
        String probe = eqtl.getRight();
        for (int f = 0; f < filesInDir.length; f++) {
            Integer id = index.get(f + "-" + snp + "-" + probe);
            if (id != null) {
                EQTL eq = allEQTLs[f][id];
                eqtls.add(eq);
            }

        }

        // meta-analyze the collected EQTLs
        double[] zscores = new double[eqtls.size()];
        int[] samplesize = new int[eqtls.size()];
        int nrSamples = 0;
        String[] datsets = new String[eqtls.size()];
        EQTL first = null;
        if (eqtls.size() >= minimalNrDatasets) {
            for (int q = 0; q < eqtls.size(); q++) {
                // if this is not the first eQTL
                // check whether we should flip the allele...
                EQTL e = eqtls.get(q);
                Boolean flipZ = false;
                if (q > 0) {
                    flipZ = BaseAnnot.flipalleles(first.getAlleles(), first.getAlleleAssessed(), e.getAlleles(), e.getAlleleAssessed());
                    if (flipZ == null) {
                        System.err.println("ERROR: alleles not compatible! " + e.getRsName() + "\t" + first.getAlleles() + "\t" + e.getAlleles());
                    }
                } else {
                    first = e;
                }

                // flip the allele if required
                if (flipZ != null) {
                    if (flipZ) {
                        zscores[q] = -e.getZscore();
                    } else {
                        zscores[q] = e.getZscore();
                    }
                }

                samplesize[q] = e.getDatasetsSamples()[0];
                nrSamples += samplesize[q];
            }

            if (nrSamples >= minimalNrSamples) {
                // calculate meta statistics
                double metaZ = ZScores.getWeightedZ(zscores, samplesize);
                double pvalue = ZScores.zToP(metaZ);

                // format
                // PValue  SNPName SNPChr  SNPChrPos       ProbeName       ProbeChr        ProbeCenterChrPos       CisTrans        SNPType AlleleAssessed  OverallZScore   DatasetsWhereSNPProbePairIsAvailableAndPassesQC DatasetsZScores DatasetsNrSamples       IncludedDatasetsMeanProbeExpression     IncludedDatasetsProbeExpressionVariance HGNCName        IncludedDatasetsCorrelationCoefficient

                String outStr =
                        pvalue + "\t"
                        + snp + "\t-\t-\t"
                        + probe + "\t-\t-\ttrans\t"
                        + eqtls.get(0).getAlleles() + "\t"
                        + eqtls.get(0).getAlleleAssessed() + "\t"
                        + metaZ + "\t"
                        + Strings.concat(datsets, Strings.comma) + "\t"
                        + Strings.concat(zscores, Strings.comma) + "\t" + Strings.concat(samplesize, Strings.comma) + "\t-\t-\t-\t-\t-\t-\t-";

                return outStr;
            } else {
                return null;
            }
        }
        return null;
    }
}
