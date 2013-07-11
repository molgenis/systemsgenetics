/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.util;

import eqtlmappingpipeline.binarymeta.meta.MetaAnalyze;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.zip.DataFormatException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class SNPAlleleCheck extends MetaAnalyze {

    @Override
    public void analyze() throws IOException, DataFormatException, Exception {
        System.out.println("");
        System.out.println("Starting analysis!");

        String[] datasets = new String[m_settings.getDatasetnames().size()];
        for (int i = 0; i < m_settings.getDatasetnames().size(); i++) {
            datasets[i] = m_settings.getDatasetnames().get(i);
        }

        if (!m_settings.getOutput().endsWith("/")) {
            m_settings.setOutput(m_settings.getOutput() + "/MetaAnalysis/");
        }

        if (!Gpio.exists(m_settings.getOutput())) {
            Gpio.createDir(m_settings.getOutput());
        }

        String[] locations = new String[m_settings.getDatasetnames().size()];
        for (int i = 0; i < locations.length; i++) {
            locations[i] = m_settings.getDatasetlocations().get(i);
        }

        ds = new BinaryResultDataset[m_settings.getDatasetlocations().size()];

        pvaluedistribution = null;
        eQTLBuffer = null;
        finalEQTLBuffer = null;
        nrInFinalBuffer = 0;


        uniqueProbes = new HashSet<String>();
        uniqueSNPs = new HashSet<String>();

        int numDatasets = ds.length;
        probes = new ArrayList<String>();

        snps = new ArrayList<String>();
        snpChr = new ArrayList<Byte>();
        snpChrPos = new ArrayList<Integer>();

        nrTotalSamples = 0;
        String[] probeName = probeTranslation.getProbes();
        probes.addAll(Arrays.asList(probeName));

        initdatasets(locations, 0, -1);
        Descriptives.lookupSqrt(nrTotalSamples);


        ProgressBar pb = new ProgressBar(snps.size());
        for (int s = 0; s < snps.size(); s++) {

            Double[][] zscorePerDataset = new Double[ds.length][probes.size()];

            for (int d = 0; d < ds.length; d++) {
                BinaryResultSNP firstSNPPassingQC = null;
                Integer snpId = snpTranslation[d][s];
                if (snpId != null) {
                    BinaryResultSNP snpObject = ds[d].getSnps()[snpId];

                    long pointer = snpObject.getzScoreIndex();
                    long nextpointer = -1;

                    if (snpId + 1 < ds[d].getSnps().length) {
                        BinaryResultSNP snpObject2 = ds[d].getSnps()[snpId + 1];
                        nextpointer = snpObject2.getzScoreIndex();
                    }

                    Float[] zscores = ds[d].getMatrix().read(pointer, nextpointer, ds[d].getNumProbes());
                    for (int p = 0; p < probes.size(); p++) {
                        Integer probeId = probeTranslationLookupTable[d][p];
                        byte probechr = probeTranslation.getProbeChr(p);
                        int probechrpos = probeTranslation.getProbeChrPos(p);
                        boolean testprobe = false;
                        if (probeId != null && testprobe) {
                            if (!zscores[probeId].isNaN()) {
                                int nrSamples = snpObject.getNumsamples();

                                double weight = Descriptives.getSqrt(nrSamples);
                                double zscore = zscores[probeId];

                                if (firstSNPPassingQC == null) {
                                    firstSNPPassingQC = snpObject;
                                } else {
                                    Boolean flipalleles = flipalleles(firstSNPPassingQC, snpObject);
                                    if (flipalleles == null) {
                                        System.err.println("ERROR! SNP alleles cannot be matched for snp\t" + snpObject.getName() + "\tin dataset\t" + datasets[d]);
                                        System.err.println("This SNP will be excluded from further research");
                                    } else if (flipalleles) {
                                        zscore = -zscore;
                                    }
                                }

                                zscorePerDataset[d][p] = zscore;

                            } else {
                                System.out.println("SNP\t"
                                        + snps.get(s) + ", "
                                        + snpChr.get(s) + ", "
                                        + snpChrPos.get(s) + "\thas not been tested for probe\t"
                                        + probes.get(p) + ", "
                                        + probeTranslation.getProbeChr(p) + ", "
                                        + probeTranslation.getProbeChrPos(p) + "\tin dataset\t" + datasets[d]);
                            }
                        }
                    }
                }
            }

            for (int d = 0; d < ds.length; d++) {
                for (int d2 = 0; d2 < ds.length; d2++) {
                    for (int p = 0; p < probes.size(); p++) {
                        if (d2 != d && zscorePerDataset[d2][p] != null && zscorePerDataset[d][p] != null) {
                            double z1 = zscorePerDataset[d][p];
                            Integer snpId = snpTranslation[d][s];
                            BinaryResultSNP snpObject1 = ds[d].getSnps()[snpId];
                            snpId = snpTranslation[d2][s];
                            BinaryResultSNP snpObject2 = ds[d2].getSnps()[snpId];

                            double z2 = zscorePerDataset[d2][p];
                            if ((z2 >= 4 && z1 <= -4) || (z1 >= 4 && z2 <= -4)) {
                                System.out.println("Discordant pair: ");
                                System.out.println(snpObject1.getName() + "\t"
                                        + BaseAnnot.toString(snpObject1.getAlleles()[0]) + "/" + BaseAnnot.toString(snpObject1.getAlleles()[1]) + "="
                                        + BaseAnnot.toString(snpObject1.getAssessedAllele()) + "\t"
                                        + snpObject1.getMinorAllele() + "/" + snpObject1.getMaf());
                                System.out.println(snpObject2.getName() + "\t"
                                        + BaseAnnot.toString(snpObject2.getAlleles()[0]) + "/" + BaseAnnot.toString(snpObject2.getAlleles()[1]) + "="
                                        + BaseAnnot.toString(snpObject2.getAssessedAllele()) + "\t"
                                        + snpObject1.getMinorAllele() + "/" + snpObject2.getMaf());
                                System.out.println("");
                            }
                        }
                    }
                }
            }
            pb.set(s);
        }
        pb.close();
    }

    // TODO: AT / GC SNPs??
    private Boolean flipalleles(BinaryResultSNP firstSNPPassingQC, BinaryResultSNP snpObject) {
        byte[] allelesfirst = firstSNPPassingQC.getAlleles();
        byte allelefirstassessed = firstSNPPassingQC.getAssessedAllele();

        byte[] allelessecond = snpObject.getAlleles();
        byte allelesecondassessed = snpObject.getAssessedAllele();

        int nridenticalalleles = 0;

        for (int i = 0; i < allelesfirst.length; i++) {
            byte allele1 = allelesfirst[i];
            for (int j = 0; j < allelessecond.length; j++) {
                if (allelessecond[j] == allele1) {
                    nridenticalalleles++;
                }
            }
        }

        if (nridenticalalleles == 2) {
            // alleles are identical. check if same allele was assessed...
            if (allelefirstassessed == allelesecondassessed) {
                return false;
            } else {
                return true;
            }
        } else {
            // try complement
            allelessecond = convertToComplementaryAlleles(allelessecond);
            allelesecondassessed = BaseAnnot.getComplement(allelesecondassessed);
            nridenticalalleles = 0;

            for (int i = 0; i < allelesfirst.length; i++) {
                byte allele1 = allelesfirst[i];
                for (int j = 0; j < allelessecond.length; j++) {
                    if (allelessecond[j] == allele1) {
                        nridenticalalleles++;
                    }
                }
            }

            if (nridenticalalleles == 2) {
                // alleles are identical. check if same allele was assessed...
                if (allelefirstassessed == allelesecondassessed) {
                    return false;
                } else {
                    return true;
                }
            }
        }
        return null;
    }

    private byte[] convertToComplementaryAlleles(byte[] allelesToCompare) {
        byte[] allelesComplementary = new byte[2];
        for (int a = 0; a < 2; a++) {
            allelesComplementary[a] = BaseAnnot.getComplement(allelesToCompare[a]);
        }
        return allelesComplementary;
    }
}
