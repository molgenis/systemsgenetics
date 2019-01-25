/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.cis;

import com.itextpdf.text.DocumentException;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLDotPlot;
import eqtlmappingpipeline.binarymeta.meta.MetaAnalyze;
import eqtlmappingpipeline.binarymeta.meta.MetaSettings;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import eqtlmappingpipeline.util.QTLFileSorter;
import java.io.EOFException;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultProbe;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CisAnalysis extends MetaAnalyze {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        // determine the sets of snps and probes to analyze
        // -- load snp + probe summary


        // create a matrix containing all the zscores

        // calculate on this matrix


        if (args.length < 1) {
            System.out.println("Usage: meta settings.xml");
        } else {
            String config = args[0];
            try {
                CisAnalysis c = new CisAnalysis(config);
                c.run();
            } catch (IOException ex) {
                Logger.getLogger(CisAnalysis.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public CisAnalysis(String args) throws IOException {
        m_settings = new MetaSettings();
        m_settings.parse(args, null, null);
        probeTranslation = new ProbeTranslation();
        probeTranslation.load(m_settings.getProbetranslationfile());
    }

    public void run() throws IOException {
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

        int permstart = 0;
        int permstop = m_settings.getNrPermutations() + 1;

        if (m_settings.getRunonlypermutation() > -1) {
            permstart = m_settings.getRunonlypermutation();
            permstop = m_settings.getRunonlypermutation() + m_settings.getNrPermutations();
        }

        for (int perm = permstart; perm < permstop; perm++) {
            ds = new BinaryResultDataset[m_settings.getDatasetlocations().size()];
            this.runCalculationRound(perm, locations, datasets, -1);
        }

        // sort the files.

        Gpio.createDir(m_settings.getOutput() + "Sorted/");

        for (int abs = 0; abs < 1; abs++) {
            for (int perm = 0; perm < m_settings.getNrPermutations() + 1; perm++) {
                QTLFileSorter sorter = new QTLFileSorter();
                String suffix = "eQTLs.txt.gz";


                //
                if (perm > 0) {
                    suffix = "PermutedEQTLsPermutationRound" + perm + ".txt.gz";

                }

                if (abs > 0) {
                    suffix += "-Absolute.txt.gz";
                }
                sorter.run(m_settings.getOutput() + suffix, null, QTLFileSorter.SORTBY.Z);
                System.out.println("Moving file.");
                System.out.println("SRC: " + m_settings.getOutput() + suffix + "_sorted.txt");
                System.out.println("DST: " + m_settings.getOutput() + "/Sorted/" + suffix);
                Gpio.moveFile(m_settings.getOutput() + suffix + "_sorted.txt", m_settings.getOutput() + "/Sorted/" + suffix);
            }
        }



        if (m_settings.getRunonlypermutation() == -1) {

            if (m_settings.getNrPermutations() > 0) {

                TextFile eQTLFile = new TextFile(m_settings.getOutput() + "/Sorted/eQTLs.txt.gz", TextFile.R);
                int nrEQTLs = eQTLFile.countLines() - 1;
                eQTLFile.close();

                FDR.calculateFDR(m_settings.getOutput() + "/Sorted/", m_settings.getNrPermutations(), nrEQTLs, m_settings.getFdrthreshold(), true, null, null, FDR.FDRMethod.ALL, true);
                EQTLDotPlot edp = new EQTLDotPlot();
                try {
                    edp.draw(m_settings.getOutput() + "/eQTLsFDR" + m_settings.getFdrthreshold() + ".txt", m_settings.getOutput() + "/DotPlot-FDR" + m_settings.getFdrthreshold() + ".pdf", EQTLDotPlot.Output.PDF); // "/eQTLsFDR" + fdrCutOff + ".txt", outputReportsDir + "/eQTLsFDR" + fdrCutOff + "DotPlot.png"
                } catch (DocumentException ex) {
                    Logger.getLogger(CisAnalysis.class.getName()).log(Level.SEVERE, null, ex);
                }
                edp = null;
            }
        }

    }

    @Override
    public void runCalculationRound(int perm, String[] locations, String[] datasets, int dToUse) throws IOException {


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

        super.initdatasets(locations, perm, dToUse);


//        System.gc();
//        System.gc();

        Descriptives.lookupSqrt(nrTotalSamples);

        System.out.println("Performing the meta-analysis now: ");
        System.out.println(nrTotalSamples + "\t total samples");


        // load SNP annotation

        // now we have loaded all essential info.. determine which snps will be tested against which probes..
        // and give each combination a unique number //
        int effectcounter = 0;
        // set up bidirectional map
        HashMap<Pair<Integer, Integer>, Integer> SNPProbeToEffectMap = new HashMap<Pair<Integer, Integer>, Integer>(14000000);
        HashMap<Integer, Pair<Integer, Integer>> EffectToSNPProbeMap = new HashMap<Integer, Pair<Integer, Integer>>(14000000);

        HashSet<String> snpsForWhichThereAreEffects = new HashSet<String>();

        System.out.println();

        boolean[] probespresentarray = new boolean[probes.size()];
        for (int p = 0; p < probes.size(); p++) {
            boolean probeIsPresent = false;
            for (int d = 0; d < ds.length; d++) {
                if (probeTranslationLookupTable[d][p] != null) {
                    probeIsPresent = true;
                }
            }
            probespresentarray[p] = probeIsPresent;
        }

        ProgressBar pb = new ProgressBar(snps.size(), "Building indexes.. please wait. ");
        HashSet<String> probesWithEffects = new HashSet<String>();
        boolean[] snpsHavingEffects = new boolean[snps.size()];
        boolean[] probeActuallyTested = new boolean[probes.size()];
        int metaprobenr = probes.size();
        int metasnpnr = snps.size();

        // determine which effects to test...
        for (int snp = 0; snp < metasnpnr; snp++) {
            byte chr = snpChr.get(snp);
            Integer pos = snpChrPos.get(snp);
            for (int probe = 0; probe < metaprobenr; probe++) {
                if (probespresentarray[probe]) {
                    byte probeChr = probeTranslation.getProbeChr(probe);
                    if (probeChr == chr) {
                        Integer probeChrPos = probeTranslation.getProbeChrPos(probe);
                        if (Math.abs(pos - probeChrPos) <= m_settings.getCisdistance()) {

                            probeActuallyTested[probe] = true;
                            snpsHavingEffects[snp] = true;


                            Pair<Integer, Integer> snpprobepair = new Pair<Integer, Integer>(snp, probe);
                            if (SNPProbeToEffectMap.get(snpprobepair) == null) {
                                SNPProbeToEffectMap.put(snpprobepair, effectcounter);
                                EffectToSNPProbeMap.put(effectcounter, snpprobepair);
                                snpsForWhichThereAreEffects.add(snps.get(snp));
                                probesWithEffects.add(probes.get(probe));

                                effectcounter++;
                            }
                        }
                    }
                }
            }
            pb.iterate();

        }
        pb.close();
        TextFile outfx = new TextFile(m_settings.getOutput() + "NumberOfEffects.txt", TextFile.W);
        System.out.println(effectcounter + " effects detected..");
        outfx.writeln(effectcounter + " effects detected..");

        outfx.writeln(snpsForWhichThereAreEffects.size() + " snps");
        System.out.println(snpsForWhichThereAreEffects.size() + " snps");

        outfx.writeln(probesWithEffects.size() + " probes");
        System.out.println(probesWithEffects.size() + " probes");
        outfx.close();

        snpsForWhichThereAreEffects = null;
        probesWithEffects = null;

        System.out.println("In total, number of effects: " + effectcounter);
        System.out.println("Now building zscore matrix");
        Float[][] allZScores = new Float[ds.length][effectcounter];

        System.out.println("Building reverse lookup table");

        // this converts a dataset snp number to a meta-snp id
        // build reverse SNP lookup table..
        Integer[][] reverseSNPLookupTable = new Integer[ds.length][];
        for (int d = 0; d < ds.length; d++) {
            reverseSNPLookupTable[d] = new Integer[ds[d].getSnps().length];
            for (int snp = 0; snp < snps.size(); snp++) {
                BinaryResultSNP snpObj = ds[d].getStringToSNP().get(snps.get(snp));
                if (snpObj != null) {
                    Integer id = snpObj.getId();
                    reverseSNPLookupTable[d][id] = snp;
                }
            }
        }

        // now load in the FCKING zscores...
        System.out.println("Loading zscores.. ");

        for (int d = 0; d < ds.length; d++) {
            // open binary file
            String filename = m_settings.getOutput() + m_settings.getDatasetnames().get(d) + "-eQTLs.dat";
            if (perm > 0) {
                filename = m_settings.getOutput() + m_settings.getDatasetnames().get(d) + "-PermutedDataPermutationRound-" + perm + ".dat";;
            }

            System.out.println("Loading file: " + filename);
            BinaryFile bf = new BinaryFile(filename, BinaryFile.R);

            boolean eof = false;
            String datasetannotation = m_settings.getDatasetannotations().get(d);
            int effectctr = 0;
            int lnctr = 0;
            while (!eof) {
                try {
                    Integer snpId = bf.readInt();
                    Integer probeId = bf.readInt();
                    Float zscore = bf.readFloat();


                    Integer metaSNPId = reverseSNPLookupTable[d][snpId];

                    BinaryResultProbe probeObj = ds[d].getProbes()[probeId];
                    Integer metaProbeId = probeTranslation.getProbeId(datasetannotation + probeObj.getName());

                    if (metaSNPId != null && metaProbeId != null) {
                        Integer effect = SNPProbeToEffectMap.get(new Pair<Integer, Integer>(metaSNPId, metaProbeId));
                        if (effect != null) {
                            allZScores[d][effect] = zscore;
                            effectctr++;
                        }
                        if (d == 0 && snps.get(metaSNPId).equals("rs4475691")) {
                            System.out.println(effect + "\t" + zscore + "\t" + snps.get(metaSNPId));
                        }
                    }

                } catch (EOFException e) {
                    eof = true;
                }
                lnctr++;
            }

            System.out.println("");
            System.out.println("Done");
            bf.close();
        }

        // all zscores loaded.. now perform the meta-analysis
        String outFileName = "eQTLs.txt.gz";
        if (perm > 0) {
            outFileName = "PermutedEQTLsPermutationRound" + perm + ".txt.gz";
        }
        System.out.println("Writing file: " + m_settings.getOutput() + outFileName);
        TextFile eQTLOutput = new TextFile(m_settings.getOutput() + outFileName, TextFile.W);
        TextFile eQTLOutputAbs = new TextFile(m_settings.getOutput() + outFileName + "-Absolute.txt.gz", TextFile.W);
        eQTLOutput.writeln(QTLTextFile.header);
        eQTLOutputAbs.writeln(QTLTextFile.header);


        pb = new ProgressBar(effectcounter, "Now performing meta-analysis...");

        for (int effect = 0; effect < effectcounter; effect++) {
            int nrDSHavingZScore = 0;
            int finalNrOfSamples = 0;
            double zsum = 0;
            double zsumabs = 0;

            Pair<Integer, Integer> snpProbePair = EffectToSNPProbeMap.get(effect);
            Integer snpId = snpProbePair.getLeft();
            Integer probeId = snpProbePair.getRight();

            BinaryResultSNP firstSNPToPassQC = null;


            BinaryResultSNP[] snpsPerDataset = new BinaryResultSNP[ds.length];
            // flip zscores
            for (int d = 0; d < ds.length; d++) {
                Float dsZScore = allZScores[d][effect];
                Integer dsSNPId = snpTranslation[d][snpId];

                if (dsZScore != null && dsSNPId != null) {

                    BinaryResultSNP dsSNPObj = ds[d].getSnps()[dsSNPId];

                    int nrsamples = dsSNPObj.getNumsamples();
                    double weight = Descriptives.getSqrt(nrsamples);
                    nrDSHavingZScore++;
                    finalNrOfSamples += nrsamples;
                    snpsPerDataset[d] = dsSNPObj;

                    if (firstSNPToPassQC == null) {
                        firstSNPToPassQC = dsSNPObj;
                    } else {
                        Boolean flipalleles = flipalleles(firstSNPToPassQC, dsSNPObj);
                        if (flipalleles == null) {
                            System.err.println("ERROR! SNP alleles cannot be matched for snp\t" + dsSNPObj.getName() + "\tin dataset\t" + d);
                            System.err.println("This SNP will be excluded from further research");
                        } else {
                            if (flipalleles) {
                                allZScores[d][effect] = -dsZScore;
                            }
                        }

                    }
                    Float finalZ = allZScores[d][effect];
                    zsum += (finalZ * weight);
                    zsumabs += Math.abs(finalZ * weight);
                }
            }

            // all zscores flipped.. run the weighted analysis..

            int nrDatasetsThatShouldHaveEffect = m_settings.getSnpDatasetPresenceThreshold();
            if (nrDatasetsThatShouldHaveEffect == 0) {
                nrDatasetsThatShouldHaveEffect = 1;
            }

            if (nrDSHavingZScore >= nrDatasetsThatShouldHaveEffect) {

                double zSumVal = zsum;

                double sqrtSample = Descriptives.getSqrt(finalNrOfSamples);
                double metaZScore = zSumVal / sqrtSample;
                double pValueOverall = Descriptives.convertZscoreToPvalue(metaZScore);

                double zSumValAbsolute = zsumabs;
                double zScoreAbs = zSumValAbsolute / sqrtSample;
                double pValueOverallAbs = Descriptives.convertZscoreToPvalue(zScoreAbs);


                for (int abs = 0; abs < 1; abs++) {
                    StringBuilder output = new StringBuilder();
                    if (abs == 0) {
                        output.append(pValueOverall).append("\t");
                    } else {
                        output.append(pValueOverallAbs).append("\t");
                    }

                    output.append(firstSNPToPassQC.getName()).append("\t");
                    output.append(snpChr.get(snpId)).append("\t");
                    output.append(snpChrPos.get(snpId)).append("\t");
                    output.append(probes.get(probeId)).append("\t");
                    Integer probeTranslationId = probeId;

                    output.append(probeTranslation.getProbeChr(probeTranslationId)).append("\t");
                    output.append(probeTranslation.getProbeChrPos(probeTranslationId)).append("\t");

                    output.append("cis").append("\t");

                    output.append(BaseAnnot.toString(firstSNPToPassQC.getAlleles()[0])).append("/").append(BaseAnnot.toString(firstSNPToPassQC.getAlleles()[1])).append("\t");

                    output.append(BaseAnnot.toString(firstSNPToPassQC.getAssessedAllele())).append("\t");

                    if (abs == 0) {
                        output.append(metaZScore).append("\t");
                    } else {
                        output.append(zScoreAbs).append("\t");
                    }

                    String[] datasetsPassingQC = new String[ds.length];
                    String[] datasetZScores = new String[ds.length];
                    String[] numSamplesPerDataset = new String[ds.length];
                    String[] emptyvalues = new String[ds.length];
                    for (int d = 0; d < ds.length; d++) {
                        emptyvalues[d] = "-";
                        if (allZScores[d][effect] == null) {
                            datasetZScores[d] = "-";
                            datasetsPassingQC[d] = "-";
                            numSamplesPerDataset[d] = "-";
                        } else {
                            datasetZScores[d] = "" + allZScores[d][effect];
                            datasetsPassingQC[d] = m_settings.getDatasetnames().get(d);
                            numSamplesPerDataset[d] = "" + snpsPerDataset[d].getNumsamples();
                        }
                    }

                    output.append(Strings.concat(datasetsPassingQC, Strings.comma)).append("\t");
                    output.append(Strings.concat(datasetZScores, Strings.comma)).append("\t");
                    output.append(Strings.concat(numSamplesPerDataset, Strings.comma)).append("\t");

                    // some empty fields
                    output.append(Strings.concat(emptyvalues, Strings.comma)).append("\t");
                    output.append(Strings.concat(emptyvalues, Strings.comma)).append("\t");

                    // hugo
                    output.append(probeTranslation.getProbeSymbol(probeTranslationId)).append("\t");

                    output.append(Strings.concat(emptyvalues, Strings.comma)).append("\t");
                    output.append(Strings.concat(emptyvalues, Strings.comma)).append("\t");
                    output.append(Strings.concat(emptyvalues, Strings.comma)).append("\t");
                    output.append(Strings.concat(emptyvalues, Strings.comma)).append("\t");

                    output.append(firstSNPToPassQC.getAssessedAllele()).append("\t");
                    if (abs == 0) {
                        eQTLOutput.writeln(output.toString());
                    } else {
                        eQTLOutputAbs.writeln(output.toString());
                    }

                }


            }
            pb.set(effect);

        }

        pb.close();
        eQTLOutput.close();
        eQTLOutputAbs.close();

        System.out.println("Done. Have a nice day.");

    }

    // TODO: AT / GC SNPs??
    public Boolean flipalleles(BinaryResultSNP firstSNPPassingQC, BinaryResultSNP snpObject) {
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

    public byte[] convertToComplementaryAlleles(byte[] allelesToCompare) {
        byte[] allelesComplementary = new byte[2];
        for (int a = 0; a < 2; a++) {
            allelesComplementary[a] = BaseAnnot.getComplement(allelesToCompare[a]);
        }
        return allelesComplementary;
    }
}
