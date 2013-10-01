/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CisEQTLProbeSNPLDCheck {

    double cr = 0.95;
    double hwe = 0.001;
    double maf = 0.01;

    //For reference only
//    public static void main(String[] args) {
//        try {
//
//            CisEQTLProbeSNPLDCheck p = new CisEQTLProbeSNPLDCheck();
//
//            int probecol = 0;
//            int probechrcol = 2;
//            int probechrposcol = 3;
//            String probeMapFile = "/target/gpfs2/gcc/groups/wijmenga/home/mjbonder/VCF/QTL-Filtering/IlluminaHT-12-V3_EQtlMappingFileAA2_4.txt";
//            
//            String snpProbeFile = "/target/gpfs2/gcc/groups/wijmenga/home/mjbonder/VCF/QTL-Filtering/All_eQTL-SNP-Probe_Combies_MJ-26-9-2013.txt";
//            String[] referenceDatasets = new String[]{
//                "/target/gpfs2/gcc/groups/wijmenga/home/mjbonder/VCF/1kg/",
//                "/target/gpfs2/gcc/groups/wijmenga/home/mjbonder/VCF/HapMap3/",
//                "/target/gpfs2/gcc/groups/wijmenga/home/mjbonder/VCF/GoNL_v4/"
//            };
//            String[] referenceDatasetNames = new String[]{
//                "1000Genomes",
//                "HapMapCEU",
//                "GoNL_V4"
//            };
//            String outputdirectory = "/target/gpfs2/gcc/groups/wijmenga/home/mjbonder/VCF/QTL-Filtering/Expression_";
//            p.runNew(probeMapFile, snpProbeFile, referenceDatasets, referenceDatasetNames, outputdirectory, 0.2, 0.2, probecol, probechrposcol, probechrcol);
//
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//    }
    public void runNew(String probeMapFile, String snpprobefile, String[] referenceDatasets, String[] referenceDatasetNames, String outputDirectory, double rSquaredThreshold, double dPrimeThreshold,
            int probecol, int probechrposcol, int probechrcol) throws IOException {
        ArrayList<Byte> probeChr = new ArrayList<Byte>();
        ArrayList<String> probeChrPos = new ArrayList<String>();
        ArrayList<String> probes = new ArrayList<String>();
        HashMap<String, Integer> probeToId = new HashMap<String, Integer>();



        TextFile tf = new TextFile(probeMapFile, TextFile.R);
        tf.readLine();
        String[] probePositions = tf.readLineElems(TextFile.tab);
        int id = 0;
        while (probePositions != null) {
            String probeStr = probePositions[probecol];
            String probechrposStr = probePositions[probechrposcol];
            String probeChrStr = probePositions[probechrcol];
            probeChr.add(ChrAnnotation.parseChr(probeChrStr));
            probeChrPos.add(probechrposStr);
            probes.add(probeStr);
            probeToId.put(probeStr, id);
            id++;
            probePositions = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(snpprobefile, TextFile.R);
        HashSet<Pair<String, String>> snpProbePairs = tf2.readAsPairs(0, 1);
        tf2.close();

        DetermineLD ldcalc = new DetermineLD();

        HashMap<Pair<String, String>, Boolean> qcResultsBoolean = new HashMap<Pair<String, String>, Boolean>();
        HashMap<Pair<String, String>, String[]> qcResultsStr = new HashMap<Pair<String, String>, String[]>();
        HashSet<Pair<String, String>> snpProbePairsTested = new HashSet<Pair<String, String>>();

        for (int d = 0; d < referenceDatasets.length; d++) {
            HashSet<Pair<String, String>> snpProbePairsTestedForThisReferenceDs = new HashSet<Pair<String, String>>();
            String referenceDatasetLocation = referenceDatasets[d];
            String referenceDatasetName = referenceDatasetNames[d];

            TriTyperGenotypeData reference = new TriTyperGenotypeData(referenceDatasetLocation);
            SNPLoader loader = reference.createSNPLoader();
            String[] snps = reference.getSNPs();
            TextFile output = new TextFile(outputDirectory + referenceDatasetName + "-QCOutput.txt", TextFile.W);
            output.writeln("SNPProbe\tFailsQC\tReason");

            HashSet<Pair<String, String>> snpProbePairsToTestForThisReferenceDataset = new HashSet<Pair<String, String>>();

            for (Pair<String, String> snpProbePair : snpProbePairs) {
                String snp = snpProbePair.getLeft();
                Integer snpId = reference.getSnpToSNPId().get(snp);

                if (snpId == null) {
                    snpProbePairsTested.add(snpProbePair);
                    snpProbePairsTestedForThisReferenceDs.add(snpProbePair);
                    output.writeln(snpProbePair.toString() + "\tUNKNOWN: SNP not present in " + referenceDatasetName);


                    String[] qcStrArr = qcResultsStr.get(snpProbePair);
                    if (qcStrArr == null) {
                        qcStrArr = new String[referenceDatasetNames.length];
                    }
                    qcStrArr[d] = "SNP not present in " + referenceDatasetName;
                    qcResultsStr.put(snpProbePair, qcStrArr);

                    Boolean b = qcResultsBoolean.get(snpProbePair);
                    if (b == null) {
                        qcResultsBoolean.put(snpProbePair, null);
                    }

                } else if (reference.getChr(snpId) == null || reference.getChr(snpId) < 1 || reference.getChr(snpId) > 23) {
                    // get the chromosome location
                    snpProbePairsTested.add(snpProbePair);
                    snpProbePairsTestedForThisReferenceDs.add(snpProbePair);
                    String[] qcStrArr = qcResultsStr.get(snpProbePair);
                    if (qcStrArr == null) {
                        qcStrArr = new String[referenceDatasetNames.length];
                    }
                    output.writeln(snpProbePair.toString() + "\tUNKNOWN: SNP maps to chromosome " + reference.getChr(snpId) + " in dataset " + referenceDatasetName);
                    qcStrArr[d] = "SNP maps to chromosome " + reference.getChr(snpId) + " in dataset " + referenceDatasetName;
                    qcResultsStr.put(snpProbePair, qcStrArr);
                    Boolean b = qcResultsBoolean.get(snpProbePair);
                    if (b == null) {
                        qcResultsBoolean.put(snpProbePair, null);
                    }
                } else {
                    snpProbePairsToTestForThisReferenceDataset.add(snpProbePair);
                }
            }

            for (byte chr = 1; chr < 23; chr++) {
                ArrayList<SortableSNP> snpsOnChr = new ArrayList<SortableSNP>();
                ArrayList<Integer> probesOnChr = new ArrayList<Integer>();

                int nrProbes = probes.size();
                for (int p = 0; p < nrProbes; p++) {
                    if (probeChr.get(p).equals(chr)) {
                        probesOnChr.add(p);
                    }
                }

                HashSet<String> snpsOnChrStrHash = new HashSet<String>();
                for (int s = 0; s < snps.length; s++) {
                    if (reference.getChr(s) == chr) {
                        snpsOnChrStrHash.add(snps[s]);
                        snpsOnChr.add(new SortableSNP(snps[s], s, chr, reference.getChrPos(s), SortableSNP.SORTBY.ID));
                    }
                }

                System.out.println(snpsOnChr.size() + "\tSNPs on Chr " + chr);
                System.out.println(probesOnChr.size() + "\tProbes on Chr " + chr);
                Collections.sort(snpsOnChr);

                HashMap<Integer, HashSet<Integer>> snpsInProbes = new HashMap<Integer, HashSet<Integer>>();

                for (int p = 0; p < probesOnChr.size(); p++) {
                    Integer probeId = probesOnChr.get(p);
                    String actualAnnotation = probeChrPos.get(probeId);
                    probePositions = actualAnnotation.split(":");
                    HashSet<Integer> snpsInProbe = new HashSet<Integer>();
                    for (String pos : probePositions) {
                        String[] probePositionElements = pos.split("-");
                        if (probePositionElements.length < 2) {
                            System.err.println("ERROR: " + pos + "\tis not parseable for probe: " + probes.get(probeId));
                        } else {
                            try {
                                Integer start = Integer.parseInt(probePositionElements[0]);
                                Integer stop = Integer.parseInt(probePositionElements[1]);
                                for (SortableSNP s : snpsOnChr) {
                                    int chrPos = s.chrpos;
                                    if (chrPos >= start && chrPos <= stop) {
                                        snpsInProbe.add(s.id);
                                    }
                                }
                            } catch (NumberFormatException e) {
                                System.err.println(pos + "\tis not parseable for probe: " + probes.get(probeId));
                            }
                        }
                    }
                    snpsInProbes.put(probeId, snpsInProbe);
                }

                for (Pair<String, String> snpProbePair : snpProbePairsToTestForThisReferenceDataset) {
                    if (snpProbePairsTestedForThisReferenceDs.contains(snpProbePair)) {
                        // do not test again
                        continue;
                    }

                    String snp = snpProbePair.getLeft();
                    Integer snpId = reference.getSnpToSNPId().get(snp);
                    if (snpsOnChrStrHash.contains(snp)) {
                        snpProbePairsTestedForThisReferenceDs.add(snpProbePair);
                        snpProbePairsTested.add(snpProbePair);
                        String probe = snpProbePair.getRight();
                        Integer probeId = probeToId.get(probe);
                        if (probeId == null) {
                            System.err.println("ERROR: no annotation loaded for probe: " + probe);
                            output.writeln(snpProbePair.toString() + "\tUNKNOWN: \tProbe annotation not loaded.");
                            String[] qcStrArr = qcResultsStr.get(snpProbePair);
                            if (qcStrArr == null) {
                                qcStrArr = new String[referenceDatasetNames.length];
                            }
                            qcStrArr[d] = "Probe annotation not loaded.";
                            qcResultsStr.put(snpProbePair, qcStrArr);
                            qcResultsBoolean.put(snpProbePair, null);
                        } else {
                            HashSet<Integer> probeSNPs = snpsInProbes.get(probeId);
                            if (probeSNPs.isEmpty()) {
                                output.writeln(snpProbePair.toString() + "\tFALSE\tNo SNPs underneath probe in " + referenceDatasetName);
                                String[] qcStrArr = qcResultsStr.get(snpProbePair);
                                if (qcStrArr == null) {
                                    qcStrArr = new String[referenceDatasetNames.length];
                                }
                                qcStrArr[d] = "No SNPs underneath probe in " + referenceDatasetName;
                                qcResultsStr.put(snpProbePair, qcStrArr);
                                qcResultsBoolean.put(snpProbePair, false);
                            } else {
                                SNP snpObj1 = reference.getSNPObject(snpId);
                                loader.loadGenotypes(snpObj1);
                                String qcStr = "";
                                Boolean failsQC = false;
                                if (snpObj1.getMAF() > maf && snpObj1.getHWEP() > hwe && snpObj1.getCR() > cr) {
                                    for (Integer otherSNP : probeSNPs) {
                                        SNP snpObj2 = reference.getSNPObject(otherSNP);
                                        loader.loadGenotypes(snpObj2);

                                        if (snpObj2.getMAF() > maf && snpObj2.getHWEP() > hwe && snpObj2.getCR() > cr) {
                                            Pair<Double, Double> ld = ldcalc.getLD(snpObj1, snpObj2, reference, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                            if (ld.getLeft() > rSquaredThreshold || ld.getRight() > dPrimeThreshold) {
                                                failsQC = true;
                                                if (qcStr.length() == 0) {
                                                    qcStr += snpObj2.getName() + "; TRUE (rsq: " + ld.getLeft() + ",D': " + ld.getRight() + ", MAF: " + snpObj2.getMAF() + ", CR: " + snpObj2.getCR() + ", HWEP: " + snpObj2.getHWEP() + ")";
                                                } else {
                                                    qcStr += "; " + snpObj2.getName() + "; TRUE (rsq: " + ld.getLeft() + ",D': " + ld.getRight() + ", MAF: " + snpObj2.getMAF() + ", CR: " + snpObj2.getCR() + ", HWEP: " + snpObj2.getHWEP() + ")";
                                                }
                                            } else {
                                                if (qcStr.length() == 0) {
                                                    qcStr += snpObj2.getName() + "; FALSE (rsq: " + ld.getLeft() + ",D': " + ld.getRight() + ", MAF: " + snpObj2.getMAF() + ", CR: " + snpObj2.getCR() + ", HWEP: " + snpObj2.getHWEP() + ")";
                                                } else {
                                                    qcStr += "; " + snpObj2.getName() + "; FALSE (rsq: " + ld.getLeft() + ",D': " + ld.getRight() + ", MAF: " + snpObj2.getMAF() + ", CR: " + snpObj2.getCR() + ", HWEP: " + snpObj2.getHWEP() + ")";
                                                }
                                            }
                                        } else {
                                            if (qcStr.length() == 0) {
                                                qcStr += snpObj2.getName() + "; UNKNOWN (MAF: " + snpObj2.getMAF() + ", CR: " + snpObj2.getCR() + ", HWEP: " + snpObj2.getHWEP() + ")";
                                            } else {
                                                qcStr += "; " + snpObj2.getName() + "; UNKNOWN (MAF: " + snpObj2.getMAF() + ", CR: " + snpObj2.getCR() + ", HWEP: " + snpObj2.getHWEP() + ")";
                                            }
                                        }
                                        snpObj2.clearGenotypes();
                                        snpObj2 = null;
                                    }
                                    output.writeln(snpProbePair.toString() + "\t" + failsQC + "\t" + qcStr);

                                } else {
                                    qcStr += "SNP does not pass QC thresholds in " + referenceDatasetName + "(MAF: " + snpObj1.getMAF() + ", CR: " + snpObj1.getCR() + ", HWEP: " + snpObj1.getHWEP() + ")";
                                    output.writeln(snpProbePair.toString() + qcStr);
                                    failsQC = null;
                                }
                                Boolean b = qcResultsBoolean.get(snpProbePair);

                                if (b == null && failsQC == null) {
                                    qcResultsBoolean.put(snpProbePair, null);
                                } else {
                                    boolean finalb = false;
                                    if (failsQC == null && b != null) {
                                        finalb = b.booleanValue();
                                    } else if (failsQC != null && b == null) {
                                        finalb = failsQC;
                                    } else {
                                        boolean bb = b.booleanValue();
                                        if (!failsQC && !bb) {
                                            finalb = false;
                                        } else {
                                            finalb = true;
                                        }
                                    }
                                    qcResultsBoolean.put(snpProbePair, finalb);
                                }

                                String[] qcStrArr = qcResultsStr.get(snpProbePair);
                                if (qcStrArr == null) {
                                    qcStrArr = new String[referenceDatasetNames.length];
                                }
                                qcStrArr[d] = qcStr;
                                qcResultsStr.put(snpProbePair, qcStrArr);
                                snpObj1.clearGenotypes();
                                snpObj1 = null;
                            }

                        }
                    }
                }
            }

            output.close();
            loader.close();

        }

        TextFile output2 = new TextFile(outputDirectory + "Merged-CombinationsNotTestedByQC.txt", TextFile.W);
        for (Pair<String, String> p : snpProbePairs) {
            if (!snpProbePairsTested.contains(p)) {
                output2.writeln(p.toString());
            }
        }
        output2.close();


        TextFile output3 = new TextFile(outputDirectory + "Merged-QCResults.txt", TextFile.W);

        String header = "SNP-Probe\tFailsQC";
        for (String refDs : referenceDatasetNames) {
            header += "\tQCOutput-" + refDs;
        }

        output3.writeln(header);
        for (Pair<String, String> p : snpProbePairs) {
            if (!qcResultsStr.containsKey(p)) {
            } else {
                String[] qcStr = qcResultsStr.get(p);
                Boolean b = qcResultsBoolean.get(p);
                if (b == null) {
                    output3.writeln(p.toString() + "\tUNKNOWN\t" + Strings.concat(qcStr, Strings.tab));
                } else {
                    output3.writeln(p.toString() + "\t" + b + "\t" + Strings.concat(qcStr, Strings.tab));
                }

            }
        }
        output3.close();


    }

    public void determineSNPProbePairsWhichMayHaveFalsePositiveEffect(String probeTranslation, String reference, String outdir) throws IOException {

        ProbeTranslation pb = new ProbeTranslation();
        pb.load(probeTranslation);


        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(reference);

        String[] snps = ds.getSNPs();

        TextFile output = new TextFile(outdir + "SNPProbeCombosWithPossibleHybArtifacts1Kg.txt", TextFile.W);
        for (byte chr = 1; chr < 23; chr++) {
            ArrayList<SortableSNP> snpsOnChr = new ArrayList<SortableSNP>();
            ArrayList<Integer> probesOnChr = new ArrayList<Integer>();

            int nrProbes = pb.getNumProbes();
            for (int p = 0; p < nrProbes; p++) {
                if (pb.getProbeChr(p) == chr) {
                    probesOnChr.add(p);
                }
            }

            for (int s = 0; s < snps.length; s++) {
                if (ds.getChr(s) == chr) {
                    snpsOnChr.add(new SortableSNP(snps[s], s, chr, ds.getChrPos(s), SortableSNP.SORTBY.ID));
                }
            }

            System.out.println(snpsOnChr.size() + "\tSNPs on Chr " + chr);
            System.out.println(probesOnChr.size() + "\tProbes on Chr " + chr);
            Collections.sort(snpsOnChr);

            HashMap<Integer, HashSet<Integer>> snpsInProbes = new HashMap<Integer, HashSet<Integer>>();

            for (int p = 0; p < probesOnChr.size(); p++) {
                String actualAnnotation = pb.getActualMappingPosition(probesOnChr.get(p));
                String[] elems = actualAnnotation.split(":");

                for (String pos : elems) {
                    String[] elems2 = pos.split("-");
                    Integer start = Integer.parseInt(elems2[0]);
                    Integer stop = Integer.parseInt(elems2[1]);

                    for (SortableSNP s : snpsOnChr) {
                        int chrPos = s.chrpos;
                        if (chrPos >= start && chrPos <= stop) {
                            HashSet<Integer> snpsInProbe = snpsInProbes.get(probesOnChr.get(p));
                            if (snpsInProbe == null) {
                                snpsInProbe = new HashSet<Integer>();
                            }
                            snpsInProbe.add(s.id);
                            snpsInProbes.put(probesOnChr.get(p), snpsInProbe);
                        }
                    }
                }
            }
            System.out.println(snpsInProbes.size() + "\tprobes with SNPs");

            SNPLoader loader = ds.createSNPLoader();
            DetermineLD ldcalc = new DetermineLD();

            ProgressBar progress = new ProgressBar(snpsOnChr.size(), "Testing chr: " + chr);

            HashSet<Integer> snpsNotPassingQC = new HashSet<Integer>();
            for (SortableSNP s : snpsOnChr) {
                // for each probe within 1Mb, check
                if (!snpsNotPassingQC.contains(s.id)) {
                    SNP snpObj1 = ds.getSNPObject(s.id);
                    loader.loadGenotypes(snpObj1);

                    if (snpObj1.getMAF() > maf && snpObj1.getHWEP() > hwe && snpObj1.getCR() > cr) {
                        for (int p = 0; p < probesOnChr.size(); p++) {
                            HashSet<Integer> snpsInProbe = snpsInProbes.get(probesOnChr.get(p));
                            if (snpsInProbe != null) {
                                boolean failsQC = false;
                                String qcStr = "";
                                for (Integer snpInProbe : snpsInProbe) {
                                    if (!snpsNotPassingQC.contains(snpInProbe)) {
                                        SNP snpObj2 = ds.getSNPObject(snpInProbe);
                                        loader.loadGenotypes(snpObj2);
                                        if (snpObj2.getMAF() > maf && snpObj2.getHWEP() > hwe && snpObj2.getCR() > cr) {
                                            Pair<Double, Double> ld = ldcalc.getLD(snpObj1, snpObj2, ds, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                            if (ld.getLeft() > 0.2 || ld.getRight() > 0.2) {
                                                failsQC = true;
                                                qcStr += "\t" + snpObj2.getName() + " (" + ld.getLeft() + ", " + ld.getRight() + ")";
                                            }
                                        } else {
                                            snpsNotPassingQC.add(snpInProbe);
                                        }
                                        snpObj2.clearGenotypes();
                                    }
                                }
                                if (failsQC) {
                                    String outputStr = s.name + "\t" + pb.getProbes()[probesOnChr.get(p)] + qcStr;
                                    output.writeln(outputStr);
                                }
                            }
                        }
                    } else {
                        snpsNotPassingQC.add(s.id);
                    }
                    snpObj1.clearGenotypes();
                }
                progress.iterate();
            }
            progress.close();
        }
        output.close();
    }
}
