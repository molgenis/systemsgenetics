/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util.HiC;

import static eqtlmappingpipeline.util.HiC.HiCQTLAnnotatorSnpBased.getNumericResolution;
import static eqtlmappingpipeline.util.HiC.HiCQTLAnnotatorSnpBased.readInQtlInformation;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.chrContacts.DesiredChrContact;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.MinimalEQTL;

/**
 *
 * @author MarcJan
 */
class HiCQTLAnnotatorBlockbased {

    public static void main(String[] args) throws IOException {
        String folderHighC = "G:\\Contacts\\";
        String resolution = "1kb"; //5kb / 1kb
        String qualityCutOff = "E30"; //0 or E30
        String normMethod = null;
//        String normMethod = "KRnorm"; //null / "KRnorm" / "SQRTVCnorm" / "VCnorm"
        double minValueQuality = 0.0;
        boolean permutationFile = true;
        boolean alternativePermutation = false;
        String probeMap = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Illumina450K_MQtlMappingFile_MJB.txt";
        String snpMap = "D:\\Werk\\UMCGundefinedUMCG\\Projects\\LL-DeepBBMRI_Methylation450K\\meQTLs\\SNPMappings\\SNPMappings.txt";

        for (int i = 35; i < 57; ++i) {
//            String QTLfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered.txt";
//            String QTLfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\eQTLProbesFDR0.05-1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-ProbeLevel.txt_LdQtlFilterd.txt";
            String QTLfile = "D:\\Werk\\UMCGundefinedUMCG\\Projects\\LL-DeepBBMRI_Methylation450K\\meQTLs\\trans-QTLs\\PermutedEQTLsPermutationRound" + i + ".head.txt";
            String proxyfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\proxiesMeQTLSnps2.txt";
//            String proxyfile = null;
            String QTLoutfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\HiC_Annot\\New_HiC_1Kb_E30\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_Perm" + i + "_HiC_LD_1kb_E30_annotated.txt";
//            String QTLoutfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\HiC\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_HiC_1kb_E30_annotated.txt";

            addAnnotationToQTLOutput(
                    QTLfile,
                    proxyfile,
                    folderHighC,
                    resolution,
                    qualityCutOff,
                    normMethod,
                    minValueQuality,
                    alternativePermutation,
                    permutationFile,
                    probeMap,
                    snpMap,
                    QTLoutfile);
        }

    }

    private static void addAnnotationToQTLOutput(String in, String inProxies, String folderHighC, String resolution, String qualityCutOff, String normMethod, double minValue, boolean alternativePermutation, boolean permutationFile, String probeMap, String snpMap, String out) throws IOException {

        HashMap<String, ArrayList<DesiredChrContact>> qtls = readInQtlTransformBlocks(in, inProxies, probeMap, snpMap, permutationFile, resolution, alternativePermutation);

        ProgressBar pb = new ProgressBar(qtls.size(), "Checking for contacts for: " + qtls.size() + " Chromosome combinations");

        TextFile outWriter = new TextFile(out, TextFile.W);

        for (Entry<String, ArrayList<DesiredChrContact>> contactsToCheck : qtls.entrySet()) {
            Collections.sort(contactsToCheck.getValue());

            String[] chrs = contactsToCheck.getKey().split("-");

            String ChrSmaller = chrs[0];
            String ChrLarger = chrs[1];

            String baseName;
            String fileToReads;
            boolean intra = false;

            if (ChrSmaller.equals(ChrLarger)) {
                baseName = folderHighC + "\\GM12878_combined_intrachromosomal\\" + resolution + "_resolution_intrachromosomal\\chr" + ChrSmaller + "\\MAPQG" + qualityCutOff;
                fileToReads = baseName + "\\chr" + ChrSmaller + "_" + resolution + ".RAWobserved";
                intra = true;
//                continue;
            } else {
                baseName = folderHighC + "\\GM12878_combined_interchromosomal\\" + resolution + "_resolution_interchromosomal\\chr" + ChrSmaller + "_chr" + ChrLarger + "\\MAPQG" + qualityCutOff;
                fileToReads = baseName + "\\chr" + ChrSmaller + "_" + ChrLarger + "_" + resolution + ".RAWobserved";
            }

//            if (normMethod == null) {
            processRawContactInformation(fileToReads, minValue, contactsToCheck.getValue(), intra);
//            } else {
//                if (intra) {
//                    processNormalizedIntraContactInformation(fileToReads, baseName, normMethod, ChrSmaller, contactsToCheck.getValue(), resolution, minValue, outWriter);
//                } else {
//                    processNormalizedInterContactInformation(fileToReads, baseName, normMethod, ChrSmaller, ChrLarger, contactsToCheck.getValue(), resolution, minValue, outWriter);
//                }
//            }
            printOutContacts(contactsToCheck.getValue(), outWriter);
            pb.iterate();
        }
        pb.close();
        outWriter.close();
    }

    //For example, here is a line from the 5kb chr1 MAPQGE30 raw observed contact matrix (GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
//40000000 40100000 59.0
    private static void processRawContactInformation(String fileToRead, double minValue, ArrayList<DesiredChrContact> contactsToCheck, boolean intra) throws IOException {

        //Check if sorted version is available
        //If not make sorted available.
        if (!Gpio.exists(fileToRead + ".sorted")) {
            if (intra) {
                umcg.genetica.io.chrContacts.SortIntraChrContacts.readNonSortedWriteSorted(fileToRead, fileToRead + ".sorted");
            } else {
                umcg.genetica.io.chrContacts.SortInterChrContacts.readNonSortedWriteSorted(fileToRead, fileToRead + ".sorted");
            }

        }

        int numberToBeMatched = 0;

        LineIterator it = FileUtils.lineIterator(new File(fileToRead + ".sorted"), "UTF-8");

        try {
            while (it.hasNext()) {
                String[] parts = StringUtils.split(it.nextLine(), '\t');

                int posChr1 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[0]);
                int posChr2 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[1]);

                while (numberToBeMatched < contactsToCheck.size()) {
                    if (posChr1 < contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        break;
                    } else if (posChr1 == contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        if (posChr2 < contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            break;
                        }
                        if (posChr2 == contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            double contact = org.apache.commons.lang.math.NumberUtils.createDouble(parts[2]);
                            if (contact >= minValue) {
                                contactsToCheck.get(numberToBeMatched).setContact();
                                numberToBeMatched++;
                            } else {
                                numberToBeMatched++;
                            }
                        } else if (posChr2 > contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            numberToBeMatched++;
                        }
                    } else if (posChr1 > contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        numberToBeMatched++;
                    }
                }
            }
        } finally {
            LineIterator.closeQuietly(it);
        }

    }

    //For example, here is a line from the 5kb chr1 MAPQGE30 raw observed contact matrix (GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
    //40000000 40100000 59.0
    //To normalize this entry using the KR normalization vector, one would divide 59.0 by the 8001st line ((40000000/5000)+1=8001) and the 8021st line ((40100000/5000)+1=8021)
    //of GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.KRnorm. The 8001st line of the KR norm file is 1.2988778370674694;
    //The 8021st line of the KR norm file is 1.6080499717941548. So the corresponding KR normalized entry for the entry above is 59.0/(1.2988778370674694*1.6080499717941548)
    //or 28.24776973966101.
    //If the KR normalization vector file is empty or all NaNs, then the KR algorithm didn’t converge on that particular matrix (likely due to sparsity of the matrix).
    private static void processNormalizedInterContactInformation(String fileToRead, String baseName, String normMethod, String chrSmaller, String chrLarger, ArrayList<DesiredChrContact> contactsToCheck, String resolution, double minValue, TextFile outWriter) throws IOException {

        //ReadIn normalization chr1
        TextFile inputNormChr1 = new TextFile(baseName + "\\chr" + chrSmaller + "_" + resolution + "." + normMethod, TextFile.R);
        ArrayList<String> normFactorSmallerChr = inputNormChr1.readAsArrayList();
        inputNormChr1.close();

//        System.out.println("Done reading norm factor 1");
        //ReadIn normalization chr2
        TextFile inputNormChr2 = new TextFile(baseName + "\\chr" + chrLarger + "_" + resolution + "." + normMethod, TextFile.R);
        ArrayList<String> normFactorLargerChr = inputNormChr2.readAsArrayList();
        inputNormChr2.close();

//        System.out.println("Done reading norm factor 2");
        if (!Gpio.exists(fileToRead + ".sorted")) {
            umcg.genetica.io.chrContacts.SortInterChrContacts.readNonSortedWriteSorted(fileToRead, fileToRead + ".sorted");
        }

        int numberToBeMatched = 0;

        LineIterator it = FileUtils.lineIterator(new File(fileToRead + ".sorted"), "UTF-8");

        try {
            while (it.hasNext()) {
                String[] parts = StringUtils.split(it.nextLine(), '\t');

                int posChr1 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[0]);
                int posChr2 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[1]);

                while (numberToBeMatched < contactsToCheck.size()) {
                    if (posChr1 < contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        break;
                    } else if (posChr1 == contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        if (posChr2 < contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            break;
                        }
                        if (posChr2 == contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            if (((posChr1 / getNumericResolution(resolution)) + 1) > normFactorSmallerChr.size()) {
                                System.out.println(baseName);
                                System.out.println("Smaller");
                                System.out.println((posChr1 / getNumericResolution(resolution) + 1));
                                System.out.println(normFactorSmallerChr.size());
                                System.exit(-1);
                            }
                            if (((posChr2 / getNumericResolution(resolution)) + 1) > normFactorLargerChr.size()) {
                                System.out.println(baseName);
                                System.out.println("Larger");
                                System.out.println((posChr2 / getNumericResolution(resolution)) + 1);
                                System.out.println(normFactorLargerChr.size());
                                System.exit(-1);
                            }
                            String factor1Base = normFactorSmallerChr.get((posChr1 / getNumericResolution(resolution)) + 1);
                            String factor2Base = normFactorLargerChr.get((posChr2 / getNumericResolution(resolution)) + 1);

                            double factor1 = 1.0;
                            double factor2 = 1.0;

                            if (NumberUtils.isNumber(factor1Base) && NumberUtils.isNumber(factor2Base)) {
                                factor1 = Double.parseDouble(factor1Base);
                                factor2 = Double.parseDouble(factor2Base);
                            } else if (NumberUtils.isNumber(factor1Base)) {
                                factor1 = Double.parseDouble(factor1Base);
                                System.out.println("Error in files.");
                                System.out.println("Base 2 is reset to 1");
                            } else if (NumberUtils.isNumber(factor2Base)) {
                                factor2 = Double.parseDouble(factor2Base);
                                System.out.println("Error in files.");
                                System.out.println("Base 1 is reset to 1");
                            }

                            double contact = org.apache.commons.lang.math.NumberUtils.createDouble(parts[2]) / (factor1 * factor2);
                            if (contact >= minValue) {
                                outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\tContact\t" + contact + "\t" + org.apache.commons.lang.math.NumberUtils.createDouble(parts[2]));
                                numberToBeMatched++;
                            } else {
                                outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-\t-");
                                numberToBeMatched++;
                            }

                        } else if (posChr2 > contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-\t-");
                            numberToBeMatched++;
                        }
                    } else if (posChr1 > contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-\t-");
                        numberToBeMatched++;
                    }
                }
            }
        } finally {
            LineIterator.closeQuietly(it);
        }

    }

    private static void processNormalizedIntraContactInformation(String fileToRead, String baseName, String normMethod, String chrSmaller, ArrayList<DesiredChrContact> contactsToCheck, String resolution, double minValue, TextFile outWriter) throws IOException {

        //ReadIn normalization chr1
        TextFile inputNormChr1 = new TextFile(baseName + "\\chr" + chrSmaller + "_" + resolution + "." + normMethod, TextFile.R);
        ArrayList<String> normFactorSmallerChr = inputNormChr1.readAsArrayList();
        inputNormChr1.close();

//        System.out.println("Done reading norm factor 1");
        if (!Gpio.exists(fileToRead + ".sorted")) {
            umcg.genetica.io.chrContacts.SortIntraChrContacts.readNonSortedWriteSorted(fileToRead, fileToRead + ".sorted");
        }

        int numberToBeMatched = 0;

        LineIterator it = FileUtils.lineIterator(new File(fileToRead + ".sorted"), "UTF-8");

        try {
            while (it.hasNext()) {
                String[] parts = StringUtils.split(it.nextLine(), '\t');

                int posChr1 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[0]);
                int posChr2 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[1]);

                while (numberToBeMatched < contactsToCheck.size()) {
                    if (posChr1 < contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        break;
                    } else if (posChr1 == contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        if (posChr2 < contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            break;
                        }
                        if (posChr2 == contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {

                            String factor1Base = normFactorSmallerChr.get((posChr1 / getNumericResolution(resolution)) + 1);
                            String factor2Base = normFactorSmallerChr.get((posChr2 / getNumericResolution(resolution)) + 1);

                            double factor1;
                            double factor2;

                            if (StringUtils.isNumeric(factor1Base) && StringUtils.isNumeric(factor2Base)) {
                                factor1 = org.apache.commons.lang.math.NumberUtils.createDouble(factor1Base);
                                factor2 = org.apache.commons.lang.math.NumberUtils.createDouble(factor2Base);

                                double contact = org.apache.commons.lang.math.NumberUtils.createDouble(parts[2]) / (factor1 * factor2);
                                if (contact >= minValue) {
                                    outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\tContact\t" + contact + "\t" + org.apache.commons.lang.math.NumberUtils.createDouble(parts[2]));
                                    numberToBeMatched++;
                                } else {
                                    outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-\t-");
                                    numberToBeMatched++;
                                }
                            } else {
                                System.out.println("Error in files.");
                                numberToBeMatched++;
                            }
                        } else if (posChr2 > contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-\t-");
                            numberToBeMatched++;
                        }
                    } else if (posChr1 > contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-\t-");
                        numberToBeMatched++;
                    }
                }
            }
        } finally {
            LineIterator.closeQuietly(it);
        }

    }

    private static HashMap<String, ArrayList<DesiredChrContact>> readInQtlTransformBlocks(String in, String inProxies, String probeMap, String snpMap, boolean permutationFile, String resolution, boolean alternativePermutation) {
        ArrayList<MinimalEQTL> qtls = null;
        boolean renameSnpProxie = false;
        System.out.println("Read in QTL info.");
        try {
            qtls = readInQtlInformation(in, inProxies, probeMap, snpMap, permutationFile, renameSnpProxie);
        } catch (IOException ex) {
            Logger.getLogger(HiCQTLAnnotatorBlockbased.class.getName()).log(Level.SEVERE, null, ex);
        }
        HashMap<String, ArrayList<DesiredChrContact>> desiredContacts;
        HashSet<String> keysAdded = new HashSet<>();
        if (!alternativePermutation) {
            System.out.println("Transforming to blocks");
            desiredContacts = new HashMap<>(qtls.size());
            for (MinimalEQTL qtl : qtls) {

                String chrProbe = String.valueOf(qtl.getProbeChr());
                String chrSnp = String.valueOf(qtl.getRsChr());
                int biny = qtl.getRsChrPos() - (qtl.getRsChrPos() % getNumericResolution(resolution));
                int binx = qtl.getProbeChrPos() - (qtl.getProbeChrPos() % getNumericResolution(resolution));
                int bin1;
                int bin2;
                String ChrSmaller;
                String ChrLarger;

                if (chrProbe.equals(chrSnp)) {
                    ChrSmaller = chrProbe;
                    ChrLarger = chrSnp;
                    if (binx < biny) {
                        bin1 = binx;
                        bin2 = biny;
                    } else {
                        bin2 = binx;
                        bin1 = biny;
                    }

                } else {
                    if (Integer.parseInt(chrProbe) < Integer.parseInt(chrSnp)) {
                        ChrSmaller = chrProbe;
                        ChrLarger = chrSnp;
                        bin1 = binx;
                        bin2 = biny;
                    } else {
                        ChrSmaller = chrSnp;
                        ChrLarger = chrProbe;
                        bin1 = biny;
                        bin2 = binx;
                    }
                }

                String key = ChrSmaller + "-" + ChrLarger;
                if (!desiredContacts.containsKey(key)) {

                    desiredContacts.put(key, new ArrayList<DesiredChrContact>());
                }

                String extendedKey = key + "-" + bin1 + "-" + bin2 + "-" + qtl.getRsName() + "-" + qtl.getProbe();

                if (!keysAdded.contains(extendedKey)) {
                    desiredContacts.get(key).add(new DesiredChrContact(bin1, bin2, 0, qtl.getRsName(), qtl.getProbe()));
                    keysAdded.add(extendedKey);
                }
            }
        } else {
            System.out.println("Alternative permutation and transforming to blocks");
            HashMap<String, Pair<String, Integer>> cpgs = new HashMap<>();
            HashMap<String, Pair<String, Integer>> snps = new HashMap<>();
            HashSet<String> trueCombinations = new HashSet<>(qtls.size());

            //Here we need to merge creation of opposites and the blocking
            System.out.println("step1");
            for (MinimalEQTL e : qtls) {
                if (!cpgs.containsKey(e.getProbe())) {
                    cpgs.put(e.getProbe(), new Pair<>(String.valueOf(e.getProbeChr()), (e.getProbeChrPos() - (e.getProbeChrPos() % getNumericResolution(resolution)))));
                }
                if (!snps.containsKey(e.getRsName())) {
                    snps.put(e.getRsName(), new Pair<>(String.valueOf(e.getRsChr()), (e.getRsChrPos() - (e.getRsChrPos() % getNumericResolution(resolution)))));
                }
                trueCombinations.add(e.getProbe() + "-" + e.getRsName());
            }
            qtls = null;

            desiredContacts = new HashMap<>(snps.size() * cpgs.size());
            System.out.println("step2");
            for (Entry<String, Pair<String, Integer>> snp : snps.entrySet()) {

                String chrSnp = String.valueOf(snp.getValue().getLeft());
                int biny = snp.getValue().getRight();

                for (Entry<String, Pair<String, Integer>> cpg : cpgs.entrySet()) {
                    if (!trueCombinations.contains(cpg.getKey() + "-" + snp.getKey())) {

                        String chrProbe = cpg.getValue().getLeft();
                        int binx = cpg.getValue().getRight();

                        int bin1;
                        int bin2;
                        String ChrSmaller;
                        String ChrLarger;

                        if (chrProbe.equals(chrSnp)) {
                            ChrSmaller = chrProbe;
                            ChrLarger = chrSnp;
                            if (binx < biny) {
                                bin1 = binx;
                                bin2 = biny;
                            } else {
                                bin2 = binx;
                                bin1 = biny;
                            }

                        } else {
                            if (Integer.parseInt(chrProbe) < Integer.parseInt(chrSnp)) {
                                ChrSmaller = chrProbe;
                                ChrLarger = chrSnp;
                                bin1 = binx;
                                bin2 = biny;
                            } else {
                                ChrSmaller = chrSnp;
                                ChrLarger = chrProbe;
                                bin1 = biny;
                                bin2 = binx;
                            }
                        }

                        String key = ChrSmaller + "-" + ChrLarger;
                        if (!desiredContacts.containsKey(key)) {
                            desiredContacts.put(key, new ArrayList<DesiredChrContact>());
                        }
                        String extendedKey = key + "-" + bin1 + "-" + bin2 + "-" + snp.getKey() + "-" + cpg.getKey();
                        if (!keysAdded.contains(extendedKey)) {
                            desiredContacts.get(key).add(new DesiredChrContact(bin1, bin2, 0, snp.getKey(), cpg.getKey()));
                            keysAdded.add(extendedKey);
                        }
                    }
                }
            }

        }
        System.out.println("Number of desiredContacts: " + keysAdded.size());
        return desiredContacts;
    }

    private static void printOutContacts(ArrayList<DesiredChrContact> contacts, TextFile outWriter) throws IOException {
//        System.out.println("Write contacts to file.");
        HashMap<String, Boolean> textToStore = new HashMap<>();

        for (DesiredChrContact c : contacts) {
            String key = c.getProbeName() + "-" + c.getSnpName();

            if (c.hasContact()) {
                textToStore.put(key, Boolean.TRUE);
            } else if (!textToStore.containsKey(key)) {
                textToStore.put(key, Boolean.FALSE);
            }
        }
        for (Entry<String, Boolean> contactInfo : textToStore.entrySet()) {
            outWriter.write(contactInfo.getKey() + "\t" + contactInfo.getValue() + "\n");
        }
    }

}
