/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util.HiC;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.chrContacts.DesiredChrContact;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.MinimalEQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author MarcJan
 */
class HiCQTLAnnotatorSnpBased {

    private static final Pattern SPLIT_TAB = Pattern.compile("\t");

    public static void main(String[] args) throws IOException {
        String folderHighC = "G:\\Contacts\\";
        String resolution = "1kb"; //5kb / 1kb
        String qualityCutOff = "E30"; //0 or E30
//        String normMethod = null;
        String normMethod = "KRnorm"; //null / "KRnorm" / "SQRTVCnorm" / "VCnorm"
        double minValueQuality = 0.0;
        boolean permutationFile = true;
        String probeMap = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Illumina450K_MQtlMappingFile_MJB.txt";
        String snpMap = "D:\\Werk\\UMCGundefinedUMCG\\Projects\\LL-DeepBBMRI_Methylation450K\\meQTLs\\SNPMappings\\SNPMappings.txt";
        System.out.println("Running");
//        for (int i = 1; i < 11; ++i) {
//            String QTLfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered.txt";
//            String QTLfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\eQTLProbesFDR0.05-1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-ProbeLevel.txt_LdQtlFilterd.txt";
        String QTLfile = "D:\\Werk\\UMCGundefinedUMCG\\Projects\\LL-DeepBBMRI_Methylation450K\\meQTLs\\trans-QTLs\\PermutedEQTLsPermutationRound57.head.txt";
//            String proxyfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\proxiesMeQTLSnps2.txt";
        String proxyfile = null;
        String QTLoutfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_AlternativePermutation_HiC_LD_Kr_1kb_E30_annotated.txt";
//            String QTLoutfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Cis_Pc22c_meQTLs\\HiC\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_HiC_1kb_E30_annotated.txt";

        addAnnotationToQTLOutput(
                QTLfile,
                proxyfile,
                folderHighC,
                resolution,
                qualityCutOff,
                normMethod,
                minValueQuality,
                permutationFile,
                probeMap,
                snpMap,
                QTLoutfile);
//        }

    }

    private static void addAnnotationToQTLOutput(String in, String inProxies, String folderHighC, String resolution, String qualityCutOff, String normMethod, double minValue, boolean permutationFile, String probeMap, String snpMap, String out) throws IOException {

        HashMap<String, ArrayList<DesiredChrContact>> qtls = readInQtlInformation2(in, inProxies, probeMap, snpMap, permutationFile, resolution);

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

            if (normMethod == null) {
                processRawContactInformation(fileToReads, minValue, contactsToCheck.getValue(), outWriter, intra);
            } else {
                if (intra) {
                    processNormalizedIntraContactInformation(fileToReads, baseName, normMethod, ChrSmaller, contactsToCheck.getValue(), resolution, minValue, outWriter);
                } else {
                    processNormalizedInterContactInformation(fileToReads, baseName, normMethod, ChrSmaller, ChrLarger, contactsToCheck.getValue(), resolution, minValue, outWriter);
                }
            }

            pb.iterate();
        }
        pb.close();
        outWriter.close();
    }

    //For example, here is a line from the 5kb chr1 MAPQGE30 raw observed contact matrix (GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
//40000000 40100000 59.0
    private static void processRawContactInformation(String fileToRead, double minValue, ArrayList<DesiredChrContact> contactsToCheck, TextFile outWriter, boolean intra) throws IOException {

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
                                outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\tContact\t" + contact);
                                numberToBeMatched++;
                            } else {
                                outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-");
                                numberToBeMatched++;
                            }
                        } else if (posChr2 > contactsToCheck.get(numberToBeMatched).getChrLocationLarger()) {
                            outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-");
                            numberToBeMatched++;
                        }
                    } else if (posChr1 > contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        outWriter.writeln(contactsToCheck.get(numberToBeMatched).getSnpName() + "\t" + contactsToCheck.get(numberToBeMatched).getProbeName() + "\t" + posChr1 + "\t" + posChr2 + "\t-\t-");
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

    public static int getNumericResolution(String resolution) {
        if (resolution.equals("1kb")) {
            return 1000;
        } else if (resolution.equals("5kb")) {
            return 5000;
        } else {
            System.out.println("\nError in resolution setting!\n");
            System.exit(-1);
        }
        return 0;
    }

    private static HashMap<String, ArrayList<DesiredChrContact>> readInQtlInformation2(String in, String inProxies, String probeMap, String snpMap, boolean permutationFile, String resolution) {
        ArrayList<MinimalEQTL> qtls = null;
        boolean renameSnpProxie = true;
        try {
            qtls = readInQtlInformation(in, inProxies, probeMap, snpMap, permutationFile, renameSnpProxie);
        } catch (IOException ex) {
            Logger.getLogger(HiCQTLAnnotatorSnpBased.class.getName()).log(Level.SEVERE, null, ex);
        }

        HashMap<String, ArrayList<DesiredChrContact>> desiredContacts = new HashMap<String, ArrayList<DesiredChrContact>>();

        for (MinimalEQTL qtl : qtls) {

            String chrProbe = String.valueOf(qtl.getProbeChr());
            String chrSnp = String.valueOf(qtl.getRsChr());
            int posChrSmaller;
            int posChrLarger;
            int bin1;
            int bin2;
            String ChrSmaller;
            String ChrLarger;

            if (chrProbe.equals(chrSnp)) {
                ChrSmaller = chrProbe;
                ChrLarger = chrSnp;

                int binx = qtl.getProbeChrPos() - (qtl.getProbeChrPos() % getNumericResolution(resolution));
                //        System.out.println("\t"+bin1);
                //Determine bin2
                int biny = qtl.getRsChrPos() - (qtl.getRsChrPos() % getNumericResolution(resolution));

                //Check if logic is consistent.
                if (binx < biny) {
                    bin1 = binx;
                    bin2 = biny;
                } else {
                    bin2 = binx;
                    bin1 = biny;
                }

            } else {

                if (Integer.parseInt(chrProbe) < Integer.parseInt(chrSnp)) {
                    posChrSmaller = qtl.getProbeChrPos();
                    posChrLarger = qtl.getRsChrPos();

                    //Determine bin1
                    //Startscounting at 0-resulution
                    bin1 = posChrSmaller - (posChrSmaller % getNumericResolution(resolution));
                    //        System.out.println("\t"+bin1);
                    //Determine bin2
                    bin2 = posChrLarger - (posChrLarger % getNumericResolution(resolution));

                    ChrSmaller = chrProbe;
                    ChrLarger = chrSnp;
                } else {
                    posChrSmaller = qtl.getRsChrPos();
                    posChrLarger = qtl.getProbeChrPos();

                    bin1 = posChrSmaller - (posChrSmaller % getNumericResolution(resolution));
                    //        System.out.println("\t"+bin1);
                    //Determine bin2
                    bin2 = posChrLarger - (posChrLarger % getNumericResolution(resolution));

                    ChrSmaller = chrSnp;
                    ChrLarger = chrProbe;
                }
            }

            String key = ChrSmaller + "-" + ChrLarger;
            if (!desiredContacts.containsKey(key)) {
                desiredContacts.put(key, new ArrayList<DesiredChrContact>());
            }
            desiredContacts.get(key).add(new DesiredChrContact(bin1, bin2, 0, qtl.getRsName(), qtl.getProbe()));
        }
        return desiredContacts;
    }

    public static ArrayList<MinimalEQTL> readInQtlInformation(String in, String inProxies, String probeMap, String snpMap, boolean permutationFile, boolean changeNameOfSnp) throws IOException {
        ArrayList<MinimalEQTL> qtls = null;
        if (!permutationFile) {
            QTLTextFile eqtlTextFile = new QTLTextFile(in, QTLTextFile.R);
            ArrayList<EQTL> qtlsBuffer = eqtlTextFile.readList();
            qtls = MinimalEQTL.convertArray(qtlsBuffer);
            eqtlTextFile.close();
        } else {

            HashMap<String, Pair<Byte, Integer>> probeLocation = readChrLocation(probeMap, 1, 3, 5, true);
            HashMap<String, Pair<Byte, Integer>> snpLocation = readChrLocation(snpMap, 2, 0, 1, true);

            TextFile textFile = new TextFile(in, TextFile.R);
            qtls = new ArrayList<>();

            String row = textFile.readLine();
            while ((row = textFile.readLine()) != null) {
                String[] parts = StringUtils.split(row, '\t');
//                System.out.println(row);

                MinimalEQTL newQtl = new MinimalEQTL();
                newQtl.setProbe(parts[2]);

                newQtl.setProbeChr(probeLocation.get(parts[2]).getLeft());
                newQtl.setProbeChrPos(probeLocation.get(parts[2]).getRight());

                newQtl.setRsName(parts[1]);
                if (snpLocation.containsKey(parts[1])) {
                    newQtl.setRsChr(snpLocation.get(parts[1]).getLeft());
                    newQtl.setRsChrPos(snpLocation.get(parts[1]).getRight());
                } else {
                    System.out.println("Error SNP: " + parts[1] + " not present in your SNP mapping file.");
                    System.exit(0);
                }
                qtls.add(newQtl);

            }

            textFile.close();
        }
        System.out.println("Number of QTLs input: " + qtls.size());
        if (inProxies != null) {
            System.out.println("Including proxies");
            qtls = includeProxyInfo(qtls, inProxies, changeNameOfSnp);
            System.out.println("Number of QTLs after including proxies: " + qtls.size());
        }

////        Write file to disk the input information.
//        QTLTextFile out = new QTLTextFile(in + ".extended.txt", QTLTextFile.W);
//        out.writeMinimal(qtls);
//        out.close();
//        System.exit(0);
        return qtls;
    }

    private static ArrayList<MinimalEQTL> includeProxyInfo(ArrayList<MinimalEQTL> qtls, String inProxies, boolean changeNameOfSnp) throws IOException {
        ArrayList<MinimalEQTL> newQtlList = new ArrayList<>();

        TextFile readProxies = new TextFile(inProxies, TextFile.R);

        String line = readProxies.readLine();
//        System.out.println(line);
        while ((line = readProxies.readLine()) != null) {
//            System.out.println(line);
            String[] lineParts = SPLIT_TAB.split(line);
            String chr = lineParts[4];
            int chrPos = Integer.parseInt(lineParts[5]);
            int chrNewPos = Integer.parseInt(lineParts[8]);
            for (MinimalEQTL e : qtls) {
                if (String.valueOf(e.getRsChr()).equals(chr) && e.getRsChrPos() == chrPos) {
                    MinimalEQTL newQtl = new MinimalEQTL();
                    newQtl.setProbe(e.getProbe());
                    newQtl.setProbeChr(e.getProbeChr());
                    newQtl.setProbeChrPos(e.getProbeChrPos());

                    if (changeNameOfSnp) {
                        newQtl.setRsName(e.getRsName() + "-" + lineParts[1]);
                    } else {
                        newQtl.setRsName(e.getRsName());
                    }
                    newQtl.setRsChr(e.getRsChr());
                    newQtl.setRsChrPos(chrNewPos);
                    newQtlList.add(newQtl);
                }
            }
        }

        for (MinimalEQTL e : qtls) {
            newQtlList.add(e);
        }

        return newQtlList;
    }

    private static HashMap<String, Pair<Byte, Integer>> readChrLocation(String file, int key, int value1, int value2, boolean skipFirstRow) throws IOException {
        HashMap<String, Pair<Byte, Integer>> locationInfo = new HashMap<>();
        TextFile textFile = new TextFile(file, TextFile.R);
        ArrayList<EQTL> qtls = new ArrayList<>();

        String row;
        if (skipFirstRow) {
            row = textFile.readLine();
        }
        while ((row = textFile.readLine()) != null) {
//            System.out.println(row);
            String[] parts = SPLIT_TAB.split(row, '\t');
            locationInfo.put(parts[key], new Pair<>(ChrAnnotation.parseChr(parts[value1]), Integer.parseInt(parts[value2])));
        }
        textFile.close();
        return locationInfo;
    }

}
