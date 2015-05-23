/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author MarcJan
 */
class HiCTransQTLAnnotator {

    private static final Pattern SPLIT_TAB = Pattern.compile("\t");

    public static void main(String[] args) throws IOException {

        String QTLfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered.txt";
        String proxyfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\proxiesMeQTLSnps.txt";
//        String proxyfile = null;
        String QTLoutfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\meQTLs\\Trans_Pc22c_CisMeQTLc_meQTLs\\RegressedOut_CisEffects_New\\eQTLsFDR0.05-ProbeLevel_BsFiltered&Filtered_HiC_LD_annotated.txt";
        String folderHighC = "F:\\Contacts\\GM12878_combined_interchromosomal\\";
        String resolution = "1kb";
        String qualityCutOff = "E30"; //0 or E30
        String normMethod = null; //null / KRnorm / SQRTVCnorm / VCnorm
        double minValueQuality = 0;

        addAnnotationToQTLOutput(
                QTLfile,
                proxyfile,
                folderHighC,
                resolution,
                qualityCutOff,
                normMethod,
                minValueQuality,
                QTLoutfile);

    }

    static void addAnnotationToQTLOutput(String in, String inProxies, String folderHighC, String resolution, String qualityCutOff, String normMethod, double minValue, String out) throws IOException {
        QTLTextFile eqtlTextFile = new QTLTextFile(in, QTLTextFile.R);

        ArrayList<EQTL> qtls = eqtlTextFile.readList();

        if (inProxies != null) {
            qtls = includeProxyInfo(qtls, inProxies);
        }
        
        System.out.println("Contacts to test: "+qtls.size());

        //Here we need to make a new Type to store the potentialy inflated files.
        TextFile outWriter = new TextFile(out, TextFile.W);
        HashMap<String, String> contactBuffer = new HashMap<String, String>();
        for (EQTL eqtl : qtls) {
            String chrProbe = String.valueOf(eqtl.getProbeChr());
            String chrSnp = String.valueOf(eqtl.getRsChr());

//            System.out.println(chrProbe+"\t"+chrSnp);
            if (chrProbe.equals(chrSnp)) {
                //Here we need to check how to normalize and treat intra-chromosomal data.
                continue;
            }

            int posChrSmaller;
            int posChrLarger;
            String ChrSmaller;
            String ChrLarger;
            int bin1;
            int bin2;
            String baseName;
            String fileToReads;

            if (Integer.parseInt(chrProbe) < Integer.parseInt(chrSnp)) {
                posChrSmaller = eqtl.getProbeChrPos();
                posChrLarger = eqtl.getRsChrPos();

                ChrSmaller = chrProbe;
                ChrLarger = chrSnp;

                //Determine bin1
                //Startscounting at 0-resulution
                bin1 = posChrSmaller - (posChrSmaller % getNumericResolution(resolution));
                //        System.out.println("\t"+bin1);
                //Determine bin2
                bin2 = posChrLarger - (posChrLarger % getNumericResolution(resolution));

                baseName = folderHighC + resolution + "_resolution_interchromosomal\\chr" + chrProbe + "_chr" + chrSnp + "\\MAPQG" + qualityCutOff;
                fileToReads = baseName + "\\chr" + chrProbe + "_" + chrSnp + "_" + resolution + ".RAWobserved";
//                    System.out.println("Reading: " + fileToReads);

            } else {
                posChrSmaller = eqtl.getRsChrPos();
                posChrLarger = eqtl.getProbeChrPos();

                ChrSmaller = chrSnp;
                ChrLarger = chrProbe;

                bin1 = posChrSmaller - (posChrSmaller % getNumericResolution(resolution));
                //        System.out.println("\t"+bin1);
                //Determine bin2
                bin2 = posChrLarger - (posChrLarger % getNumericResolution(resolution));

                baseName = folderHighC + resolution + "_resolution_interchromosomal\\chr" + chrSnp + "_chr" + chrProbe + "\\MAPQG" + qualityCutOff;
                fileToReads = baseName + "\\chr" + chrSnp + "_" + chrProbe + "_" + resolution + ".RAWobserved";
//                    System.out.println("Reading: " + fileToReads);

            }

            if (normMethod == null) {
                if (contactBuffer.containsKey(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2)) {
                    outWriter.writeln(eqtl.getRsName() + "\t" + eqtl.getProbe() + "\t" + contactBuffer.get(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2));
                } else {
                    if (readRawInterContactInformation(fileToReads, minValue, bin1, bin2)) {
                        outWriter.writeln(eqtl.getRsName() + "\t" + eqtl.getProbe() + "\tContact");
                        contactBuffer.put(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2, "Contact");
                    } else {
                        outWriter.writeln(eqtl.getRsName() + "\t" + eqtl.getProbe() + "\t-");
                        contactBuffer.put(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2, "-");
                    }
                }
            } else {
                if (contactBuffer.containsKey(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2)) {
                    outWriter.writeln(eqtl.getRsName() + "\t" + eqtl.getProbe() + "\t" + contactBuffer.get(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2));
                } else {
                    if (readNormalizedInterContactInformation(fileToReads, baseName, normMethod, ChrSmaller, ChrLarger, posChrSmaller, posChrLarger, resolution, minValue)) {
                        outWriter.writeln(eqtl.getRsName() + "\t" + eqtl.getProbe() + "\tContact");
                        contactBuffer.put(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2, "Contact");
                    } else {
                        outWriter.writeln(eqtl.getRsName() + "\t" + eqtl.getProbe() + "\t-");
                        contactBuffer.put(ChrSmaller + "_" + ChrLarger + "_" + bin1 + "_" + bin2, "-");
                    }
                }
            }

        }
        outWriter.close();
    }

    private static ArrayList<EQTL> includeProxyInfo(ArrayList<EQTL> qtls, String inProxies) throws IOException {
        ArrayList<EQTL> newQtlList = new ArrayList<EQTL>();

        TextFile readProxies = new TextFile(inProxies, TextFile.R);

        String line = readProxies.readLine();
//        System.out.println(line);
        while ((line = readProxies.readLine()) != null) {
//            System.out.println(line);
            String[] lineParts = SPLIT_TAB.split(line);
            String chr = lineParts[4];
            int chrPos = Integer.parseInt(lineParts[5]);
            int chrNewPos = Integer.parseInt(lineParts[8]);
            for (EQTL e : qtls) {
                if (String.valueOf(e.getRsChr()).equals(chr) && e.getRsChrPos() == chrPos) {
                    EQTL newQtl = new EQTL();
                    newQtl.setProbe(e.getProbe());
                    newQtl.setProbeChr(e.getProbeChr());
                    newQtl.setProbeChrPos(e.getProbeChrPos());

                    newQtl.setRsName(e.getRsName() + "-" + lineParts[1]);
                    newQtl.setRsChr(e.getRsChr());
                    newQtl.setRsChrPos(chrNewPos);
                    newQtlList.add(newQtl);
                }
            }
        }

        for (EQTL e : qtls) {
            newQtlList.add(e);
        }

        return newQtlList;
    }

    //For example, here is a line from the 5kb chr1 MAPQGE30 raw observed contact matrix (GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
//40000000 40100000 59.0
    private static boolean readRawInterContactInformation(String fileToReads, double minValue, int bin1, int bin2) throws IOException {
//        System.out.println("\t\t"+fileToReads);
//        System.out.println("\t"+bin2);
        //See if bin1 and bin2 are in the file.
        boolean contactFound = false;

        //Check if sorted version is available
        //If not make sorted available.
        if (!Gpio.exists(fileToReads + ".sorted")) {
            umcg.genetica.io.chrContacts.SortInterChrContacts.readNonSortedWriteSorted(fileToReads, fileToReads + ".sorted");
        }

        BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(fileToReads + ".sorted"), "UTF-8"));

        String row;

        while ((row = input.readLine()) != null) {
            String[] parts = StringUtils.split(row, '\t');
//            System.out.println(row);
            int posChr1 = Integer.parseInt(parts[0]);
            if (posChr1 == bin1) {
                int posChr2 = Integer.parseInt(parts[1]);
                if (posChr2 == bin2) {
                    double contact = Double.parseDouble(parts[2]);
                    if (contact >= minValue) {
                        contactFound = true;
                    }
                    break;
                } else if (posChr2 > bin2) {
                    break;
                }
            } else if (posChr1 > bin1) {
                break;
            }
        }
        input.close();
        return contactFound;
    }

    //For example, here is a line from the 5kb chr1 MAPQGE30 raw observed contact matrix (GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
    //40000000 40100000 59.0
    //To normalize this entry using the KR normalization vector, one would divide 59.0 by the 8001st line ((40000000/5000)+1=8001) and the 8021st line ((40100000/5000)+1=8021) 
    //of GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.KRnorm. The 8001st line of the KR norm file is 1.2988778370674694;
    //The 8021st line of the KR norm file is 1.6080499717941548. So the corresponding KR normalized entry for the entry above is 59.0/(1.2988778370674694*1.6080499717941548)
    //or 28.24776973966101. 
    //If the KR normalization vector file is empty or all NaNs, then the KR algorithm didn’t converge on that particular matrix (likely due to sparsity of the matrix).
    private static boolean readNormalizedInterContactInformation(String fileToRead, String baseName, String normMethod, String chrSmaller, String chrLarger, int posChrSmaller, int posChrLarger, String resolution, double minValue) throws IOException {
        //Determine bin1
        //Starts counting at 0-resulution
        int bin1 = posChrSmaller - (posChrSmaller % getNumericResolution(resolution));

        //Determine bin2
        int bin2 = posChrLarger - (posChrLarger % getNumericResolution(resolution));

        //ReadIn normalization chr1
        TextFile inputNormChr1 = new TextFile(baseName + "\\chr" + chrSmaller + "_" + resolution + "." + normMethod, TextFile.R);
        ArrayList<String> normFactorSmallerChr = inputNormChr1.readAsArrayList();
        inputNormChr1.close();

        //ReadIn normalization chr2
        TextFile inputNormChr2 = new TextFile(baseName + "\\chr" + chrLarger + "_" + resolution + "." + normMethod, TextFile.R);
        ArrayList<String> normFactorLargerChr = inputNormChr2.readAsArrayList();

        inputNormChr2.close();

        if (!Gpio.exists(fileToRead + ".sorted")) {
            umcg.genetica.io.chrContacts.SortInterChrContacts.readNonSortedWriteSorted(fileToRead, fileToRead + ".sorted");
        }

        BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(fileToRead + ".sorted"), "UTF-8"));

        String row;

        //See if bin1 and bin2 are in the file.
        boolean contactFound = false;

        while ((row = input.readLine()) != null) {
            String[] parts = StringUtils.split(row, '\t');

            int posChr1 = Integer.parseInt(parts[0]);
            if (posChr1 == bin1) {
                int posChr2 = Integer.parseInt(parts[1]);
                if (posChr2 == bin2) {

                    String factor1Base = normFactorSmallerChr.get((posChr1 / getNumericResolution(resolution)) + 1);
                    String factor2Base = normFactorLargerChr.get((posChr2 / getNumericResolution(resolution)) + 1);

                    double factor1;
                    double factor2;

                    if (StringUtils.isNumeric(factor1Base) && StringUtils.isNumeric(factor2Base)) {
                        factor1 = Double.parseDouble(factor1Base);
                        factor2 = Double.parseDouble(factor2Base);

                        double contact = Double.parseDouble(parts[2]) / (factor1 * factor2);
                        if (contact >= minValue) {
                            contactFound = true;
                        }
                        break;
                    }

                } else if (posChr2 > bin2) {
                    break;
                }
            } else if (posChr1 > bin1) {
                break;
            }
        }
        input.close();
        return contactFound;
    }

    private static int getNumericResolution(String resolution) {
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
}
