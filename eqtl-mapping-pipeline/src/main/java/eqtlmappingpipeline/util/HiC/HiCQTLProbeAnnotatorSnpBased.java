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
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.chrContacts.DesiredChrContact;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.MinimalEQTL;

/**
 *
 * @author MarcJan
 */
class HiCQTLProbeAnnotatorSnpBased {

    public static void main(String[] args) throws IOException {
        String folderHighC = "G:\\Contacts\\";
        String resolution = "1kb"; //5kb / 1kb
        String qualityCutOff = "E30"; //0 or E30
//        String normMethod = null;
        String normMethod = null; //null / "KRnorm" / "SQRTVCnorm" / "VCnorm"
        double minValueQuality = 0.0;
        boolean permutationFile = false;
        String probeMap = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Illumina450K_MQtlMappingFile_MJB.txt";
        String snpMap = "D:\\Werk\\UMCGundefinedUMCG\\Projects\\LL-DeepBBMRI_Methylation450K\\meQTLs\\SNPMappings\\SNPMappings.txt";
        System.out.println("Running");

        String QTLfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\ArrayInQTLFormatForAnnot.txt";

        String QTLoutfile = "D:\\WebFolders\\OwnCloud\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\ArrayInQTLFormatForAnnot_ProbeBased_HiC_1kb_E30_annotated.txt";

        addAnnotationToQTLOutput(
                QTLfile,
                folderHighC,
                resolution,
                qualityCutOff,
                normMethod,
                minValueQuality,
                permutationFile,
                probeMap,
                snpMap,
                QTLoutfile);

    }

    private static void addAnnotationToQTLOutput(String in, String folderHighC, String resolution, String qualityCutOff, String normMethod, double minValue, boolean permutationFile, String probeMap, String snpMap, String out) throws IOException {
        int[] chrLocations = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
        HashMap<String, ArrayList<DesiredChrContact>> qtls = readInQtlProbeInformation(in, probeMap, snpMap, permutationFile, resolution);

        ProgressBar pb = new ProgressBar(qtls.size(), "Checking for contacts for: " + qtls.size() + " Chromosome combinations");

        TextFile outWriter = new TextFile(out, TextFile.W);

        for (Entry<String, ArrayList<DesiredChrContact>> probeChr : qtls.entrySet()) {
            Collections.sort(probeChr.getValue());
            ArrayList<DesiredChrContact> copy = new ArrayList<>();
            for (DesiredChrContact d : probeChr.getValue()) {
                copy.add(d.clone());
            }

            for (int otherChr : chrLocations) {

                String ChrSmaller = probeChr.getKey();
                String ChrLarger = String.valueOf(otherChr);
                boolean probeIsSmaller = true;

                if (otherChr < Integer.parseInt(ChrSmaller)) {
                    ChrSmaller = String.valueOf(otherChr);
                    ChrLarger = probeChr.getKey();
                    probeIsSmaller = false;
                }

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
                    if (intra) {
                        processRawContactInformation(fileToReads, minValue, copy, intra, probeIsSmaller);
                    } else {
                        processRawContactInformation(fileToReads, minValue, probeChr.getValue(), intra, probeIsSmaller);
                    }
                } else {
                    throw new UnsupportedOperationException("Not supported yet.");
                }
            }

            writeContactInformation(outWriter, probeChr.getValue(), copy);
            pb.iterate();
        }
        pb.close();
        outWriter.close();
    }

    //For example, here is a line from the 5kb chr1 MAPQGE30 raw observed contact matrix (GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
//40000000 40100000 59.0
    private static void processRawContactInformation(String fileToRead, double minValue, ArrayList<DesiredChrContact> contactsToCheck, boolean intra, boolean probeIsSmaller) throws IOException {

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

                int posChr1;

                if (probeIsSmaller) {
                    posChr1 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[0]);
                } else {
                    posChr1 = org.apache.commons.lang.math.NumberUtils.createInteger(parts[1]);
                }

                while (numberToBeMatched < contactsToCheck.size()) {
                    if (posChr1 < contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        break;
                    } else if (posChr1 == contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        double contact = org.apache.commons.lang.math.NumberUtils.createDouble(parts[2]);
                        if (contact >= minValue) {
                            contactsToCheck.get(numberToBeMatched).setNormalizedContactValue(contactsToCheck.get(numberToBeMatched).getNormalizedContactValue() + 1);
                        }
                        numberToBeMatched++;

                    } else if (posChr1 > contactsToCheck.get(numberToBeMatched).getChrLocationSmaller()) {
                        numberToBeMatched++;
                    }
                }
            }
        } finally {
            LineIterator.closeQuietly(it);
        }

    }

    private static HashMap<String, ArrayList<DesiredChrContact>> readInQtlProbeInformation(String in, String probeMap, String snpMap, boolean permutationFile, String resolution) {
        ArrayList<MinimalEQTL> qtls = null;
        boolean renameSnpProxie = true;
        try {
            qtls = HiCQTLAnnotatorSnpBased.readInQtlInformation(in, null, probeMap, snpMap, permutationFile, renameSnpProxie);
        } catch (IOException ex) {
            Logger.getLogger(HiCQTLProbeAnnotatorSnpBased.class.getName()).log(Level.SEVERE, null, ex);
        }

        HashMap<String, ArrayList<DesiredChrContact>> desiredContacts = new HashMap<>();
        HashSet<String> probes = new HashSet<>();
        for (MinimalEQTL qtl : qtls) {

            int bin1 = qtl.getProbeChrPos() - (qtl.getProbeChrPos() % HiCQTLAnnotatorSnpBased.getNumericResolution(resolution));
            String ChrSmaller = String.valueOf(qtl.getProbeChr());

            String key = ChrSmaller;
            if (!desiredContacts.containsKey(key)) {
                desiredContacts.put(key, new ArrayList<DesiredChrContact>());
            }
            if (!probes.contains(qtl.getProbe())) {
                desiredContacts.get(key).add(new DesiredChrContact(bin1, 0, 0, null, qtl.getProbe()));
                probes.add(qtl.getProbe());
            }
        }
        return desiredContacts;
    }

    private static void writeContactInformation(TextFile outWriter, ArrayList<DesiredChrContact> valueIntraChr, ArrayList<DesiredChrContact> valueInterChr) throws IOException {
        int i = 0;
        for (DesiredChrContact d : valueIntraChr) {
            String contactIntra = "-";
            if (d.getNormalizedContactValue() > 0) {
                contactIntra = "contact";
            }
            String contactInter = "-";
            if (valueInterChr.get(i).getNormalizedContactValue() > 0) {
                contactInter = "contact";
            }
            outWriter.writeln(d.getProbeName() + "\t" + contactIntra + "\t" + d.getNormalizedContactValue() + "\t" + contactInter + "\t" + valueInterChr.get(i).getNormalizedContactValue());
            ++i;
        }

    }

}
