/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.fasta.ReferenceGenomeFasta;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class WriteExtendedFastaOfInterestRegion {

    public static void main(String[] args) {
//            String referenceFasta = args[0];
//            String fileWithPositions = args[1];
//            String outputFile = args[2];

        String referenceFasta = "E:\\OnlineFolders\\BitSync\\UMCG\\ProbeMapping\\human_g1k_v37.fasta";
        
//        String fileWithPositions = "E:\\OnlineFolders\\BitSync\\UMCG\\ProbeMapping\\Info\\V71_2\\RNA_Seq\\IlluminaMatchedRNASeqProbes.gff";
//        String outputFile = "E:\\OnlineFolders\\BitSync\\UMCG\\ProbeMapping\\Fasta_Urmo_1.txt";
//        boolean bedFile = false;
        
        String fileWithPositions = "E:\\OnlineFolders\\AeroFS\\RP3_BIOS_Methylation\\Annotations\\Illumina450K_MQtlMappingFile_MJB.txt";
        String outputFile = "E:\\OnlineFolders\\BitSync\\UMCG\\Fastatmp.txt";
        boolean bedFile = true;
        
        int maxScore = 5;
        int minDiff = 30;
        
        ReferenceGenomeFasta refGen = null;
        HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> interestStrings = null;

        try {
            if(!bedFile){
                System.out.print("Reading reference GTF ...... ");
                interestStrings = readFileWithPositions(fileWithPositions);
            } else {
                System.out.println("Reading reference BED .....");
                interestStrings = readFileWithPositionsBed(fileWithPositions);
            }
        } catch (Exception ex) {
            Logger.getLogger(WriteExtendedFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println(" done");
        try {
            System.out.print("Reading reference FASTA ......");
            refGen = new ReferenceGenomeFasta(new File(referenceFasta), ReferenceGenomeFasta.HUMAN_NORMAL_CHR);
        } catch (Exception ex) {
            Logger.getLogger(WriteExtendedFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println(" done");
        if (refGen != null && interestStrings != null) {
            try {
                writeFastaIncDirectionGuessing(refGen, interestStrings, outputFile, maxScore, minDiff);
            } catch (IOException ex) {
                Logger.getLogger(WriteExtendedFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    private static HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> readFileWithPositions(String fileWithPositions) throws IOException {
        HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> chromosomalProbeInformation = new HashMap<>();
        TextFile gffFileReader = new TextFile(fileWithPositions, TextFile.R);

        int counter = 0;
        String str;
        while ((str = gffFileReader.readLine()) != null) {
            String[] parts = str.split("\t");
//            System.out.println(str);

            if (!chromosomalProbeInformation.containsKey(parts[0])) {
                chromosomalProbeInformation.put(parts[0], new HashMap<String, Triple<String, Integer, Integer>>());
            }

            chromosomalProbeInformation.get(parts[0]).put(parts[8], new Triple<String, Integer, Integer>(parts[6], Integer.parseInt(parts[3]), Integer.parseInt(parts[4])));
            counter++;
        }
        System.out.println("Number of entries read in: " + counter);
        gffFileReader.close();
        return chromosomalProbeInformation;
    }
    
    private static HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> readFileWithPositionsBed(String fileWithPositions) throws IOException {
        HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> chromosomalProbeInformation = new HashMap<>();
        TextFile gffFileReader = new TextFile(fileWithPositions, TextFile.R);

        int counter = 0;
        String str = gffFileReader.readLine();
        while ((str = gffFileReader.readLine()) != null) {
            String[] parts = str.split("\t");
//            System.out.println(str);

            if (!chromosomalProbeInformation.containsKey(parts[3])) {
                chromosomalProbeInformation.put(parts[3], new HashMap<String, Triple<String, Integer, Integer>>());
            }
//            //\"; IlluminaProbeSeq \" + parts[7]
            chromosomalProbeInformation.get(parts[3]).put(parts[1]+"\"; IlluminaProbeSeq \"" + parts[7], new Triple<String, Integer, Integer>(".", Integer.parseInt(parts[4]), Integer.parseInt(parts[5])));
            counter++;
        }
        System.out.println("Number of entries read in: " + counter);
        gffFileReader.close();
        return chromosomalProbeInformation;
    }

    private static void writeFastaIncDirectionGuessing(ReferenceGenomeFasta refGen, HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> interestStrings, String outputFile, int maxScore, int minDiff) throws IOException {
        TextFile writer = new TextFile(outputFile, TextFile.W);
        

        for (Entry<String, HashMap<String, Triple<String, Integer, Integer>>> probesPerChr : interestStrings.entrySet()) {
//            System.out.println(probesPerChr.getKey());

            for (Entry<String, Triple<String, Integer, Integer>> probe : probesPerChr.getValue().entrySet()) {
//                System.out.println("\t"+probe.getKey() + "      "+ probe.getValue().getLeft());
                try {
                    String[] partsKey = probe.getKey().split("\";");
                    StringBuilder s = new StringBuilder();
                    s.append(partsKey[0]).append("\t");
                    String probeSeq = refGen.getNucleotides(probesPerChr.getKey(), probe.getValue().getMiddle(), probe.getValue().getRight()).toString();

                    if (probe.getKey().contains(probeSeq)) {
                        s.append("+").append("\n");
//                        s.append(probeSeq).append("\n");
//                        System.out.println(s.toString());
                        writer.write(s.toString());
                    } else if (probe.getKey().contains(getFullComplement(reverse(probeSeq)))) {
                        s.append("-").append("\n");
//                        s.append(getFullComplement(reverse(probeSeq))).append("\n");
//                        System.out.println(s.toString());
                        writer.write(s.toString());
                    } else {
                        
                        partsKey[1] = partsKey[1].replace(" IlluminaProbeSeq \"", "");
                        if (partsKey[1].length() == probeSeq.length()) {
//                        System.out.println(partsKey[1]);
                            int editDistancePos = calcEditDistance(partsKey[1], probeSeq);
                            int editDistanceNeg = calcEditDistance(partsKey[1], getFullComplement(reverse(probeSeq)).toString());

                            if (editDistancePos < 15 || editDistanceNeg < 15 ||  Math.abs(editDistancePos-editDistanceNeg)>10) {
                                if (editDistancePos > editDistanceNeg) {
                                    s.append("-").append("\n");
//                                    s.append(getFullComplement(reverse(probeSeq))).append("\n");
                                    //                        System.out.println(s.toString());
                                    writer.write(s.toString());
                                } else {
                                    s.append("+").append("\n");
//                                    s.append(probeSeq).append("\n");
                                    //                        System.out.println(s.toString());
                                    writer.write(s.toString());
                                }
                            } else {
                                s.append(partsKey[1]).append("\t?\n");
                                s.append("").append("\n");
                                System.out.println(s.toString());
                                System.out.println("editDistance: "+  editDistancePos + "\t"+ probeSeq);
                                System.out.println("editDistance: "+  editDistanceNeg + "\t" + getFullComplement(reverse(probeSeq)));
                            }

                        } 
//                        else {
//                        s.append("?").append("\n");
//                        s.append("").append("\n");
//                        System.out.println(s.toString());
//
//                        System.out.println(probeSeq);
//                        System.out.println(getFullComplement(reverse(probeSeq)));
//                        }
                    }

                } catch (Exception ex) {
                    Logger.getLogger(WriteExtendedFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        writer.close();
    }

    public static String reverse(String str) {
        return new StringBuilder(str).reverse().toString();
    }

    public static StringBuilder getFullComplement(String longstr) {
        StringBuilder buffer = new StringBuilder(longstr.length());
        for (int i = 0; i < longstr.length(); i++) {
            buffer.append(getComplement(longstr.charAt(i)));
        }
        return buffer;
    }

    public static char getComplement(char x) {
        if (x == 'A') {
            return 'T';
        }
        if (x == 'T') {
            return 'A';
        }
        if (x == 'U') {
            return 'A';
        }
        if (x == 'C') {
            return 'G';
        }
        if (x == 'G') {
            return 'C';
        }
        // inversion genotypes
        if (x == 'N') {
            return 'N';
        }
        if (x == 'I') {
            return 'I';
        }
        return '0';
    }

    private static int calcEditDistance(String targetSeq, String probeSeq) {
        int editDistance = 0;

        for (int i = 0; i < targetSeq.length(); i++) {
            if (targetSeq.charAt(i) != probeSeq.charAt(i)) {
                editDistance++;
            }
        }

        return editDistance;
    }
}
