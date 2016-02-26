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
public class WriteFastaOfInterestRegion {

    public static void main(String[] args) {
            String referenceFasta = args[0];
            String fileWithPositions = args[1];
            String outputFile = args[2];

//        String referenceFasta = "E:\\OnlineFolders\\BitSync\\UMCG\\ProbeMapping\\human_g1k_v37.fasta";
//        String fileWithPositions = "E:\\OnlineFolders\\BitSync\\UMCG\\ProbeMapping\\Info\\V71_2\\RNA_Seq\\IlluminaMatchedRNASeqProbes.gff";
//        String outputFile = "E:\\OnlineFolders\\BitSync\\UMCG\\ProbeMapping\\Fasta_Urmo_1.txt";

        ReferenceGenomeFasta refGen = null;
        HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> interestStrings = null;

        try {
            System.out.println("Reading reference GTF.");
            interestStrings = readFileWithPositions(fileWithPositions);
        } catch (Exception ex) {
            Logger.getLogger(WriteFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
        }

        try {
            System.out.println("Reading reference FASTA.");
            refGen = new ReferenceGenomeFasta(new File(referenceFasta), ReferenceGenomeFasta.HUMAN_NORMAL_CHR);
        } catch (Exception ex) {
            Logger.getLogger(WriteFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
        }

        if (refGen != null && interestStrings != null) {
            try {
                writeFasta(refGen, interestStrings, outputFile);
            } catch (IOException ex) {
                Logger.getLogger(WriteFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
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
            System.out.println(str);

            if (!chromosomalProbeInformation.containsKey(parts[0])) {
                chromosomalProbeInformation.put(parts[0], new HashMap<String, Triple<String, Integer, Integer>>());
            }
            int start = Integer.parseInt(parts[3]);
            int stop = Integer.parseInt(parts[4]);
            if(stop<start){
               int tmp = stop;
               stop = start;
               start = tmp;
            }
            chromosomalProbeInformation.get(parts[0]).put(parts[8]+"|_"+parts[6], new Triple<String, Integer, Integer>(parts[6], start, stop));
            counter++;
        }
        System.out.println("Number of entries read in: " + counter);
        gffFileReader.close();
        return chromosomalProbeInformation;
    }

    private static void writeFasta(ReferenceGenomeFasta refGen, HashMap<String, HashMap<String, Triple<String, Integer, Integer>>> interestStrings, String outputFile) throws IOException {
        TextFile writer = new TextFile(outputFile, TextFile.W);

        for (Entry<String, HashMap<String, Triple<String, Integer, Integer>>> probesPerChr : interestStrings.entrySet()) {
//            System.out.println(probesPerChr.getKey());

            for (Entry<String, Triple<String, Integer, Integer>> probe : probesPerChr.getValue().entrySet()) {
//                System.out.println("\t"+probe.getKey() + "      "+ probe.getValue().getLeft());
                try {
                    StringBuilder s = new StringBuilder();
                    s.append(">").append(probe.getKey()).append("\n");
                    String probeSeq = refGen.getNucleotides(probesPerChr.getKey(), probe.getValue().getMiddle(), probe.getValue().getRight()).toString();

                    if (probe.getKey().endsWith("|_+")) {
                        s.append(probeSeq).append("\n");
//                        System.out.println(s.toString());
                        writer.write(s.toString());
                    } else if (probe.getKey().endsWith("|_-")) {
                        s.append(WriteExtendedFastaOfInterestRegion.getFullComplement(WriteExtendedFastaOfInterestRegion.reverse(probeSeq))).append("\n");
//                        System.out.println(s.toString());
                        writer.write(s.toString());
                    } else {
                        System.out.println("Error with the directionality of the probe: "+ probe.getKey());
                    }

                } catch (Exception ex) {
                    Logger.getLogger(WriteFastaOfInterestRegion.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        writer.close();
    }
}
