/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MarcJan
 */
public class ConverteMappingAndSnpFile {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\t");
    private static Pattern SPLIT_ON_COLON = Pattern.compile(":");

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

        String snpMappings = "E:\\LLD_HRC_IMPUTATION\\SNPMappings.txt";
        String snps = "E:\\LLD_HRC_IMPUTATION\\SNPs.txt";

        String bedFolder = "E:\\dbSNP_147\\";

        String snpsOut = "E:\\LLD_HRC_IMPUTATION\\SNPs2.txt";
        String snpMappingsOut = "E:\\LLD_HRC_IMPUTATION\\SNPMappings2.txt";

        System.out.println("Processing original SNP file: ");
        LinkedHashSet<String> OrderingEst = readSnpOrderFile(snps);
        System.out.println("Done found: " + OrderingEst.size() + " SNPs");

        System.out.println("Processing original mapping file: ");
        HashMap<String, String> mappingEst = readSnpMappings(snpMappings, false, (OrderingEst.size()));
        System.out.println("Done found: " + mappingEst.size() + " mappings");

        System.out.println("Processing dbSNP bed folder: ");
        HashMap<String, String> mappingUmcg = readSnpMappingsBed(bedFolder, mappingEst, mappingEst.size());
        System.out.println("Done found: " + mappingUmcg.size() + " mappings");

        System.out.println("Remapping original data.");
        remapPositionsAndSnps(OrderingEst, mappingEst, mappingUmcg, snpsOut, snpMappingsOut);

    }

    private static HashMap<String, String> readSnpMappings(String file1, boolean byPosition, int initialSize) {
        HashMap<String, String> snpMapping = new HashMap<String, String>((int) Math.ceil(initialSize / 0.75));

        try {
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(file1)));
            String str;
            //Pair<String, String> t = new Pair<String, String>();
            while ((str = in.readLine()) != null) {
                String[] parts = SPLIT_ON_TAB.split(str);
                if (byPosition) {
                    if (!snpMapping.containsKey((parts[0] + ":" + parts[1]))) {
                        snpMapping.put((parts[0] + ":" + parts[1]), parts[2]);
                    }
                } else if (!snpMapping.containsKey((parts[2]))) {
                    snpMapping.put(parts[2], (parts[0] + ":" + parts[1]));
                }
            }
            in.close();

        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }
        return (snpMapping);
    }

    private static HashMap<String, String> readSnpMappingsBed(String file1, HashMap<String, String> inputMap, int initialSize) {
        HashMap<String, String> snpMapping = new HashMap<String, String>((int) Math.ceil(initialSize / 0.75));
        
        HashSet<String> interestCombinations = new HashSet<>();
        
        for(Entry<String, String> snp : inputMap.entrySet()){
            interestCombinations.add(snp.getValue());
        }

        File folderIn = new File(file1);
        for (File f : folderIn.listFiles()) {
            String file = f.getAbsolutePath();
            if (file.endsWith(".bed") || file.endsWith(".bed.gz") ) {
                
                String chr = f.getName().split("_")[2].replace(".gz", "").replace(".bed", "");
//                System.out.println(chr);
                try {
                    TextFile in = new TextFile(file, TextFile.R);
                    String str= in.readLine();
                    while ((str = in.readLine()) != null) {
                        String[] parts = SPLIT_ON_TAB.split(str);
//                        parts[0] = parts[0].replace("chr", "");
                        if(interestCombinations.contains(chr+":"+parts[2])){
                            snpMapping.put(chr+":"+parts[2], parts[3]);
                        }
                    }
                    in.close();
                } catch (IOException e) {
                    System.out.println(e.getMessage());
                    System.exit(-1);
                }
            }
        }
        return (snpMapping);
    }

    private static LinkedHashSet<String> readSnpOrderFile(String file1) {
        LinkedHashSet<String> snpOrdering = new LinkedHashSet<String>();
        try {
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(file1)));
            String str;
            int i = 0;
            while ((str = in.readLine()) != null) {
                if (!snpOrdering.contains((str))) {
                    snpOrdering.add(str);
                } else {
                    System.out.println("\tDuplicate entry: " + str);
                }
            }

            in.close();

        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }
        return (snpOrdering);
    }

    private static void remapPositionsAndSnps(LinkedHashSet<String> OrderingEst, HashMap<String, String> mappingEst, HashMap<String, String> mappingUmcg, String outFile, String outFile2) {
        try {
            TextFile out = new TextFile(outFile, TextFile.W);
            TextFile out2 = new TextFile(outFile2, TextFile.W);

            for (String entry : OrderingEst) {
//                System.out.println(entry);
                String snpPos = mappingEst.get(entry);
                String[] t = SPLIT_ON_COLON.split(snpPos);
                if(mappingUmcg.containsKey(snpPos)){
                    out.writeln(mappingUmcg.get(snpPos));
                    out2.writeln(t[0] + "\t" + t[1] + "\t" + mappingUmcg.get(snpPos));
                }else {
                    System.out.println("Problem: "+ entry);
                    out.writeln(entry);
                    out2.writeln(t[0] + "\t" + t[1] + "\t" + entry);
                }
            }
            out.close();
            out2.close();
        } catch (IOException ex) {
            Logger.getLogger(ConverteMappingAndSnpFile.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
