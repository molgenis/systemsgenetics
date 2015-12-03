/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.probemapping;

import gnu.trove.map.hash.THashMap;
import java.awt.TextField;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class reading {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\t");
    private static Pattern SPLIT_ON_SEMICOLON = Pattern.compile(";");
    private static Pattern SPLIT_ON_SEMICOLON2 = Pattern.compile("; ");
    private static Pattern SPLIT_ON_SPACE = Pattern.compile(" ");
    protected static final String ENCODING = "ISO-8859-1";

    /**
     * Read sam files
     * @param folderIn
     * @param fileExtention
     * @return 
     */
    public static HashMap<String, ArrayList<ArrayList<String>>> readInMultipleSamFiles(String folderIn, String fileExtention) {
        HashMap<String, ArrayList<ArrayList<String>>> readAlignments = new HashMap<String, ArrayList<ArrayList<String>>>();

        File file = new File(folderIn);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(fileExtention)) {
                vecFiles.add(files[f]);
            }
        }

        for (int f = 0; f < vecFiles.size(); f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            String id = vecFiles.get(f).toString();
            id = id.replace(currentFile.getParent(), "");

            if (id.contains(".chromosome.")) {
                id = id.split(".chromosome.")[1];
                id = id.replace(".sam", "");
                id = "Chr" + id;
            } else if (id.contains(".nchr.")) {
                id = "Non Chromosomal region";
            } else if (id.contains(".ncrna.")) {
                id = "Non Coding RNA Transcripts";
            } else if (id.contains(".cdna.")) {
                id = "Transcripts";
            }

            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                
                String str;

                while ((str = in.readLine()) != null) {
                    //System.out.println(str);
                    if (!str.startsWith("@")) {
                        String[] parts = SPLIT_ON_TAB.split(str);
                        if (!(parts[1].equals("4") || parts[1].equals("20"))) {
                            ArrayList<ArrayList<String>> lines = new ArrayList<ArrayList<String>>(7);
                            ArrayList<String> line = new ArrayList<String>();

                            line.add(id);
                            line.add(parts[2]);

                            if (parts[1].equals("16")) {
                                line.add("-");
                            } else {
                                line.add("+");
                            }

                            line.add(parts[3]);

                            line.add(String.valueOf(Integer.parseInt(parts[3]) + (parts[9].length() - 1)));
                            line.add(parts[5]);

                            for (int i = 6; i < parts.length; ++i) {
                                if (parts[i].startsWith("NM:")) {
                                    line.add(parts[i].replace("NM:i:", ""));
                                }
                                if (parts[i].startsWith("XA:Z:")) {
                                    String[] parts2 = SPLIT_ON_SEMICOLON.split(parts[i].replace("XA:Z:", ""));

                                    for (String piece : parts2) {
                                        ArrayList<String> t = new ArrayList<String>(7);
                                        t.add(id);
                                        String[] parts3 = piece.split(",");
                                        for (int j = 0; j < parts3.length; ++j) {
                                            if (j == 1) {
                                                if (parts3[j].startsWith("+")) {
                                                    t.add("+");
                                                    t.add(parts3[j].replace("+", ""));
                                                    t.add(String.valueOf(Integer.parseInt(parts3[j].replace("+", "")) + (parts[9].length() - 1)));
                                                } else {
                                                    t.add("-");
                                                    t.add(parts3[j].replace("-", ""));
                                                    t.add(String.valueOf(Integer.parseInt(parts3[j].replace("-", "")) + (parts[9].length() - 1)));
                                                }

                                            }
                                            t.add(parts3[j]);
                                        }

                                        lines.add(t);
                                    }
                                }
                            }

                            lines.add(line);

                            if (readAlignments.containsKey(parts[0])) {
                                readAlignments.get(parts[0]).addAll(lines);
                            } else {
                                readAlignments.put(parts[0], lines);
                            }
                        }
                    }

                }

            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        return (readAlignments);
    }

    /**
     * Read sam files
     * if there are more than maxDiff differences in the read the row is skipped
     * @param folderIn
     * @param fileExtention
     * @param maxDiff
     * @return 
     */
    public static HashMap<String, ArrayList<ArrayList<String>>> readInMultipleSamFiles2(String folderIn, String fileExtention, int maxDiff) {
        HashMap<String, ArrayList<ArrayList<String>>> readAlignments = new HashMap<String, ArrayList<ArrayList<String>>>();

        File file = new File(folderIn);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(fileExtention)) {
                vecFiles.add(files[f]);
            }
        }

        for (int f = 0; f < vecFiles.size(); f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            String id = vecFiles.get(f).toString();
            id = id.replace(currentFile.getParent(), "");

            if (id.contains(".chromosome.")) {
                id = id.split(".chromosome.")[1];
                id = id.replace(".sam", "");
                id = "Chr" + id;
            } else if (id.contains(".nchr.")) {
                id = "Non Chromosomal region";
            } else if (id.contains(".ncrna.")) {
                id = "Non Coding RNA Transcripts";
            } else if (id.contains(".cdna.")) {
                id = "Transcripts";
            }

            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str;

                while ((str = in.readLine()) != null) {
                    //System.out.println(str);
                    if (!str.startsWith("@")) {
                        String[] parts = SPLIT_ON_TAB.split(str);
                        if (!(parts[1].equals("4") || parts[1].equals("20"))) {
                            ArrayList<ArrayList<String>> lines = new ArrayList<ArrayList<String>>(7);
                            ArrayList<String> line = new ArrayList<String>();

                            line.add(id);
                            line.add(parts[2]);

                            if (parts[1].equals("16")) {
                                line.add("-");
                            } else {
                                line.add("+");
                            }

                            line.add(parts[3]);

                            line.add(String.valueOf(Integer.parseInt(parts[3]) + (parts[9].length() - 1)));
                            line.add(parts[5]);

                            for (int i = 6; i < parts.length; ++i) {
                                if (parts[i].startsWith("NM:")) {
                                    line.add(parts[i].replace("NM:i:", ""));
                                }
                                if (parts[i].startsWith("XA:Z:")) {
                                    String[] parts2 = SPLIT_ON_SEMICOLON.split(parts[i].replace("XA:Z:", ""));

                                    for (String piece : parts2) {
                                        ArrayList<String> t = new ArrayList<String>(7);
                                        t.add(id);
                                        String[] parts3 = piece.split(",");
                                        for (int j = 0; j < parts3.length; ++j) {
                                            if (j == 1) {
                                                if (parts3[j].startsWith("+")) {
                                                    t.add("+");
                                                    t.add(parts3[j].replace("+", ""));
                                                    t.add(String.valueOf(Integer.parseInt(parts3[j].replace("+", "")) + (parts[9].length() - 1)));
                                                } else {
                                                    t.add("-");
                                                    t.add(parts3[j].replace("-", ""));
                                                    t.add(String.valueOf(Integer.parseInt(parts3[j].replace("-", "")) + (parts[9].length() - 1)));
                                                }

                                            }
                                            t.add(parts3[j]);
                                        }
                                        if (Integer.parseInt(t.get(t.size() - 1)) <= maxDiff) {
                                            lines.add(t);
                                        }
                                    }
                                }
                            }
                            if (Integer.parseInt(line.get(line.size() - 1)) <= maxDiff) {
                                lines.add(line);
                            }


                            if (readAlignments.containsKey(parts[0])) {
                                readAlignments.get(parts[0]).addAll(lines);
                            } else {
                                readAlignments.put(parts[0], lines);
                            }
                        }
                    }

                }

            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        return (readAlignments);
    }

    /**
     * Read sam files
     * if there are more than maxDiff differences in the read the row is skipped
     * Degenerate bases are not counted as a mismatch (removed from edit distance)
     * 
     * @param folderIn
     * @param fileExtention
     * @param maxDiff
     * @return 
     */
    public static HashMap<String, ArrayList<ArrayList<String>>> readInMultipleSamFiles2DG(String folderIn, String fileExtention, int maxDiff) {
        HashMap<String, ArrayList<ArrayList<String>>> readAlignments = new HashMap<String, ArrayList<ArrayList<String>>>();

        File file = new File(folderIn);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(fileExtention)) {
                vecFiles.add(files[f]);
            }
        }

        for (int f = 0; f < vecFiles.size(); f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            String id = vecFiles.get(f).toString();
            id = id.replace(currentFile.getParent(), "");

            if (id.contains(".chromosome.")) {
                id = id.split(".chromosome.")[1];
                id = id.replace(".sam", "");
                id = "Chr" + id;
            } else if (id.contains(".nchr.")) {
                id = "Non Chromosomal region";
            } else if (id.contains(".ncrna.")) {
                id = "Non Coding RNA Transcripts";
            } else if (id.contains(".cdna.")) {
                id = "Transcripts";
            }

            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str;

                while ((str = in.readLine()) != null) {
                    //System.out.println(str);
                    if (!str.startsWith("@")) {
                        String[] parts = SPLIT_ON_TAB.split(str);
                        if (!(parts[1].equals("4") || parts[1].equals("20"))) {
                            ArrayList<ArrayList<String>> lines = new ArrayList<ArrayList<String>>(7);
                            ArrayList<String> line = new ArrayList<String>();

                            line.add(id);
                            line.add(parts[2]);

                            if (parts[1].equals("16")) {
                                line.add("-");
                            } else {
                                line.add("+");
                            }

                            line.add(parts[3]);

                            line.add(String.valueOf(Integer.parseInt(parts[3]) + (parts[9].length() - 1)));
                            //System.out.println(parts[9]);
                            int nrNs = getNrNs(parts[9]);
                            line.add(parts[5]);

                            for (int i = 10; i < parts.length; ++i) {
                                if (parts[i].startsWith("NM:")) {
                                    int maxDif = Integer.parseInt(parts[i].replace("NM:i:", ""));
                                    line.add(String.valueOf(maxDif - nrNs));
                                }
                                if (parts[i].startsWith("XA:Z:")) {
                                    System.out.println("Skrewed");
                                }
                            }
                            if (Integer.parseInt(line.get(line.size() - 1)) <= maxDiff) {
                                lines.add(line);
                            }


                            if (readAlignments.containsKey(parts[0])) {
                                readAlignments.get(parts[0]).addAll(lines);
                            } else {
                                readAlignments.put(parts[0], lines);
                            }
                        }
                    }

                }

            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        return (readAlignments);
    }

    /**
     * Return number N's in a sequence.
     * @param string
     * @return 
     */
    private static int getNrNs(String string) {
        char[] characterArray = string.toCharArray();

        int numberN = 0;

        for (int i = 0; i < characterArray.length; ++i) {
            if (characterArray[i] == 'n' || characterArray[i] == 'N') {
                numberN++;
            }
        }

        return (numberN);
    }
    
    /**
     * Read annotation file.
     * General read in
     * @param annotationFile
     * @param storingId
     * @param sizeMap
     * @return 
     */
    public static THashMap<String, THashMap<String, String>> readAnnotationFile(String annotationFile, int storingId, int sizeMap) {
        THashMap<String, THashMap<String, String>> probeInfo = new THashMap<String, THashMap<String, String>>((int) Math.ceil(sizeMap / 0.75));
        int entryId = 0;
        try {
            TextFile in = new TextFile(annotationFile, TextFile.R);
            String str = "";

            str = in.readLine();
            String[] header = SPLIT_ON_TAB.split(str);


            while ((str = in.readLine()) != null) {
                String[] strParts = SPLIT_ON_TAB.split(str);
                THashMap<String, String> t = new THashMap<String, String>((int) Math.ceil(header.length / 0.75));
                for (int i = 0; i < strParts.length; ++i) {
                    if (i != storingId) {
                        t.put(header[i], strParts[i]);
                    }
                }
                if (storingId == -1) {
                    probeInfo.put(String.valueOf(entryId), t);
                    entryId++;
                } else if (storingId == -2) {
                    probeInfo.put(strParts[0]+"-"+strParts[1]+"-"+strParts[28], t);
                    entryId++;
                }else {
                    probeInfo.put(strParts[storingId], t);
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probeInfo);
    }

    /**
     * Read annotation file
     * Supply int key and int for value
     * 
     * @param annotationFile
     * @param firstRowAsHeader
     * @param key
     * @param val
     * @param sizeMap
     * @return 
     */
    public static HashMap<String, String> readAnnotationFileHash(String annotationFile, boolean firstRowAsHeader, int key, int val, int sizeMap) {
        HashMap<String, String> probeInfo = new HashMap<String, String>((int) Math.ceil(sizeMap / 0.75));

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile)), ENCODING), 8096);
            String str = "";
            if (firstRowAsHeader) {
                str = in.readLine();
            }

            while ((str = in.readLine()) != null) {
                String[] strParts = SPLIT_ON_TAB.split(str);
                if(val==-1){
                    probeInfo.put(strParts[key], str);
                } else {
                    probeInfo.put(strParts[key], strParts[val]);
                }
                
            }

        } catch (IOException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probeInfo);
    }
    
    /**
     * Read annotation file
     * Supply int key and int for value
     * 
     * @param annotationFile
     * @param firstRowAsHeader
     * @param key
     * @param val
     * @param sizeMap
     * @return 
     */
    public static HashMap<String, Triple<Integer, Integer, Integer>> readAnnotationFileHashMap(String annotationFile, boolean firstRowAsHeader, int key, int val1, int val2, int val3, int sizeMap) {
        HashMap<String, Triple<Integer, Integer, Integer>> probeInfo = new HashMap<String, Triple<Integer, Integer, Integer>>((int) Math.ceil(sizeMap / 0.75));

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile)), ENCODING), 8096);
            String str = "";
            if (firstRowAsHeader) {
                str = in.readLine();
            }

            while ((str = in.readLine()) != null) {
                String[] strParts = SPLIT_ON_TAB.split(str);
                Triple<Integer, Integer, Integer> tmp;
                if(strParts[val2].contains(":")){
                    strParts[val2] = strParts[val2].split(":")[0];
                } 
                if(strParts[val3].contains(":")){
                    strParts[val3] = strParts[val3].split(":")[1];
                } 
                if(strParts[val1].equals("Y")){
                    strParts[val1] = "24";
                } else if(strParts[val1].equals("X")){
                    strParts[val1] = "23";
                } else {
                    
                }
                tmp = new Triple<Integer, Integer, Integer>(Integer.parseInt(strParts[val1]),Integer.parseInt(strParts[val2]),Integer.parseInt(strParts[val3]));
                probeInfo.put(strParts[key], tmp);
            }

        } catch (IOException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probeInfo);
    }
    
    /**
     * Read annotation file
     * Supply int key and int for value
     * 
     * @param annotationFile
     * @param firstRowAsHeader
     * @param key
     * @param val
     * @param sizeMap
     * @return 
     */
    public static ArrayList< Triple<Integer, Integer, Integer>> readAnnotationFileArrayList(String annotationFile, int col1, int col2, int col3, boolean firstRowAsHeader, int sizeMap) {
        ArrayList< Triple<Integer, Integer, Integer>> probeInfo = new ArrayList< Triple<Integer, Integer, Integer>>((int) Math.ceil(sizeMap / 0.75));

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile)), ENCODING), 8096);
            String str;
            if (firstRowAsHeader) {
                str = in.readLine();
            }

            while ((str = in.readLine()) != null) {
                String[] strParts = SPLIT_ON_TAB.split(str);
                Triple<Integer, Integer, Integer> tmp;
                if(strParts[col1].contains("chr")){
                    strParts[col1] = strParts[col1].replace("chr", "");
                } 
                if(strParts[col1].equalsIgnoreCase("Y")){
                    strParts[col1] = "24";
                } else if(strParts[col1].equalsIgnoreCase("X")){
                    strParts[col1] = "23";
                }
                if(strParts[col1].length()<=2){
                    tmp = new Triple<Integer, Integer, Integer>(Integer.parseInt(strParts[col1]),Integer.parseInt(strParts[col2]),Integer.parseInt(strParts[col3]));
                    probeInfo.add(tmp);
                }
                
            }

        } catch (IOException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probeInfo);
    }

    /**
     * read GTF file (gene code information)
     * @param annotationFile
     * @param sizeMap
     * @return 
     */
    public static HashMap<String, String> readGTFAnnotationFileHash(String annotationFile, int sizeMap) {
        HashMap<String, String> probeInfo = new HashMap<String, String>((int) Math.ceil(sizeMap / 0.75));

        try {
            TextFile in = new TextFile(annotationFile, TextFile.R);
            String str = "";

            while ((str = in.readLine()) != null) {
                String[] strParts = SPLIT_ON_TAB.split(str);
                if (strParts.length == 9) {

                    String[] strParts2 = SPLIT_ON_SEMICOLON2.split(strParts[8]);

                    HashMap<String, String> tmpHash = new HashMap<String, String>();
                    for (String tmp : strParts2) {
                        tmp = tmp.replaceAll("\"", "");
//                        System.out.println(tmp);
                        String[] tmpPart = SPLIT_ON_SPACE.split(tmp);
                        if(tmpPart.length==2){
                            tmpHash.put(tmpPart[0], tmpPart[1]);
                        }
                    }

                    String tmp = tmpHash.get("gene_id");

                    tmpHash.put("gene_id", tmp);
                    tmp = tmp.split("\\.")[0];
                    tmpHash.put("gene_id", tmp);

                    tmp = tmpHash.get("transcript_id");
                    tmp = tmp.split("\\.")[0];
                    tmpHash.put("transcript_id", tmp);

                    if (probeInfo.containsKey(tmpHash.get("gene_id"))) {
                        if (!(probeInfo.get(tmpHash.get("gene_id")).equals(tmpHash.get("gene_name")))) {
                            System.out.println(tmpHash.get("gene_id") + "\t" + probeInfo.get(tmpHash.get("gene_id")) + "\t" + tmpHash.get("gene_name"));
                            System.exit(0);
                        }
                    } else {
                        //System.out.println(tmpHash.get("gene_id")+"\t"+tmpHash.get("gene_name"));
                        probeInfo.put(tmpHash.get("gene_id"), tmpHash.get("gene_name"));
                    }

                    if (probeInfo.containsKey(tmpHash.get("transcript_id"))) {
                        if (!(probeInfo.get(tmpHash.get("transcript_id")).equals(tmpHash.get("gene_name")))) {
                            System.out.println(tmpHash.get("transcript_id") + "\t" + probeInfo.get(tmpHash.get("transcript_id")) + "\t" + tmpHash.get("gene_name"));
                            System.exit(0);
                        }
                    } else {
                        //System.out.println(tmpHash.get("gene_id")+"\t"+tmpHash.get("gene_name"));
                        probeInfo.put(tmpHash.get("transcript_id"), tmpHash.get("gene_name"));
                    }
                }
//                probeInfo.put(strParts[key], strParts[val]);
            }

        } catch (IOException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probeInfo);
    }

    /**
     * Read SNP information files
     * @param annotationFileFolder
     * @param minMaf
     * @param firstRowAsHeader
     * @return 
     */
    public static HashMap<String, HashSet<Integer>> readMultipleSNPAnnotationFilesSmall(String annotationFileFolder, double minMaf, boolean firstRowAsHeader) {
        HashMap<String, HashSet<Integer>> snpInfo = new HashMap<String, HashSet<Integer>>();
        File file = new File(annotationFileFolder);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }

        for (int f = 0; f < vecFiles.size(); f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str = "";
                if (firstRowAsHeader) {
                    str = in.readLine();
                }

                String[] header = SPLIT_ON_TAB.split(str);

                while ((str = in.readLine()) != null) {

                    String[] strParts = SPLIT_ON_TAB.split(str);

                    if (!strParts[5].isEmpty()) {
                        //System.out.println(strParts[5]);
                        if (Double.parseDouble(strParts[5]) > minMaf) {
                            if (snpInfo.containsKey(strParts[2])) {
                                snpInfo.get(strParts[2]).add(Integer.parseInt(strParts[3]));

                            } else {
                                HashSet<Integer> locations = new HashSet<Integer>();
                                snpInfo.put(strParts[2], locations);
                            }
                        }
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        return (snpInfo);
    }

    /**
     * Read SNP information files
     * @param annotationFileFolder
     * @param minMaf
     * @param firstRowAsHeader
     * @return 
     */
    public static HashMap<String, LinkedHashMap<Integer, String>> readMultipleSNPAnnotationFiles(String annotationFileFolder, double minMaf, boolean firstRowAsHeader) {
        HashMap<String, LinkedHashMap<Integer, String>> snpInfo = new HashMap<String, LinkedHashMap<Integer, String>>(25);
        File file = new File(annotationFileFolder);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }

        for (int f = 0; f < vecFiles.size(); f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            String currectChr = "";
            HashMap<Integer, String> locations = new HashMap<Integer, String>();
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str = "";
                if (firstRowAsHeader) {
                    str = in.readLine();
                }

                String[] header = SPLIT_ON_TAB.split(str);

                while ((str = in.readLine()) != null) {

                    String[] strParts = SPLIT_ON_TAB.split(str);
                    
                    StringBuilder keys = new StringBuilder(strParts[0]);
                    if(strParts.length>8 && !strParts[8].isEmpty()){
                        keys.append(";").append(strParts[8]);
                    }
                    
                    if (!strParts[5].isEmpty()) {
                        //System.out.println(strParts[5]);
                        if (Double.parseDouble(strParts[5]) >= minMaf) {
                            if (locations.size()>0) {
                                if(locations.containsKey(Integer.parseInt(strParts[3])) && strParts.length>8){
                                    StringBuilder newKeys = new StringBuilder(locations.get(Integer.parseInt(strParts[3])));
                                    newKeys.append(";").append(strParts[8]);
                                    locations.put(Integer.parseInt(strParts[3]), newKeys.toString());
                                } else{
                                    locations.put(Integer.parseInt(strParts[3]), keys.toString());
                                }
                            } else {
                                locations.put(Integer.parseInt(strParts[3]), keys.toString());
                                currectChr = strParts[2];
                            }
                        }
                    }
                }
                
                ArrayList<Integer> keys = new ArrayList<Integer>();
                keys.addAll(locations.keySet());
                Collections.sort(keys);

                LinkedHashMap<Integer, String> locations2 = new LinkedHashMap<Integer, String>((int)(Math.round((double)locations.size() / 0.75)));
                for(Integer i : keys){
                    locations2.put(i, locations.get(i));
                }
                locations = null;
                snpInfo.put(currectChr, locations2);
                
            } catch (IOException e) {
                e.printStackTrace();
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        
        return (snpInfo);
    }
    
    /**
     * Read SNP information files
     * @param annotationFileFolder
     * @param minMaf
     * @param firstRowAsHeader
     * @return 
     */
    public static HashMap<String, HashMap<Integer, String>> readMultipleSNPAnnotationFiles2(String annotationFileFolder, double minMaf, boolean firstRowAsHeader) {
        HashMap<String, HashMap<Integer, String>> snpInfo = new HashMap<String, HashMap<Integer, String>>(25);
        File file = new File(annotationFileFolder);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }

        for (int f = 0; f < vecFiles.size(); f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            String currectChr = "";
            HashMap<Integer, String> locations = new HashMap<Integer, String>();
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str = "";
                if (firstRowAsHeader) {
                    str = in.readLine();
                }

                String[] header = SPLIT_ON_TAB.split(str);

                while ((str = in.readLine()) != null) {

                    String[] strParts = SPLIT_ON_TAB.split(str);
                    
                    StringBuilder keys = new StringBuilder(strParts[0]);
                    if(strParts.length>8 && !strParts[8].isEmpty()){
                        keys.append(";").append(strParts[8]);
                    }
                    
                    if (!strParts[5].isEmpty()) {
                        //System.out.println(strParts[5]);
                        if (Double.parseDouble(strParts[5]) >= minMaf) {
                            if (locations.size()>0) {
                                if(locations.containsKey(Integer.parseInt(strParts[3])) && strParts.length>8){
                                    StringBuilder newKeys = new StringBuilder(locations.get(Integer.parseInt(strParts[3])));
                                    newKeys.append(";").append(strParts[8]);
                                    locations.put(Integer.parseInt(strParts[3]), newKeys.toString());
                                } else{
                                    locations.put(Integer.parseInt(strParts[3]), keys.toString());
                                }
                            } else {
                                locations.put(Integer.parseInt(strParts[3]), keys.toString());
                                currectChr = strParts[2];
                            }
                        }
                    }
                }

                snpInfo.put(currectChr, locations);
                
            } catch (IOException e) {
                e.printStackTrace();
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        
        return (snpInfo);
    }
    
    /**
     * Read SNP information files
     * @param annotationFileFolder
     * @param minMaf
     * @param firstRowAsHeader
     * @return 
     */
    public static HashMap<String, TreeMap<Integer, String>> readMultipleSNPAnnotationFiles3(String annotationFileFolder, double minMaf, boolean firstRowAsHeader) {
        HashMap<String, TreeMap<Integer, String>> snpInfo = new HashMap<String, TreeMap<Integer, String>>(25);
        File file = new File(annotationFileFolder);
        File[] files = file.listFiles();
        ArrayList<File> vecFiles = new ArrayList<File>();
        for (int f = 0; f < files.length; f++) {
            if (files[f].getAbsolutePath().endsWith(".txt")) {
                vecFiles.add(files[f]);
            }
        }

        for (int f = 0; f < vecFiles.size(); f++) {
            File currentFile = vecFiles.get(f);
            System.out.println("Processing:\t" + f + "\t" + currentFile.getAbsolutePath());
            String currectChr = "";
            TreeMap<Integer, String> locations = new TreeMap<Integer, String>();
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(currentFile), ENCODING), 8096);
                String str = "";
                if (firstRowAsHeader) {
                    str = in.readLine();
                }

                String[] header = SPLIT_ON_TAB.split(str);

                while ((str = in.readLine()) != null) {

                    String[] strParts = SPLIT_ON_TAB.split(str);
                    
                    StringBuilder keys = new StringBuilder(strParts[0]);
                    if(strParts.length>8 && !strParts[8].isEmpty()){
                        keys.append(";").append(strParts[8]);
                    }
                    
                    if (!strParts[5].isEmpty()) {
                        //System.out.println(strParts[5]);
                        if (Double.parseDouble(strParts[5]) >= minMaf) {
                            if (locations.size()>0) {
                                if(locations.containsKey(Integer.parseInt(strParts[3])) && strParts.length>8){
                                    StringBuilder newKeys = new StringBuilder(locations.get(Integer.parseInt(strParts[3])));
                                    newKeys.append(";").append(strParts[8]);
                                    locations.put(Integer.parseInt(strParts[3]), newKeys.toString());
                                } else{
                                    locations.put(Integer.parseInt(strParts[3]), keys.toString());
                                }
                            } else {
                                locations.put(Integer.parseInt(strParts[3]), keys.toString());
                                currectChr = strParts[2];
                            }
                        }
                    }
                }
                snpInfo.put(currectChr, locations);
                
            } catch (IOException e) {
                e.printStackTrace();
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        
        return (snpInfo);
    }
    
    /**
     * Read one file into HashSet<String>
     *
     * @param fileWithAnnotation
     * @return Sample annotation
     */
    public static HashSet<String> readFilterHash(String probeFilteringFiles) {
        HashSet<String> probesToBeRemoved = new HashSet<String>();

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(probeFilteringFiles)), ENCODING), 8096);
            String str;

            while ((str = in.readLine()) != null) {
                probesToBeRemoved.add(str);
            }

            in.close();

        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probesToBeRemoved);
    }
    
    /**
     * Read one file into HashSet<String>
     *
     * @param fileWithAnnotation
     * @return Sample annotation
     */
    public static ArrayList<String> readListToArrayList(String probeFilteringFiles) {
        ArrayList<String> probesToBeRemoved = new ArrayList<String>();

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(probeFilteringFiles)), ENCODING), 8096);
            String str;

            while ((str = in.readLine()) != null) {
                probesToBeRemoved.add(str);
            }

            in.close();

        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probesToBeRemoved);
    }
    
    /**
     * Read multiple file into HashSet<String>
     *
     * @param fileWithAnnotation
     * @return Sample annotation
     */
    public static HashSet<String> readFilterHash2(String[] probeFilteringFiles) {
        ArrayList<HashSet<String>> probesToBeRemoved = new ArrayList<HashSet<String>>();
        for (String s : probeFilteringFiles) {
            try {
                BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(s)), ENCODING), 8096);
                String str;
                HashSet<String> tmpPprobesToBeRemoved = new HashSet<String>();
                while ((str = in.readLine()) != null) {
                    tmpPprobesToBeRemoved.add(str);
                }

                in.close();
                probesToBeRemoved.add(tmpPprobesToBeRemoved);
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(-1);
            }
        }
        HashSet<String> finalSet = probesToBeRemoved.get(0);
        for (int i = 1; i < probesToBeRemoved.size(); ++i) {
            finalSet.retainAll(probesToBeRemoved.get(i));
        }

        return (finalSet);
    }
    
    /**
     * Read double matrix file restricting to given rows Eigenvector file / pc
     * file / probe matrix
     *
     * @param eigenVectorFile
     * @return
     */
    public static DoubleMatrixDataset<String, String> readDoubleMatrixFile(String eigenVectorFile, Set<String> rowsToInclude) {

        DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDataset<String, String>();
        try {
            if (rowsToInclude == null) {
                tmp = new DoubleMatrixDataset<String, String>(eigenVectorFile);
            } else {
                tmp = new DoubleMatrixDataset<String, String>(eigenVectorFile, null, rowsToInclude);
            }
        } catch (IOException ex) {
            Logger.getLogger(reading.class.getName()).log(Level.SEVERE, null, ex);
        }

        return (tmp);
    }

    public static HashMap<String, Pair<String, Double>> readMetaAnalysisResults(String metaAnalysisScores, boolean firstRowAsHeader, int key, int val1, int val2, int sizeMap) {
        HashMap<String, Pair<String, Double>> probeInfo = new HashMap<String, Pair<String, Double>>((int) Math.ceil(sizeMap / 0.75));

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(metaAnalysisScores)), ENCODING), 8096);
            String str = "";
            if (firstRowAsHeader) {
                str = in.readLine();
            }

            while ((str = in.readLine()) != null) {
                String[] strParts = SPLIT_ON_TAB.split(str);
                Pair<String, Double> tmp = new Pair<String, Double>(strParts[val1],Double.parseDouble(strParts[val2]));
                probeInfo.put(strParts[key], tmp);
            }

        } catch (IOException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probeInfo);
    }

    public static HashMap<Integer, HashMap<String, Pair<Integer, Integer>>> readAnnotationFileHashMap2(String annotationFile, boolean firstRowAsHeader, int key, int val1, int val2, int val3, int sizeMap) {
        HashMap<Integer, HashMap<String, Pair<Integer, Integer>>> probeInfo = new HashMap<Integer, HashMap<String, Pair<Integer, Integer>>>((int) Math.ceil(sizeMap / 0.75));

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile)), ENCODING), 8096);
            String str = "";
            if (firstRowAsHeader) {
                str = in.readLine();
            }

            while ((str = in.readLine()) != null) {
                String[] strParts = SPLIT_ON_TAB.split(str);
                Pair<Integer, Integer> tmp;
                if(strParts[val2].contains(":")){
                    strParts[val2] = strParts[val2].split(":")[0];
                } 
                if(strParts[val3].contains(":")){
                    strParts[val3] = strParts[val3].split(":")[1];
                } 
                if(strParts[val1].equals("Y")){
                    strParts[val1] = "24";
                } else if(strParts[val1].equals("X")){
                    strParts[val1] = "23";
                } else {
                    
                }
                
                int chr = Integer.parseInt(strParts[val1]);
                tmp = new Pair<Integer, Integer>(Integer.parseInt(strParts[val2]),Integer.parseInt(strParts[val3]));
                
                if(probeInfo.containsKey(chr)){
                    probeInfo.get(chr).put(strParts[key], tmp);
                } else {
                    HashMap<String, Pair<Integer, Integer>> tmpje = new HashMap<String, Pair<Integer, Integer>>();
                    tmpje.put(strParts[key], tmp);
                    probeInfo.put(chr, tmpje);
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            System.exit(-1);
        }

        return (probeInfo);
    }
}
