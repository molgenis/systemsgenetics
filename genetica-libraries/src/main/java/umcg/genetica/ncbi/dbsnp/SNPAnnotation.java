/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.ncbi.dbsnp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class SNPAnnotation implements java.io.Serializable {

    public HashMap<Integer, Integer> rsToChrPos = new HashMap<Integer, Integer>();
    private HashSet<String> contigs;
    private HashMap<String, String> contigChr;
    private HashMap<String, String> contigOrientation;
    public HashMap<String, ArrayList<Integer>> chrToSNP = new HashMap<String, ArrayList<Integer>>();
    public HashMap<String, HashSet<Integer>> chrToUniquePositions = new HashMap<String, HashSet<Integer>>();
    Integer[][] rsToRs = null;

    public SNPAnnotation(String contigInfoLoc, String snpInfoLoc) throws IOException {
        loadContigInfo(contigInfoLoc);
        loadSNPAnnotation(snpInfoLoc);
    }

    public SNPAnnotation(String contigInfoLoc, String snpInfoLoc, String ref) throws IOException {
        loadContigInfo(contigInfoLoc, ref);
        loadSNPAnnotation(snpInfoLoc);
    }

    private void loadContigInfo(String location) throws IOException {
        System.out.println("Loading contigs from: " + location);
        contigs = new HashSet<String>();
        contigChr = new HashMap<String, String>();
        contigOrientation = new HashMap<String, String>();

        TextFile in = new TextFile(location, TextFile.R);

        String[] elems = in.readLineElemsReturnObjects(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            contigs.add(elems[0]);
            contigChr.put(elems[0], elems[5]);
            System.out.println(elems[0] + "\t" + elems[5]);
            contigOrientation.put(elems[0], elems[8]);
            elems = in.readLineElemsReturnObjects(TextFile.tab);
        }
        in.close();

    }

    private void loadContigInfo(String location, String buildString) throws IOException {
        System.out.println("Loading contigs from: " + location);
        contigs = new HashSet<String>();
        contigChr = new HashMap<String, String>();
        contigOrientation = new HashMap<String, String>();

        TextFile in = new TextFile(location, TextFile.R);

        String[] elems = in.readLineElemsReturnObjects(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            if (elems[12].equals(buildString)) {
//		if(elems[1].equals()){
//
//		}


                contigs.add(elems[0]);
                contigChr.put(elems[0], elems[5]);
                System.out.println(elems[0] + "\t" + elems[5]);
                contigOrientation.put(elems[0], elems[8]);
            }


            elems = in.readLineElemsReturnObjects(TextFile.tab);
        }
        in.close();

    }

    private void loadSNPAnnotation(String loc) throws IOException {
        System.out.println("Loading SNP Annotation from: " + loc);

        int linecounter = 0;
        TextFile in = new TextFile(loc, TextFile.R);


        int ctr = 0;
        String[] elems = in.readLineElemsReturnObjects(TextFile.tab);
        while (elems != null) {

            String contig = elems[2];
            if (contigs.contains(contig)) {
                String chrStr = contigChr.get(contig);
                Integer chrPos = null;
                try {
                    chrPos = Integer.parseInt(elems[10]);
                } catch (Exception e) {
//		    System.out.println("Cannot parse chr pos: "+elems[10]+"\t"+Strings.concat(elems, Strings.tab));
                }

                try {

                    Integer RSnum = Integer.parseInt(elems[1]);
                    ArrayList<Integer> chrSNPs = chrToSNP.get(chrStr);
                    HashSet<Integer> uniquePositions = chrToUniquePositions.get(chrStr);

                    if (uniquePositions == null) {
                        uniquePositions = new HashSet<Integer>();
                    }

                    if (chrSNPs == null) {
                        chrSNPs = new ArrayList<Integer>();
                    }

                    chrSNPs.add(RSnum);

                    if (rsToChrPos.get(RSnum) != null) {
                        rsToChrPos.put(RSnum, -1);
                    } else {

                        rsToChrPos.put(RSnum, chrPos);
                    }

                    if (chrPos != null) {
                        uniquePositions.add(chrPos);
                    }

                    chrToUniquePositions.put(chrStr, uniquePositions);
                    chrToSNP.put(chrStr, chrSNPs);
                } catch (java.lang.ArrayIndexOutOfBoundsException u) {
                    System.out.println(Strings.concat(elems, Strings.tab));
                }
            }

//		if (linecounter % 100000 == 0) {
//		    System.out.println(linecounter + " lines parsed");
//		}
            linecounter++;
            elems = in.readLineElemsReturnObjects(TextFile.tab);
        }
        in.close();

        int numTotalPositions = 0;
        Set< Entry<String, HashSet<Integer>>> it = chrToUniquePositions.entrySet();
        for (Entry<String, HashSet<Integer>> e : it) {
            if (e.getValue() != null) {
                numTotalPositions += e.getValue().size();
            }
        }

        System.out.println(numTotalPositions + " positions loaded");

    }

    public HashSet<String> getContigs() {
        return contigs;
    }
}
