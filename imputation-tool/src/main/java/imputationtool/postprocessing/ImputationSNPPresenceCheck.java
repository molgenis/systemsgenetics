/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputationtool.postprocessing;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ImputationSNPPresenceCheck {

    /**
     * @param args the command line arguments
     */
    public void run(String[] args) throws IOException {
        // TODO code application logic here


        if (args.length < 5) {
            System.out.println("snps snpmappings pedmapreferenceloc pedmapdatasetloc dosagefileloc");
        } else {
            ImputationSNPPresenceCheck m = new ImputationSNPPresenceCheck();
            // load the list of SNPs to check
            System.out.println("Loading SNPs: " + args[1]);
            m.loadSNPs(args[0]);

            System.out.println("Loading SNP mappings: " + args[1]);
            // load the snp mappings
            m.loadSNPMappings(args[1]);

            System.out.println("Checking reference exclusion files: " + args[2]);
            // check reference exclusion files
            m.checkReferenceExclusionFiles(args[2], true);

            // check reference map files
            System.out.println("Checking reference map files: " + args[2]);
            m.checkReferenceMapFiles(args[2], true);

//	    System.out.println("Checking dataset exclusion files: "+args[3]);
//	    // check reference exclusion files
//	    m.checkReferenceExclusionFiles(args[2], false);
//
//	    // check reference map files
//	    System.out.println("Checking dataset map files: "+args[3]);
//	    m.checkReferenceMapFiles(args[3], false);

//	    // check dosage files
//	    m.checkDosageFiles(args[4]);

            // output
            m.output();
        }

        System.exit(0);

    }
    private HashMap<String, Integer> snpMappings;
    private ArrayList<String> snps;
    private HashMap<Integer, HashSet<String>> snpsPerChr;
    private HashMap<String, String> reasonExcluded;

    public ImputationSNPPresenceCheck() {
        snps = new ArrayList<String>();
        snpMappings = new HashMap<String, Integer>();
        snpsPerChr = new HashMap<Integer, HashSet<String>>();
        reasonExcluded = new HashMap<String, String>();
    }

    private void loadSNPs(String string) throws IOException {

        TextFile in = new TextFile(string, TextFile.R);
        String line = "";

        while ((line = in.readLine()) != null) {
            snps.add(line);
        }
        in.close();


    }

    private void loadSNPMappings(String string) throws IOException {

        TextFile in = new TextFile(string, TextFile.R);
        String line = "";

        while ((line = in.readLine()) != null) {
            String[] elems = line.split("\t");
            String snp = elems[2];
            Integer chr = null;
            try {
                chr = Integer.parseInt(elems[0]);

            } catch (NumberFormatException e) {
            }

            if (chr != null) {
                snpMappings.put(snp, chr);
            }

        }
        in.close();


        for (String s : snps) {
            Integer chr = snpMappings.get(s);
            HashSet<String> v = snpsPerChr.get(chr);
            if (v == null) {
                v = new HashSet<String>();

            }

            v.add(s);
            snpsPerChr.put(chr, v);

        }

    }

    private void checkReferenceExclusionFiles(String string, boolean reference) throws IOException {

        for (int chr = 0; chr < 23; chr++) {
            HashSet<String> v = snpsPerChr.get(chr);
            if (v != null) {
                TextFile in = new TextFile(string + "chr" + chr + ".excludedsnps.txt", TextFile.R);

                String line = "";
                while ((line = in.readLine()) != null) {
                    String[] elems = line.split("\t");
                    String snp = elems[0];
                    if (v.contains(snp)) {
                        String reason = reasonExcluded.get(snp);
                        if (reason == null) {
                            reason = "";
                        }
                        if (reference) {
                            reason += "\tReference: " + elems[2];
                        } else {
                            reason += "\tDataset: " + elems[2];
                        }
                        reasonExcluded.put(snp, reason);
                    }
                }
                in.close();
            }
        }
    }

    private void checkReferenceMapFiles(String string, boolean reference) throws IOException {

        for (int chr = 0; chr < 23; chr++) {
            HashSet<String> v = snpsPerChr.get(chr);
            if (v != null) {
                TextFile in = new TextFile(string + "chr" + chr + ".map", TextFile.R);

                String line = "";
                while ((line = in.readLine()) != null) {
                    String[] elems = line.split("\t");
                    String snp = elems[1];
                    if (v.contains(snp)) {
                        String reason = reasonExcluded.get(snp);
                        if (reason == null) {
                            reason = "";
                        }
                        if (reference) {
                            reason += "\tReference: in map file";
                        } else {
                            reason += "\tDataset: in map file";
                        }

                        reasonExcluded.put(snp, reason);
                    }

                }
                in.close();

            }
        }
    }

    private void checkDosageFiles(String string) throws IOException {
        for (int chr = 0; chr < 23; chr++) {
            HashSet<String> v = snpsPerChr.get(chr);
            if (v != null) {
                TextFile in = new TextFile(string + "ImputedGenotypeDosageFormatPLINK-Chr" + chr + ".dose", TextFile.R);

                String line = "";
                while ((line = in.readLine()) != null) {
                    String[] elems = line.split("\t");
                    String snp = elems[0];
                    if (v.contains(snp)) {
                        String reason = "" + reasonExcluded.get(snp);
                        reason += "\tIsIncludedInMapFile";
                        reasonExcluded.put(snp, reason);
                    }

                }
                in.close();

            }
        }


    }

    private void output() {
        for (String s : snps) {
            String reason = reasonExcluded.get(s);
            System.out.println(s + "\t" + snpMappings.get(s) + reason);
        }
    }
}
