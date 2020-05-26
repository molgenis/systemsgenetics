/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import umcg.genetica.containers.SortableSNP;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class ProbeSNPMapper {

    /**
     * @param args the command line arguments
     */
    public void mapprobes(String snpmappings, String probeannotation, String probelist) {
        // read conditional file 

        try {
            // read snp position file
            TextFile tf = new TextFile(snpmappings, TextFile.R);
            HashMap<Byte, ArrayList<SortableSNP>> snpsPerChr = new HashMap<Byte, ArrayList<SortableSNP>>();
            String[] elems = tf.readLineElemsReturnReference(TextFile.tab);

            while (elems != null) {
                Byte chr = ChrAnnotation.parseChr(elems[0]);
                Integer chrPos = Integer.parseInt(elems[1]);
                String name = elems[2];
                if (chr > 0) {
                    ArrayList<SortableSNP> snps = snpsPerChr.get(chr);
                    if (snps == null) {
                        snps = new ArrayList<SortableSNP>();
                    }

                    snps.add(new SortableSNP(name, 0, chr, chrPos, SortableSNP.SORTBY.CHRPOS));

                    snpsPerChr.put(chr, snps);
                }

                elems = tf.readLineElemsReturnReference(TextFile.tab);
            }

            tf.close();

            HashMap<String, String> probeToChrPos = new HashMap<String, String>();
            for (byte chr = 0; chr < 26; chr++) {
                ArrayList<SortableSNP> snps = snpsPerChr.get(chr);
                if (snps != null) {
                    java.util.Collections.sort(snps);
                    snpsPerChr.put(chr, snps);
                    System.out.println("chr:\t" + chr + "\tsnps:\t" + snps.size());
                }

            }

            // 
            tf = new TextFile(probeannotation, TextFile.R);
            elems = tf.readLineElemsReturnReference(TextFile.tab);
            while (elems != null) {
                if (!elems[0].equals("-")) {
                    probeToChrPos.put(elems[1], elems[3] + ";" + elems[4] + ":" + elems[5]);
                }
                elems = tf.readLineElemsReturnReference(TextFile.tab);
            }
            tf.close();

            tf = new TextFile(probelist, TextFile.R);
            elems = tf.readLineElemsReturnReference(TextFile.tab);
            while (elems != null) {
                String probe = elems[0];
                String chrpos = probeToChrPos.get(probe);
                if (chrpos != null) {


                    String[] chrposelems = chrpos.split(";");
                    byte chr = ChrAnnotation.parseChr(chrposelems[0]);
                    if (chr == -1) {
                        System.out.println("THERE IS NO ANNOTATION FOR PROBE: " + probe + "\t" + chrpos);
                    } else {
                        chrposelems = chrposelems[1].split(":");

                        Integer chrstart = Integer.parseInt(chrposelems[0]);
                        Integer chrend = Integer.parseInt(chrposelems[1]);

                        if (chrstart.equals(chrend)) {
                            System.out.println("ERROR: probe chr start == probe chr end");
                        } else {
                            ArrayList<SortableSNP> snps = snpsPerChr.get(chr);
                            String snpsInProbe = "";
                            for (int i = 0; i < snps.size(); i++) {
                                if (snps.get(i).chrpos >= chrstart && snps.get(i).chrpos <= chrend) {
                                    snpsInProbe += snps.get(i).name + ",";
                                }
                            }
                            System.out.println(probe + "\t" + chr + "\t" + chrstart + "\t" + chrend + "\t" + snpsInProbe);
                        }

                    }
                } else {
                    System.out.println("THERE IS NO ANNOTATION FOR PROBE: " + probe);
                }
                elems = tf.readLineElemsReturnReference(TextFile.tab);
            }
            tf.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(0);
        }
    }
}
