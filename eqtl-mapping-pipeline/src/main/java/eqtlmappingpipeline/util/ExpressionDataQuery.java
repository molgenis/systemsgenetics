/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.ProbeAnnotation;
import umcg.genetica.io.trityper.TriTyperSNPMappings;

/**
 *
 * @author harmjan
 */
public class ExpressionDataQuery {

    public static void getProbesAroundSNPs(String probeAnnotationFile, String snplistfile, String snpAnnotationFile, int eqtlwindow, String outfile) throws IOException {
        HashSet<String> querySNPs = null;

        TextFile snpFile = new TextFile(snplistfile, TextFile.R);
        querySNPs = (HashSet<String>) snpFile.readAsSet(0, TextFile.tab);
        snpFile.close();
        System.out.println(querySNPs.size() + " SNPs loaded " + snplistfile);

        TriTyperSNPMappings snpMappings = new TriTyperSNPMappings(snpAnnotationFile, querySNPs);



        ProbeAnnotation pb = new ProbeAnnotation(probeAnnotationFile);

        HashSet<String> selectedProbes = new HashSet<String>();
        HashSet<String> snpsWithProbes = new HashSet<String>();
        int nrProbes = pb.getProbes().length;
        System.out.println(nrProbes + " probes loaded.");
        for (String snp : querySNPs) {
            Integer snpId = snpMappings.getSNPId(snp);
            if (snpId != null) {
                // iterate probes...
                short snpchr = (short) snpMappings.getChr(snpId);
                int snpchrpos = snpMappings.getChrPos(snpId);

                System.out.println(snp + "\t" + snpchr + "\t" + snpchrpos);

                for (int p = 0; p < nrProbes; p++) {
                    short chr = pb.getChr()[p];
                    int chrpos = (pb.getChrStart()[p] + pb.getChrEnd()[p]) / 2;

//                    System.out.println("\t"+p+"\t"+pb.getProbes()[p]+"\t"+chr+"\t"+chrpos);
                    if (chr == snpchr) {
//                        System.out.println("Probe: "+p+" has same chr "+chr);

                        int distance = Math.abs(snpchrpos - chrpos);
                        if (distance <= eqtlwindow) {
                            selectedProbes.add(pb.getProbes()[p]);
                            snpsWithProbes.add(snp);
                        }
                    }

                }
            } else {
                System.out.println(snp + " not present in dataset");
            }
        }

        System.out.println(selectedProbes.size() + " probes selected for " + snpsWithProbes.size() + " SNPs");

        TextFile out = new TextFile(outfile, TextFile.W);
        out.writeList(Arrays.asList(selectedProbes.toArray(new String[0])));
        out.close();


    }
}
