/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import java.io.IOException;
import java.util.zip.DataFormatException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultProbe;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class BinaryResultFileConverter {

    public void export(String indir, String prefix, String outdir, int nrPermsToExport) throws IOException {
        System.out.println("About to export binary result data to text");
        System.out.println("In:  " + indir);
        System.out.println("Pre: " + prefix);
        System.out.println("Out: " + outdir);

        System.out.println("");
        if (!outdir.endsWith("/")) {
            outdir += "/";
        }

        Gpio.createDir(outdir);

        for (int perm = 0; perm < nrPermsToExport; perm++) {
            BinaryResultDataset ds = new BinaryResultDataset(indir, prefix, perm);

            BinaryResultSNP[] snps = ds.getSnps();
            BinaryResultProbe[] probes = ds.getProbes();

            int snpCtr = 0;

            String desc = "RealData";
            if (perm > 0) {
                desc = "Permutation" + perm;
            }

            TextFile outProbes = new TextFile(outdir + "Probes." + desc + ".txt.gz", TextFile.W);

            outProbes.writeln("Name\tChr\tChrStart\tChrEnd\tMidpoint\tAnnotation");
            for (int probeId = 0; probeId < probes.length; probeId++) {
                StringBuilder p = new StringBuilder();
                BinaryResultProbe probe = probes[probeId];
                p.append(probe.getName()).append("\t");
                p.append(probe.getChr()).append("\t");
                p.append(probe.getStart()).append("\t");
                p.append(probe.getStop()).append("\t");
                p.append(probe.getMidpoint()).append("\t");
                p.append(probe.getAnnotation());
                outProbes.writeln(p.toString());
            }
            outProbes.close();
            
            ProgressBar pb = new ProgressBar(snps.length, "Exporting " + desc);
            TextFile outZ = new TextFile(outdir + "ZScores." + desc + ".txt.gz", TextFile.W);
            outZ.writeln("Z\tSNP\tProbe\tType");
            TextFile outSNP = new TextFile(outdir + "SNPs." + desc + ".txt.gz", TextFile.W);
            outSNP.writeln("SNP\tChr\tChrPos\tMAF\tCR\tHWEP\tAllele1/Allele2\tMinorAllele/AlleleAssessed\tNumSamples");
            for (BinaryResultSNP snp : snps) {
                try {
                    Float[] zscores = ds.readSNPZScores(snp);
                    for (int probe = 0; probe < zscores.length; probe++) {
                        if (zscores[probe] != null) {
                            StringBuilder s = new StringBuilder();
                            StringBuilder z = new StringBuilder();

                            z.append(zscores[probe]).append("\t");
                            z.append(snp.getName()).append("\t");
                            z.append(probes[probe].getName()).append("\t");

                            if (!snp.getChr().equals(probes[probe].getChr())) {
                                z.append("Trans");
                            } else {
                                int distance = Math.abs(snp.getChrpos() - probes[probe].getMidpoint());
                                if (distance > 250000) {
                                    z.append("Trans");
                                } else {
                                    z.append("Cis;" + distance);
                                }
                            }

                            s.append(snp.getName()).append("\t");
                            s.append(snp.getChr()).append("\t");
                            s.append(snp.getChrpos()).append("\t");
                            s.append(snp.getMaf()).append("\t");
                            s.append(snp.getCr()).append("\t");
                            s.append(snp.getHwe()).append("\t");
                            s.append(BaseAnnot.toString(snp.getAlleles()[0])).append("/").append(BaseAnnot.toString(snp.getAlleles()[1])).append("\t");
                            s.append(BaseAnnot.toString(snp.getMinorAllele())).append("/").append(BaseAnnot.toString(snp.getAssessedAllele())).append("\t");
                            s.append(snp.getNumsamples()).append("\t");

                            outZ.writeln(z.toString());
                            outSNP.writeln(s.toString());
                        }
                    }
                } catch (DataFormatException e) {

                    System.err.println("There is an error with your binary data!");
                    e.printStackTrace();
                }
                snpCtr++;

                pb.set(snpCtr);
            }

            outSNP.close();


            pb.close();
            outZ.close();
            ds.close();
        }

    }
}
