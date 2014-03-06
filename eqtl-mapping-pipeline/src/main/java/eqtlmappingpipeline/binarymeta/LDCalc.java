/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class LDCalc {

    public void LDCalc(String snpfile, String datasetLoc, String out) throws IOException, Exception {

        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(datasetLoc);



        HashSet<String> uniqueSNPs = new HashSet<String>();

        TextFile tf = new TextFile(snpfile, TextFile.R);
//	tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            uniqueSNPs.add(elems[0].trim());
            elems = tf.readLineElems(TextFile.tab);
        }


        String[] snpsToQuery = new String[uniqueSNPs.size()];
        snpsToQuery = uniqueSNPs.toArray(snpsToQuery);

        SNPLoader loader = ds.createSNPLoader();
        DetermineLD ld = new DetermineLD();
        TextFile outfile = new TextFile(out, TextFile.W);
        outfile.writeln("snpA\tsnpB\tr2\tsnpAChr\tsnpBChr\tDistance");

        for (int i = 0; i < snpsToQuery.length; i++) {
            Integer id = ds.getSnpToSNPId().get(snpsToQuery[i]);


            if (id != -9) {

                SNP snp1 = ds.getSNPObject(id);
                loader.loadGenotypes(snp1);
                loader.loadDosage(snp1);

                if (snp1.getChr() > 0) {


                    for (int j = i + 1; j < snpsToQuery.length; j++) {

                        Integer id2 = ds.getSnpToSNPId().get(snpsToQuery[j]);


                        if (id2 != -9) {
                            SNP snp2 = ds.getSNPObject(id2);


                            if (snp1.getChr() == snp2.getChr()) {
                                loader.loadGenotypes(snp2);
                                loader.loadDosage(snp2);
                                double r2 = ld.getRSquared(snp1, snp2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

                                if (r2 > 0.5) {
                                    String outln = snp1.getName() + "\t" + snp2.getName() + "\t" + r2 + "\t" + snp1.getChr() + "\t" + snp2.getChr() + "\t" + Math.abs(snp1.getChrPos() - snp2.getChrPos());
                                    System.out.println(outln);

                                    outfile.writeln(outln);
                                }
                                snp2.clearGenotypes();

                            }

                        }

                    }
                }

                snp1.clearGenotypes();
            }


        }

        outfile.close();

    }
}
