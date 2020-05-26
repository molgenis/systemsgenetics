/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class InversionCallsToTriTyper {

    public void convert(String inFile, String outDir) throws IOException {

        if (inFile == null || inFile.trim().length() == 0 || !Gpio.exists(inFile)) {
            throw new IOException("Could not find file: " + inFile);
        }

        if (outDir == null || outDir.trim().length() == 0) {
            throw new IOException("Invalid directory name: " + inFile);
        }

        outDir = Gpio.formatAsDirectory(outDir);
        Gpio.createDir(outDir);

        TextFile tf = new TextFile(inFile, TextFile.R);

        // header
        String[] inversions = tf.readLineElems(TextFile.tab);

        // now determine the number of individuals..
        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<String> individuals = new ArrayList<String>();
        while (elems != null) {
            if (elems.length >= inversions.length) {
                individuals.add(elems[0]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        System.out.println("Detected " + inversions.length + " inversions for " + individuals.size() + " individuals.");

        TextFile snpOut = new TextFile(outDir + "SNPs.txt", TextFile.W);
        TextFile snpOutMapping = new TextFile(outDir + "SNPMappings.txt", TextFile.W);
        for (int i = 0; i < inversions.length; i++) {
            snpOut.writeln(inversions[i]);
            snpOutMapping.writeln("1\t1\t" + inversions[i]);
        }
        snpOutMapping.close();
        snpOut.close();

        TextFile indOut = new TextFile(outDir + "Individuals.txt", TextFile.W);
        TextFile indOutPhe = new TextFile(outDir + "PhenotypeInformation.txt", TextFile.W);
        for (int i = 0; i < individuals.size(); i++) {
            indOut.writeln(individuals.get(i));
            indOutPhe.writeln(individuals.get(i) + "\tunknown\tinclude\tunknown");
        }
        indOutPhe.close();
        indOut.close();


        WGAFileMatrixGenotype matrix = new WGAFileMatrixGenotype(inversions.length, individuals.size(), new File(outDir + "GenotypeMatrix.dat"), false);

        // parse the file again
        tf.open();
        tf.readLine();
        elems = tf.readLineElems(TextFile.tab);

        int indCtr = 0;
        while (elems != null) {
            if (elems.length >= inversions.length) {
                for (int i = 1; i < elems.length; i++) {
                    String genotype = elems[i];
                    String[] genotypeElems = genotype.split("/");

                    if (genotypeElems[0].equals("NI")) {
                        genotypeElems[0] = "N";
                    }
                    if (genotypeElems[1].equals("NI")) {
                        genotypeElems[1] = "N";
                    }

                    byte allele1 = BaseAnnot.toByte(genotypeElems[0]);
                    byte allele2 = BaseAnnot.toByte(genotypeElems[1]);

                    matrix.setAllele1(i - 1, indCtr, allele1);
                    matrix.setAllele2(i - 1, indCtr, allele2);
                }
                indCtr++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        matrix.close();
    }
}
