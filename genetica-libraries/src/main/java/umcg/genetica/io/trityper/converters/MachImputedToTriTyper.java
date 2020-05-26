/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;

/**
 *
 * @author lude
 */
public class MachImputedToTriTyper {

    public MachImputedToTriTyper(String dir, String outputDir) throws IOException {
        Gpio.createDir(outputDir);
        HashMap<String, Integer> hashSNPs = new HashMap<String, Integer>();
        ArrayList<String> vecSNPs = new ArrayList<String>();

        //Load all imputed SNPs:
        System.out.println("Determining number of unique SNPs:");
        TextFile outSNPs = new TextFile(outputDir + "/SNPs.txt", TextFile.W);

        String[] fileList = Gpio.getListOfFiles(dir, ".mlinfo");
        TextFile in = null;
        String[] str = null;
        for (String filename : fileList) {
            System.out.println("Processing file:\t" + filename);

            in = new TextFile(filename, TextFile.R);
            str = in.readLineElemsReturnReference(TextFile.tab);

            while (str != null) {
                String snp = new String(str[0].getBytes());
                if (!hashSNPs.containsKey(snp)) {
                    hashSNPs.put(snp, vecSNPs.size());
                    vecSNPs.add(snp);
                }
                outSNPs.write(snp + "\n");
                str = in.readLineElemsReturnReference(TextFile.tab);
            }
            in.close();
            in = null;
        }
        outSNPs.close();

        //Load all individuals:
        System.out.println("Determining unique number of individuals:");
        HashMap<String, Integer> hashInds = new HashMap<String, Integer>();
        ArrayList<String> vecInds = new ArrayList<String>();
        fileList = Gpio.getListOfFiles(dir, ".mlgeno");
        for (String filename : fileList) {
            in = new TextFile(filename, TextFile.R);
            str = in.readLineElemsReturnReference(TextFile.space);
            while (str != null) {
                System.out.println(str[0]);
                if (!hashInds.containsKey(str[0])) {
                    hashInds.put(str[0], vecInds.size());
                    vecInds.add(str[0]);
                }
                str = in.readLineElemsReturnReference(TextFile.space);
            }
        }

        fileList = Gpio.getListOfFiles(dir, ".mldose");
        for (String filename : fileList) {
            in = new TextFile(filename, TextFile.R);
            str = in.readLineElemsReturnReference(TextFile.space);
            while (str != null) {
                System.out.println(str[0]);
                if (!hashInds.containsKey(str[0])) {
                    hashInds.put(str[0], vecInds.size());
                    vecInds.add(str[0]);
                }
                str = in.readLineElemsReturnReference(TextFile.space);
            }
            in.close();
        }

        TextFile outInds = new TextFile(outputDir + "/Individuals.txt", TextFile.W);
        for (String ind : vecInds) {
            String individual = ind.split(">")[1];
            outInds.write(individual + "\n");
        }
        outInds.close();

        System.out.println("Number of unique SNPs:\t" + vecSNPs.size());
        System.out.println("Number of unique individuals:\t" + vecInds.size());

        //Process genotypes:
        System.out.println("Importing genotypes:");
        File fileGenotypeMatrix = new File(outputDir + "/GenotypeMatrix.dat");
        WGAFileMatrixGenotype matrixGenotype = new WGAFileMatrixGenotype(vecSNPs.size(), vecInds.size(), fileGenotypeMatrix, false);

        fileList = Gpio.getListOfFiles(dir, ".mlinfo");

        for (String filename : fileList) {
            System.out.println("Processing " + filename);
            int[] snpIDArray = new int[1000000];
            in = new TextFile(filename, TextFile.R);
            str = in.readLineElemsReturnReference(TextFile.tab);
            int counter = 0;
            while (str != null) {
                String snp = new String(str[0].getBytes());
                snpIDArray[counter] = hashSNPs.get(snp);
                counter++;
                str = in.readLineElemsReturnReference(TextFile.tab);
            }
            in.close();

            String genofilename = filename.substring(0, filename.length() - 7);
            genofilename += ".mlgeno";
            in = new TextFile(genofilename, TextFile.R);
            str = in.readLineElemsReturnReference(TextFile.space);

            while (str != null) {
                int indID = hashInds.get(str[0]);
                for (int c = 2; c < str.length; c++) {
                    int snpID = snpIDArray[c - 2]; // A C 
                    byte[] allele1 = new byte[1];
                    allele1[0] = str[c].getBytes()[0];
                    byte[] allele2 = new byte[1];
                    allele2[0] = str[c].getBytes()[2];
                    matrixGenotype.setAllele1(snpID, indID, allele1);
                    matrixGenotype.setAllele2(snpID, indID, allele2);
                }
                str = in.readLineElemsReturnReference(TextFile.space);
            }
            in.close();
        }

        matrixGenotype.close();

        // Process dosage values
        File fileImputedDosageMatrix = new File(outputDir + "/ImputedDosageMatrix.dat");
        WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(vecSNPs.size(), vecInds.size(), fileImputedDosageMatrix, false);

        fileList = Gpio.getListOfFiles(dir, ".mlinfo");
        System.out.println("Importing dosage information:");
        for (String filename : fileList) {
            System.out.println("Processing " + filename);
            int[] snpIDArray = new int[1000000];
            in = new TextFile(filename, TextFile.R);
            str = in.readLineElemsReturnReference(TextFile.tab);
            int counter = 0;
            while (str != null) {
                String snp = new String(str[0].getBytes());
                snpIDArray[counter] = hashSNPs.get(snp);
                counter++;
                str = in.readLineElemsReturnReference(TextFile.tab);
            }
            in.close();

            String genofilename = filename.substring(0, filename.length() - 7);
            genofilename += ".mldose";
            in = new TextFile(genofilename, TextFile.R);
            str = in.readLineElemsReturnReference(TextFile.space);

            while (str != null) {

                int indID = hashInds.get(str[0]);
                for (int c = 2; c < str.length; c++) {
                    int snpID = snpIDArray[c - 2];
                    int dosageInt = (int) Math.round(Double.parseDouble(str[c]) * 100d);
                    if (dosageInt < 0 || dosageInt > 200) {
                        System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpID + "\t" + str[c]);
                    }
                    byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                    byte[] dosage = new byte[1];
                    dosage[0] = (byte) dosageByte;
                    matrixImputedDosage.setDosage(snpID, indID, dosage);
                }
                System.out.println("Sample\t" + str[0] + "\thas been processed.");
            }

            in.close();

        }

        matrixImputedDosage.close();

        System.out.println("\n\n\n");
        System.out.println("MACH Imputed data has been imported.");
        System.out.println("\n\n\n");

    }
}
