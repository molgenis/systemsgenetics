/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class MachImputedTransposedDosageToTriTyper {

    public MachImputedTransposedDosageToTriTyper(String location, String outputloc, String headerfile) throws IOException {
        System.out.println("Starting the conversion of MACH imputed dosage");

        Gpio.createDir(location);
        Gpio.createDir(outputloc);

        if (!Gpio.exists(location)) {
            throw new IOException("Could not create directory: " + location);
        }

        if (!Gpio.exists(outputloc)) {
            throw new IOException("Could not create directory: " + outputloc);
        }

        if (headerfile == null || !Gpio.exists(headerfile)) {
            throw new IOException("Could not load headerfile: " + headerfile);
        }

        int numindividuals = 0;
        TextFile individuals = new TextFile(outputloc + "Individuals.txt", TextFile.W);
        TextFile phenotype = new TextFile(outputloc + "PhenotypeInformation.txt", TextFile.W);

            TextFile t = new TextFile(headerfile, TextFile.R);
            String[] elems = t.readLineElemsReturnReference(TextFile.tab);
            for (int i = 7; i < elems.length; i++) {
                individuals.write(elems[i] + "\n");
                phenotype.write(elems[i] + "\tcontrol\tinclude\tunknown\n");
                numindividuals++;
            }
            t.close();

        individuals.close();
        phenotype.close();

        String[] fileList = Gpio.getListOfFiles(location);
        if (fileList == null || fileList.length == 0) {
            System.out.println("Error: could not find any text files or gipped files at your --in location: " + location);
            System.exit(-1);
        }
        boolean[] gzipped = new boolean[fileList.length];
        int f = 0;
        for (String file : fileList) {
            System.out.println("Detected file: " + file);
            if (file.endsWith(".gz")) {
                gzipped[f] = true;
            }
            f++;
        }

        // col 10 >> echte data
        // SNP	CHR	POS	Al1	Al2	Freq1	AVG_MACH_QUAL   RS3-3555

        TextFile snps = new TextFile(outputloc + "SNPs.txt", TextFile.W);
        TextFile snpmappings = new TextFile(outputloc + "SNPMappings.txt", TextFile.W);

        BufferedOutputStream matrixImputedDosage = new BufferedOutputStream(new FileOutputStream(outputloc + "ImputedDosageMatrix.dat"));
        DataOutputStream outDosage = new DataOutputStream(matrixImputedDosage);

        BufferedOutputStream matrixGenotype = new BufferedOutputStream(new FileOutputStream(outputloc + "GenotypeMatrix.dat"));
        DataOutputStream outGenotype = new DataOutputStream(matrixGenotype);

        int numsnps = 0;

        boolean headerparsed = false;
        f = 0;
        for (String file : fileList) {
            System.out.println("Pre-Processing file: " + file);
            TextFile g = new TextFile(location + file, false);

            elems = g.readLineElemsReturnReference(Strings.whitespace);

            try {
                boolean missingvalues = false;
                while (elems != null) {
                    String snp = elems[0];
                    String chr = elems[1];
                    String pos = elems[2];

                    snps.write(snp + "\n");
                    snpmappings.write(chr + "\t" + pos + "\t" + snp + "\tFreq:" + elems[5] + ";QUAL:" + elems[6] + "\n");


                    byte[] dosage = new byte[numindividuals];
                    byte[] b_allele1 = new byte[numindividuals];
                    byte[] b_allele2 = new byte[numindividuals];
                    double dos = 0;

                    String allele1 = elems[3];
                    String allele2 = elems[4];

                    for (int i = 7; i < elems.length; i++) {
                        int ind = i - 7;
                        boolean isMissingGenotype = false;
                        Double d_dosage = null;
                        byte dosageByte = -1;

                        try {
                            d_dosage = Double.parseDouble(elems[i]);
                            isMissingGenotype = false;
                        } catch (Exception e) {
                            System.err.println("Could not parse " + elems[i] + " as double in file " + file + ". It seems there is an error in the format of your genotype input! - We will handle these values as missing genotypes!");
                            missingvalues = true;
                            isMissingGenotype = true;
                        }

                        if (isMissingGenotype) {
                            b_allele1[ind] = 0;
                            b_allele2[ind] = 0;
                            dosageByte = -1;
                        } else {
                            dos = d_dosage;
                            if (d_dosage < 0.5) {
                                b_allele1[ind] = allele2.getBytes()[0];
                                b_allele2[ind] = allele2.getBytes()[0];

                            } else if (d_dosage > 1.5) {
                                b_allele1[ind] = allele1.getBytes()[0];
                                b_allele2[ind] = allele1.getBytes()[0];

                            } else {
                                b_allele1[ind] = allele1.getBytes()[0];
                                b_allele2[ind] = allele2.getBytes()[0];
                            }
                            int dosageInt = (int) Math.round(d_dosage * 100d);
                            if (dosageInt < 0 || dosageInt > 200) {
                                System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snp + "\t" + d_dosage);
                            }

                            dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                        }

                        //}

                        dosage[ind] = dosageByte;


                    }

                    // write output
                    outDosage.write(dosage);
                    outGenotype.write(b_allele1);
                    outGenotype.write(b_allele2);

                    elems = g.readLineElemsReturnReference(Strings.whitespace);
                    numsnps++;
                }

                if (missingvalues) {
                    System.err.println("WARNING: your file " + location + file + " has missing values for a number of SNPs.");
                }
            } catch (Exception e) {
                e.printStackTrace();
                System.err.println("There is an error in your genotype input data in file: " + file + ". Was it formatted correctly?");
                System.exit(0);
            }



            g.close();

            f++;

        }

        snps.close();
        snpmappings.close();


        System.out.println("Found " + numindividuals + " samples and " + numsnps + " snps");

        outGenotype.close();
        outDosage.close();

        matrixImputedDosage.close();
        matrixGenotype.close();
    }
}
