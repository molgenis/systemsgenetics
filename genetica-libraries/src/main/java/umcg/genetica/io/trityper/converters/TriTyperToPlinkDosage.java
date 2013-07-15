/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class TriTyperToPlinkDosage {

    public void outputDosageInformation(String inDir, String outputDir, String famFileLoc, boolean splitperchromosome) throws IOException {

        String snpfile = "/Data/Projects/NoortjeFesten/2011-08-16-LifeLines-CommonTraitGenetics.txt";

        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(inDir);

        TextFile out = null;
        TextFile[] outPC = null;

        if (splitperchromosome) {
            System.out.println("Splitting per chromosome");
            outPC = new TextFile[23];
            for (int i = 1; i < 23; i++) {
                outPC[i] = new TextFile(outputDir + "/ImputedGenotypeDosageFormatPLINK-Chr" + i + ".dose", TextFile.W);
            }
        } else {
            System.out.println("Generating a single dosage file");
            out = new TextFile(outputDir + "/ImputedGenotypeDosageFormatPLINK.dose", TextFile.W);
        }

        HashSet<String> snpsToInclude = null;
        if (snpfile != null) {
            snpsToInclude = new HashSet<String>();

            BufferedReader in = new BufferedReader(new FileReader(new File(snpfile)));
            String line = "";
            while ((line = in.readLine()) != null) {
                snpsToInclude.add(line.trim());
            }
            in.close();
        }


        int nrInds = ds.getIndividuals().length;
        String[] individuals = ds.getIndividuals();

        int nrSamplesIncluded = 0;
        for (int ind = 0; ind < nrInds; ind++) {
            if (ds.getIsIncluded()[ind]) {
                nrSamplesIncluded++;
            }
        }
        int[] samplesToOutput = new int[nrSamplesIncluded];
        int sampleItr = 0;
        for (int ind = 0; ind < nrInds; ind++) {
            if (ds.getIsIncluded()[ind]) {
                samplesToOutput[sampleItr] = ind;
                sampleItr++;
            }
        }

        HashMap<String, String> sampleToFamId = new HashMap<String, String>();
        if (famFileLoc != null) {
            TextFile famFile = new TextFile(famFileLoc, TextFile.R);


            String line = "";

            while ((line = famFile.readLine()) != null) {
                String[] elems = line.split(" ");
                if (elems.length > 1) {
                    sampleToFamId.put(elems[1], elems[0]);
                }
            }

            famFile.close();
        }

        String header = "SNP\tA1\tA2";
        for (int i = 0; i < nrSamplesIncluded; i++) {
            String sample = individuals[samplesToOutput[i]];
            String familyId = sampleToFamId.get(sample);
            if (familyId == null) {
                familyId = "1";
            }

            header += "\t" + familyId + "\t" + sample;
        }



        if (splitperchromosome) {

            for (int chr = 1; chr < 23; chr++) {
                outPC[chr].write(header + "\n");
            }
        } else {
            out.write(header + "\n");
        }



        long counter = 0;

        String[] SNPs = ds.getSNPs();
        int numSNPs = SNPs.length;

        SNPLoader loader = ds.createSNPLoader();

        java.text.DecimalFormat df = new java.text.DecimalFormat("0.00", new java.text.DecimalFormatSymbols(java.util.Locale.US));
        for (int snpID = 0; snpID < numSNPs; snpID++) {
            //for (int snpID = 0; snpID < 100000; snpID++) {

            SNP snp = ds.getSNPObject(snpID);
            loader.loadGenotypes(snp);
            loader.loadDosage(snp);

            if (snpsToInclude == null || snpsToInclude.contains(snp.getName())) {


                //System.out.println(snpDataObject.rsName + "\t" + snpDataObject.hweControlsAndCasesExactP + "\t" + snpDataObject.mafOverall);

                if (snp.getMAF() > 0) {

                    boolean takeComplement = false;
                    for (int i = 0; i < nrSamplesIncluded; i++) {
                        int ind = samplesToOutput[i];

                        double dosagevalue = snp.getDosageValues()[ind];
                        if (snp.getGenotypes()[ind] == 0 && dosagevalue > 1) {
                            takeComplement = true;
                            break;
                        }
                        if (snp.getGenotypes()[ind] == 2 && dosagevalue < 1) {
                            takeComplement = true;
                            break;
                        }
                    }


                    StringBuilder sb = new StringBuilder(snp.getName() + "\t");
                    byte[] alleles = snp.getAlleles();
                    if (takeComplement) {
                        sb.append((char) alleles[1]);
                        sb.append("\t");
                        sb.append((char) alleles[0]);
                    } else {
                        sb.append((char) alleles[0]);
                        sb.append("\t");
                        sb.append((char) alleles[1]);
                    }


                    for (int i = 0; i < nrSamplesIncluded; i++) {
                        int ind = samplesToOutput[i];
                        double dosage = 2 - snp.getDosageValues()[ind];
                        String dosageString = df.format(dosage);
                        sb.append("\t");
                        sb.append(dosageString);
                    }


                    sb.append("\n");
                    if (splitperchromosome) {

                        int chr = ds.getChr(snpID);
                        if (chr > 0 && chr < 23) {
                            outPC[chr].write(sb.toString());
                        }
                    } else {
                        out.write(sb.toString());
                    }

                    counter++;
                    if (counter % 1000 == 0) {
                        System.out.println("SNPs outputted:\t" + counter);
                    }

                }
            }
        }


        if (splitperchromosome) {
            for (int chr = 1; chr < 23; chr++) {
                outPC[chr].close();
            }
        } else {
            out.close();
        }





        System.exit(0);


    }

    public void splitDosageInformationPerChromosome(String inputFile, String beagledir, String template, String batchname, String outputDir) throws IOException {

        HashMap<String, Integer> hashSNPs = new HashMap<String, Integer>();
        for (int chr = 1; chr <= 22; chr++) {
            String templatecopy = new String(template);
            templatecopy = templatecopy.replace("BATCH", batchname);
            templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

            String fileName = beagledir + "/" + templatecopy + ".r2";
            TextFile in = new TextFile(fileName, TextFile.R);
            String str = "";
            while ((str = in.readLine()) != null) {
                String data[] = str.split("\t");
                hashSNPs.put(data[0], chr);
            }
            in.close();

        }


        TextFile out[] = new TextFile[23];

        for (int chr = 1; chr <= 22; chr++) {
            out[chr] = new TextFile(outputDir + "/ImputedGenotypeDosageFormatPLINKChr" + chr + ".dose", TextFile.W);
        }

        TextFile in = new TextFile(inputFile, TextFile.R);
        String str = in.readLine();
        for (int chr = 1; chr <= 22; chr++) {
            out[chr].write(str + "\n");
        }
        int counter = 0;
        while ((str = in.readLine()) != null) {
            String data[] = str.split("\t");
            int chr = hashSNPs.get(data[0]);
            out[chr].write(str + "\n");
            System.out.println(counter + "\t" + chr);
            counter++;
        }
        for (int chr = 1; chr <= 22; chr++) {

            out[chr].close();
        }

    }
}
