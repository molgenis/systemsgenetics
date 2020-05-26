/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
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

        String snpfile = null;

        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(inDir);

        HashSet<String> snpsToInclude = null;
        if (snpfile != null) {
            snpsToInclude = new HashSet<String>();

            TextFile in = new TextFile(snpfile, TextFile.R);
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

        int startchr = -1;
        int stopchr = 0;
        if (splitperchromosome) {
            startchr = 1;
            stopchr = 23;
        }

        long counter = 0;

        String[] SNPs = ds.getSNPs();
        int numSNPs = SNPs.length;

        SNPLoader loader = ds.createSNPLoader();
        java.text.DecimalFormat df = new java.text.DecimalFormat("0.00", new java.text.DecimalFormatSymbols(java.util.Locale.US));

        int buffersize = 1000;

        for (int chr = startchr; chr < stopchr; chr++) {

            String fileName = outputDir + "/ImputedGenotypeDosageFormatPLINK.dose.gz";
            if (splitperchromosome) {
                fileName = outputDir + "/ImputedGenotypeDosageFormatPLINK-Chr" + chr + ".dose.gz";
            }
            TextFile out = new TextFile(fileName, TextFile.W);
            String header = "SNP\tA1\tA2";
            for (int i = 0; i < nrSamplesIncluded; i++) {
                String sample = individuals[samplesToOutput[i]];
                String familyId = sampleToFamId.get(sample);
                if (familyId == null) {
                    familyId = "1";
                }

                header += "\t" + familyId + "\t" + sample;
            }
            out.write(header + "\n");

            ArrayList<Integer> snpsToFinallyInclude = new ArrayList<Integer>();;
            if (splitperchromosome) {
                for (int snp = 0; snp < SNPs.length; snp++) {
                    String snpName = SNPs[snp];
                    if (snpsToInclude == null || snpsToInclude.contains(snpName)) {
                        Byte chrName = ds.getChr(snp);
                        if (chrName != null && chrName == chr) {
                            snpsToFinallyInclude.add(snp);
                        }
                    }
                }
                System.out.println("Exporting " + snpsToFinallyInclude.size() + " SNPs for chromosome: " + chr);
            } else {
                for (int snp = 0; snp < SNPs.length; snp++) {
                    String snpName = SNPs[snp];
                    if (snpsToInclude == null || snpsToInclude.contains(snpName)) {
                        snpsToFinallyInclude.add(snp);
                    }
                }
                System.out.println("Exporting " + snpsToFinallyInclude.size() + " SNPs");
            }

            int bufferCounter = 0;
            int snpCounter = 0;
            StringBuilder sb = new StringBuilder();
            ProgressBar pb = new ProgressBar(snpsToFinallyInclude.size(), "Exporting " + snpsToFinallyInclude.size() + " SNPs");
            while (snpCounter < snpsToFinallyInclude.size()) {
                Integer SNPIdToExport = snpsToFinallyInclude.get(snpCounter);
                SNP snp = ds.getSNPObject(SNPIdToExport);
                loader.loadGenotypes(snp);

                if (snp.getMAF() > 0) {

                    loader.loadDosage(snp);
                    boolean takeComplement = false;
                    byte[] genotypes = snp.getGenotypes();
                    double[] dosageValues = snp.getDosageValues();
                    for (int i = 0; i < nrSamplesIncluded; i++) {
                        int ind = samplesToOutput[i];
                        byte genotype  = genotypes[ind];
                        double dosagevalue = dosageValues[ind];
                        if (genotype == 0 && dosagevalue > 1) {
                            takeComplement = true;
                            break;
                        }
                        if (genotype == 2 && dosagevalue < 1) {
                            takeComplement = true;
                            break;
                        }
                    }

                    sb.append(snp.getName()).append("\t");
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
                        double dosage = 2 - dosageValues[ind];
                        sb.append("\t");
                        sb.append(df.format(dosage));
                    }
                    sb.append("\n");

                }
                bufferCounter++;
                snpCounter++;

                snp.clearGenotypes();
                if (snpCounter % buffersize == 0) {
                    out.write(sb.toString());
                    sb = new StringBuilder();
                    bufferCounter = 0;
//                    System.out.println(snpCounter + " SNPs exported");
                }
                pb.iterate();
            }
            pb.close();

            if (bufferCounter != 0) {
                out.write(sb.toString());
            }


            out.close();
        }
        loader.close();


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
