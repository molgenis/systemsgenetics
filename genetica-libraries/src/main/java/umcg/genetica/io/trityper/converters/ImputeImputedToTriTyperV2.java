/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Set;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;

/**
 *
 * @author harmjan, Patrick Deelen
 */
public class ImputeImputedToTriTyperV2 {

    public void importImputedDataWithProbabilityInformationImpute(String inputDir, String outputDir, int nrSamples, String sampleListFile, String listOfSamplesToIncludeFile) throws Exception {

        System.out.println("This converter will process files containing the following patterns: chr#, chr_# and chr-#");
        ArrayList<String> snps = new ArrayList<String>();
        ArrayList<String> snpMappings = new ArrayList<String>();

        ArrayList<String> allSamples = null;
        int[] sampleToId = null;

        int nrSamplesToInclude = 0;

        if (sampleListFile != null) {
            try {
                TextFile tf = new TextFile(sampleListFile, TextFile.R);
                allSamples = tf.readAsArrayList();
                tf.close();
                System.out.println(sampleListFile + "\tcontains identifiers for " + allSamples.size() + " individuals");
                if (allSamples.size() != nrSamples) {
                    System.err.println("The number of samples you specified in " + sampleListFile + " does not correspond to the number of samples you have specified using --nrSamples.");
                    System.exit(-1);
                }
            } catch (IOException e) {
                System.err.println("Error: could not find or read file: " + sampleListFile);
            }
        }

        if (listOfSamplesToIncludeFile != null) {
            if (sampleListFile == null || allSamples == null) {
                System.err.println("ERROR: you can only specify which samples to include when you also provide a sample list file using --samples");
                System.exit(-1);
            } else {
                try {
                    Set<String> samplesToInclude = null;
                    TextFile tf = new TextFile(listOfSamplesToIncludeFile, TextFile.R);
                    samplesToInclude = tf.readAsSet(0, TextFile.tab);
                    tf.close();

                    if (samplesToInclude.size() > allSamples.size()) {
                        System.out.println("WARNING: the number of samples in " + listOfSamplesToIncludeFile + " is larger than the actual number of samples in the imputed data (according to your specification).");
                    }

                    System.out.println("About to include: " + samplesToInclude.size() + "\t out of " + allSamples.size());

                    // give each sample a new ID..
                    int ctr = 0;
                    sampleToId = new int[allSamples.size()];
                    for (int i = 0; i < allSamples.size(); i++) {
                        String sample = allSamples.get(i);
                        if (samplesToInclude.contains(sample)) {
                            sampleToId[i] = ctr;
                            ctr++;
                        } else {
                            sampleToId[i] = -1;
                        }
                    }
                    nrSamplesToInclude = ctr;

                    if (nrSamplesToInclude == 0) {
                        System.err.println("ERROR: none of the samples will be included. Check the sample identifier overlap between " + sampleListFile + " and " + listOfSamplesToIncludeFile);
                        System.exit(-1);
                    }

                    System.out.println("Total number of samples that will be eventually included: " + nrSamplesToInclude);

                } catch (IOException e) {
                    System.err.println("Error: could not find or read file: " + listOfSamplesToIncludeFile);
                }
            }
        } else {
            nrSamplesToInclude = nrSamples;
        }


        long nrSNPsAvailable = 0;
        boolean proceed = true;
        for (int chr = 1; chr <= 22; chr++) {
            // make a file list of batches for this chr....
            String[] fileList = makeFileList(inputDir, chr);
            System.out.println("Found " + fileList.length + " files for chr " + chr);
            for (int f = 0; f < fileList.length; f++) {
                String fileName = inputDir + "/" + fileList[f];
                System.out.println("Processing file:\t" + fileName);
                try {
                    TextFile in = new TextFile(fileName, TextFile.R);
                    String str;
                    while ((str = in.readLine()) != null) {
                        while (str.contains("  ")) {
                            str = str.replace("  ", " ");
                        }
                        String data[] = str.split(" ");
                        String snp = new String(data[1].getBytes());
                        String snpPos = new String(data[2].getBytes());

                        String snpMapping = chr + "\t" + snpPos + "\t" + snp;
                        snpMappings.add(snpMapping);

                        snps.add(snp);
                        nrSNPsAvailable++;
                    }
                    System.out.println("Number of SNPs parsed so far:\t" + nrSNPsAvailable);
                    System.out.println("");
                    in.close();
                } catch (java.io.EOFException e) {
                    System.err.println("Error parsing GZipped file. Possibly this file is corrupted: " + fileName);
                    proceed = false;
                } catch (IOException e) {
                    System.err.println("Error parsing file:\t" + e.getMessage());
                    proceed = false;
                }
            }
            System.out.println("");
        }

        if (!proceed) {
            System.err.println("Your data contains errors.. Will not proceed.");
            System.exit(-1);
        }

        if (snps.isEmpty()) {
            System.err.println("No snps found to convert.");
        } else {

            System.out.println("Found a total of " + snps.size() + " snps.");
            System.out.println("");

            System.out.println("Writing SNP mappings to file:");
            try {
                TextFile outSNP = new TextFile(outputDir + "SNPMappings.txt", TextFile.W);
                for (int snp = 0; snp < snpMappings.size(); snp++) {
                    outSNP.write(((String) snpMappings.get(snp)) + "\n");
                    if (snp % 2000 == 1999) {
                        System.out.print(".");
                    }
                }
                System.out.println("");
                outSNP.close();
            } catch (IOException e) {
                System.err.println("Error writing SNPs.txt file:\t" + outputDir + "SNPMappings.txt");
                System.exit(-1);
            }

            System.out.println("Writing marker definition to file:");
            try {
                TextFile outSNP = new TextFile(outputDir + "SNPs.txt", TextFile.W);
                for (int snp = 0; snp < snps.size(); snp++) {
                    outSNP.write(((String) snps.get(snp)) + "\n");
                    if (snp % 2000 == 1999) {
                        System.out.print(".");
                    }
                }
                System.out.println("");
                outSNP.close();
            } catch (IOException e) {
                System.err.println("Error writing SNPs.txt file:\t" + e.getMessage());
                System.exit(-1);
            }

            if (allSamples != null) {
                System.out.println("Writing samples to file.");

                try {
                    TextFile tfInd = new TextFile(outputDir + "Individuals.txt", TextFile.W);
                    TextFile tfPheno = new TextFile(outputDir + "PhenotypeInformation.txt", TextFile.W);
                    for (int i = 0; i < allSamples.size(); i++) {
                        if (sampleToId == null || sampleToId[i] != -1) {
                            tfInd.writeln(allSamples.get(i));
                            tfPheno.writeln(allSamples.get(i) + "\tunknown\tinclude\tunknown");
                        }
                    }
                    tfInd.close();
                    tfPheno.close();
                } catch (IOException e) {
                    System.err.println("Error writing Individuals.txt and/or PhenotypeInformation.txt files:\t" + e.getMessage());
                    System.exit(-1);
                }


            }

            int nrSNPs = (int) nrSNPsAvailable;
            WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamplesToInclude, new File(outputDir + "GenotypeMatrix.dat"), false);
            WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamplesToInclude, new File(outputDir + "/ImputedDosageMatrix.dat"), false);

            int snpIndex = 0;
            for (int chr = 1; chr <= 22; chr++) {
                // make a file list of batches for this chr....
                String[] fileList = makeFileList(inputDir, chr);


                for (int f = 0; f < fileList.length; f++) {
                    String fileName = inputDir + "/" + fileList[f];
                    //String fileName = inputDir + "/chr" + chr + ".probs";
                    System.out.println("Processing file:\t" + fileName);
                    String snp = "";
                    int lnctr = 0;
                    try {
                        TextFile in = new TextFile(fileName, TextFile.R);
                        String str;
                        while ((str = in.readLine()) != null) {
                            while (str.contains("  ")) {
                                str = str.replace("  ", " ");
                            }
                            lnctr++;
                            String data[] = str.split(" ");
                            snp = new String(data[0].getBytes());
                            byte[] allele1 = new byte[nrSamplesToInclude];
                            byte[] allele2 = new byte[nrSamplesToInclude];
                            byte[] alleles = new byte[2];
                            alleles[0] = data[3].getBytes()[0];
                            alleles[1] = data[4].getBytes()[0];
                            byte[] dosage = new byte[nrSamplesToInclude];

                            if (data.length != (nrSamples * 3) + 5) {
                                throw new Exception("Expected " + nrSamples + "samples. Found: " + (data.length - 5) / 3f + " samples");
                            }


                            for (int sample = 0; sample < nrSamples; sample++) {

                                if (sampleToId == null || sampleToId[sample] != -1) {

                                    int sampleId = sample;
                                    if (sampleToId != null) {
                                        sampleId = sampleToId[sample];
                                    }

                                    double dosageValue = Double.parseDouble(data[sample * 3 + 5 + 1]) + 2 * Double.parseDouble(data[sample * 3 + 5 + 2]);

                                    int dosageInt = (int) Math.round(dosageValue * 100d);
                                    byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                                    if (dosageInt < 0 || dosageInt > 200) {
                                        System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpIndex + "\t" + data[sample * 3 + 5] + "\t" + data[sample * 3 + 5 + 1] + "\t" + data[sample * 3 + 5 + 2]);
                                    } else {
                                        dosage[sampleId] = (byte) dosageByte;
                                    }
                                    if (dosageValue < 0.5) {
                                        allele1[sampleId] = alleles[0];
                                        allele2[sampleId] = alleles[0];
                                    } else {
                                        if (dosageValue > 1.5) {
                                            allele1[sampleId] = alleles[1];
                                            allele2[sampleId] = alleles[1];
                                        } else {
                                            allele1[sampleId] = alleles[0];
                                            allele2[sampleId] = alleles[1];
                                        }
                                    }
                                }


                            }
                            fileMatrixGenotype.setAllele1(snpIndex, 0, allele1);
                            fileMatrixGenotype.setAllele2(snpIndex, 0, allele2);
                            matrixImputedDosage.setDosage(snpIndex, 0, dosage);
                            snpIndex++;
                        }
                        in.close();
                    } catch (IOException e) {
                        System.err.println("Error parsing file:\t" + fileName);
                        System.exit(-1);
                    } catch (NumberFormatException e) {
                        System.err.println("Error parsing dosage value for SNP " + snp + " on line: " + lnctr);
                        System.exit(-1);
                    }
                }
            }

            fileMatrixGenotype.close();
            matrixImputedDosage.close();
        }

        System.exit(0);
    }

    private String[] makeFileList(String inputDir, int chr) {

        File dir = new File(inputDir);
        String[] files = dir.list();

        ArrayList<String> filelist = new ArrayList<String>();

        for (int i = 0; i < files.length; i++) {
            if (files[i].toLowerCase().endsWith(".txt.gz") || files[i].toLowerCase().endsWith(".txt") || files[i].toLowerCase().endsWith(".gz")) {
                if (files[i].toLowerCase().contains("chr" + chr) || files[i].toLowerCase().contains("chr_" + chr) || files[i].toLowerCase().contains("chr-" + chr)) {
                    filelist.add(files[i]);
                }
            }
        }

        String[] list = new String[filelist.size()];
        int i = 0;
        for (String f : filelist) {
            list[i] = f;
            i++;
        }

        java.util.Arrays.sort(list);
        return list;
    }
}
