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
import java.util.ArrayList;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;

/**
 *
 * @author harmjan
 */
public class ImputeImputedToTriTyper {

    public void importImputedDataWithProbabilityInformationImpute(String inputDir, String outputDir) throws IOException {


        ArrayList<String> vectorSNP = new ArrayList<String>();
        ArrayList<String> vectorSNPMappings = new ArrayList<String>();
        ArrayList<String> vectorInd = new ArrayList<String>();
        long nrSNPsAvailable = 0;
        for (int chr = 1; chr <= 22; chr++) {

            // make a file list of batches for this chr....
            String[] fileList = makeFileList(inputDir, chr);


            for (int f = 0; f < fileList.length; f++) {
                String fileName = inputDir + "/" + fileList[f];
                System.out.println("Processing file:\t" + fileName);
                try {
                    BufferedReader in = new BufferedReader(new FileReader(new File(fileName)));
                    String str = "";
                    while ((str = in.readLine()) != null) {
                        while (str.contains("  ")) {
                            str = str.replace("  ", " ");
                        }
                        String data[] = str.split(" ");
                        String snp = new String(data[1].getBytes());
                        String snpPos = new String(data[2].getBytes());
                        while (snpPos.length() < 9) {
                            snpPos = "0" + snpPos;
                        }
                        String snpMapping = chr + "\t" + snpPos + "\t" + snp;
                        vectorSNPMappings.add(snpMapping);

                        vectorSNP.add(snp);
                        nrSNPsAvailable++;
                    }
                    System.out.println("Number of SNPs parsed so far:\t" + nrSNPsAvailable);
                    System.out.println("");
                    in.close();
                } catch (IOException e) {
                    System.out.println("Error parsing file:\t" + e.getMessage());
                    e.printStackTrace();
                }
            }


        }

        System.out.println("\nWriting SNP mappings to file:");

        BufferedWriter outSNP = new BufferedWriter(new FileWriter(outputDir + "SNPMappings.txt"));
        for (int snp = 0; snp < vectorSNPMappings.size(); snp++) {
            outSNP.write(((String) vectorSNPMappings.get(snp)) + "\n");
            outSNP.flush();
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        outSNP.close();


        System.out.println("\nWriting marker definition to file:");

        BufferedWriter outSNPFile = new BufferedWriter(new FileWriter(outputDir + "SNPs.txt"));
        for (int snp = 0; snp < vectorSNP.size(); snp++) {
            outSNPFile.write(((String) vectorSNP.get(snp)) + "\n");
            outSNP.flush();
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        outSNP.close();


        int nrSNPs = (int) nrSNPsAvailable;
        int nrSamples = vectorInd.size();
        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, new File(outputDir + "GenotypeMatrix.dat"), false);
        WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamples, new File(outputDir + "/ImputedDosageMatrix.dat"), false);

        int snpIndex = 0;
        for (int chr = 1; chr <= 22; chr++) {
            // make a file list of batches for this chr....
            String[] fileList = makeFileList(inputDir, chr);


            for (int f = 0; f < fileList.length; f++) {
                String fileName = inputDir + "/" + fileList[f];
                //String fileName = inputDir + "/chr" + chr + ".probs";
                System.out.println("Processing file:\t" + fileName);
                try {
                    TextFile in = new TextFile(fileName, TextFile.R);
                    String str = "";
                    while ((str = in.readLine()) != null) {
                        while (str.contains("  ")) {
                            str = str.replace("  ", " ");
                        }
                        String data[] = str.split(" ");
                        String snp = new String(data[0].getBytes());
                        byte[] allele1 = new byte[nrSamples];
                        byte[] allele2 = new byte[nrSamples];
                        byte[] alleles = new byte[2];
                        alleles[0] = data[3].getBytes()[0];
                        alleles[1] = data[4].getBytes()[0];
                        byte[] dosage = new byte[nrSamples];
                        for (int sample = 0; sample < nrSamples; sample++) {
                            double dosageValue = Double.parseDouble(data[sample * 3 + 5 + 1]) + 2 * Double.parseDouble(data[sample * 3 + 5 + 2]);
                            int dosageInt = (int) Math.round(dosageValue * 100d);
                            byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                            if (dosageInt < 0 || dosageInt > 200) {
                                System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpIndex + "\t" + data[sample * 3 + 5] + "\t" + data[sample * 3 + 5 + 1] + "\t" + data[sample * 3 + 5 + 2]);
                            } else {
                                dosage[sample] = (byte) dosageByte;
                            }
                            if (dosageValue < 0.5) {
                                allele1[sample] = alleles[0];
                                allele2[sample] = alleles[0];
                            } else {
                                if (dosageValue > 1.5) {
                                    allele1[sample] = alleles[1];
                                    allele2[sample] = alleles[1];
                                } else {
                                    allele1[sample] = alleles[0];
                                    allele2[sample] = alleles[1];
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
                    System.out.println("Error parsing file:\t" + e.getMessage());
                    e.printStackTrace();
                }
            }
        }

        fileMatrixGenotype.close();
        matrixImputedDosage.close();

        System.exit(0);
    }

    private String[] makeFileList(String inputDir, int chr) {

        File dir = new File(inputDir);
        String[] files = dir.list();

        ArrayList< String> filelist = new ArrayList<String>();

        for (int i = 0; i < files.length; i++) {
            if (files[i].startsWith("chr_" + chr + "_")) {
                filelist.add(files[i]);
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
