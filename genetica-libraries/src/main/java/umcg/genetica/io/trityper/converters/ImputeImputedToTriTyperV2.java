/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;

/**
 *
 * @author harmjan, Patrick Deelen
 */
public class ImputeImputedToTriTyperV2 {

    public void importImputedDataWithProbabilityInformationImpute(String inputDir, String outputDir, int nrSamples) throws Exception {

        System.out.println("This version will procces files complying to this pattern: chr# and chr_#");


        ArrayList<String> Snps = new ArrayList<String>();
        ArrayList<String> SnpMappings = new ArrayList<String>();

        long nrSNPsAvailable = 0;
        for (int chr = 1; chr <= 22; chr++) {

            // make a file list of batches for this chr....
            String[] fileList = makeFileList(inputDir, chr);


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
//						while (snpPos.length() < 9) {
//							snpPos = "0" + snpPos;
//						}
                        String snpMapping = chr + "\t" + snpPos + "\t" + snp;
                        SnpMappings.add(snpMapping);

                        Snps.add(snp);
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
        try {
            TextFile outSNP = new TextFile(outputDir + "SNPMappings.txt", TextFile.W);
            for (int snp = 0; snp < SnpMappings.size(); snp++) {
                outSNP.write(((String) SnpMappings.get(snp)) + "\n");
                if (snp % 2000 == 1999) {
                    System.out.print(".");
                }
            }
            System.out.println("");
            outSNP.close();
        } catch (IOException e) {
            System.out.println("Error writing SNPs.txt file:\t" + e.getMessage());
            for (int ex = 0; ex < e.getStackTrace().length; ex++) {
                System.out.println(e.getStackTrace()[ex].getClassName() + "\t" + e.getStackTrace()[ex].getMethodName() + "\t" + e.getStackTrace()[ex].getLineNumber());
            }
        }

        System.out.println("\nWriting marker definition to file:");
        try {
            TextFile outSNP = new TextFile(outputDir + "SNPs.txt", TextFile.W);
            for (int snp = 0; snp < Snps.size(); snp++) {
                outSNP.write(((String) Snps.get(snp)) + "\n");
                if (snp % 2000 == 1999) {
                    System.out.print(".");
                }
            }
            System.out.println("");
            outSNP.close();
        } catch (IOException e) {
            System.out.println("Error writing SNPs.txt file:\t" + e.getMessage());
            for (int ex = 0; ex < e.getStackTrace().length; ex++) {
                System.out.println(e.getStackTrace()[ex].getClassName() + "\t" + e.getStackTrace()[ex].getMethodName() + "\t" + e.getStackTrace()[ex].getLineNumber());
            }
        }

        int nrSNPs = (int) nrSNPsAvailable;
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
                    String str;
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

                        if (data.length != (nrSamples * 3) + 5) {
                            throw new Exception("Expected " + nrSamples + "samples. Found: " + (data.length - 5) / 3f + " samples");
                        }


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

//        System.out.println(inputDir+"\t"+chr);
        File dir = new File(inputDir);
        String[] files = dir.list();

        ArrayList<String> filelist = new ArrayList<String>();

        for (int i = 0; i < files.length; i++) {
            //if(files[i].equals("chr_"+chr+"_")){
//            System.out.println(files[i]);
            String[] split = files[i].split("\\.");
			if(files[i].equals("chr" + chr)){
				filelist.add(files[i]);
			}
			else if(split.length > 1){
//                System.out.println(split[0]);
                if(split[0].equals("chr_" + chr) || split[0].equals("chr" + chr)){
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
