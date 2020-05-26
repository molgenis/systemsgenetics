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
 * @author harmjan
 */
public class BeagleImputedToTriTyper {

    private HashMap<String, String> familyData;
    private boolean familyDataLoaded;
    private ArrayList<String> vArrayListInd;
    private HashMap<String, Integer> vhashInd;
    private ArrayList<String> arrayListInd;
    private ArrayList<String> arrayListSNP;
    private HashMap<String, Double> SNPR2;
    private HashMap<String, Integer> SNPR2Present;
    private ArrayList<String> arrayListSNPMappings;
    private HashMap<String, Integer> hashInd;
    private HashMap<String, Integer> hashSNP;

    private String[] getBatches(int numBatches) {

        String[] batches = new String[numBatches];
        String[] alphabet = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"};

        String firstletter = "a";
        int alphacounter = 0;
        int betacounter = 0;
        for (int i = 0; i < numBatches; i++) {
            if (i % 26 == 0) {
                firstletter = alphabet[alphacounter];
                alphacounter++;
                betacounter = 0;
            }

            batches[i] = firstletter + alphabet[betacounter];
            betacounter++;
        }
        return batches;
    }

    public boolean importImputedDataWithDosageInformationBeagleBatches(String inputDir, String template, int numBatches, String outputDir, int chrStart, int chrEnd) throws IOException {

        String[] batchNames = getBatches(numBatches);

//	int chrStart = 1;
//	int chrEnd = 22;
        // String[] batchNames = {"1", "2", "3", "4", "5"};
        int nrBatches = batchNames.length;
        int[] batchSampleIndex = new int[nrBatches + 1];
        int[] batchNrSamples = new int[nrBatches + 1];

        // int[] chromosomes = {2,9,10,16,17,18,20};

        hashSNP = new HashMap<String, Integer>();
        arrayListSNP = new ArrayList<String>();
        arrayListSNPMappings = new ArrayList<String>();
        arrayListInd = new ArrayList<String>();
        hashInd = new HashMap<String, Integer>();
        long nrSNPsAvailable = 0;

        // goddaf1-10QC-imput_Chr1-HM2-0.goddaf1-10QC-imput_Chr1-0-BEAGLE.gprobs
        // COPD-Asia_Chr22-HM2-g.COPD-Asia_Chr22-g-BEAGLE.gprobs.gz
        boolean allFilesAvailable = true;
        int filesnotfound = 0;
        for (int batch = 0; batch < nrBatches; batch++) {
            batchSampleIndex[batch] = arrayListInd.size();
            int nrSamplesThisBatch = 0;
            for (int chr = chrStart; chr <= chrEnd; chr++) {
                // for(int chr: chromosomes){
                String templatecopy = new String(template);
                templatecopy = templatecopy.replace("BATCH", batchNames[batch]);
                templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

                // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
                String fileName = inputDir + "/" + templatecopy + ".gprobs.gz";
                String r2FileName = inputDir + "/" + templatecopy + ".r2";


                if (!Gpio.canRead(fileName)) {
                    System.out.println("Cannot open file:\t" + fileName);
                    allFilesAvailable = false;
                    filesnotfound++;
                }

                if (!Gpio.canRead(r2FileName)) {
                    System.out.println("Cannot open file:\t" + r2FileName);
                    allFilesAvailable = false;
                    filesnotfound++;
                }

            }
        }

        if (!allFilesAvailable) {
            System.out.println("Not all imputed dosage files are available!!! (" + filesnotfound + " out of " + (nrBatches * 22) + "). Exiting...");
            return false;
        }

        boolean preprocessingCompleted = false;

        String snpfile = outputDir + "SNPs.txt";
        String indfile = outputDir + "Individuals.txt";
        if (Gpio.canRead(snpfile) && Gpio.canRead(indfile)) {
            preprocessingCompleted = true;
        }
        if (preprocessingCompleted) {

            System.out.println("Preprocessing has already been completed, not necessary to conduct preprocessing again (SNPs.txt and Individuals.txt already exist):");
            System.out.println("Parsing SNPs.txt:");

            TextFile in = new TextFile(snpfile, TextFile.R);
            String str = null;
            while ((str = in.readLine()) != null) {
                String snp = str.trim();
                hashSNP.put(snp, arrayListSNP.size());
                arrayListSNP.add(snp);
                nrSNPsAvailable++;
            }
            in.close();


            for (int batch = 0; batch < nrBatches; batch++) {
                batchSampleIndex[batch] = arrayListInd.size();
                System.out.println("BatchSampleIndex:\t" + batch + "\t" + batchSampleIndex[batch]);
                int nrSamplesThisBatch = 0;
                for (int chr = chrStart; chr <= chrEnd; chr++) {

                    String templatecopy = new String(template);
                    templatecopy = templatecopy.replace("BATCH", batchNames[batch]);
                    templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

                    String fileName = inputDir + "/" + templatecopy + ".gprobs.gz";
                    String r2FileName = inputDir + "/" + templatecopy + ".r2";

                    System.out.print("Processing file:\t" + fileName);

                    in = new TextFile(fileName, TextFile.R);
                    str = in.readLine();
                    String data[] = str.split(" ");
                    for (int c = 3; c < data.length; c++) {
                        if (!hashInd.containsKey(data[c])) {
                            hashInd.put(data[c], arrayListInd.size());
                            arrayListInd.add(data[c]);
                            nrSamplesThisBatch++;
                        }
                    }
                    System.out.println("Number of individuals parsed so far:\t" + arrayListInd.size());
                    in.close();
                }
                System.out.println("Number of unique samples for batch:\t" + batch + "\t" + nrSamplesThisBatch);
                batchNrSamples[batch] = nrSamplesThisBatch;
            }

        } else {
            SNPR2 = new HashMap<String, Double>();
            SNPR2Present = new HashMap<String, Integer>();

            for (int batch = 0; batch < nrBatches; batch++) {
                batchSampleIndex[batch] = arrayListInd.size();
                int nrSamplesThisBatch = 0;
                // for(int chr: chromosomes){
                for (int chr = chrStart; chr <= chrEnd; chr++) {
                    String templatecopy = new String(template);
                    templatecopy = templatecopy.replace("BATCH", batchNames[batch]);
                    templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

                    // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
                    String fileName = inputDir + "/" + templatecopy + ".gprobs.gz";
                    String r2FileName = inputDir + "/" + templatecopy + ".r2";

                    System.out.print("Processing file:\t" + fileName);
                    // read batch r2 scores
                    TextFile r2reader = new TextFile(r2FileName, TextFile.R);
                    String ln = "";

                    while ((ln = r2reader.readLine()) != null) {
                        String[] elems = ln.split("\t");
                        if (elems.length == 2) {
                            try {
                                Double r2Val = Double.parseDouble(elems[1]);
                                Double r2Prev = SNPR2.get(elems[0]);
                                if (r2Prev == null) {
                                    SNPR2.put(elems[0], r2Val);
                                    SNPR2Present.put(elems[0], 1);
                                } else {
                                    SNPR2.put(elems[0], r2Prev + r2Val);
                                    SNPR2Present.put(elems[0], SNPR2Present.get(elems[0]) + 1);
                                }
                            } catch (NumberFormatException e) {
                                if (elems[1].toLowerCase().equals("nan")) {
                                }
                            }
                        }
                    }

                    r2reader.close();


                    TextFile in = new TextFile(fileName, TextFile.R);

                    String str = in.readLine();
                    String data[] = str.split(" ");
                    for (int c = 3; c < data.length; c++) {
                        if (!hashInd.containsKey(data[c])) {
                            hashInd.put(data[c], arrayListInd.size());
                            arrayListInd.add(data[c]);
                            nrSamplesThisBatch++;
                        }
                    }

                    int line = 0;
                    int prevLine = -1;
                    while ((str = in.readLine()) != null) {
                        //			while (str.contains("  ")) {
                        //			    str = str.replace("  ", " ");
                        //			}
                        data = str.split(" ");
                        String snp = new String(data[0].getBytes());
                        if (!hashSNP.containsKey(snp)) {
                            String snpPos = "1";
                            String snpMapping = chr + "\t" + snpPos + "\t" + snp;
                            arrayListSNPMappings.add(snpMapping);
                            hashSNP.put(snp, arrayListSNP.size());
                            arrayListSNP.add(snp);
                            nrSNPsAvailable++;
                        }
                        if (line % 10000 == 0 && line > prevLine) {
                            System.out.print(".");
                            prevLine = line;
                        }
                        line++;
                    }
                    System.out.println("");
                    System.out.println("Number of SNPs parsed so far:\t" + nrSNPsAvailable);
                    System.out.println("Number of individuals parsed so far:\t" + arrayListInd.size());
                    System.out.println("");
                    in.close();
                }
                System.out.println("Number of unique samples for batch:\t" + batch + "\t" + nrSamplesThisBatch);
                batchNrSamples[batch] = nrSamplesThisBatch;
            }


        }
        
        writeIndividuals(outputDir);
        writeSNPs(outputDir);
        
        int nrSNPs = (int) nrSNPsAvailable;
        int nrSamples = arrayListInd.size();
        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, new File(outputDir + "GenotypeMatrix.dat"), false);
        WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamples, new File(outputDir + "/ImputedDosageMatrix.dat"), false);



        //for (int chr = chrStart; chr <= chrEnd; chr++) {
        for (int chr = chrEnd; chr >= chrStart; chr--) {
            
            TextFile[] in = new TextFile[nrBatches];
            for (int batch = 0; batch < nrBatches; batch++) {
                String templatecopy = new String(template);
                templatecopy = templatecopy.replace("BATCH", batchNames[batch]);
                templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);
                String fileName = inputDir + "/" + templatecopy + ".gprobs.gz";
                System.out.println("Processing file:\t" + fileName);
                in[batch] = new TextFile(fileName, TextFile.R);
                in[batch].readLine();
            }
            System.out.print("Progress:");
            int line = 0;
            while (1 == 1) {

                String[] str = new String[nrBatches];
                String[] snp = new String[nrBatches];
                byte[] alleles = new byte[2];

                for (int batch = 0; batch < nrBatches; batch++) {
                    str[batch] = in[batch].readLine();
                }
                if (str[0] == null) {
                    break;
                }

                for (int batch = 0; batch < nrBatches; batch++) {
                    while (str[batch].contains("  ")) {
                        str[batch] = str[batch].replace("  ", " ");
                    }
                    String data[] = str[batch].split(" ");
                    snp[batch] = new String(data[0].getBytes());
                    alleles[0] = data[1].getBytes()[0];
                    alleles[1] = data[2].getBytes()[0];
                }


                //Check whether SNP names are identical:
                String firstSNP = snp[0];
                for (int batch = 0; batch < nrBatches; batch++) {
                    if (!firstSNP.equals(snp[batch])) {
                        System.out.println("Error! Format of different batches are not identical, different SNP names are found for the same readline!:\t" + firstSNP + "\t" + snp[batch]);
                        System.exit(-1);
                    }
                }

                //SNP names are identical, proceed:
                Integer snpIndex = hashSNP.get(firstSNP);

                byte[] allele1 = new byte[nrSamples];
                byte[] allele2 = new byte[nrSamples];
                byte[] dosage = new byte[nrSamples];

                for (int batch = 0; batch < nrBatches; batch++) {
                    String data[] = str[batch].split(" ");
                    for (int sample = 0; sample < batchNrSamples[batch]; sample++) {
                        int sampleIndex = batchSampleIndex[batch] + sample;
                        double dosageValue = Double.parseDouble(data[sample * 3 + 4]) * 1d + Double.parseDouble(data[sample * 3 + 5]) * 2d;
                        int dosageInt = (int) Math.round(dosageValue * 100d);
                        byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                        if (dosageInt < 0 || dosageInt > 200) {
                            System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpIndex + "\t" + data[sample * 3 + 3] + "-" + data[sample * 3 + 4] + "-" + data[sample * 3 + 5]);
                        } else {
                            dosage[sampleIndex] = (byte) dosageByte;
                        }
                        if (dosageValue < 0.5) {
                            allele1[sampleIndex] = alleles[0];
                            allele2[sampleIndex] = alleles[0];
                        } else {
                            if (dosageValue > 1.5) {
                                allele1[sampleIndex] = alleles[1];
                                allele2[sampleIndex] = alleles[1];
                            } else {
                                allele1[sampleIndex] = alleles[0];
                                allele2[sampleIndex] = alleles[1];
                            }
                        }
                    }
                }
                fileMatrixGenotype.setAllele1(snpIndex, 0, allele1);
                fileMatrixGenotype.setAllele2(snpIndex, 0, allele2);
                matrixImputedDosage.setDosage(snpIndex, 0, dosage);
                if (line % 10000 == 0) {
                    System.out.print(".");
                }
                line++;

            }
            System.out.println(" Completed.");

            for (int batch = 0; batch < nrBatches; batch++) {
                in[batch].close();
            }
        }
        fileMatrixGenotype.close();
        matrixImputedDosage.close();

        return true;
    }

    public void importImputedDataWithDosageInformationBeagle(String inputLocation, String template, String extension, String outputDir) throws IOException {

        boolean allFilesAvailable = true;
        int filesnotfound = 0;
        int chrStart = 1;
        int chrEnd = 22;
        for (int chr = chrStart; chr <= chrEnd; chr++) {
            // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
            String templatecopy = new String(template);
            templatecopy = templatecopy.replace("BATCH", "aa");
            templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

            // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
            String fileName = inputLocation + "/" + templatecopy + "." + extension;
            // String fileName = inputLocation + "/" + descriptor + "-chr" + chr +"."+  extension;
            if (!Gpio.canRead(fileName)) {
                System.out.println("Cannot open file:\t" + fileName);
                allFilesAvailable = false;
                filesnotfound++;
            }
        }


        if (!allFilesAvailable) {
            System.out.println("Not all imputed dosage files are available!!! (" + filesnotfound + " out of 22 ). Exiting...");
            System.exit(-1);
        }

        hashSNP = new HashMap<String, Integer>();
        arrayListSNP = new ArrayList<String>();
        arrayListSNPMappings = new ArrayList<String>();
        arrayListInd = new ArrayList<String>();
        hashInd = new HashMap<String, Integer>();
        long nrSNPsAvailable = 0;

        for (int chr = chrStart; chr <= chrEnd; chr++) {
            // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
            String templatecopy = new String(template);
            templatecopy = templatecopy.replace("BATCH", "aa");
            templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

            // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
            String fileName = inputLocation + "/" + templatecopy + "." + extension;
            System.out.print("Processing file:\t" + fileName);


            TextFile in = new TextFile(fileName, TextFile.R);

            String str;
            String data[];
            if (!familyDataLoaded) {
                str = in.readLine();
                // System.out.println(str);
                data = str.split(" ");

                for (int c = 3; c < data.length; c++) {
                    if (!hashInd.containsKey(data[c])) {
                        hashInd.put(data[c], arrayListInd.size());
                        arrayListInd.add(data[c]);
                        System.out.println("Found new individual:" + data[c]);

                        // System.out.println(data[c]);
                    }
                }
            } else {
                hashInd = vhashInd;
                arrayListInd = vArrayListInd;
            }
            int prevLine = -1;
            int line = 0;
            int numValues = 0;
            while ((str = in.readLine()) != null) {
                while (str.contains("  ")) {
                    str = str.replace("  ", " ");
                }
                data = str.split(" ");

                numValues = (data.length - 3) / 3;

                String snp = new String(data[0].getBytes());
                if (!hashSNP.containsKey(snp)) {
                    String snpPos = "1";
                    String snpMapping = chr + "\t" + snpPos + "\t" + snp;
                    arrayListSNPMappings.add(snpMapping);
                    hashSNP.put(snp, arrayListSNP.size());
                    arrayListSNP.add(snp);
                    nrSNPsAvailable++;
                }

                if (line % 10000 == 0 && line > prevLine) {
                    System.out.print(".");
                    prevLine = line;
                }
                line++;
            }
            System.out.println("");
            System.out.println("Number of SNPs parsed so far:\t" + nrSNPsAvailable + " for " + numValues + " samples");

            System.out.println("");
            in.close();
        }

        System.out.println("Number of individuals parsed:\t" + arrayListInd.size());

        writeIndividuals(outputDir);
        writeSNPs(outputDir);
        
        int nrSNPs = (int) nrSNPsAvailable;
        int nrSamples = arrayListInd.size();
        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, new File(outputDir + "GenotypeMatrix.dat"), false);
        WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamples, new File(outputDir + "/ImputedDosageMatrix.dat"), false);
        for (int chr = chrStart; chr <= chrEnd; chr++) {
            String templatecopy = new String(template);
            templatecopy = templatecopy.replace("BATCH", "aa");
            templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

            // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
            String fileName = inputLocation + "/" + templatecopy + "." + extension;

            System.out.print("Processing file:\t" + fileName);

            TextFile in = new TextFile(fileName, TextFile.R);

            String str = in.readLine();

            int prevLine = -1;
            int line = 0;
            while ((str = in.readLine()) != null) {
                while (str.contains("  ")) {
                    str = str.replace("  ", " ");
                }
                String[] data = str.split(" ");
                String snp = new String(data[0].getBytes());
                int snpIndex = ((Integer) hashSNP.get(snp)).intValue();
                byte[] allele1 = new byte[nrSamples];
                byte[] allele2 = new byte[nrSamples];
                byte[] alleles = new byte[2];
                alleles[0] = data[1].getBytes()[0];
                alleles[1] = data[2].getBytes()[0];
                byte[] dosage = new byte[nrSamples];
                for (int sample = 0; sample < nrSamples; sample++) {
                    // AB BB
                    double dosageValue = Double.parseDouble(data[sample * 3 + 4]) * 1d + Double.parseDouble(data[sample * 3 + 5]) * 2d;

                    int dosageInt = (int) Math.round(dosageValue * 100d);
                    byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                    if (dosageInt < 0 || dosageInt > 200) {
                        System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpIndex + "\t" + data[sample * 3 + 3] + "-" + data[sample * 3 + 4] + "-" + data[sample * 3 + 5]);
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
                    fileMatrixGenotype.setAllele1(snpIndex, sample, allele1);
                    fileMatrixGenotype.setAllele2(snpIndex, sample, allele2);
                    matrixImputedDosage.setDosage(snpIndex, sample, dosage);
                }

                if (line % 10000 == 0 && line > prevLine) {
                    System.out.print(".");
                    prevLine = line;
                }
                line++;

            }
            in.close();

            System.out.println("");
        }

        fileMatrixGenotype.close();
        matrixImputedDosage.close();

    }

    public void loadFamFile(String file) {
        vArrayListInd = new ArrayList();
        vhashInd = new HashMap();
        System.out.println("Loading FAM file:\t" + file);
        try {
            TextFile in = new TextFile(file, TextFile.R);

            String line = "";
            familyData = new HashMap<String, String>();

            while ((line = in.readLine()) != null) {
                String[] elems = line.split(" ");

                if (elems.length >= 6) {
                    String sample = elems[1];

                    String fid = elems[0];
                    String pid = elems[2];
                    String mid = elems[3];
                    String sex = elems[4];
                    String phe = elems[5];

                    if (sex.equals("1")) {
                        sex = "male";
                    } else if (sex.equals("2")) {
                        sex = "female";
                    } else {
                        sex = "unknown";
                    }

                    if (phe.equals("1")) {
                        phe = "control";
                    } else if (phe.equals("2")) {
                        phe = "case";
                    } else {
                        phe = "unknown";
                    }

                    familyData.put(sample, sample + "\t" + phe + "\tinclude\t" + sex + "\t" + fid + "\t" + pid + "\t" + mid + "\n");

                    if (!vhashInd.containsKey(sample)) {
                        vhashInd.put(sample, vArrayListInd.size());
                        vArrayListInd.add(sample);
                    }

                }

            }

            in.close();

            familyDataLoaded = true;

        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }



    }

    private void writeSNPs(String outputDir) throws IOException {
        System.out.println("\nWriting SNP mappings to file:");

        TextFile outSNP = new TextFile(outputDir + "SNPMappings.txt", TextFile.W);
        for (int snp = 0; snp < arrayListSNPMappings.size(); snp++) {
            outSNP.write(arrayListSNPMappings.get(snp) + "\n");
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        outSNP.close();

        outSNP = new TextFile(outputDir + "SNPR2Scores.txt", TextFile.W);
        for (int snp = 0; snp < arrayListSNP.size(); snp++) {
            String rsName = arrayListSNP.get(snp);
            Double r2val = SNPR2.get(rsName);
            if (r2val == null) {
                r2val = Double.NaN;
            } else {
                r2val /= SNPR2Present.get(rsName);
            }
            outSNP.write(rsName + "\t" + r2val + "\t" + SNPR2Present.get(rsName) + "\n");
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        outSNP.close();

        System.out.println("\nWriting marker definition to file:");
        outSNP = new TextFile(outputDir + "SNPs.txt", TextFile.W);
        for (int snp = 0; snp < arrayListSNP.size(); snp++) {
            outSNP.write(arrayListSNP.get(snp) + "\n");
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        outSNP.close();
    }

    private void writeIndividuals(String outputDir) throws IOException {
        System.out.println("\nWriting individuals to file:");

        TextFile outInd = new TextFile(outputDir + "Individuals.txt", TextFile.W);
        TextFile outPhe = new TextFile(outputDir + "PhenotypeInformation.txt", TextFile.W);
        for (int ind = 0; ind < arrayListInd.size(); ind++) {
            outInd.write(arrayListInd.get(ind) + "\n");
            if (familyDataLoaded) {
                if (familyData.get(arrayListInd.get(ind)) != null) {
                    outPhe.write(familyData.get(arrayListInd.get(ind)));
                } else {
                    outPhe.write(arrayListInd.get(ind) + "\tcontrol\tinclude\tunknown" + "\n");
                }

            } else {
                outPhe.write(arrayListInd.get(ind) + "\tcontrol\tinclude\tfemale" + "\n");
            }
            if (ind % 5 == 4) {
                System.out.print(".");
            }
        }

        System.out.println("");
        outInd.close();
        outPhe.close();
    }
}
