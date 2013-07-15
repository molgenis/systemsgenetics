/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collections;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class FinalReportToTriTyper {

    final static double twoDividedByPI = 2.0d / Math.PI;

    /**
     * Starts parsing a genotype report file, which can be in different formats.
     * This import program can accomodate many formats, and uses buffering to
     * achieve sufficient performance.
     */
    public FinalReportToTriTyper(String inputFile, String outputDirString, boolean isIlluminaFinalReportFile, String delimiter, String decimalSeparator) throws IOException {

        //Check whether we can write to the output directory:
        File outputDir = new File(outputDirString);
        if (!outputDir.isDirectory()) {
            System.out.println("Your output directory does not exist!");
            System.exit(-1);
        }

        //ArrayLists and hashes for determining file size of final report file:
        HashMap<String, Integer> hashInd = new HashMap<String, Integer>();
        ArrayList<String> vecInd = new ArrayList<String>();
        HashMap<String, Integer> hashSNP = new HashMap<String, Integer>();
        ArrayList<String> vecSNP = new ArrayList<String>();

        //First parse file, determine what the amount of unique samples and SNPs is.
        System.out.println("");
        System.out.println("TriTyperImporter V1.0, 2008, Lude Franke, University Medical Centre Utrecht, Lude@ludesign.nl");
        System.out.println("");
        System.out.println("Processing file:\t" + inputFile);
        System.out.println("Inventorizing input file, determining number of unique SNPs and samples:");
        int columnSample = -1;
        int columnSNP = -1;
        int columnAllele1 = -1;
        int columnAllele2 = -1;
        int columnTheta = -1;
        int columnR = -1;
        int columnX = -1;
        int columnY = -1;
        boolean rawDataAvailable = false;

        //Try to open the input file:
        if (!Gpio.canRead(inputFile)) {
            System.out.println("");
            System.out.println("Cannot open file:\t" + inputFile);
            System.out.println("Are you sure it is located at this place?");
            System.exit(-1);
        }
        TextFile in = new TextFile(inputFile, TextFile.R);

        //If this is an Illumina Final Report file, first process the irrelevant header:
        String str = null;
        if (isIlluminaFinalReportFile) {
            int countIlluminaFinalReport = 0;
            while ((str = in.readLine()) != null) {
                String[] data = str.split(delimiter);
                if (data[0].trim().equals("[Data]")) {
                    break;
                }
                countIlluminaFinalReport++;
                if (countIlluminaFinalReport > 100) {
                    System.out.println("\nError: You have defined that this file is a final report file, which it does not seem to be as a row with the word [Data] cannot be found!");
                    System.exit(-1);
                }
            }
        }

        //Now parse the column identifiers:
        str = in.readLine();

        //Check whether we actually are dealing with a Final Report file, user might have forgotten to instruct this:
        if (str.toLowerCase().startsWith("[header]")) {
            while ((str = in.readLine()) != null) {
                String[] data = str.split(delimiter);
                if (data[0].trim().equals("[Data]")) {
                    break;
                }
            }
            str = in.readLine();
        }

        String[] data = str.split(delimiter);
        if (data.length <= 1) {
            System.out.println("");
            System.out.println("Error parsing input file! The file cannot be delimited!");
            String delimiterDescription = "tab";
            if (delimiter.equals(" ")) {
                delimiterDescription = "space";
            }
            if (delimiter.equals(",")) {
                delimiterDescription = "comma";
            }
            if (delimiter.equals(";")) {
                delimiterDescription = "semicolon";
            }
            System.out.println("Are you sure it is " + delimiterDescription + " delimited ?");
            System.exit(-1);
        }

        for (int d = 0; d < data.length; d++) {
            String column = data[d].trim().toLowerCase();
            if (column.equals("sample id")) {
                columnSample = d;
            }
            if (column.equals("snp name")) {
                columnSNP = d;
            }
            if (column.contains("allele1")) {
                columnAllele1 = d;
            }
            if (column.contains("allele 1")) {
                columnAllele1 = d;
            }
            if (column.contains("allele2")) {
                columnAllele2 = d;
            }
            if (column.contains("allele 2")) {
                columnAllele2 = d;
            }
            if (column.equals("r")) {
                columnR = d;
            }
            if (column.equals("theta")) {
                columnTheta = d;
            }
        }
        if (columnSample == -1) {
            System.out.println("\nError: Within the header of this file the sample id column (Sample ID) cannot be found!");
            System.exit(-1);
        }
        if (columnAllele1 == -1) {
            System.out.println("\nError: Within the header of this file the allele 1 column (Allele1) cannot be found!");
            System.exit(-1);
        }
        if (columnAllele2 == -1) {
            System.out.println("\nError: Within the header of this file the allele 2 column (Allele2) cannot be found!");
            System.exit(-1);
        }
        if (columnSNP == -1) {
            System.out.println("\nError: Within the header of this file the SNP name column (SNP Name) cannot be found!");
            System.exit(-1);
        }
        rawDataAvailable = true;
        for (int d = 0; d < data.length; d++) {
            String column = data[d].trim().toLowerCase();
            if (column.equals("x")) {
                columnX = d;
            }
            if (column.equals("y")) {
                columnY = d;
            }
        }

        if ((columnR == -1 || columnTheta == -1) && (columnX == -1 || columnY == -1)) {
            System.out.println("Within the header of this file no raw intensity data is present (either R and Theta, or X and Y). Only imputation of triallelic SNPs will be possible");
            rawDataAvailable = false;
        }

        System.out.println("");

        boolean fileAlreadyInventorized = false;
        if ((new File(outputDirString + "Individuals.txt")).exists() && (new File(outputDirString + "SNPs.txt")).exists()) {
            fileAlreadyInventorized = true;
        }

        if (!fileAlreadyInventorized) {

            //Start processing this file
            String previousSNP = null;
            String previousInd = null;
            long linesProcessed = 0;
            while ((str = in.readLine()) != null) {
                // System.out.println(str);
                data = str.split(delimiter);



                if (data.length <= 1) {
                    System.out.println("\nError parsing input file! The file cannot be delimited!");
                    String delimiterDescription = "tab";
                    if (delimiter.equals(" ")) {
                        delimiterDescription = "space";
                    }
                    if (delimiter.equals(",")) {
                        delimiterDescription = "comma";
                    }
                    if (delimiter.equals(";")) {
                        delimiterDescription = "semicolon";
                    }
                    System.out.println("Are you sure it is " + delimiterDescription + " delimited ?");

                }
                if (data.length <= columnSNP || data.length <= columnSample) {
                    System.out.println("\nError: For record entry " + (linesProcessed + 1) + " the SNP or sample cannot be parsed! Record: " + str);
                    System.exit(-1);
                }
                String snp = data[columnSNP];
                String ind = data[columnSample];

                if (!snp.equals(previousSNP) && !hashSNP.containsKey(snp)) {
                    hashSNP.put(snp, vecSNP.size());
                    vecSNP.add(snp);
                }

                if (!ind.equals(previousInd) && !hashInd.containsKey(ind)) {
                    hashInd.put(ind, vecInd.size());
                    vecInd.add(ind);
                }

                previousSNP = snp;
                previousInd = ind;

                linesProcessed++;
                if (linesProcessed % 500000 == 0) {
                    System.out.println(linesProcessed + "\tLines processed. Number of unique SNPs read so far:\t" + vecSNP.size() + "\tNumber of unique Individuals read so far:\t" + vecInd.size());
                }
            }
            System.out.println(linesProcessed + "\tLines processed. Number of unique SNPs read:\t" + vecSNP.size() + "\tNumber of unique Individuals read:\t" + vecInd.size());
            in.close();

            //Check whether SNPMappings.txt is available. This will improve processing speed considerably in subsequent operations:
            String fileSNPMappings = new String(outputDirString + "SNPMappings.txt");
            if (!Gpio.canRead(fileSNPMappings)) {
                System.out.println("\nNon critical warning: SNPMappings.txt can not be found in the output directory. Data will not be stored in optimized way, which will negatively affect the speed of TriTyper.\n");
            } else {

                System.out.println("\nLoading SNP mappings from file:\t" + fileSNPMappings);
                TextFile inSNP = new TextFile(fileSNPMappings, TextFile.R);
                String str2;
                ArrayList<String> vectorTemp = new ArrayList<String>();
                boolean needsSorting = false;
                while ((str2 = inSNP.readLine()) != null) {
                    data = str2.split("\t");
                    if (hashSNP.containsKey(data[2])) {
                        if (data[1].length() != 9) {
                            needsSorting = true;
                            while (data[1].length() < 9) {
                                data[1] = "0" + data[1];
                            }
                        }
                        vectorTemp.add(data[0] + "\t" + data[1] + "\t" + data[2]);
                    }
                }
                inSNP.close();
                if (needsSorting) {
                    System.out.println("Sorting SNPs on chromosome and physical position that are present in SNP mappings file:");
                    Collections.sort(vectorTemp);
                }

                HashMap<String, Integer> hashSNPMappings = new HashMap<String, Integer>();
                ArrayList<String> vecSNPMappings = new ArrayList<String>();
                for (int snp = 0; snp < vectorTemp.size(); snp++) {
                    String snpString = vectorTemp.get(snp);
                    hashSNPMappings.put(snpString.split("\t")[2], vecSNPMappings.size());
                    vecSNPMappings.add(snpString.split("\t")[2]);
                }
                System.out.println("Number of SNPs with available physical mappings:\t" + vecSNPMappings.size());

                //Now sort the processed SNPs and arrange them, according to what is known:
                boolean[] snpMappingsUsed = new boolean[vecSNPMappings.size()];
                ArrayList vecSNPCopy = new ArrayList();
                for (int snp = 0; snp < vecSNP.size(); snp++) {
                    String rsName = vecSNP.get(snp);
                    if (hashSNPMappings.containsKey(rsName)) {
                        snpMappingsUsed[hashSNPMappings.get(rsName)] = true;
                    }
                }
                ArrayList<String> vecSNPNew = new ArrayList<String>();
                HashMap<String, Integer> hashSNPNew = new HashMap<String, Integer>();
                for (int snp = 0; snp < vecSNPMappings.size(); snp++) {
                    if (snpMappingsUsed[snp]) {
                        String rsName = vecSNPMappings.get(snp);
                        hashSNPNew.put(rsName, vecSNPNew.size());
                        vecSNPNew.add(rsName);
                    }
                }

                //Now add the SNPs for which no mapping is available. These will be imported, but we cannot do anything with them:
                ArrayList<String> snpsWithoutMapping = new ArrayList<String>();
                for (int snp = 0; snp < vecSNP.size(); snp++) {
                    String rsName = vecSNP.get(snp);
                    if (!hashSNPNew.containsKey(rsName)) {
                        hashSNPNew.put(rsName, vecSNPNew.size());
                        vecSNPNew.add(rsName);
                        snpsWithoutMapping.add(rsName);
                    }
                }
                if (snpsWithoutMapping.size() > 0) {
                    System.out.println("Non critical warning: No physical mapping is available for SNPs:");
                    for (int s = 0; s < snpsWithoutMapping.size(); s++) {
                        System.out.println(snpsWithoutMapping.get(s));
                    }
                    System.out.println("");
                }

                //Replace the SNP hashmap and vector.
                vecSNP.clear();
                hashSNP.clear();
                vecSNP = vecSNPNew;
                hashSNP = hashSNPNew;
            }

            //Write individuals to file:

            System.out.print("Writing individuals to file:\t");
            TextFile outInd = new TextFile(outputDirString + "Individuals.txt", TextFile.W);
            for (int ind = 0; ind < vecInd.size(); ind++) {
                String individual = ((String) vecInd.get(ind));
                outInd.write(individual + "\n");
            }
            outInd.close();
            System.out.println("OK");



            System.out.print("Writing SNPs to file:\t");
            TextFile outSNP = new TextFile(outputDirString + "SNPs.txt", TextFile.W);
            for (int snp = 0; snp < vecSNP.size(); snp++) {
                outSNP.write(((String) vecSNP.get(snp)) + "\n");

            }
            outSNP.close();
            System.out.println("OK");


        } else {

            //Load individuals from file:
            vecInd.clear();
            hashInd.clear();

            TextFile inInd = new TextFile(outputDirString + "Individuals.txt", TextFile.R);
            while ((str = inInd.readLine()) != null) {
                hashInd.put(str, vecInd.size());
                vecInd.add(str);
            }
            inInd.close();

            //Load SNPs from file:
            vecSNP.clear();
            hashSNP.clear();

            TextFile inSNP = new TextFile(outputDirString + "SNPs.txt", TextFile.R);
            while ((str = inSNP.readLine()) != null) {
                hashSNP.put(str, vecSNP.size());
                vecSNP.add(str);
            }
            inSNP.close();
        }



        int nrInds = vecInd.size();
        int nrSNPs = vecSNP.size();

        //We now have inventorized the file and have generated the SNPs.txt and Individuals.txt files.
        //Now try to determine the order of genotypes, so we can chose a buffering technique. If no order can be found we do not buffer, but importing will be extremely slow.
        boolean fileOrderPerSampleAllSNPs = false;
        boolean fileOrderPerSNPAllSamples = false;

        //Try to open the input file:
        in = new TextFile(inputFile, TextFile.R);

        //If this is an Illumina Final Report file, first process the irrelevant header:
        str = null;
        if (isIlluminaFinalReportFile) {
            while ((str = in.readLine()) != null) {
                data = str.split(delimiter);
                if (data[0].trim().equals("[Data]")) {
                    break;
                }
            }
        }

        //Now parse the column identifiers:
        str = in.readLine();

        //Check whether we actually are dealing with a Final Report file, user might have forgotten to instruct this:
        if (str.toLowerCase().startsWith("[header]")) {
            while ((str = in.readLine()) != null) {
                data = str.split(delimiter);
                if (data[0].trim().equals("[Data]")) {
                    break;
                }
            }
            str = in.readLine();
        }

        data = str.split(delimiter);
        int previousIndID = -1;
        int previousSNPID = -1;
        while ((str = in.readLine()) != null) {
            if (str.indexOf("\"") != -1) {
                str.replaceAll("\"", "");
            }
            if (str.indexOf("\'") != -1) {
                str.replaceAll("\'", "");
            }
            data = str.split(delimiter);
            String snp = data[columnSNP];
            String ind = data[columnSample];
            int snpID = hashSNP.get(snp);
            int indID = hashInd.get(ind);
            if (previousIndID != -1 && previousSNPID != -1) {
                if (snpID == previousSNPID && indID != previousIndID) {
                    fileOrderPerSNPAllSamples = true;
                    System.out.println("Based on the import file, TriTyper Importer assumes that the order of the file is such that for each SNP all samples are underneath each other in the import file. This assumptions increases importing performance.");
                }
                if (snpID != previousSNPID && indID == previousIndID) {
                    fileOrderPerSampleAllSNPs = true;
                    System.out.println("Based on the import file, TriTyper Importer assumes that the order of the file is such that for each sample all SNPs are underneath each other in the import file. This assumptions increases importing performance.");
                }
                break;
            }
            previousIndID = indID;
            previousSNPID = snpID;
        }




        System.out.print("Initializing binary data files:\t");
        RandomAccessFile file = new RandomAccessFile(outputDirString + "GenotypeMatrix.dat", "rw");
        RandomAccessFile fileRawData = null;
        if (rawDataAvailable) {
            fileRawData = new RandomAccessFile(outputDirString + "RawDataMatrix.dat", "rw");
        }
        System.out.println("OK");

        //Fill files with zeros:
        long size = (long) vecSNP.size() * (long) vecInd.size();
        long sizeGenotypeMatrix = size * 2;
        long sizeRawDataMatrix = size * 3;

        //Set size of files:
        file.setLength(0);
        if (rawDataAvailable) {
            fileRawData.setLength(0);
        }

        System.out.print("Making binary files zero:\t");
        //Quickly fill using buffers:
        file.seek(0);
        if (rawDataAvailable) {
            fileRawData.seek(0);
        }
        byte[] emptyString = new byte[10000];
        for (int s = 0; s < 10000; s++) {
            emptyString[s] = 0;
        }
        for (long a = 0; a < size / 10000; a++) {
            file.write(emptyString);
            file.write(emptyString);
            if (rawDataAvailable) {
                fileRawData.write(emptyString);
                fileRawData.write(emptyString);
                fileRawData.write(emptyString);
            }
        }

        //Fill rest with bytes:
        long rest = size % 10000;
        for (int a = 0; a < rest; a++) {
            byte emptyByte = 0;
            file.write(emptyByte);
            file.write(emptyByte);
            if (rawDataAvailable) {
                fileRawData.write(emptyByte);
                fileRawData.write(emptyByte);
                fileRawData.write(emptyByte);
            }
        }
        System.out.println("OK");

        System.out.println("Processing input file:");

        //Seek to beginning of file:
        file.seek(0);
        if (rawDataAvailable) {
            fileRawData.seek(0);
        }

        //Try to open the input file:

        in = new TextFile(inputFile, TextFile.R);

        //If this is an Illumina Final Report file, first process the irrelevant header:
        str = null;
        if (isIlluminaFinalReportFile) {
            while ((str = in.readLine()) != null) {
                data = str.split(delimiter);
                if (data[0].trim().equals("[Data]")) {
                    break;
                }
            }
        }

        //Now parse the column identifiers:
        str = in.readLine();

        //Check whether we actually are dealing with a Final Report file, user might have forgotten to instruct this:
        if (str.toLowerCase().startsWith("[header]")) {
            while ((str = in.readLine()) != null) {
                data = str.split(delimiter);
                if (data[0].trim().equals("[Data]")) {
                    break;
                }
            }
            str = in.readLine();
        }

        data = str.split(delimiter);

        //If the file has such an order that for each sample all SNPs are underneath each other, we use a buffering approach:
        byte[][] bufferAllele1 = null;
        byte[][] bufferAllele2 = null;
        byte[][] bufferR = null;
        byte[][] bufferTheta = null;
        int bufferFirstInd = 0;
        int bufferCurrentPos = 0;
        if (fileOrderPerSampleAllSNPs) {
            bufferAllele1 = new byte[nrSNPs][100];
            bufferAllele2 = new byte[nrSNPs][100];
            bufferR = new byte[nrSNPs][100];
            bufferTheta = new byte[nrSNPs][100];
        }
        if (fileOrderPerSNPAllSamples) {
            bufferAllele1 = new byte[1][nrInds];
            bufferAllele2 = new byte[1][nrInds];
            bufferR = new byte[1][nrInds];
            bufferTheta = new byte[1][nrInds];
        }

        //Start processing this file
        long linesProcessed = 0;
        previousIndID = -1;
        previousSNPID = -1;
        boolean warningGivenOnABGenotypeDefinition = false;
        while ((str = in.readLine()) != null) {

            //Remove quotes, if they exist:
            if (str.indexOf("\"") != -1) {
                str.replaceAll("\"", "");
            }
            if (str.indexOf("\'") != -1) {
                str.replaceAll("\'", "");
            }

            //Get individual values:
            data = str.split(delimiter);
            String snp = data[columnSNP];
            String ind = data[columnSample];
            double r = 0;
            double theta = 0;
            if (rawDataAvailable) {
                if (columnR != -1 && columnTheta != -1) {
                    if (data.length <= columnR || data.length <= columnTheta) {
                        System.out.println("\nError: For record entry " + (linesProcessed + 1) + " R or Theta values cannot be parsed! Record: " + str);
                        System.out.println("Can it be there are some entries in the file that do not have R or Theta value information?");
                        System.exit(-1);
                    }
                    String rString = data[columnR];
                    String thetaString = data[columnTheta];
                    if (!decimalSeparator.equals(".")) {
                        thetaString = thetaString.replaceAll(decimalSeparator, ".");
                        rString = rString.replaceAll(decimalSeparator, ".");
                    }
                    //Parse R value:
                    try {
                        r = Double.parseDouble(rString);
                    } catch (Exception e) {
                        System.out.println("\nError parsing R value: '" + rString + "'. Are you sure it has been saved in the correct locale?");
                        System.out.println("This method assumes R values have a decimal separator that is a dot.");
                        System.out.println("E.g. if you export a final report from within BeadStudio, using a Dutch Windows");
                        System.out.println("locale, the eventual final report file uses a comma as decimal separator.");
                        System.out.println("In that case use option '-decimalseparatoriscomma'");
                        System.exit(-1);
                    }
                    //Parse Theta value:
                    try {
                        theta = Double.parseDouble(thetaString);
                    } catch (Exception e) {
                        System.out.println("\nError parsing theta value: '" + thetaString + "'. Are you sure it has been saved in the correct locale?");
                        System.out.println("This method assumes theta values have a decimal separator that is a dot.");
                        System.out.println("E.g. if you export a final report from within BeadStudio, using a Dutch Windows");
                        System.out.println("locale, the eventual final report file uses a comma as decimal separator.");
                        System.out.println("In that case use option '-decimalseparatoriscomma'");
                        System.exit(-1);
                    }
                } else {
                    if (data.length <= columnX || data.length <= columnY) {
                        System.out.println("\nError: For record entry " + (linesProcessed + 1) + " X or Y intensities cannot be parsed! Record: " + str);
                        System.out.println("Can it be there are some entries in the file that do not have X or Y intensity information?");
                        System.exit(-1);
                    }
                    String xString = data[columnX];
                    String yString = data[columnY];
                    if (!decimalSeparator.equals(".")) {
                        xString = xString.replaceAll(decimalSeparator, ".");
                        yString = yString.replaceAll(decimalSeparator, ".");
                    }
                    double x = 0;
                    double y = 0;
                    try {
                        x = Double.parseDouble(xString);
                    } catch (Exception e) {
                        System.out.println("\nError parsing X value: '" + xString + "'. Are you sure it has been saved in the correct locale?");
                        System.out.println("This method assumes X values have a decimal separator that is a dot.");
                        System.out.println("E.g. if you export a final report from within BeadStudio, using a Dutch Windows");
                        System.out.println("locale, the eventual final report file uses a comma as decimal separator.");
                        System.out.println("In that case use option '-decimalseparatoriscomma'");
                        System.exit(-1);
                    }
                    try {
                        y = Double.parseDouble(yString);
                    } catch (Exception e) {
                        System.out.println("\nError parsing Y value: '" + yString + "'. Are you sure it has been saved in the correct locale?");
                        System.out.println("This method assumes Y values have a decimal separator that is a dot.");
                        System.out.println("E.g. if you export a final report from within BeadStudio, using a Dutch Windows");
                        System.out.println("locale, the eventual final report file uses a comma as decimal separator.");
                        System.out.println("In that case use option '-decimalseparatoriscomma'");
                        System.exit(-1);
                    }
                    //r = Math.sqrt(x * x + y * y);
                    r = x + y;
                    theta = 1;
                    if (x > 0) {
                        theta = twoDividedByPI * Math.atan2(y, x);
                    }
                }

            }
            byte rByte = (byte) (Byte.MIN_VALUE + (Math.min(255d, r * 50d)));
            byte thetaByte = (byte) (Byte.MIN_VALUE + (theta * 200d));

            //Inspect genotype calls, these either should be A, C, G or T, - will become 0:
            byte allele1 = data[columnAllele1].getBytes()[0];
            byte allele2 = data[columnAllele2].getBytes()[0];
            if (allele1 == 45) {
                allele1 = 0;
            }
            if (allele2 == 45) {
                allele2 = 0;
            }
            if (allele1 == 66) {
                allele1 = 67;
                if (!warningGivenOnABGenotypeDefinition) {
                    warningGivenOnABGenotypeDefinition = true;
                    System.out.println("\n\n\nWarning! The input genotype report file contains alleles that have been coded as B! These will be changed to C, please take this into account!!!\n\n\n");
                }
            }
            if (allele2 == 66) {
                allele2 = 67;
                if (!warningGivenOnABGenotypeDefinition) {
                    warningGivenOnABGenotypeDefinition = true;
                    System.out.println("\n\n\nWarning! The input genotype report file contains alleles that have been coded as B! These will be changed to C, please take this into account!!!\n\n\n");
                }
            }

            //Write data:
            int snpID = ((Integer) hashSNP.get(snp)).intValue();
            int indID = ((Integer) hashInd.get(ind)).intValue();

            if (fileOrderPerSampleAllSNPs || fileOrderPerSNPAllSamples) {
                if (fileOrderPerSampleAllSNPs) {
                    if (indID != previousIndID && previousIndID != -1) {
                        bufferCurrentPos++;
                    }
                    if (bufferCurrentPos == 100) {
                        //Flush buffer, hundred samples have just been processed
                        System.out.println("100 samples have been processed, flushing buffers:");
                        for (int s = 0; s < nrSNPs; s++) {
                            file.seek((long) s * (long) nrInds * 2 + (long) bufferFirstInd);
                            file.write(bufferAllele1[s]);
                            file.seek((long) s * (long) nrInds * 2 + (long) nrInds + (long) bufferFirstInd);
                            file.write(bufferAllele2[s]);
                            if (rawDataAvailable) {
                                fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) vecInd.size() + (long) bufferFirstInd);
                                fileRawData.write(bufferR[s]);
                                fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) 2 * vecInd.size() + (long) bufferFirstInd);
                                fileRawData.write(bufferTheta[s]);
                            }
                        }
                        bufferAllele1 = new byte[nrSNPs][100];
                        bufferAllele2 = new byte[nrSNPs][100];
                        bufferR = new byte[nrSNPs][100];
                        bufferTheta = new byte[nrSNPs][100];
                        bufferCurrentPos = 0;
                        bufferFirstInd = indID;
                    }
                    bufferAllele1[snpID][bufferCurrentPos] = allele1;
                    bufferAllele2[snpID][bufferCurrentPos] = allele2;
                    bufferR[snpID][bufferCurrentPos] = rByte;
                    bufferTheta[snpID][bufferCurrentPos] = thetaByte;
                } else {
                    if (snpID != previousSNPID && previousSNPID != -1) {
                        int s = previousSNPID;
                        file.seek((long) s * (long) nrInds * 2);
                        file.write(bufferAllele1[0]);
                        file.seek((long) s * (long) nrInds * 2 + (long) nrInds);
                        file.write(bufferAllele2[0]);
                        if (rawDataAvailable) {
                            fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) vecInd.size());
                            fileRawData.write(bufferR[0]);
                            fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) 2 * vecInd.size());
                            fileRawData.write(bufferTheta[0]);
                        }
                        bufferAllele1 = new byte[1][nrInds];
                        bufferAllele2 = new byte[1][nrInds];
                        bufferR = new byte[1][nrInds];
                        bufferTheta = new byte[1][nrInds];
                    }
                    bufferAllele1[0][indID] = allele1;
                    bufferAllele2[0][indID] = allele2;
                    bufferR[0][indID] = rByte;
                    bufferTheta[0][indID] = thetaByte;
                }
            } else {
                file.seek((long) snpID * (long) nrInds * 2 + (long) indID);
                file.write(allele1);
                file.seek((long) snpID * (long) nrInds * 2 + (long) nrInds + (long) indID);
                file.write(allele2);
                if (rawDataAvailable) {
                    fileRawData.seek((long) snpID * (long) vecInd.size() * 3 + (long) vecInd.size() + (long) indID);
                    fileRawData.write(rByte);
                    fileRawData.seek((long) snpID * (long) vecInd.size() * 3 + (long) 2 * vecInd.size() + (long) indID);
                    fileRawData.write(thetaByte);
                }
            }

            linesProcessed++;
            if (linesProcessed % 500000 == 0) {
                System.out.println(linesProcessed + "\tLines processed");
            }

            previousIndID = indID;
            previousSNPID = snpID;
        }

        if (fileOrderPerSampleAllSNPs || fileOrderPerSNPAllSamples) {
            if (fileOrderPerSampleAllSNPs) {
                //Flush remaining buffer:
                System.out.println("Flushing remaining buffer (" + (bufferCurrentPos + 1) + " samples):");
                for (int s = 0; s < nrSNPs; s++) {
                    byte[] bufferAllele1Subset = new byte[bufferCurrentPos + 1];
                    byte[] bufferAllele2Subset = new byte[bufferCurrentPos + 1];
                    byte[] bufferRSubset = new byte[bufferCurrentPos + 1];
                    byte[] bufferThetaSubset = new byte[bufferCurrentPos + 1];
                    for (int i = 0; i <= bufferCurrentPos; i++) {
                        bufferAllele1Subset[i] = bufferAllele1[s][i];
                        bufferAllele2Subset[i] = bufferAllele2[s][i];
                        bufferRSubset[i] = bufferR[s][i];
                        bufferThetaSubset[i] = bufferTheta[s][i];
                    }
                    file.seek((long) s * (long) nrInds * 2 + (long) bufferFirstInd);
                    file.write(bufferAllele1Subset);
                    file.seek((long) s * (long) nrInds * 2 + (long) nrInds + (long) bufferFirstInd);
                    file.write(bufferAllele2Subset);
                    if (rawDataAvailable) {
                        fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) vecInd.size() + (long) bufferFirstInd);
                        fileRawData.write(bufferRSubset);
                        fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) 2 * vecInd.size() + (long) bufferFirstInd);
                        fileRawData.write(bufferThetaSubset);
                    }
                }
            } else {
                //Flush remaining buffer:
                int s = previousSNPID;
                file.seek((long) s * (long) nrInds * 2);
                file.write(bufferAllele1[0]);
                file.seek((long) s * (long) nrInds * 2 + (long) nrInds);
                file.write(bufferAllele2[0]);
                if (rawDataAvailable) {
                    fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) vecInd.size());
                    fileRawData.write(bufferR[0]);
                    fileRawData.seek((long) s * (long) vecInd.size() * 3 + (long) 2 * vecInd.size());
                    fileRawData.write(bufferTheta[0]);
                }
            }
        }

        System.out.println(linesProcessed + "\tLines processed");

        //Close files:
        in.close();
        file.close();
        if (rawDataAvailable) {
            fileRawData.close();
        }



        //Output final remarks:
        System.out.println("Import of data has completed successfully!");
        System.out.println("");
        System.out.println("Please ensure you include a valid PhenotypeInformation.txt and SNPMappings.txt in the output directory.");
        System.out.println("These two additional files are required in order for TriTyper to function correctly.");
    }

    /**
     * @param file, the command line arguments
     */
    public static void main(String[] args) {

        if (Runtime.getRuntime().maxMemory() < 500 * 1024 * 1024) {
            System.out.println("Error! You have not dedicated at least 512 Mb to TriTyper Importer. Please use the -Xmx512m option to invoke TriTyper Importer!");
            System.exit(-1);
        }

        /*
         args = new String[4];
         args[0] = "-finalreportformat";
         args[1] = "-commadelimited";
         args[2] = "/Volumes/iPod/DMG/Data/CoeliacSasha/CD_Sasha_230107_FinalReport.csv";
         args[3] = "/Users/ludefranke/Documents/DMG/Manuscripts/CNVDeletionModel/TriTyper/";
         */

        /*
         args = new String[3];
         args[0] = "-tabdelimited";
         args[1] = "/Volumes/Backup/DMG/Data/Illumina_HumanCytoSNP-12v1.0/HumanCytoSNP-12_89CEU_FinalReport.txt";
         args[2] = "/Volumes/Backup/DMG/Data/Illumina_HumanCytoSNP-12v1.0/TriTyper/";
         */

//        args = new String[2];
//        args[0] = "/Data/GeneticalGenomicsDatasets/RotterdamStudy/GenotypeFinalReports/RS3_660array_1003samples_January2011_FinalReport.txt";
//        args[1] = "/Data/GeneticalGenomicsDatasets/RotterdamStudy/GenotypeFinalReports/TriTyper-660/";

        //args = new String[2];
        //args[0] = "/Volumes/Data/Lude/Data/CardioChip/EpicNL/EPIC_NL_111209_FinalReport_RTheta.txt";
        //args[1] = "/Volumes/Data/Lude/Data/CardioChip/EpicNL/";

        /*
         args = new String[2];
         args[0] = "/Volumes/Data/Lude/Data/GeneticalGenomics/Copd/copd_UMCG_april_2009_ALL_FinalReport.txt";
         args[1] = "/Volumes/Data/Lude/Data/GeneticalGenomics/Copd/TriTyper";
         */

        //args = new String[2];
        //args[0] = "/Volumes/Data/Lude/Data/MarcelWolfs/Marcel_Wolfs_290609_FinalReport.txt";
        //args[0] = "/Volumes/Data/Lude/Data/MarcelWolfs/FinalReports/Marcel Wolfs 190809_FinalReport.txt";
        //args[0] = "/Volumes/Data/Lude/Data/GeneticalGenomics/COPD/CNVTriTyper/FinalReport.txt";
        //args[0] = "/Volumes/Expression Backup 1/Lifelines/plaat1-2-3-4-5_100909_FinalReport_lude.txt";
        //args[0] = "/Volumes/Data/Lude/Data/CoeliacGWASII/UKCases/CNVTriTyper/FinalReport.txt";
        //args[0] = "/Volumes/Data/Lude/Data/MarcelWolfsFinal/OmniData/project Marcel_FinalReport.txt";

        //args[1] = "/Volumes/Data/Lude/Data/MarcelWolfs/TriTyper/";
        //args[1] = "/Volumes/Data/Lude/Data/MarcelWolfs/FinalReports/TriTyper/";
        //args[1] = "/Volumes/Data/Lude/Data/GeneticalGenomics/COPD/CNVTriTyper/";
        //args[1] = "/Volumes/Data/Lude/Data/CoeliacGWASII/UKCases/CNVTriTyper/";
        //args[0] = "/Volumes/Data/Lude/Data/RichardSinke/hanny_trijnie_irene_12062009_FinalReport.txt";
        //args[0] = "/Volumes/Data/Lude/Data/RichardSinke/4 sample_opnieuw180609_FinalReport.txt";
        //args[1] = "/Volumes/Data/Lude/Data/RichardSinke/4SampleOpnieuw/";
        //args[1] = "/Volumes/Data/Lude/Data/LifeLinesDNAQuantityAnalysis/";
        //args[1] = "/Volumes/Data/Lude/Data/MarcelWolfsFinal/OmniData/";


        if (args.length <= 1) {
            System.out.println("TriTyperImporter V1.0, 2008, Lude Franke, University Medical Centre Utrecht, Lude@ludesign.nl");
            System.out.println("\nUsage:\njava -jar trityperimporter.jar [-options] inputfile outputdir");
            System.out.println("\nThe default action is to process the input file, assuming it is tab delimited.");
            System.out.println("This file has to start with a header, describing which columns represent");
            System.out.println("the sample ID (Sample ID), SNP name (SNP Name), allele 1 genotype (Allele 1),");
            System.out.println("allele 2 genotype (Allele 2). If raw intensity data is available either the");
            System.out.println("R (R) and theta values (Theta) values should be included or X and Y columns.");
            System.out.println("");
            System.out.println("Files will be generated in 'outputdir' that can subsequently be used by TriTyper.");
            System.out.println("");
            System.out.println("It is recommended to use an Illumina final report file, as it is very easy to generate.");
            System.out.println("To use a final report file, please use the option: \"-finalreportformat\".");
            System.out.println("");
            System.out.println("To improve processing speed of TriTyper it is recommended to already have the");
            System.out.println("SNPMappings.txt file available in the output directory. This will sort the SNPs");
            System.out.println("per chromosome, which will considerably improve the speed of triallelic SNP discovery.");
            System.out.println("");
            System.out.println("Options:\t");
            System.out.println("-finalreportformat         Input file has been generated as a final report from within BeadStudio.");
            System.out.println("-tabdelimited              Input file is tab delimited (standard)");
            System.out.println("-commadelimited            Input file is comma delimited");
            System.out.println("-spacedelimited            Input file is space delimited");
            System.out.println("-semicolondelimited        Input file is semicolon delimited");
            System.out.println("-decimalseparatoriscomma   Within the input file decimals are separated with comma's, instead of dots (Applies e.g. to Dutch locale).");
            System.exit(0);
        }

        boolean isIlluminaFinalReportFile = false;
        String delimiter = "\t";
        String decimalSeparator = ".";
        for (int a = 0; a < args.length; a++) {
            String argument = args[a].trim();
            if (argument.startsWith("-commadelimited")) {
                delimiter = ",";
            }
            if (argument.startsWith("-tabdelimited")) {
                delimiter = "\t";
            }
            if (argument.startsWith("-spacedelimited")) {
                delimiter = " ";
            }
            if (argument.startsWith("-semicolondelimited")) {
                delimiter = ";";
            }
            if (argument.startsWith("-finalreportformat")) {
                isIlluminaFinalReportFile = true;
            }
            if (argument.startsWith("-decimalseparatoriscomma")) {
                decimalSeparator = ",";
            }
        }
        String inputFile = args[args.length - 2].trim();
        String outputDir = args[args.length - 1].trim();
        if (!outputDir.endsWith("/")) {
            outputDir += "/";
        }
        try {
            new FinalReportToTriTyper(inputFile, outputDir, isIlluminaFinalReportFile, delimiter, decimalSeparator);
        } catch (IOException ex) {
            Logger.getLogger(FinalReportToTriTyper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
