/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;

/**
 *
 * @author harmjan
 */
public class PedAndMapToTriTyper {

    public void importPEDFile(String[] mapFile, String mapFileDelimiter, int mapFileChrColumn, int mapFileChrPosColumn, int mapFileSNPColumn, String[] pedFile, String pedFileDelimiter, String outputDir, String[] casesToInclude) throws IOException {

        System.out.println("PED File importer v1.0");

        boolean ignoreFirstLinePEDFile = false;
        int genotypeColumnOffset = 6;
        if (mapFile == null) {
            ignoreFirstLinePEDFile = true;
            System.out.println("Assuming the SNPs that are present within the PED file(s) are described in the first line of the PED file:");

            HashMap<String, Integer> hashUniqueSNPs = new HashMap<String, Integer>();
            ArrayList<String> vecUniqueSNPs = new ArrayList<String>();
            System.out.println("\nProcessing first line of PED file(s):");
            mapFile = new String[pedFile.length];

            for (int pf = 0; pf < pedFile.length; pf++) {
                mapFile[pf] = outputDir + "MAPFile.txt";
                TextFile in = new TextFile(pedFile[pf], TextFile.R);
                String str2 = in.readLine();
				
				//this does not work when splitting on \\s
                //str2 = str2.replaceAll("\t", pedFileDelimiter);
				
                String data[] = str2.split(pedFileDelimiter);
                for (int c = 6; c < data.length; c++) {
                    if (!hashUniqueSNPs.containsKey(data[c])) {
                        hashUniqueSNPs.put(data[c], new Integer(vecUniqueSNPs.size()));
                        vecUniqueSNPs.add(data[c]);
                    }
                }
            }
            System.out.println("Number of unique SNPs parsed from file:\t" + vecUniqueSNPs.size());
            TextFile mapOut = new TextFile(outputDir + "MAPFile.txt", TextFile.W);
            for (int s = 0; s < vecUniqueSNPs.size(); s++) {
                mapOut.write(vecUniqueSNPs.get(s) + "\n");
            }
            mapOut.close();
            mapFileChrColumn = 0;
            mapFileChrPosColumn = 0;
            mapFileSNPColumn = 0;

            System.out.println("Please make sure you update SNPMappings.txt!!! This file is currently not correct, as no physical mapping information is available for each of the SNPs!");
        }

        //Check whether amount of map and ped files is identical:
        if (mapFile.length != pedFile.length) {
            System.out.println("Critical Error!\nThe number of entered map files and ped files is not identical!");
            System.exit(-1);
        }

        //Add slash if directory does not end with slash:
        if (!outputDir.endsWith("/")) {
            outputDir += "/";
        }

        //Check for existence of directories and files:
        File fileOutputDir = new File(outputDir);
        if (!fileOutputDir.isDirectory()) {
            System.out.println("Critical Error!\nThe output directory you have provided does not exist!:\t" + fileOutputDir.getAbsolutePath());
            System.exit(-1);
        }
        if (!fileOutputDir.canWrite()) {
            System.out.println("Critical Error!\nCannot write to output directory:\t" + outputDir);
            System.exit(-1);
        }

        for (int mf = 0; mf < mapFile.length; mf++) {
            if (!(new File(mapFile[mf])).canRead()) {
                System.out.println("Critical Error!\nCannot read from file:\t" + mapFile[mf]);
                System.exit(-1);
            }
            if (!(new File(pedFile[mf])).canRead()) {
                System.out.println("Critical Error!\nCannot read from file:\t" + pedFile[mf]);
                System.exit(-1);
            }
        }

        //Load markers:
        HashMap<String, Integer> hashSNP = new HashMap<String, Integer>();
        ArrayList<String> vectorSNP = new ArrayList<String>();
        ArrayList<String> vectorSNPMappings = new ArrayList<String>();
        int[] mapFileNrSNPs = new int[mapFile.length];
        HashMap[] hashMapFileSNPIndex = new HashMap[mapFile.length];
        for (int mf = 0; mf < mapFile.length; mf++) {
            hashMapFileSNPIndex[mf] = new HashMap();
        }

        System.out.println("\nDetermining number of makers per .map file:");

        for (int mf = 0; mf < mapFile.length; mf++) {
            TextFile inSNP = new TextFile(mapFile[mf], TextFile.R);
            String str2;
            while ((str2 = inSNP.readLine()) != null) {
                str2 = new String(str2.getBytes());
                String data[] = str2.split(mapFileDelimiter);
                if (str2.trim().length() > 0) {
                    //Store SNP:
                    hashSNP.put(data[mapFileSNPColumn], new Integer(vectorSNP.size()));
                    hashMapFileSNPIndex[mf].put(mapFileNrSNPs[mf], vectorSNP.size());
                    vectorSNP.add(data[mapFileSNPColumn]);
                    mapFileNrSNPs[mf]++;
                    //Store SNP Mapping:
                    String snpMapping = data[mapFileChrColumn] + "\t" + data[mapFileChrPosColumn] + "\t" + data[mapFileSNPColumn];
                    vectorSNPMappings.add(snpMapping);
                }
            }
            System.out.println("Number of markers in map file '" + mapFile[mf] + "':\t" + mapFileNrSNPs[mf]);
            inSNP.close();
        }
        System.out.println("Number of markers:\t" + vectorSNP.size() + "\t" + hashSNP.size());

        int nrSNPs = vectorSNP.size();

        //Load individuals:
        HashMap<String, Integer> hashInd = new HashMap<String, Integer>();
        ArrayList<String> vectorInd = new ArrayList<String>();

        System.out.println("\nProcessing individuals");
        TextFile out = new TextFile(outputDir + "PhenotypeInformation.txt", TextFile.W);

        for (int pf = 0; pf < pedFile.length; pf++) {
            System.out.println("Parsing: " + pedFile[pf]);

            TextFile in = new TextFile(pedFile[pf], TextFile.R);

            String str2 = "";
            if (ignoreFirstLinePEDFile) {
                str2 = in.readLine();
            }

            Pattern p = Pattern.compile(pedFileDelimiter);
            Pattern idmatcher = Pattern.compile("^([a-zA-Z_0-9\\S]*)\\s([\\Sa-zA-Z_0-9]*)\\s([a-zA-Z_0-9\\S]*)\\s([a-zA-Z_0-9\\S]*)\\s(-{0,1}\\d*\\.{0,1}\\d+)\\s(-{0,1}\\d*\\.{0,1}\\d+)\\s");
            while ((str2 = in.readLine()) != null) {

                Matcher m = idmatcher.matcher(str2);


                if (m.find()) {
                    String group = m.group(0);

                    String famid = m.group(1);
                    String samid = m.group(2);
                    String fatid = m.group(3);
                    String momid = m.group(4);
                    String sexid = m.group(5);
                    String affid = m.group(6);

                    String individual = famid + "-" + samid;
                    String sex = "male";
                    if (sexid.equals("2")) {
                        sex = "female";
                    }
                    String affectionStatus = "unknown";
                    String include = "include";

                    if (affid.equals("0")) {
                        include = "include";
                    }
                    if (affid.equals("1")) {
                        affectionStatus = "control";
                        include = "include";
                    }

                    for (int c = 0; c < casesToInclude.length; c++) {
                        if (affid.equals(casesToInclude[c])) {
                            affectionStatus = "case";
                            include = "include";
                        }
                    }
                    

                    if (!hashInd.containsKey(individual)) {
                        out.write(individual + "\t" + affectionStatus + "\t" + include + "\t" + sex + "\n");
                        hashInd.put(individual, new Integer(vectorInd.size()));
                        vectorInd.add(individual);
                    }


                } else {
                    System.out.println("Line does not match PED pattern: check the format in your file!");
                    System.out.println(str2);
                    System.exit(-1);
                }
            }
            System.out.println("Number of total individuals parsed:\t" + vectorInd.size());
            in.close();
        }
        out.close();

        System.out.println("Number of unique individuals in all datasets:\t" + vectorInd.size());

        int numSamples = vectorInd.size();

        //Write individuals file:
//        System.out.println("\nWriting individuals to file:");

        ProgressBar pb = new ProgressBar(vectorInd.size(), "Writing individuals to file: " + outputDir + "Individuals.txt");
        TextFile outInd = new TextFile(outputDir + "Individuals.txt", TextFile.W);
        for (int ind = 0; ind < vectorInd.size(); ind++) {
            outInd.write(vectorInd.get(ind) + "\n");

            pb.iterate();
        }
        pb.close();

        outInd.close();


        //Write individuals and SNPs file:
        pb = new ProgressBar(vectorSNP.size(), "Writing marker definition to file: " + outputDir + "SNPs.txt");
        TextFile outSNP = new TextFile(outputDir + "SNPs.txt", TextFile.W);
        for (int snp = 0; snp < vectorSNP.size(); snp++) {
            outSNP.write(vectorSNP.get(snp) + "\n");

            pb.iterate();
        }
        pb.close();
        outSNP.close();


        //Write individuals and SNPs file:
        pb = new ProgressBar(vectorSNP.size(), "Writing marker mapping definition to file: " + outputDir + "SNPMappings.txt");
        outSNP = new TextFile(outputDir + "SNPMappings.txt", TextFile.W);
        for (int snp = 0; snp < vectorSNP.size(); snp++) {
            outSNP.write(vectorSNPMappings.get(snp) + "\n");

            pb.iterate();
        }
        pb.close();
        outSNP.close();


        //Now write genotype data:
        int nrSamples = vectorInd.size();
        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, new File(outputDir + "GenotypeMatrix.dat"), false);

        System.out.println("\nLoading genotype data from PED file and writing into TriTyper format:");
        for (int pf = 0; pf < pedFile.length; pf++) {
            TextFile in = new TextFile(pedFile[pf], TextFile.R);
            System.out.println();
            String str2 = "";
            if (ignoreFirstLinePEDFile) {
                str2 = in.readLine();
            }
            byte conversionAllele[] = new byte[5];
            conversionAllele[0] = 0;
            conversionAllele[1] = 65;
            conversionAllele[2] = 67;
            conversionAllele[3] = 71;
            conversionAllele[4] = 84;

            boolean endOfFile = false;
            int individualCounterTotal = 0;

            int sampleBufferSize = 250;
            if (sampleBufferSize > vectorInd.size()) {
                sampleBufferSize = vectorInd.size();
            }
            pb = new ProgressBar(vectorInd.size(), "Parsing file:\t" + pedFile[pf]);
            while (!endOfFile) {
                int individualCounter = 0;
                byte[][] dataAllele1 = new byte[sampleBufferSize][mapFileNrSNPs[pf]];
                byte[][] dataAllele2 = new byte[sampleBufferSize][mapFileNrSNPs[pf]];


                for (int i = 0; i < sampleBufferSize; i++) {
                    // System.out.println("Loading individual:\t" + (individualCounter + individualCounterTotal));
                    str2 = in.readLine();
                    if (str2 != null) {
						
						//this does not work when splitting on \\s
//                        if (!pedFileDelimiter.equals("\t")) {
//                            while (str2.contains("\t")) {
//                                str2 = str2.replace("\t", pedFileDelimiter);
//                            }
//                        }

                        //str2 = str2.replaceAll("\t", pedFileDelimiter);
                        String data[] = str2.split(pedFileDelimiter);
//                        String individual = data[0] + "-" + data[1];

                        for (int snp = 0; snp < mapFileNrSNPs[pf]; snp++) {
                            byte allele1 = data[genotypeColumnOffset + snp * 2].getBytes()[0];
                            byte allele2 = data[genotypeColumnOffset + 1 + snp * 2].getBytes()[0];
                            if (allele1 == 48) {
                                allele1 = 0;
                            }
                            if (allele2 == 48) {
                                allele2 = 0;
                            }
                            if (allele1 == 73 || allele1 == 68) {
                                allele1 = 0;
                                //System.out.println("Strange D/I genotype for sample:\t" + individual + "\tSNP:\t" + snp + "\t" + (String) vectorSNPMappings.get(snp));
                            }
                            if (allele2 == 73 || allele2 == 68) {
                                allele2 = 0;
                                //System.out.println("Strange D/I genotype for sample:\t" + individual + "\tSNP:\t" + snp + "\t" + (String) vectorSNPMappings.get(snp));
                            }

                            if (allele1 == 78) {
                                allele1 = 0;
                            }

                            if (allele2 == 78) {
                                allele2 = 0;
                            }
                            if (allele1 != 0 && allele1 != 65 && allele1 != 67 && allele1 != 71 && allele1 != 84) {
                                //Assume PED file is coded in 0, 1, 2, 3, 4 format, convert to proper alleles:
                                allele1 = conversionAllele[Integer.parseInt(data[genotypeColumnOffset + snp * 2])];
                                allele2 = conversionAllele[Integer.parseInt(data[genotypeColumnOffset + 1 + snp * 2])];
                            }
                            dataAllele1[individualCounter][snp] = allele1;
                            dataAllele2[individualCounter][snp] = allele2;
                        }
                    } else {
                        endOfFile = true;
                        break;
                    }
                    individualCounter++;
                    //


                }

                //Write data to file:
                for (int snp = 0; snp < mapFileNrSNPs[pf]; snp++) {
                    byte[] allele1 = new byte[individualCounter];
                    byte[] allele2 = new byte[individualCounter];
                    for (int b = 0; b < individualCounter; b++) {
                        allele1[b] = dataAllele1[b][snp];
                        allele2[b] = dataAllele2[b][snp];
                    }
                    int snpIndex = ((Integer) hashMapFileSNPIndex[pf].get(snp)).intValue();
                    fileMatrixGenotype.setAllele1(snpIndex, individualCounterTotal, allele1);
                    fileMatrixGenotype.setAllele2(snpIndex, individualCounterTotal, allele2);

                }

                individualCounterTotal += individualCounter;
                pb.set(individualCounterTotal);
            }
//            int individualCounter = individualCounterTotal;
//            System.out.println("Number of individuals:\t" + individualCounter);
            in.close();
            pb.close();
        }


        fileMatrixGenotype.close();

        System.out.println("");
                
        System.out.println(
                  
                  "Import of PED files has completed.\n"
                + "----------------------------------\n"
                + "|          Please note:          |\n"
                + "----------------------------------\n"
                + "Your TriTyper files have been placed in the folder: "+outputDir + "\n"
                + "Please check whether the PhenotypeInformation.txt file in this folder reflects the data for your samples.\n"
                + "Have a nice day.\n");

    }

    public void importPEDFile(String dataLocation, String mapdelimiter, int chrcol, int chrpos, int snpcol, String peddelimiter, String outputLocation, String[] casesToInclude) throws IOException {
        String[] pedFile = Gpio.getListOfFiles(dataLocation, "ped");
        String[] mapFile = Gpio.getListOfFiles(dataLocation, "map");

        boolean filesok = true;
        if (pedFile.length == 0) {
            System.out.println("Error!: Directory does not contain any PED files");
            filesok = false;
        }

        if (mapFile.length == 0) {
            System.out.println("Error!: Directory does not contain any MAP files");
            filesok = false;
        }

        if (!filesok) {
            System.exit(-1);
        } else {
            importPEDFile(mapFile, mapdelimiter, chrcol, chrpos, snpcol, pedFile, peddelimiter, outputLocation, casesToInclude);
        }
    }

    public void importPEDFile(String dataLocation, String outputLocation) throws IOException {
        String[] pedFile = Gpio.getListOfFiles(dataLocation, "ped");
        String[] mapFile = Gpio.getListOfFiles(dataLocation, "map");


        boolean filesok = true;
        if (pedFile.length == 0) {
            System.out.println("Error!: Directory does not contain any PED files");
            filesok = false;
        }

        if (mapFile.length == 0) {
            System.out.println("Error!: Directory does not contain any MAP files");
            filesok = false;
        }

        if (!filesok) {
            System.exit(-1);
        } else {
            String[] casesToInclude = {"2"};
            // try to find the file splitter
            String pedSplit = determineFileSplitter(pedFile[0]);

            String mapSplit = determineFileSplitter(mapFile[0]);
            if (pedSplit == null || mapSplit == null) {
                System.err.println("ERROR: could not split your ped or map file by whitespace");
                System.exit(-1);
            } else {

                if (!outputLocation.endsWith("/")) {
                    outputLocation += "/";
                }
                Gpio.createDir(outputLocation);

                importPEDFile(mapFile, mapSplit, 0, 3, 1, pedFile, pedSplit, outputLocation, casesToInclude);
            }

        }

    }

    private String determineFileSplitter(String file) throws IOException {
//        TextFile tf = new TextFile(file, TextFile.R);
//        String firstLine = tf.readLine(); // skip the header
//        String ln = tf.readLine();
//
//        if (firstLine == null) {
//            throw new IOException("Is your file " + file + " empty?");
//        } else {
//            if (ln == null) {
//                throw new IOException("Error trying to split second line");
//            } else {
//                String[] elems = ln.split("\t");
//                if (elems.length > 1) {
//                    System.out.println(file + "\tis split by tab. Assuming possible other files of this type are as well.");
//                    return "\t";
//                }
//                elems = ln.split(" ");
//                if (elems.length > 1) {
//                    System.out.println(file + "\tis split by space. Assuming possible other files of this type are as well.");
//                    return " ";
//                }
//                elems = ln.split("\\s");
//                if (elems.length > 1) {
//                    System.out.println(file + "\tis split by whitespace. Assuming possible other files of this type are as well.");
//                    return "\\s";
//                }
//
//                tf.close();
//            }
		
		
		//ped and map do not enforce concistent use of tab or space within one fill. Splitting should always be performed with \\s
		return "\\s";
		
    }
//    private String[] getListOfFiles(String dataLocation, String extension) {
//        File loc = new File(dataLocation);
//        String[] fileList = loc.list();
//
//        ArrayList<String> al = new ArrayList<String>();
//
//        for (String fileloc : fileList) {
//            File file = new File(fileloc);
//            if (!file.isDirectory()) {
//                try {
//                    String fileName = file.getName();
//                    // System.out.println(fileloc+"-"+fileName);
//                    int mid = fileName.lastIndexOf(".");
//                    String fname = fileName.substring(0, mid);
//                    String ext = fileName.substring(mid + 1, fileName.length());
//
//
//                    // System.out.println(ext);
//                    if (ext.toLowerCase().equals(extension)) {
//                        al.add(loc.getAbsolutePath() + "/" + fileName);
//                    }
//                } catch (StringIndexOutOfBoundsException e) {
//                }
//            }
//        }
//
//        String[] output = new String[al.size()];
//        int i = 0;
//        for (String item : al) {
//            output[i] = al.get(i);
//            i++;
//        }
//
//        loc = null;
//
//        return output;
//    }
}
