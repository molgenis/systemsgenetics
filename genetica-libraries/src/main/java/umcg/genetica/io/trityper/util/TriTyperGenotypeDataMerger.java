/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;

/**
 *
 * @author harmjan
 */
public class TriTyperGenotypeDataMerger {

    public static void main(String[] args) {

        try {
            TriTyperGenotypeDataMerger merger = new TriTyperGenotypeDataMerger();

            String[] datasets = new String[2];


            datasets[0] = "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverCyto\\";
            datasets[1] = "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverCyto2\\";


//            for (int i = 0; i < 23; i++) {
//                int chr = i + 1;
//                if (chr == 23) {
//                    datasets[i] = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/ChrX/";
//                } else {
//                    datasets[i] = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Chr" + chr + "/";
//                }
//            }

            String outdir = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/SatVatLiverMuscle/LiverOmni/CytoAndOmniSampleMerge/";
//            merger.combinePrioritizerDatasetsMergeCommonSNPs(datasets[0],datasets[1], outdir, null);
            merger.mergeDatasetsOnCommonSamples(datasets, outdir);
            //merger.checkMerge(datasets, outdir);

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    // this merges two datasets on a set of common snps where both datasets contain different individuals
    public void combinePrioritizerDatasetsMergeCommonSNPs(String baseDir1, String baseDir2, String outputDir, String snps) throws IOException {
        //BROKEN!!!!
        System.out.println("TriTyper Dataset Combiner");
        System.out.println("\n\n");
        System.out.println("Starting to combine dataset '" + baseDir1 + "' and dataset '" + baseDir2 + "'");
        System.out.println("The output directory will be placed in '" + outputDir + "'");
        System.out.println("");

        outputDir = Gpio.formatAsDirectory(outputDir);
        Gpio.createDir(outputDir);
        
        HashSet<String> hashSNPsConfine = new HashSet<String>();
        if (snps != null) {
            System.out.println("Loading snp file from " + snps);
            TextFile snpfile = new TextFile(snps, TextFile.R);
            String[] snpstoquery = snpfile.readAsArray();
            hashSNPsConfine.addAll(Arrays.asList(snpstoquery));
            snpfile.close();
            System.out.println("Will merge at most " + hashSNPsConfine.size() + " snps.");
        }

        System.out.println("\nLoading data from dataset 1:");
        TriTyperGenotypeData genotypeDataset1 = new TriTyperGenotypeData();
        genotypeDataset1.load(baseDir1);


        System.out.println("");
        System.out.println("\nLoading data from dataset 2:");
        TriTyperGenotypeData genotypeDataset2 = new TriTyperGenotypeData();
        genotypeDataset2.load(baseDir2);

        System.out.println("\n\n");
        ArrayList<String> vectorSNP = new ArrayList<String>();


        String[] snps1 = genotypeDataset1.getSNPs();

        for (int snpID1 = 0; snpID1 < snps1.length; snpID1++) {

            String rsName = snps1[snpID1];
            Integer snp2Id = genotypeDataset2.getSnpToSNPId().get(rsName);
            if (snp2Id != -9) {
                if (hashSNPsConfine.isEmpty()) {

                    vectorSNP.add(rsName);
                } else {
                    if (hashSNPsConfine.contains(rsName)) {
                        vectorSNP.add(rsName);
                    }
                }
            }
        }
        System.out.println("Number of unique SNPs or probes that are present in both datasets, and will be included in combined dataset:\t" + vectorSNP.size());
        System.out.println("\n\n");

        System.out.println("\nCombining phenotype information files:");
        HashMap<String, Integer> hashInd = new HashMap<String, Integer>();
        ArrayList<String> vectorInd = new ArrayList<String>();


        TextFile phenotypeInformationOut = new TextFile(outputDir + "PhenotypeInformation.txt", TextFile.W);

        String[] inds1 = genotypeDataset1.getIndividuals();
        String[] inds2 = genotypeDataset2.getIndividuals();

        for (int ind = 0; ind < inds1.length; ind++) {
            String individual = inds1[ind];
            String sex = "male";

            if (genotypeDataset1.getIsFemale()[ind] == null) {
                System.out.println(individual + " is missing phenotype status. Exiting!");
                System.exit(0);
            }

            if (genotypeDataset1.getIsFemale()[ind]) {
                sex = "female";
            }


            String affectionStatus = "unknown";


            if (genotypeDataset1.getIsCase()[ind] == null){
               affectionStatus = "unknown";
            } else if (genotypeDataset1.getIsCase()[ind]) {
                affectionStatus = "case";
            } else {
                affectionStatus = "control";
            }

            String include = "include";

            if (!genotypeDataset1.getIsIncluded()[ind]) {
                include = "exclude";
            }


            phenotypeInformationOut.write(individual + "\t" + affectionStatus + "\t" + include + "\t" + sex + "\n");
            hashInd.put(individual, new Integer(vectorInd.size()));
            vectorInd.add(individual);
        }

        for (int ind = 0; ind < inds2.length; ind++) {
            String individual = inds2[ind];
            String sex = "male";

            if (genotypeDataset2.getIsFemale()[ind] == null) {
                System.out.println(individual + " is missing phenotype status. Exiting!");
                System.exit(0);
            }
            if (genotypeDataset2.getIsFemale()[ind]) {
                sex = "female";
            }
            String affectionStatus = "unknown";
            if (genotypeDataset2.getIsCase()[ind]) {
                affectionStatus = "case";
            } else {
                affectionStatus = "control";
            }

            String include = "include";

            if (!genotypeDataset2.getIsIncluded()[ind]) {
                include = "exclude";
            }
            phenotypeInformationOut.write(individual + "\t" + affectionStatus + "\t" + include + "\t" + sex + "\n");
            hashInd.put(individual, new Integer(vectorInd.size()));
            vectorInd.add(individual);
        }

        System.out.println("Total number of individuals:\t" + vectorInd.size());

        phenotypeInformationOut.close();


        int numSamples = vectorInd.size();

        //Write individuals file:
        System.out.println("\nWriting combined individuals to file:");

        TextFile outInd = new TextFile(outputDir + "Individuals.txt", TextFile.W);
        for (int ind = 0; ind < vectorInd.size(); ind++) {
            outInd.write(((String) vectorInd.get(ind)) + "\n");

            if (ind % 5 == 4) {
                System.out.print(".");
            }
        }
        System.out.println("");
        outInd.close();


        //Write individuals and SNPs file:
        System.out.println("\nWriting unique SNP / probe definition to file:");

        TextFile outSNP = new TextFile(outputDir + "SNPs.txt", TextFile.W);
        for (int snp = 0; snp < vectorSNP.size(); snp++) {
            outSNP.write(vectorSNP.get(snp) + "\n");
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        outSNP.close();


        byte[] complementAllele = new byte[256];
        complementAllele[84] = 65;
        complementAllele[65] = 84;
        complementAllele[67] = 71;
        complementAllele[71] = 67;

        //Open the genotypes file:
        System.out.println("\nCombining genotype and raw data from both analyses:");


        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(vectorSNP.size(), numSamples, new File(outputDir + "GenotypeMatrix.dat"), false);

        SNPLoader loader1 = genotypeDataset1.createSNPLoader();
        SNPLoader loader2 = genotypeDataset2.createSNPLoader();

        for (int x = 0; x < vectorSNP.size(); x++) {

            if (x % 10000 == 0) {
                System.out.println("Number of unique SNPs parsed so far:\t" + x);
            }

            String rsName = vectorSNP.get(x);
            int snpID1 = genotypeDataset1.getSnpToSNPId().get(rsName);
            int snpID2 = genotypeDataset2.getSnpToSNPId().get(rsName);
            //System.out.println(rsName + "\t" + snpID1 + "\t" + snpID2);
            SNP snpDataObject1 = genotypeDataset1.getSNPObject(snpID1);
            loader1.loadGenotypes(snpDataObject1);

            SNP snpDataObject2 = genotypeDataset2.getSNPObject(snpID2);
            loader2.loadGenotypes(snpDataObject2);

            //Determining minor allele dataset 1:
            byte[] alleles = new byte[4];
            int alleleItr = 0;
            int alleleCount[][] = new int[2][3];
            for (int individualID = 0; individualID < genotypeDataset1.getIndividuals().length; individualID++) {
                byte allele1Byte = snpDataObject1.getAllele1()[individualID];
                byte allele2Byte = snpDataObject1.getAllele2()[individualID];
                if (allele1Byte != 0 && allele2Byte != 0 && allele1Byte != 48 && allele2Byte != 48) {
                    int allelecode1 = -1;
                    int allelecode2 = -1;
                    for (int i = 0; i < 2; i++) {
                        if (alleles[i] == allele1Byte) {
                            allelecode1 = i;
                        }
                    }
                    if (allelecode1 == -1) {
                        alleles[alleleItr] = allele1Byte;
                        allelecode1 = alleleItr;
                        alleleItr++;
                    }

                    for (int i = 0; i < 2; i++) {
                        if (alleles[i] == allele2Byte) {
                            allelecode2 = i;
                        }
                    }
                    if (allelecode2 == -1) {
                        alleles[alleleItr] = allele2Byte;
                        allelecode2 = alleleItr;
                        alleleItr++;
                    }
                    alleleCount[0][allelecode2]++;
                    alleleCount[0][allelecode1]++;
                }
            }
            //System.out.println(rsName + "\tDataset 1, 1\t" + alleles[0] + "\t2\t" + alleles[1] + "\t3\t" + alleles[2] + "\t4\t" + alleles[3]);

            //Determine whether there are alleles in dataset 2 that are different from alleles in dataset 1:
            boolean allelesDifferent = false;
            int oldAlleleItr = alleleItr;
            //alleleItr = 0; alleles[0] = 0; alleles[1] = 0;
            for (int individualID = 0; individualID < genotypeDataset2.getIndividuals().length; individualID++) {
                byte allele1Byte = snpDataObject2.getAllele1()[individualID];
                byte allele2Byte = snpDataObject2.getAllele2()[individualID];
                if (allele1Byte != 0 && allele2Byte != 0 && allele1Byte != 48 && allele2Byte != 48) {
                    int allelecode1 = -1;
                    int allelecode2 = -1;
                    for (int i = 0; i < alleleItr; i++) {
                        if (alleles[i] == allele1Byte) {
                            allelecode1 = i;
                        }
                    }
                    if (allelecode1 == -1) {
                        alleles[alleleItr] = allele1Byte;
                        allelecode1 = alleleItr;
                        alleleItr++;
                        if (alleleItr > 2) {
                            allelesDifferent = true;
                        }
                    }
                    for (int i = 0; i < alleleItr; i++) {
                        if (alleles[i] == allele2Byte) {
                            allelecode2 = i;
                        }
                    }
                    if (allelecode2 == -1) {
                        alleles[alleleItr] = allele2Byte;
                        allelecode2 = alleleItr;
                        alleleItr++;
                        if (alleleItr > 2) {
                            allelesDifferent = true;
                        }
                    }
                }
            }
            alleleItr = oldAlleleItr;
            //System.out.println(rsName + "\tDataset 2, 1\t" + alleles[0] + "\t2\t" + alleles[1] + "\t3\t" + alleles[2] + "\t4\t" + alleles[3]);

            //Determining minor allele dataset 2:
            boolean alleleDifferenceError = false;
            //if (allelesDifferent) System.out.println(rsName + "\tAlleles are different, probably TOP / BOT issue, taking complement alleles.");
            for (int individualID = 0; individualID < genotypeDataset2.getIndividuals().length; individualID++) {
                byte allele1Byte = snpDataObject2.getAllele1()[individualID];
                byte allele2Byte = snpDataObject2.getAllele2()[individualID];
                //System.out.println(rsName + "\t" + individualID + "\t" + allele1Byte[0] + "\t" + allele2Byte[0] + "\t" + alleleItr + "\t" + alleles[0] + "\t" + alleles[1] + "\t" + alleles[2]);
                if (allele1Byte != 0 && allele2Byte != 0 && allele1Byte != 48 && allele2Byte != 48) {
                    if (allelesDifferent) {
                        allele1Byte = complementAllele[allele1Byte];
                        allele2Byte = complementAllele[allele2Byte];
                    }
                    int allelecode1 = -1;
                    int allelecode2 = -1;
                    for (int i = 0; i < 2; i++) {
                        if (alleles[i] == allele1Byte) {
                            allelecode1 = i;
                        }
                    }
                    if (allelecode1 == -1) {
                        alleles[alleleItr] = allele1Byte;
                        allelecode1 = alleleItr;
                        alleleItr++;
                        if (alleleItr > 2) {
                            alleleDifferenceError = true;
                            break;
                        }
                    }

                    for (int i = 0; i < 2; i++) {
                        if (alleles[i] == allele2Byte) {
                            allelecode2 = i;
                        }
                    }
                    if (allelecode2 == -1) {
                        alleles[alleleItr] = allele2Byte;
                        allelecode2 = alleleItr;
                        alleleItr++;
                        if (alleleItr > 2) {
                            alleleDifferenceError = true;
                            break;
                        }
                    }
                    alleleCount[1][allelecode2]++;
                    alleleCount[1][allelecode1]++;
                }
            }
            double maf1 = (double) alleleCount[0][0] / (double) (alleleCount[0][0] + alleleCount[0][1]);
            double maf2 = (double) alleleCount[1][0] / (double) (alleleCount[1][0] + alleleCount[1][1]);

            if (alleleDifferenceError) {
                System.out.println("\nError! SNP\t" + rsName + "\thas more than two different alleles in the two datasets, excluding this SNP!!!");
                System.out.println("Please ensure that you used the same genotype allele naming convention within BeadStudio.");
                System.out.println("It is highly recommended to only use datasets that have been generated based on the Illumina 'TOP Allele' naming convention.");
                alleleCount[0][0] = 0;
                alleleCount[1][0] = 0;
                alleleCount[0][1] = 0;
                alleleCount[1][1] = 0;
            } else {
                if (alleleCount[0][0] > alleleCount[0][1]) {
                    if (alleleCount[1][0] < alleleCount[1][1]) {
                        if (maf1 / maf2 > 1.5 || maf2 / maf1 > 1.5) {
                            System.out.println("Warning: " + rsName + " has quite different allele frequency between dataset 1 (Allele freq. = " + maf1 + ") and dataset 2 (Same allele freq. = " + maf2 + ")");
                        }
                    }
                }
                byte[] allele1 = new byte[numSamples];
                byte[] allele2 = new byte[numSamples];
//		byte[] gcScores = new byte[numSamples];
//		byte[] rValues = new byte[numSamples];
//		byte[] thetaValues = new byte[numSamples];
                for (int y = 0; y < numSamples; y++) {
                    if (y < genotypeDataset1.getIndividuals().length) {
                        allele1[y] = snpDataObject1.getAllele1()[y];
                        allele2[y] = snpDataObject1.getAllele2()[y];
//			if (fileMatrixRawData != null) {
//			    gcScores[y] = snpDataObject1.getGcScores()[y];
//			    rValues[y] = snpDataObject1.getRValues()[y];
//			    thetaValues[y] = snpDataObject1.getThetaValues()[y];
//			}
                    } else {
                        allele1[y] = snpDataObject2.getAllele1()[y - genotypeDataset1.getIndividuals().length];
                        allele2[y] = snpDataObject2.getAllele2()[y - genotypeDataset1.getIndividuals().length];
                        if (allelesDifferent) {
                            allele1[y] = complementAllele[allele1[y]];
                            allele2[y] = complementAllele[allele2[y]];
                        }
//			if (fileMatrixRawData != null) {
//			    gcScores[y] = snpDataObject2.getGcScores()[y - genotypeDataset1.getNumIndividuals()];
//			    rValues[y] = snpDataObject2.getRValues()[y - genotypeDataset1.getNumIndividuals()];
//			    thetaValues[y] = snpDataObject2.getThetaValues()[y - genotypeDataset1.getNumIndividuals()];
//			}
                    }
                }
                fileMatrixGenotype.setAllele1(x, 0, allele1);
                fileMatrixGenotype.setAllele2(x, 0, allele2);
//		if (fileMatrixRawData != null) {
//		    fileMatrixRawData.setGCScore(x, 0, gcScores);
//		    fileMatrixRawData.setR(x, 0, rValues);
//		    fileMatrixRawData.setTheta(x, 0, thetaValues);
//		}
            }
        }
        fileMatrixGenotype.close();
//	if (fileMatrixRawData != null) {
//	    fileMatrixRawData.close();
//	}
        System.out.println("Number of unique SNPs parsed in total:t" + vectorSNP.size());

        System.out.println("\nCombining of TriTyper datasets has finished. Please ensure the errors and warnings that might have been observed are acceptable.");
        System.out.println("\nPlease ensure you copy a SNPMappings.txt file to the outputfolder holding the combined dataset.");
        System.out.println("After copying this file, you should be able to run TriTyper on these two combined datasets.");
    }

    // this merges non-common SNPs for the same samples across different (for example imputed) datasets.
    public void mergeDatasetsOnCommonSamples(String[] datasetLocations, String outputDir) throws IOException {
        if (datasetLocations.length < 2) {
            throw new IllegalArgumentException("Error: Nothing to combine, only one dataset presented to the program.");
        }
        if (outputDir == null) {
            throw new IllegalArgumentException("Error: No outputdir selected.");
        }
        if (!outputDir.endsWith("/")) {
            outputDir += "/";
        }
        Gpio.createDir(outputDir);

        TextFile log = new TextFile(outputDir + "log.txt", TextFile.W);
        // check the input folders.
        // creating a new GenotypeDataObject implicitly performs these checks.
        TriTyperGenotypeData[] ds = new TriTyperGenotypeData[datasetLocations.length];
        for (int d = 0; d < ds.length; d++) {
            log.writeln("Loading\t: " + datasetLocations[d]);
            System.out.println("Loading\t: " + datasetLocations[d]);
            ds[d] = new TriTyperGenotypeData();
            ds[d].load(datasetLocations[d]);
            log.writeln("SNPs\t: " + ds[d].getSNPs().length);
            log.writeln("Samples\t: " + ds[d].getIndividuals().length);
            log.writeln();
        }

        // check whether some of the datasets have duplicate samples
        HashSet<String> duplicateIndividualsOverAllDatasets = new HashSet<String>();
        for (int d = 0; d < ds.length; d++) {
            HashSet<String> visitedIndividuals = new HashSet<String>();
            String[] individualsInDatasetD = ds[d].getIndividuals();
            for (int indId = 0; indId < individualsInDatasetD.length; indId++) {
                String s = individualsInDatasetD[indId];
                if (visitedIndividuals.contains(s)) {
                    duplicateIndividualsOverAllDatasets.add(s);
                    log.writeln("Dataset\t" + d + "\tcontains duplicate sample which will be excluded:\t" + s);
                    System.out.println("Dataset\t" + d + "\tcontains duplicate sample which will be excluded:\t" + s);
                }
                visitedIndividuals.add(s);
            }
        }

        // now check whether there are shared individuals over all datasets.
        HashMap<String, Integer> individualCounterAcrossDatasets = new HashMap<String, Integer>();
        HashSet<String> uniqueIndividualsInAllDatasets = new HashSet<String>();
        for (int d = 0; d < ds.length; d++) {
            String[] individualsInDatasetD = ds[d].getIndividuals();
            for (int indId = 0; indId < individualsInDatasetD.length; indId++) {
                String ind = individualsInDatasetD[indId];

                // check whether the individual is included and not a duplicate in any dataset
                if (ds[d].getIsIncluded()[indId] && !duplicateIndividualsOverAllDatasets.contains(ind)) {
                    Integer ctr = individualCounterAcrossDatasets.get(ind);
                    if (ctr == null) {
                        ctr = 0;
                    }
                    ctr++;
                    individualCounterAcrossDatasets.put(ind, ctr);
                    uniqueIndividualsInAllDatasets.add(ind);
                }
            }
        }
        
        duplicateIndividualsOverAllDatasets = null;
        String[] availableIndividualsForAllDatasets = uniqueIndividualsInAllDatasets.toArray(new String[0]);

        // we need a map to convert individualnames to an index
        HashMap<String, Integer> includedIndividualsToId = new HashMap<String, Integer>();
        ArrayList<String> includedIndividuals = new ArrayList<String>();
        int indId = 0;
        for (String s : availableIndividualsForAllDatasets) {
            Integer ctr = individualCounterAcrossDatasets.get(s);
            if (ctr != null && ctr > 1 && ctr < ds.length + 1) {
                includedIndividualsToId.put(s, indId);
                includedIndividuals.add(s);
                indId++;
            } else {
                log.writeln("Individual\t" + s + "\tis excluded because it is either a duplicate or because it is present " + ctr + " times, while expected num == " + ds.length);
                System.out.println("Individual\t" + s + "\tis excluded because it is either a duplicate or because it is present " + ctr + " times, while expected num == " + ds.length);
            }
        }

        availableIndividualsForAllDatasets = null;
        individualCounterAcrossDatasets = null;
        uniqueIndividualsInAllDatasets = null;

        log.writeln("Number of samples that are shared by all datasets:\t" + includedIndividuals.size());
        System.out.println("Number of samples that are shared by all datasets:\t" + includedIndividuals.size());
        if (includedIndividuals.isEmpty()) {
            log.writeln("Nothing to merge, since no included dataset shares at least 1 individual");
            System.out.println("Nothing to merge, since no included dataset shares at least 1 individual");
            log.close();
            System.exit(0);
        }

        // now check whether any of the SNPs is present in more than one dataset (comparison to first dataset)
        // we want to exclude these duplicates, or at least warn the user about their presence..

        // now check for unique SNPs accross datasets

        HashSet<String> visitedSNPsAcrossDatasets = new HashSet<String>();
        HashSet<String> snpsDuplicateAcrossDatasets = new HashSet<String>();

        TextFile duplicateOut = new TextFile(outputDir + "DuplicateSNPWithinDatasets.txt", TextFile.W);
        duplicateOut.writeln("dataset\tsnp");
        String[][] snpsDuplicatePerDataset = new String[ds.length][0];
        for (int d = 0; d < ds.length; d++) {
            String[] snpsInDataset = ds[d].getSNPs();

            HashSet<String> snpsVisitedInThisDataset = new HashSet<String>();
            HashSet<String> duplicateSNPsInThisDataset = new HashSet<String>();
            for (String s : snpsInDataset) {
                if (!visitedSNPsAcrossDatasets.contains(s)) {
                    visitedSNPsAcrossDatasets.add(s);
                } else {
                    snpsDuplicateAcrossDatasets.add(s);
                }
                if (!snpsVisitedInThisDataset.contains(s)) {
                    snpsVisitedInThisDataset.add(s);
                } else {
                    duplicateOut.writeln(d + "\t" + s);
                    duplicateSNPsInThisDataset.add(s);
                }
            }
            snpsDuplicatePerDataset[d] = duplicateSNPsInThisDataset.toArray(new String[0]);
            System.out.println("Dataset\t" + d + "\thas\t" + duplicateSNPsInThisDataset.size() + "\tduplicate SNPs");
        }
        duplicateOut.close();

        for (int d = 0; d < ds.length; d++) {
            if (snpsDuplicatePerDataset[d].length > 0) {
                log.writeln("Dataset " + d + " has duplicate SNPs. Please remove them before you continue");
                System.out.println("Dataset " + d + " has duplicate SNPs. Please remove them before you continue");
                log.close();
                System.exit(-1);
            }
        }

        log.writeln("Duplicate SNPs accross datasets: " + snpsDuplicateAcrossDatasets.size());
        System.out.println("Duplicate SNPs accross datasets: " + snpsDuplicateAcrossDatasets.size());

        duplicateOut = new TextFile(outputDir + "DuplicateSNPAccrossDatasets.txt", TextFile.W);


        String header = "SNP";
        for (int d = 0; d < ds.length; d++) {
            header += "\t" + d + " Alleles\t" + d + " MinorAllele\t"+d+" MAF\t"+d+" HWEP\t"+d+" CR";
        }
        duplicateOut.writeln(header);

        SNPLoader[] loaders = new SNPLoader[ds.length];
        for (int d = 0; d < ds.length; d++) {
            loaders[d] = ds[d].createSNPLoader();
        }

        TextFile genotypeComp = new TextFile(outputDir + "GenotypeComparisonOfDuplicateSNPs.txt.gz", TextFile.W);

        genotypeComp.writeln("SNP\td1\tCR\tMAF\tHWEP\tAlleles\tMinorAllele\tflipallele1\td2\tCR\tMAF\tHWEP\tAlleles\tMinorAllele\tflipallele2\tcalled\tdifferent\tidentical");
        // make a lookup table for shared individuals
        int[][] indLookup = new int[ds.length][includedIndividualsToId.size()];
        for (int d = 0; d < ds.length; d++) {
            String[] inds = ds[d].getIndividuals();
            indLookup[d] = new int[inds.length];

            for (int i = 0; i < includedIndividuals.size(); i++) {
                indLookup[d][i] = ds[d].getIndividualToId().get(includedIndividuals.get(i));
            }
        }

        int nrWithInCompatibleAlleles = 0;

        HashMap<String, Integer> duplicateSNPsSelectFromThisDataset = new HashMap<String, Integer>();
        HashSet<String> duplicateSNPsThatShouldBeIncluded = new HashSet<String>();
        for (String s : snpsDuplicateAcrossDatasets) {
            SNP[] snps = new SNP[ds.length];
            StringBuilder output = new StringBuilder();
            for (int d = 0; d < ds.length; d++) {
                Integer snpId = ds[d].getSnpToSNPId().get(s);
                if (snpId != -9) {
                    snps[d] = ds[d].getSNPObject(snpId);
                    loaders[d].loadGenotypes(snps[d]);
                    output.append("\t").append(BaseAnnot.toString(snps[d].getAlleles()[0])).append("/")
                            .append(BaseAnnot.toString(snps[d].getAlleles()[1])).append("\t")
                            .append(BaseAnnot.toString(snps[d].getMinorAllele())).append("\t")
                            .append(snps[d].getMAF()).append("\t")
                            .append(snps[d].getHWEP()).append("\t")
                            .append(snps[d].getCR());
                }
            }
            
            duplicateOut.writeln(s + output.toString());


            // make a lookup index for the shaared individuals.

            // now see whether the genotypes correspond..
            Boolean[] flipAlleles = CompareAllelicDirections.compare(snps);
            if (flipAlleles == null) {
                String out = s;
                for (int d = 0; d < ds.length; d++) {
                    if (snps[d] != null) {
                        SNP snpObj = snps[d];
                        out += "\t" + d + "\t"
                                + BaseAnnot.toString(snpObj.getAlleles()[0]) + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                                + "\t" + BaseAnnot.toString(snpObj.getMinorAllele());
                    } else {
                        out += "\tNotPresentIn" + d;
                    }

                }
                nrWithInCompatibleAlleles++;
                genotypeComp.writeln(out);
            } else {
                short[] genotypes1 = new short[includedIndividualsToId.size()];
                short[] genotypes2 = new short[includedIndividualsToId.size()];
                Integer[] nrDifferentPerDataset = new Integer[ds.length];
                for (int d1 = 0; d1 < ds.length; d1++) {


                    SNP snp1 = snps[d1];
                    if (snp1 != null) {
                        byte[] gt1 = snp1.getGenotypes();
                        for (int i = 0; i < includedIndividuals.size(); i++) {
                            short gt = gt1[indLookup[d1][i]];
                            if (flipAlleles[d1]) {
                                if (gt == 0) {
                                    gt = 2;
                                } else if (gt == 2) {
                                    gt = 0;
                                }
                            }
                            genotypes1[i] = gt;
                        }
                        for (int d2 = d1 + 1; d2 < ds.length; d2++) {
                            SNP snp2 = snps[d2];
                            if (snp2 != null) {
                                byte[] gt2 = snp2.getGenotypes();
                                for (int i = 0; i < includedIndividuals.size(); i++) {
                                    short gt = gt2[indLookup[d2][i]];
                                    if (flipAlleles[d2]) {
                                        if (gt == 0) {
                                            gt = 2;
                                        } else if (gt == 2) {
                                            gt = 0;
                                        }
                                    }
                                    genotypes2[i] = gt;
                                }

                                // everything initialized. Now continue to comparison
                                int identical = 0;
                                int different = 0;
                                int called = 0;
                                for (int g = 0; g < genotypes1.length; g++) {
                                    short g1 = genotypes1[g];
                                    short g2 = genotypes2[g];
                                    if (g1 < 0 || g2 < 0) {
                                        // one or both is uncalled
                                    } else {
                                        called++;
                                        if (g1 == g2) {
                                            identical++;
                                        } else {
                                            different++;
                                        }
                                    }
                                }
                                StringBuilder snpStats1 = new StringBuilder();
                                snpStats1.append(snp1.getCR()).append("\t").append(snp1.getMAF()).append("\t").append(snp1.getHWEP()).append("\t").append(BaseAnnot.toString(snp1.getAlleles()[0])).append("/").append(BaseAnnot.toString(snp1.getAlleles()[1])).append("\t").append(BaseAnnot.toString(snp1.getMinorAllele()));
                                StringBuilder snpStats2 = new StringBuilder();
                                snpStats2.append(snp2.getCR()).append("\t").append(snp2.getMAF()).append("\t").append(snp2.getHWEP()).append("\t").append(BaseAnnot.toString(snp2.getAlleles()[0])).append("/").append(BaseAnnot.toString(snp2.getAlleles()[1])).append("\t").append(BaseAnnot.toString(snp2.getMinorAllele()));
                                genotypeComp.writeln(s + "\t" + d1 + "\t" + snpStats1.toString() + "\t" + flipAlleles[d1] + "\t" + d2 + "\t" + snpStats2.toString() + "\t" + flipAlleles[d2] + "\t" + called + "\t" + different + "\t" + identical);
                                nrDifferentPerDataset[d1] = different;
                            }
                        }
                    }
                }
                int compsMade = 0;
                double diffSum = 0;
                for (int d = 0; d < ds.length; d++) {
                    if (nrDifferentPerDataset[d] != null) {
                        diffSum += nrDifferentPerDataset[d];
                        compsMade++;
                    }
                }
                double ratio = diffSum / compsMade;
                if (ratio >= 1d) {
                    // exclude SNP
                } else {
                    // otherwise take the SNP with the best callrate
                    int maxCRDs = -1;
                    double maxCR = 0;
                    for (int d = 0; d < ds.length; d++) {
                        SNP snpObj = snps[d];
                        if (snpObj.getCR() > maxCR) {
                            maxCR = snpObj.getCR();
                            maxCRDs = d;
                        }
                    }
                    duplicateSNPsThatShouldBeIncluded.add(s);
                    duplicateSNPsSelectFromThisDataset.put(s, maxCRDs);
                }

            }
            for (int d = 0; d < ds.length; d++) {
                snps[d].clearGenotypes();
            }
        }
        
        System.out.println("Nr of SNPs with incompatible alleles: " + nrWithInCompatibleAlleles);
        genotypeComp.close();

        duplicateOut.close();

        HashSet<String> duplicateSNPs = new HashSet<String>();
        duplicateSNPs.addAll(Arrays.asList(snpsDuplicatePerDataset[0]));
        duplicateSNPs.addAll(Arrays.asList(snpsDuplicatePerDataset[1]));
        duplicateSNPs.addAll(snpsDuplicateAcrossDatasets);
        
        snpsDuplicatePerDataset = null;
        visitedSNPsAcrossDatasets = null;
        log.writeln("Total number of duplicate SNPs: " + duplicateSNPs.size());
        System.out.println("Total number of duplicate SNPs: " + duplicateSNPs.size());

        String[] listOfDuplicateSNPs = duplicateSNPs.toArray(new String[0]);
        log.writeList(Arrays.asList(listOfDuplicateSNPs));

        ArrayList<String> uniqueSNPs = new ArrayList<String>();
        boolean[] datasethasUniqueSNPs = new boolean[ds.length];
        for (int d = 0; d < ds.length; d++) {
            String[] snpsInDataset = ds[d].getSNPs();
            for (String s : snpsInDataset) {
                if (!duplicateSNPs.contains(s)) {
                    uniqueSNPs.add(s);
                    datasethasUniqueSNPs[d] = true;
                }
            }
        }

        TextFile snpsout = new TextFile(outputDir + "UniqueSNPsAccrossDatasets.txt", TextFile.W);
        String[] listOfSNPs = uniqueSNPs.toArray(new String[0]);
        snpsout.writeList(Arrays.asList(listOfSNPs));
        snpsout.close();

        log.writeln("Unique SNPs: " + uniqueSNPs.size() + "\tselected duplicate SNPs: " + duplicateSNPsSelectFromThisDataset.size() + "\ttotal: " + (uniqueSNPs.size() + duplicateSNPsSelectFromThisDataset.size()));
        System.out.println("Unique SNPs: " + uniqueSNPs.size() + "\tselected duplicate SNPs: " + duplicateSNPsSelectFromThisDataset.size() + "\ttotal: " + (uniqueSNPs.size() + duplicateSNPsSelectFromThisDataset.size()));

        if (uniqueSNPs.isEmpty()) {
            log.writeln("No unique SNPs detected. Nothing to merge");
            System.out.println("No unique SNPs detected. Nothing to merge");
            log.close();
            System.exit(0);
        }



        TextFile outInds = new TextFile(outputDir + "Individuals.txt", TextFile.W);
        TextFile outPheno = new TextFile(outputDir + "PhenotypeInformation.txt", TextFile.W);
        // write the individual files.
        for (String ind : includedIndividuals) {
            outInds.writeln(ind);
            for (int d = 0; d < ds.length; d++) {
                Integer id = ds[d].getIndividualId(ind);
                if (id != null) {
                    String gender = "male";
                    if (ds[d].getIsFemale()[id] == null) {
                        gender = "unknown";
                    } else if (ds[d].getIsFemale()[id]) {
                        gender = "female";
                    }

                    String status = "control";
                    if (ds[d].getIsCase()[id] == null) {
                        status = "unknown";
                    } else if (ds[d].getIsCase()[id]) {
                        status = "case";
                    }
                    outPheno.writeln(ind + "\t" + status + "\tinclude\t" + gender);
                    break; // stop looking for this ind.
                }
            }
        }
        
        outInds.close();
        outPheno.close();
        
        log.writeln();
        log.writeln("Final size of matrix will be: " + includedIndividuals.size() + " x " + (uniqueSNPs.size() + duplicateSNPsSelectFromThisDataset.size()));
        System.out.println("Final size of matrix will be: " + includedIndividuals.size() + " x " + (uniqueSNPs.size() + duplicateSNPsSelectFromThisDataset.size()));

        uniqueSNPs.addAll(duplicateSNPsThatShouldBeIncluded);
        
        // remove GenotypeDataMatrix.dat if it is already there
        
        File fileGenotypeMatrix = new File(outputDir + "GenotypeMatrix.dat");
        if (fileGenotypeMatrix.exists()){
            fileGenotypeMatrix.delete();
            fileGenotypeMatrix = new File(outputDir + "GenotypeMatrix.dat");
        }
        
        WGAFileMatrixGenotype newDS = new WGAFileMatrixGenotype(uniqueSNPs.size(), includedIndividuals.size(), fileGenotypeMatrix, false);
        //System.exit(0);
        TextFile snpout = new TextFile(outputDir + "SNPs.txt", TextFile.W);
        TextFile snpMappingsOut = new TextFile(outputDir + "SNPMappings.txt", TextFile.W);

        int nrSnps = uniqueSNPs.size();
        int nrSamplesIncluded = includedIndividuals.size();
        HashSet<String> snpsProcessed = new HashSet<String>();
        ProgressBar pb = new ProgressBar(nrSnps, "Parsing SNPs");
        int[] dataFromSet = new int[ds.length];
        for (int d = 0; d < ds.length; d++) {
            dataFromSet[d] = 0;
            TriTyperGenotypeData data = ds[d];
            if (datasethasUniqueSNPs[d]) {
                // create a quick lookup table for the individuals.
                String[] inds = data.getIndividuals();
                int nrSamplesInDataset = inds.length;
                Integer[] sampleToNewSampleId = new Integer[nrSamplesInDataset];
                for (int i = 0; i < inds.length; i++) {
                    //System.out.println(inds[i]+"\t"+i+"\t"+includedIndividualsToId.get(inds[i]));
                    sampleToNewSampleId[i] = includedIndividualsToId.get(inds[i]);
                }

                System.out.println("Now processing dataset " + d + " (" + datasetLocations[d] + ")");
                SNPLoader loader = loaders[d];
                for (int s = 0; s < nrSnps; ++s) {

                    String snpname = uniqueSNPs.get(s);

                    // check whether this snp is a duplicate accross datasets
                    // then check whether this SNP should be included from this dataset.

                    Integer snpId = null;
                    if (duplicateSNPsThatShouldBeIncluded.contains(snpname)) {
                        Integer dsId = duplicateSNPsSelectFromThisDataset.get(snpname);
                        if (dsId != null && dsId.equals(d)) {
                            snpId = data.getSnpToSNPId().get(snpname);
                        }
                    } else {
                        snpId = data.getSnpToSNPId().get(snpname);
                    }

                    if (snpId != null) {
                        dataFromSet[d]++;
                        snpsProcessed.add(snpname);
                        SNP snpObj = data.getSNPObject(snpId);
                        snpout.writeln(snpname);
                        snpMappingsOut.writeln(ChrAnnotation.parseByte(snpObj.getChr()) + "\t" + snpObj.getChrPos() + "\t" + snpname);

                        loader.loadGenotypes(snpObj);
                        byte[] allele1 = snpObj.getAllele1();
                        byte[] allele2 = snpObj.getAllele2();
                        byte[] outputallele1 = new byte[nrSamplesIncluded];
                        byte[] outputallele2 = new byte[nrSamplesIncluded];
                        for (int i = 0; i < nrSamplesInDataset; i++) {
                            Integer newId = sampleToNewSampleId[i];
                            if (newId != null) {
                                outputallele1[newId] = allele1[i];
                                outputallele2[newId] = allele2[i];
                            }
                        }

                        newDS.setAlleles(s, outputallele1, outputallele2);

                        snpObj.clearGenotypes();

                        pb.iterate();
                    }
                }
                loader.close();
            }



        }
        
        pb.close();
        
        for (int d = 0; d < ds.length; d++) {
            System.out.println("Percentage of SNPs form dataset "+d+": " + (dataFromSet[d]/(double)uniqueSNPs.size())*100 );
        }
        
        if (snpsProcessed.size() != uniqueSNPs.size()) {
            log.writeln("ERROR: nr of processed SNPs unequal to nr of unique SNPs. Found: " + snpsProcessed.size() + "\tExpected: " + uniqueSNPs.size());
            System.out.println("ERROR: nr of processed SNPs unequal to nr of unique SNPs. Found: " + snpsProcessed.size() + "\tExpected: " + uniqueSNPs.size());
        } else {
            log.writeln("Everything seems to be ok. Have a nice day.");
            System.out.println("Everything seems to be ok. Have a nice day.");
        }
        snpout.close();
        snpMappingsOut.close();
        
        newDS.close();
        log.close();
    }

    public void checkMerge(String[] datasetLocations, String outputDir) throws IOException {

        System.out.println("Checking MERGE");
        HashSet<String> uniquesnps = new HashSet<String>();
        TextFile snpsIn = new TextFile(outputDir + "UniqueSNPsAccrossDatasets.txt", TextFile.R);
        uniquesnps.addAll(snpsIn.readAsArrayList());
        snpsIn.close();

        TextFile log = new TextFile(outputDir + "mergecheck.txt", TextFile.W);
        log.writeln("Checking: " + outputDir);
        System.out.println("Checking: " + outputDir);

        TriTyperGenotypeData output = new TriTyperGenotypeData();
        output.load(outputDir);
        SNPLoader outputLoader = output.createSNPLoader();
        String[] indsInOutput = output.getIndividuals();
        for (int d = 0; d < datasetLocations.length; d++) {
            log.writeln("Now parsing: " + datasetLocations[d]);
            System.out.println("Now parsing: " + datasetLocations[d]);
            TriTyperGenotypeData input = new TriTyperGenotypeData();
            input.load(datasetLocations[d]);
            Integer[] fromInputToOutput = new Integer[indsInOutput.length];
            for (int ind = 0; ind < indsInOutput.length; ind++) {
                fromInputToOutput[ind] = input.getIndividualId(indsInOutput[ind]);
            }

            // now check the merger
            SNPLoader inputLoader = input.createSNPLoader();
            String[] snpsInInput = input.getSNPs();
            Integer[] snpMap = new Integer[snpsInInput.length];
            for (int s = 0; s < snpMap.length; s++) {
                String snpName = snpsInInput[s];
                if (uniquesnps.contains(snpName)) {
                    Integer snpIdInOutput = output.getSnpToSNPId().get(snpsInInput[s]);
                    if (snpIdInOutput != -9) {
                        SNP snpObjInOutput = output.getSNPObject(snpIdInOutput);
                        SNP snpObjInInput = input.getSNPObject(s);
                        SNP[] snpObjs = new SNP[]{snpObjInInput, snpObjInOutput};

                        if (snpObjInInput == null || snpObjInOutput == null) {
                            System.out.println("WARNING: snp present in input but not in output:\t" + snpMap[s] + "\t" + s + "\t" + snpIdInOutput);
                        } else {
                            inputLoader.loadGenotypes(snpObjInInput);
                            outputLoader.loadGenotypes(snpObjInOutput);
                            Boolean[] flipAlleles = CompareAllelicDirections.compare(snpObjs);

                            if (flipAlleles == null) {
                                String outStr = snpsInInput[s];
                                for (int ds = 0; ds < snpObjs.length; ds++) {
                                    String dsName = datasetLocations[d];

                                    if (ds == 1) {
                                        dsName = "output";
                                    }

                                    if (snpObjs[ds] != null) {

                                        SNP snpObj = snpObjs[ds];
                                        outStr += "\t" + dsName + "\t"
                                                + BaseAnnot.toString(snpObj.getAlleles()[0]) + "/" + BaseAnnot.toString(snpObj.getAlleles()[1])
                                                + "\t" + BaseAnnot.toString(snpObj.getMinorAllele());
                                    } else {
                                        outStr += "\tNotPresentIn" + dsName;
                                    }

                                }
//                            log.writeln(outStr);
                            } else {



                                // compare alleles
                                short[] genotypes1 = new short[fromInputToOutput.length];
                                short[] genotypes2 = new short[fromInputToOutput.length];
                                Integer[] nrDifferentPerDataset = new Integer[2];
                                int d1 = 0;
                                int d2 = 1;
                                SNP snp1 = snpObjs[d1];

                                byte[] gt1 = snp1.getGenotypes();
                                for (int i = 0; i < fromInputToOutput.length; i++) {
                                    short gt = gt1[fromInputToOutput[i]];
                                    if (flipAlleles[d1]) {
                                        if (gt == 0) {
                                            gt = 2;
                                        } else if (gt == 2) {
                                            gt = 0;
                                        }
                                    }
                                    genotypes1[i] = gt;
                                }

                                SNP snp2 = snpObjs[d2];

                                byte[] gt2 = snp2.getGenotypes();
                                for (int i = 0; i < fromInputToOutput.length; i++) {
                                    short gt = gt2[i];
                                    if (flipAlleles[d2]) {
                                        if (gt == 0) {
                                            gt = 2;
                                        } else if (gt == 2) {
                                            gt = 0;
                                        }
                                    }
                                    genotypes2[i] = gt;
                                }

                                // everything initialized. Now continue to comparison
                                int identical = 0;
                                int different = 0;
                                int called = 0;
                                for (int g = 0; g < genotypes1.length; g++) {
                                    short g1 = genotypes1[g];
                                    short g2 = genotypes2[g];
                                    if (g1 < 0 || g2 < 0) {
                                        // one or both is uncalled
                                    } else {
                                        called++;
                                        if (g1 == g2) {
                                            identical++;
                                        } else {
                                            different++;
                                        }
                                    }
                                }
                                String snpStats1 = snp1.getName() + "\t" + snp1.getCR() + "\t" + snp1.getMAF() + "\t" + snp1.getHWEP() + "\t" + BaseAnnot.toString(snp1.getAlleles()[0]) + "/" + BaseAnnot.toString(snp1.getAlleles()[1]) + "\t" + BaseAnnot.toString(snp1.getMinorAllele());
                                String snpStats2 = snp2.getName() + "\t" + snp2.getCR() + "\t" + snp2.getMAF() + "\t" + snp2.getHWEP() + "\t" + BaseAnnot.toString(snp2.getAlleles()[0]) + "/" + BaseAnnot.toString(snp2.getAlleles()[1]) + "\t" + BaseAnnot.toString(snp2.getMinorAllele());
                                log.writeln(snpMap[s] + "\t" + datasetLocations[d] + "\t" + snpStats1 + "\t" + flipAlleles[d1] + "\toutput\t" + snpStats2 + "\t" + flipAlleles[d2] + "\t" + called + "\t" + different + "\t" + identical);
//                            System.out.println(snpMap[s] + "\t" + datasetLocations[d] + "\t" + snpStats1 + "\t" + flipAlleles[d1] + "\toutput\t" + snpStats2 + "\t" + flipAlleles[d2] + "\t" + called + "\t" + different + "\t" + identical);
                                nrDifferentPerDataset[d1] = different;





                            } // end if/else flipalleles == null
                            for (int q = 0; q < snpObjs.length; q++) {
                                snpObjs[q].clearGenotypes();
                            }
                        }


                    } // else snpIdInOutput == null
                }



            } // end for int s = 0; s < snpMap.length; s++
            inputLoader.close();
            input = null;

        }
        
        
        outputLoader.close();
        output = null;
        log.close();
    }
}
