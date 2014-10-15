/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputationtool.postprocessing;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Pattern;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class TriTyperDatasetCorrelator {

    int nrThreads = 1;
    private TriTyperGenotypeData ggDataset1;
    private TriTyperGenotypeData ggDataset2;
    private TextFile log;
    private HashSet<String> confineSNPList;
    private HashMap<String, ArrayList<Double>> beagleR2;
    private int[] beaglecorrelationfreqdistribution;

    public TriTyperDatasetCorrelator(String dataset1, String dataset1Name, String dataset2, String dataset2Name) throws IOException {
        nrThreads = 2;
        System.out.println("Running with " + nrThreads + " threads.");

        ggDataset1 = new TriTyperGenotypeData();
        ggDataset1.load(dataset1);

        ggDataset2 = new TriTyperGenotypeData();
        ggDataset2.load(dataset2);

    }

    public TriTyperDatasetCorrelator(String dataset1, String dataset1Name, String dataset2, String dataset2Name, String beagleInput, String template, Integer numBatches) throws IOException {
        nrThreads = 2;
        System.out.println("Running with " + nrThreads + " threads.");

        loadBeagleR2(beagleInput, template, numBatches);

        ggDataset1 = new TriTyperGenotypeData();
        ggDataset1.load(dataset1);

        ggDataset2 = new TriTyperGenotypeData();
        ggDataset2.load(dataset2);

    }

    public double[] getGenotypes(SNP snp, TriTyperGenotypeData gg, int[] inds) {
        double[] tmpGenotypes = new double[inds.length];

        double[] dosage = snp.getDosageValues();
        byte[] genotypes = snp.getGenotypes();

        if (dosage == null) {

            int g = 0;
            for (int i = 0; i < inds.length; i++) {
                int indId = inds[i];
                tmpGenotypes[g] = (double) genotypes[indId];
                g++;
            }
//	    System.out.println("");
        } else {
            int g = 0;
            for (int i = 0; i < inds.length; i++) {
                int indId = inds[i];
                tmpGenotypes[g] = dosage[indId];

                // System.out.print(tmpGenotypes[g] + " ");
                g++;
            }
//	    System.out.println("");
        }

//	for(int g=0; g<tmpGenotypes.length; g++ ){
//	    System.out.print(tmpGenotypes[g]+" ");
//	}
        // System.out.println("");

        return tmpGenotypes;
    }

    public void run(String outputLocation) throws IOException {

        //Complement allele information:
        byte[] complementAllele = new byte[256];
        complementAllele[0] = 0; // Nothing
        complementAllele[84] = 65; //T > A
        complementAllele[65] = 84; //A > T
        complementAllele[48] = 48; //- > -
        complementAllele[67] = 71; //C > G
        complementAllele[71] = 67; //G > C
        int[] alleleIndex = new int[256];
        for (int a = 0; a < 256; a++) {
            alleleIndex[a] = 4;
        }
        alleleIndex[65] = 0;
        alleleIndex[67] = 1;
        alleleIndex[71] = 2;
        alleleIndex[84] = 3;

        log = new TextFile(outputLocation + "/correlationOutput.txt", TextFile.W);

        determineUniqueSNPS();

        int numInd1 = ggDataset1.getIndividuals().length;
        int numInd2 = ggDataset2.getIndividuals().length;

        int[] inds1 = new int[numInd1];
        int[] inds2 = new int[numInd1];

        HashMap<Integer, Integer> ind1ToInd2 = new HashMap<Integer, Integer>();

        // determine shared individuals
        int numSharedIndividuals = 0;
        for (int i = 0; i < numInd1; i++) {
            String indName1 = ggDataset1.getIndividuals()[i];
            if (ggDataset1.getIsIncluded()[i]) {
                Integer indId2 = ggDataset2.getIndividualId(indName1);
                if (indId2 != -9) {
                    if (ggDataset2.getIsIncluded()[indId2]) {


                        inds1[i] = i;
                        inds2[i] = indId2;

                        numSharedIndividuals++;
                    } else {
                        inds1[i] = -1;
                        inds2[i] = -1;
                    }
                } else {
                    inds1[i] = -1;
                    inds2[i] = -1;
                }
            } else {
                inds1[i] = -1;
                inds2[i] = -1;
            }
        }

        System.out.println("Shared samples: " + numSharedIndividuals);
        if (numSharedIndividuals == 0) {
            System.exit(-1);
        }

        int[] inds1Final = new int[numSharedIndividuals];
        int[] inds2Final = new int[numSharedIndividuals];

        int counter = 0;
        for (int i = 0; i < numInd1; i++) {
            if (inds1[i] != -1 && inds2[i] != -1) {
                System.out.println(counter + "\t" + inds1[i] + "\t" + inds2[i]);
                inds1Final[counter] = inds1[i];
                inds2Final[counter] = inds2[i];
                counter++;
            }
        }

        counter = 0;
        int below80 = 0;
        int belowSquared80 = 0;
        int flipped = 0;

        String[] snps = ggDataset1.getSNPs();

        int[] correlationfreqdistribution = new int[11];

        int[] correlationdifferencedistribution = new int[11];

        int prevInt = 0;
        int corNo = 0;
        ScatterPlot s = null;
        if (beagleR2 != null) {

            // create image
            s = new ScatterPlot(1000);



        }

        System.out.println("Mapping SNPs");
        SNPLoader loader1 = ggDataset1.createSNPLoader();
        SNPLoader loader2 = ggDataset2.createSNPLoader();

        for (int snp1 = 0; snp1 < snps.length; snp1++) {
            Integer snp2 = ggDataset2.getSnpToSNPId().get(snps[snp1]);
            if (snp2 != -9) {
                boolean takeComplement = false;
                String excludeReason = "";
                boolean exclude = false;

                SNP snp1Object = ggDataset1.getSNPObject(snp1);
                SNP snp2Object = ggDataset2.getSNPObject(snp2);

                loader1.loadGenotypes(snp1Object);
                loader2.loadGenotypes(snp2Object);

                double[] genotypes1 = getGenotypes(snp1Object, ggDataset1, inds1Final);
                double[] genotypes2 = getGenotypes(snp2Object, ggDataset2, inds2Final);

                // search for missing genotypes
                int missingGenotypes = 0;
                for (int g = 0; g < genotypes1.length; g++) {
                    // System.out.print(genotypes1[g]+"\t"+genotypes2[g]+"\n");
                    if (genotypes1[g] == -1 || genotypes2[g] == -1) {
                        genotypes1[g] = -1;
                        genotypes2[g] = -1;
                        missingGenotypes++;
                    }
                }
                // System.out.println("");

                if ((double) missingGenotypes / genotypes1.length > 0.1) {
                    exclude = true;
                    excludeReason = "SNP has low callrate (> 10%): Missing: " + missingGenotypes + " / " + genotypes1.length;
                }

                double maf1 = snp1Object.getMAF();
                double maf2 = snp2Object.getMAF();

                if (maf1 < 0.05 || maf2 < 0.05) {
                    exclude = true;
                    excludeReason += "\tMAF < 0.05: " + maf1 + "\t" + maf2;
                }

                //Check whether the physical mapping of the two SNPs is identical:
                if (ggDataset1.getChr(snp1) != null) {
                    int chr1 = ggDataset1.getChr(snp1);
                    int pos1 = ggDataset1.getChrPos(snp1);

                    int chr2 = ggDataset2.getChr(snp2);
                    int pos2 = ggDataset2.getChrPos(snp2);
                    if (chr1 != chr2 || pos1 != pos2) {
                        // do stuff
                        exclude = true;
                        excludeReason += "\tSNPs map to different positions";
                    }
                }

                //SNP has been typed both in dataset and dataset2, do we need to take complementary alleles?
                byte[] allelesbytes = snp1Object.getAlleles();
                String alleles2 = new String(snp2Object.getAlleles());
                String alleles1 = null;
                try {
                    alleles1 = new String(snp1Object.getAlleles(), "UTF-8");
                } catch (Exception e) {
                }

                boolean allelesOk = true;

                if (allelesbytes[0] == 00 || allelesbytes[1] == 00) {
                    exclude = true;
                }


                if (alleles1 == null) {
                    exclude = true;
                    excludeReason += " SNPs has null alleles";
                }


                // hashSNPAllelesHapMap.put(rsName, allelesHapMap);

                boolean strandForward = true;
                int[] alleleIndex1 = new int[5];

                for (int ind = 0; ind < ggDataset2.getIndividuals().length; ind++) {
                    if (ggDataset2.getIsIncluded()[ind]) {
                        byte[] snpallele1 = snp2Object.getAllele1();
                        byte[] snpallele2 = snp2Object.getAllele2();
                        byte allele1Byte = snpallele1[ind];
                        alleleIndex1[alleleIndex[allele1Byte]]++;
                        byte allele2Byte = snpallele2[ind];
                        alleleIndex1[alleleIndex[allele2Byte]]++;
                    }
                }
                int[] alleleIndex2 = new int[5];
                for (int ind = 0; ind < ggDataset1.getIndividuals().length; ind++) {
                    if (ggDataset1.getIsIncluded()[ind]) {
                        byte[] snpallele1 = snp1Object.getAllele1();
                        byte[] snpallele2 = snp1Object.getAllele2();

                        byte allele1Byte = snpallele1[ind];
                        alleleIndex2[alleleIndex[allele1Byte]]++;
                        byte allele2Byte = snpallele2[ind];
                        alleleIndex2[alleleIndex[allele2Byte]]++;
                    }
                }

                double[] alleleIndexFreq1 = new double[4];
                double[] alleleIndexFreq2 = new double[4];
                int itr = 0;

                boolean issueResolved = false;

                while (!issueResolved) {
                    //Take complement alleles when necessary:
                    if (!strandForward) {
                        int[] alleleIndex1Copy = new int[4];
                        System.arraycopy(alleleIndex1, 0, alleleIndex1Copy, 0, 4);
                        alleleIndex1[0] = alleleIndex1Copy[3];
                        alleleIndex1[1] = alleleIndex1Copy[2];
                        alleleIndex1[2] = alleleIndex1Copy[1];
                        alleleIndex1[3] = alleleIndex1Copy[0];
                    }
                    //Determine total number of called alleles:
                    int totalCalled1 = 0;
                    int totalCalled2 = 0;
                    for (int a = 0; a < 4; a++) {
                        totalCalled1 += alleleIndex1[a];
                        totalCalled2 += alleleIndex2[a];
                    }
                    //Calculate allele freq:
                    for (int a = 0; a < 4; a++) {
                        alleleIndexFreq1[a] = (double) alleleIndex1[a] / (double) totalCalled1;
                        alleleIndexFreq2[a] = (double) alleleIndex2[a] / (double) totalCalled2;
                    }
                    //Check whether alleles are identical:
                    int nrDifferentAllelesPresent = 0;
                    for (int a = 0; a < 4; a++) {
                        if (alleleIndexFreq1[a] > 0 || alleleIndexFreq2[a] > 0) {
                            nrDifferentAllelesPresent++;
                        }
                    }
                    if (nrDifferentAllelesPresent > 2) {
                        strandForward = !strandForward;
                    } else {
                        //break;
                        issueResolved = true;
                    }
                    itr++;

                    if (itr >= 2) {
                        if (!issueResolved) {
                            //Taking complementary alleles does not resolve anything:
                            exclude = true;
                            issueResolved = true;
                            excludeReason += "\tIncompatibleAlleles:Dataset=" + alleles2 + ",HapMap=" + alleles1;
                        }
                    }
                }

                takeComplement = !strandForward;

                //Check whether allele freq is comparable:
                boolean concordant = true;
                for (int a = 0; a < 4; a++) {
                    if (alleleIndexFreq1[a] > 0 && alleleIndexFreq2[a] > 0) {
                        if (alleleIndexFreq1[a] > 0.5 && alleleIndexFreq2[a] < 0.5) {
                            concordant = false;
                        }
                        if (alleleIndexFreq1[a] < 0.5 && alleleIndexFreq2[a] > 0.5) {
                            concordant = false;
                        }
                    }
                }

                //If SNP is AT or CG SNP, it can be we had to take complementary allele:
                byte[] snpAlleles = snp2Object.getAlleles();
                if (snpAlleles[0] + snpAlleles[1] == 65 + 84 || snpAlleles[0] + snpAlleles[1] == 67 + 71) {
                    if (!concordant) {
                        takeComplement = !takeComplement;
                        concordant = true;
                    }
                }



                if (exclude && takeComplement) {
                    System.out.println(snp1Object.getName() + "\t" + excludeReason);
                } else if (!exclude) {
                    double[] finalGenotypes1 = null;
                    double[] finalGenotypes2 = null;
                    if (missingGenotypes > 0) {

                        finalGenotypes1 = new double[genotypes1.length - missingGenotypes];
                        finalGenotypes2 = new double[genotypes1.length - missingGenotypes];
                        int ctr = 0;
                        for (int i = 0; i < genotypes1.length; i++) {
                            if (genotypes1[i] != -1) {
                                finalGenotypes1[ctr] = genotypes1[i];

                                finalGenotypes2[ctr] = genotypes2[i];

                                ctr++;
                            }
                        }
                    } else {
                        finalGenotypes1 = genotypes1;
                        finalGenotypes2 = genotypes2;
                    }


                    // System.out.println(snps[snp1]);
                    double correlation = JSci.maths.ArrayMath.correlation(finalGenotypes1, finalGenotypes2);



                    // System.out.println(correlation);

                    // System.exit(0);
                    double absCor = Math.abs(correlation);
                    double absCorSquared = absCor * absCor;

                    if (beagleR2 != null) {
                        ArrayList<Double> r2s = beagleR2.get(snp1Object.getName());
                        double m = 0;
                        for (int i = 0; i < r2s.size(); i++) {
                            m += r2s.get(i);
                        }
                        m /= r2s.size();
                        s.plot(absCorSquared, m);

                    }

                    int binNumber = (int) (absCorSquared * 10d);
                    correlationfreqdistribution[binNumber]++;
                    if (absCor <= 0.80) {
                        below80++;
                    }

                    if (absCorSquared <= 0.80) {
                        belowSquared80++;
                    }

                    if (correlation < 0) {
                        flipped++;
                    }

                    snp1Object = null;
                    snp2Object = null;

                    counter++;

                    //corrs[corNo] = correlation * correlation;
                    corNo++;

                    if (counter % 10000 == 0 && counter > prevInt) {
                        System.out.println(counter + " SNPS processed\t" + below80 + "\t" + (100 * (double) below80 / counter) + "% R <= 0.80.\ts" + (100 * (double) belowSquared80 / counter) + "% R2 <= 0.80.\t" + flipped + " - " + (100 * (double) flipped / counter) + " flipped");
                        prevInt = counter;
                    }
                }
                snp1Object.clearGenotypes();
                snp2Object.clearGenotypes();
            }
        }

        if (beagleR2 != null) {
            s.draw(outputLocation + "/CorrelationVsBeagleR2.png");
        }

        System.out.println(counter + " SNPS processed\t" + below80 + "\t" + (100 * (double) below80 / counter) + "% R <= 0.80.\ts" + (100 * (double) belowSquared80 / counter) + "% R2 <= 0.80.\t" + flipped + " - " + (100 * (double) flipped / counter) + " flipped");
        // log.log("Mean r-squared: "+JSci.maths.ArrayMath.mean(corrs)+" variance: "+JSci.maths.ArrayMath.variance(corrs));
        prevInt = counter;
        log.writeln("Distribution of correlations (counts, R2):");
        int total = 0;
        for (int u = 0; u < correlationfreqdistribution.length; u++) {

            total += correlationfreqdistribution[u];
            log.writeln(u + "\t" + correlationfreqdistribution[u] + "\t" + total);
        }

        double sum = 0d;
        log.writeln("Distribution of correlations (frequency, R2):");
        for (int u = 0; u < correlationfreqdistribution.length; u++) {
            double freq = (double) correlationfreqdistribution[u] / total;
            sum += freq;
            log.writeln(u + "\t" + freq + "\t" + sum);
        }

        if (beaglecorrelationfreqdistribution != null) {
            log.writeln("Beagle Distribution of correlations (R2):");
            total = 0;
            for (int u = 0; u < beaglecorrelationfreqdistribution.length; u++) {
                total += beaglecorrelationfreqdistribution[u];
            }

            sum = 0d;
            for (int u = 0; u < beaglecorrelationfreqdistribution.length; u++) {
                double freq = (double) beaglecorrelationfreqdistribution[u] / total;
                sum += freq;
                log.writeln(u + "\t" + freq + "\t" + sum);
            }
        }

        log.close();
    }

    private void loadBeagleR2(String inputDir, String template, Integer numBatches) throws IOException {

        String[] batchNames = getBatches(numBatches);

        beagleR2 = new HashMap<String, ArrayList<Double>>();
        boolean allFilesAvailable = true;
        for (int b = 0; b < numBatches; b++) {

            String batchName = batchNames[b];

            for (int chr = 1; chr <= 22; chr++) {
                // for(int chr: chromosomes){
                String templatecopy = new String(template);
                templatecopy = templatecopy.replace("BATCH", batchName);
                templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

                // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
                String fileName = inputDir + "/" + templatecopy + ".r2";
                if (!Gpio.canRead(fileName)) {
                    System.out.println("Cannot open file:\t" + fileName);
                    allFilesAvailable = false;
                }
            }
        }

        if (allFilesAvailable) {

            beaglecorrelationfreqdistribution = new int[11];
            for (int b = 0; b < numBatches; b++) {

                String batchName = batchNames[b];

                for (int chr = 1; chr <= 22; chr++) {

                    String templatecopy = new String(template);
                    templatecopy = templatecopy.replace("BATCH", batchName);
                    templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);
                    String fileName = inputDir + "/" + templatecopy + ".r2";

                    TextFile in = new TextFile(fileName, TextFile.R);
                    String line = "";
                    Pattern tab = Pattern.compile("\t");
                    while ((line = in.readLine()) != null) {
                        String[] elems = tab.split(line);
                        Double r2 = null;
                        if (elems.length == 2) {
                            try {
                                r2 = Double.parseDouble(elems[1]);
                                if (r2.isNaN()) {
                                    r2 = 0d;

                                }
                                int binNumber = (int) (r2 * 10d);
                                beaglecorrelationfreqdistribution[binNumber]++;

                                ArrayList<Double> r2s = beagleR2.get(elems[0]);
                                if (r2s == null) {
                                    r2s = new ArrayList<Double>();
                                }
                                r2s.add(r2);
                                beagleR2.put(elems[0], r2s);
                            } catch (NumberFormatException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                    in.close();

                }
            }
        }
    }

    private void getRandomSNPs(int num, String output) throws IOException {
        String[] snps1 = ggDataset1.getSNPs();
        String[] snps2 = ggDataset2.getSNPs();

        int counter = 0;

        ArrayList<String> snps = new ArrayList<String>();

        for (int i = 0; i < snps1.length; i++) {
            Integer snp2Id = ggDataset2.getSnpToSNPId().get(snps1[i]);
            if (snp2Id != -9) {
                snps.add(snps1[i]);
            }
        }


        TextFile log2 = new TextFile(output, TextFile.W);
        for (int i = 0; i < num; i++) {
            log2.writeln(snps.remove((int) Math.round(Math.random() * snps.size())));
        }
        log2.close();

    }

    private void determineUniqueSNPS() throws IOException {
        String[] snps1 = ggDataset1.getSNPs();
        String[] snps2 = ggDataset2.getSNPs();

        int counter = 0;

        if (confineSNPList != null) {
            int absent1 = 0;
            int absent2 = 0;
            Iterator<String> i = confineSNPList.iterator();
            while (i.hasNext()) {
                if (ggDataset1.getSnpToSNPId().get(i.next()) == -9) {
                    absent1++;
                }
                if (ggDataset2.getSnpToSNPId().get(i.next()) == -9) {
                    absent2++;
                }
            }

            log.writeln(absent1 + " of the " + confineSNPList.size() + " SNPs in ConfineSNPList are not present in " + ggDataset1.getGenotypeFileName() + ", and " + absent2 + " are not present in " + ggDataset2.getGenotypeFileName());

        } else {

            for (int i = 0; i < snps1.length; i++) {
                Integer snp2Id = ggDataset2.getSnpToSNPId().get(snps1[i]);
                if (snp2Id == -9) {
                    counter++;
                }
            }
            log.writeln(counter + " of the " + snps1.length + " SNPs in " + ggDataset1.getGenotypeFileName() + " are not present in " + ggDataset2.getGenotypeFileName() + ": " + (((double) counter / snps1.length) * 100) + "%");

            counter = 0;
            for (int i = 0; i < snps2.length; i++) {
                Integer snp1Id = ggDataset1.getSnpToSNPId().get(snps2[i]);
                if (snp1Id == -9) {
                    counter++;
                }

            }

            log.writeln(counter + " of the " + snps2.length + " SNPs in " + ggDataset2.getGenotypeFileName() + " are not present in " + ggDataset1.getGenotypeFileName() + ": " + (((double) counter / snps2.length) * 100) + "%");
        }


    }

    private String[] getBatches(int numBatches) {

        String[] batches = new String[numBatches];
        String[] alphabet = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"};

        String firstletter = "a";
        int alphacounter = 1;
        int betacounter = 0;
        for (int i = 0; i < numBatches; i++) {
            if (i == 25) {
                firstletter = alphabet[alphacounter];
                alphacounter++;
                betacounter = 0;
            }

            batches[i] = firstletter + alphabet[betacounter];
            betacounter++;
        }
        return batches;
    }

    public void confineToSNPs(String snpList) throws IOException {
        confineSNPList = new HashSet<String>();

        TextFile tf = new TextFile(snpList, TextFile.R);
        confineSNPList.addAll(tf.readAsArrayList());
        tf.close();
    }

    private HashMap<Integer, Integer> determineSampleIDMap() {
        String[] ds1Samples = ggDataset1.getIndividuals();
        String[] ds2Samples = ggDataset2.getIndividuals();
        HashMap<Integer, Integer> samplemap = new HashMap<Integer, Integer>();
        int numShared = 0;
        for (int i = 0; i < ds1Samples.length; i++) {
            String sample1 = ds1Samples[i];
            for (int j = 0; j < ds2Samples.length; j++) {
                String sample2 = ds2Samples[j];
                if (sample1.equals(sample2)) {
                    samplemap.put(i, j);
                    numShared++;
                }
            }
        }

        System.out.println("Shared samples between datasets:" + numShared);

        if (numShared > 0) {
            return samplemap;
        }
        return null;
    }
}
