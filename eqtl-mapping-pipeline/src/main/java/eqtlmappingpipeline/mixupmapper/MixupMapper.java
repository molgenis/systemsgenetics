/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.mixupmapper;

import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;
import umcg.genetica.util.RankArray;

/**
 *
 * @author harmjan
 */
public class MixupMapper extends MetaQTL3 {

    private TriTyperGenotypeData genotypeData;
    private DoubleMatrixDataset<String, String> traitData;
    private HashMap<String, String> genotypeToTrait;
    private HashMap<String, String> traitToGenotype;

    public void run(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
            String ingt, String inexp, String inexpplatform, String inexpannot,
            String gte, String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads, Integer maxNrResults,
            String regressouteqtls, String snpprobecombofile, String inputeQTLs, boolean testall) throws IOException, Exception {


        String initialOutputdir = out;
        if (inputeQTLs == null || !Gpio.exists(inputeQTLs)) {
            if (out != null) {
                // check outdir for input
                initialOutputdir = Gpio.formatAsDirectory(initialOutputdir);
                out = initialOutputdir + "Cis-eQTLs/";
                System.out.println("Looking for file: " + out + "eQTLProbesFDR0.05.txt");
                if (Gpio.exists(out + "eQTLProbesFDR0.05.txt")) {
                    inputeQTLs = out + "eQTLProbesFDR0.05.txt";
                }
            }
        }

        initialOutputdir = Gpio.formatAsDirectory(initialOutputdir);
        if (inputeQTLs == null) {
            System.out.println("Could not find eQTL file. Will therefore perform eQTL mapping first.");

            out = initialOutputdir + "Cis-eQTLs/";
            initialize(xmlSettingsFile, texttoreplace, texttoreplacewith, null, null, ingt, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, maxNrResults, regressouteqtls, snpprobecombofile, true, true, null, 0.05d, null);
            mapEQTLs();
            inputeQTLs = m_settings.outputReportsDir + "eQTLProbesFDR0.05.txt";

            // clear memory
            m_gg = null;
        }
        
        if(!Gpio.exists(inputeQTLs)){
            System.err.println("Something went wrong during eQTL mapping: most probably the FDR calculations failed.\nCheck the files in: "+m_settings.outputReportsDir);
            System.exit(-1);
        }

        System.out.println("Using: " + inputeQTLs + " as input for MixupMapper");

        loadGenotypeData(ingt);
        loadGeneExpressionData(inexp);
        linkSamples(gte);
        out = initialOutputdir + "MixupMapper/";
        Gpio.createDir(out);
        mapMixups(inputeQTLs, testall, out, false);

    }

    private void loadGenotypeData(String ingt) throws IOException {
        genotypeData = new TriTyperGenotypeData(ingt);
    }

    private void loadGeneExpressionData(String inexp) throws IOException {


        TextFile tf = new TextFile(inexp, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        tf.close();
        if (elems[0].trim().toLowerCase().equals("probe") && elems[3].trim().toLowerCase().equals("chr") && elems[4].trim().toLowerCase().equals("chrstart")) {
            System.out.println("File is in old TriTyper format..");
            traitData = new DoubleMatrixDataset<String, String>(inexp, null, null, 9);
        } else {
            traitData = new DoubleMatrixDataset<String, String>(inexp);
        }



        // rank the data...
        traitData.rawData = rankAllExpressionData(traitData.rawData);
    }

    /**
     * Ranks all the expressiondata available in rawData
     */
    private double[][] rankAllExpressionData(double[][] matrix) {

        RankArray r = new RankArray();
        boolean rankWithTies = false;


        for (int p = 0; p < matrix.length; ++p) {
            double[] probeData = matrix[p];
            probeData = r.rank(probeData, rankWithTies);
            matrix[p] = probeData;
        }
        return matrix;
    }

    private void linkSamples(String gte) throws IOException {
        genotypeToTrait = new HashMap<String, String>();
        traitToGenotype = new HashMap<String, String>();

        String[] genotypeSamples = genotypeData.getIndividuals();
        HashSet<String> gtInds = new HashSet<String>();
        HashSet<String> tInds = new HashSet<String>();
        gtInds.addAll(Arrays.asList(genotypeSamples));
        tInds.addAll(traitData.colObjects);

        if (gte != null && Gpio.exists(gte)) {
            System.out.println("Loading genotype to trait link file: " + gte);


            TextFile tf = new TextFile(gte, TextFile.R);
            String[] gtelems = tf.readLineElems(Strings.whitespace);
            while (gtelems != null) {
                if (gtelems.length >= 2) {
                    String g = gtelems[0];

                    Integer gId = genotypeData.getIndividualId(g);
                    if (gId != -9 && genotypeData.getIsIncluded()[gId]) {
                        String t = gtelems[1];

                        if (tInds.contains(t)) {
                            if (genotypeToTrait.containsKey(g)) {
                                System.out.println("WARNING: " + g + "\tgenotype sample has already been linked to trait sample\t" + genotypeToTrait.get(g));
                            } else {
                                genotypeToTrait.put(g, t);
                            }
                            if (traitToGenotype.containsKey(t)) {
                                System.out.println("WARNING: " + t + "\ttrait sample has already been linked to genotype sample\t" + traitToGenotype.get(t));
                            } else {
                                traitToGenotype.put(t, g);
                            }
                        }
                    }


                }
                gtelems = tf.readLineElems(Strings.whitespace);
            }
            tf.close();
            System.out.println(traitToGenotype.size() + " combinations of traits and genotypes samples loaded.");
            System.out.println(genotypeToTrait.size() + " combinations of genotypes and traits samples loaded.");

        } else if (gte != null && !Gpio.exists(gte)) {
            System.out.println("ERROR: genotype to trait link file: " + gte + " does not exist.");
        } else {
            for (String ind : genotypeSamples) {
                Integer gId = genotypeData.getIndividualId(ind);
                if (gId != -9 && genotypeData.getIsIncluded()[gId]) {
                    if (traitData.hashCols.get(ind) != null) {
                        if (genotypeToTrait.get(ind) == null) {
                            genotypeToTrait.put(ind, ind);
                            traitToGenotype.put(ind, ind);
                        } else {
                            System.out.println("WARNING: " + ind + "\tgenotype sample has already been linked to trait sample\t" + genotypeToTrait.get(ind) + ". Your dataset most probably contains duplicate identifiers.");
                        }

                    }
                }
            }

            System.out.println("Linking genotypes and trait data yielded: " + traitToGenotype.size() + " links");
        }

        if (traitToGenotype.isEmpty()) {
            System.out.println("ERROR: no links between genotype and trait samples found.");
            System.exit(0);
        }
    }

    private void mapMixups(String inputeQTLs, boolean testAll, String outDir, boolean leavehalveout) throws IOException {


        // load eQTLs, check whether SNPs and Probes are in dataset.
        ArrayList<Pair<String, String>> eQTLs = new ArrayList<Pair<String, String>>();
        TextFile tf = new TextFile(inputeQTLs, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length >= 4) {
                String snp = elems[1];
                String probe = elems[4];

                if (traitData.hashRows.containsKey(probe) && genotypeData.getSnpToSNPId().get(snp) != -9) {
                    Pair<String, String> p = new Pair<String, String>(snp, probe);
                    eQTLs.add(p);
                } else {
                    System.out.println("QTL not present in datasets: " + snp + " - " + probe);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        if (eQTLs.isEmpty()) {
            System.out.println("ERROR: no eQTLs detected.");
            System.exit(-1);
        }


        HashMap<String, Integer> genotypeToRowIndex = new HashMap<String, Integer>();
        HashMap<String, Integer> traitToColIndex = new HashMap<String, Integer>();

        String[] gtInds = genotypeData.getIndividuals();
        String[] trInds = traitData.colObjects.toArray(new String[0]);

        int nrGenotypes = 0;
        int nrTraits = 0;

        // create a diagonal for easy reference
        if (testAll) {
            // give all the genotype samples a new index (counting from 0)
            System.out.println("All vs all comparison");
            int ctr = 0;
            ArrayList<String> samplesToOrder = new ArrayList<String>();
            for (int i = 0; i < gtInds.length; i++) {
                String sampleName = gtInds[i];
                if (genotypeData.getIsIncluded()[i] && genotypeToTrait.get(sampleName) != null) {
                    genotypeToRowIndex.put(sampleName, ctr);
                    ctr++;
                } else if (genotypeData.getIsIncluded()[i]) {
                    samplesToOrder.add(sampleName);
                }
            }

            for (String s : samplesToOrder) {
                genotypeToRowIndex.put(s, ctr);
                ctr++;
            }
            nrGenotypes = ctr;

            // now give the linked trait samples the same index.
            ArrayList<Integer> samplesWithoutGenotypes = new ArrayList<Integer>();

            for (int i = 0; i < trInds.length; i++) {
                String traitSample = trInds[i];
                String genotypeSample = traitToGenotype.get(traitSample);
                if (genotypeSample == null) {
                    // we need a counter of some sort to give these samples an index as well...
                    samplesWithoutGenotypes.add(i);
                } else {
                    Integer genotypeSampleIndex = genotypeToRowIndex.get(genotypeSample);
                    traitToColIndex.put(trInds[i], genotypeSampleIndex);

//                    System.out.println(trInds[i] + "\t" + genotypeSampleIndex + "\t" + genotypeSample);
                    nrTraits++;
                }
            }

            System.out.println(samplesWithoutGenotypes.size() + " trait samples could not be linked to a genotype");

            for (Integer sample : samplesWithoutGenotypes) {
                traitToColIndex.put(trInds[sample], nrTraits);
                nrTraits++;
            }

        } else {
            // test only linked samples
            int ctr = 0;
            System.out.println("Only comparing samples that have been linked to each other");
            for (int i = 0; i < gtInds.length; i++) {
                if (genotypeData.getIsIncluded()[i]) {
                    String gtInd = gtInds[i];
                    if (genotypeToTrait.get(gtInd) != null) {
                        genotypeToRowIndex.put(gtInds[i], ctr);
                        traitToColIndex.put(genotypeToTrait.get(gtInd), ctr);
                        ctr++;
                    }
                }
            }
            nrTraits = ctr;
            nrGenotypes = ctr;
        }

        System.out.println(nrTraits + " trait samples will be tested");
        System.out.println(nrGenotypes + " genotype samples will be tested");

        ////////////////////////////////
        // perform the actual testing //
        ////////////////////////////////

        // initialize matrices
        double[][] comparisonMatrix = new double[nrGenotypes][nrTraits];
        double[][] comparisonMatrixNrTested = new double[nrGenotypes][nrTraits];
        SNPLoader loader = genotypeData.createSNPLoader();

        // test all eQTLs
        int numTested = 0;
        int numNotTested = 0;
        System.out.println("Using " + eQTLs.size() + " eQTLs.");
        for (Pair<String, String> eqtl : eQTLs) {
            String snp = eqtl.getLeft();
            String probe = eqtl.getRight();

            if (!leavehalveout || (leavehalveout && Math.random() <= 0.5)) {

                Integer probeId = traitData.hashRows.get(probe);
                Integer snpId = genotypeData.getSnpToSNPId().get(snp);

                if (probeId == null || snpId == -9) {
                    // there is no such eQTL in the dataset!
                    System.out.println("Null trait or SNP ID:" + snpId + " (" + snp + ")\t" + probeId + " (" + snp + ")");
                    numNotTested++;
                } else {

                    // determine genotype mean and SD for eQTL.
                    double sdAA = -1;
                    double sdAB = -1;
                    double sdBB = -1;
                    double meanAA = -1;
                    double meanAB = -1;
                    double meanBB = -1;

                    SNP loadedSNP = genotypeData.getSNPObject(snpId);
                    loader.loadGenotypes(loadedSNP);
                    int[] genotypes = new int[genotypeToRowIndex.size()];

                    int numAA = 0;
                    int numAB = 0;
                    int numBB = 0;

                    for (int i = 0; i < gtInds.length; i++) {
                        String genotypeSample = gtInds[i];
                        Integer genotypeSampleIndex = genotypeToRowIndex.get(genotypeSample);

                        if (genotypeSampleIndex != null) {
                            int gt = loadedSNP.getGenotypes()[i];
                            genotypes[genotypeSampleIndex] = gt;

//                            if (genotypeSampleIndex > 75) {
//                                System.out.println(loadedSNP.getName() + "\t" + genotypeSample + "\t" + genotypeSampleIndex + "\t" + gt);
//                            }

                            if (gt == 0) {
                                numAA++;
                            } else if (gt == 2) {
                                numBB++;
                            } else if (gt == 1) {
                                numAB++;
                            }
                        }
                    }

                    if (numAA >= 3 && numAB >= 3 && numBB >= 3) {
                        double[] aa = new double[numAA];
                        double[] ab = new double[numAB];
                        double[] bb = new double[numBB];

                        int aaCTR = 0;
                        int abCTR = 0;
                        int bbCTR = 0;

                        for (int exp = 0; exp < trInds.length; exp++) {
                            String traitSample = trInds[exp];
                            String linkedGenotype = traitToGenotype.get(traitSample);
                            Integer linkedGenotypeIndex = genotypeToRowIndex.get(linkedGenotype);
                            // use only linked samples to recreate the eQTL
                            if (linkedGenotype != null) {
                                double expValue = traitData.rawData[probeId][exp];
                                int gt = genotypes[linkedGenotypeIndex];
                                if (gt != -1) {
                                    if (gt == 0) {
                                        aa[aaCTR] = expValue;
                                        aaCTR++;
                                    } else if (gt == 2) {
                                        bb[bbCTR] = expValue;
                                        bbCTR++;
                                    } else if (gt == 1) {
                                        ab[abCTR] = expValue;
                                        abCTR++;
                                    }
                                }
                            } else {
//                                System.err.println("No linked sample for expression sample: " + traitSample + "\t" + linkedGenotype);
                            }

                        }

                        sdAA = JSci.maths.ArrayMath.standardDeviation(aa);
                        sdBB = JSci.maths.ArrayMath.standardDeviation(bb);
                        sdAB = JSci.maths.ArrayMath.standardDeviation(ab);

                        meanAA = JSci.maths.ArrayMath.mean(aa);
                        meanBB = JSci.maths.ArrayMath.mean(bb);
                        meanAB = JSci.maths.ArrayMath.mean(ab);

                        if (sdAA > 0 && sdAB > 0 && sdBB > 0) {
                            for (int exp = 0; exp < trInds.length; exp++) {

                                String traitSample = trInds[exp];
                                Integer traitIndex = traitToColIndex.get(traitSample);
                                for (int gen = 0; gen < gtInds.length; gen++) {

                                    String genotypeSample = gtInds[gen];
                                    Integer genotypeIndex = genotypeToRowIndex.get(genotypeSample);
                                    if (traitIndex != null && genotypeIndex != null) {

                                        double expression = traitData.rawData[probeId][exp];

                                        int gt = genotypes[genotypeIndex];
                                        if (gt != -1) {
                                            double z = 0;
                                            if (gt == 0) {
                                                z = Math.abs(expression - meanAA) / sdAA;
                                            } else if (gt == 1) {
                                                z = Math.abs(expression - meanAB) / sdAB;
                                            } else {
                                                z = Math.abs(expression - meanBB) / sdBB;
                                            }

//                                        System.out.println(genotypeIndex + " - " + comparisonMatrixNrTested.length + "\t" + traitIndex + " - " + comparisonMatrixNrTested[0].length);
                                            if (!Double.isNaN(z) && z != 0) {
                                                comparisonMatrixNrTested[genotypeIndex][traitIndex]++;
                                                comparisonMatrix[genotypeIndex][traitIndex] += z;
                                            }
                                        }
                                    }
                                }
                            }
                            numTested++;
                        } else {
//                            System.out.println("Standard deviation is zero for one of the genotype groups: AA: " + sdAA + "\tAB: " + sdAB + "\tBB: " + sdBB);
                            numNotTested++;
                        }
                    } else {
                        numNotTested++;
                        // System.out.println("Minor allele frequency too low:\t" + snp + "\t" + probe + "\t" + numAA + "\t" + numAB + "\t" + numBB);
                    }
                    loadedSNP.clearGenotypes();
                }
            }

        }
        loader.close();

        System.out.println("Number QTLs tested: " + numTested + "");
        System.out.println("Number QTLs not tested: " + numNotTested + "");

        if (numTested == 0) {
            System.err.println("An error has occurred: none of the eQTLs was used during the MixupMapper test");
            System.exit(-1);
        }

        // scores have been calculated.. now visualize, and output...
        // scale result
        String[] gtRowNames = new String[nrGenotypes];
        String[] trColNames = new String[nrTraits];

        for (int exp = 0; exp < trInds.length; exp++) {
            String traitSample = trInds[exp];
            Integer traitIndex = traitToColIndex.get(traitSample);
            if (traitIndex != null) {
                trColNames[traitIndex] = traitSample;
            }
        }
        for (int gen = 0; gen < gtInds.length; gen++) {
            String genotypeSample = gtInds[gen];
            Integer gtIndex = genotypeToRowIndex.get(genotypeSample);
            if (gtIndex != null) {
                gtRowNames[gtIndex] = genotypeSample;
            }
        }

        for (int row = 0; row < comparisonMatrix.length; row++) {
            for (int col = 0; col < comparisonMatrix[row].length; col++) {
                comparisonMatrix[row][col] /= comparisonMatrixNrTested[row][col];
            }
        }

        // scale over columns
        for (int col = 0; col < nrTraits; col++) {
            double[] colData = new double[nrGenotypes];
            for (int row = 0; row < nrGenotypes; row++) {
                colData[row] = comparisonMatrix[row][col];
            }
            double mean = JSci.maths.ArrayMath.mean(colData);
            double sd = JSci.maths.ArrayMath.standardDeviation(colData);
            for (int row = 0; row < nrGenotypes; row++) {
                comparisonMatrix[row][col] = (comparisonMatrix[row][col] - mean) / sd;
            }
        }

        // scale over rows
        for (int row = 0; row < nrGenotypes; row++) {
            double[] rowData = new double[nrTraits];
            System.arraycopy(comparisonMatrix[row], 0, rowData, 0, nrTraits);

            double mean = JSci.maths.ArrayMath.mean(rowData);
            double sd = JSci.maths.ArrayMath.standardDeviation(rowData);
            for (int col = 0; col < nrTraits; col++) {
                comparisonMatrix[row][col] = (comparisonMatrix[row][col] - mean) / sd;
            }
        }

        // store..
        DoubleMatrixDataset<String, String> scoringMatrix = new DoubleMatrixDataset<String, String>(comparisonMatrix, Arrays.asList(gtRowNames), Arrays.asList(trColNames));
        scoringMatrix.save(outDir + "MixupMapperScores.txt");

        // match genotype samples
        TextFile matchedGTOut = new TextFile(outDir + "BestMatchPerGenotype.txt", TextFile.W);
        matchedGTOut.writeln("Genotype\tOriginalLinkedTrait\tOriginalLinkedTraitScore\tBestMatchingTrait\tBestMatchingTraitScore\tMixup");
        for (int row = 0; row < nrGenotypes; row++) {
            String gtSample = gtRowNames[row];
            double minValue = Double.MAX_VALUE;
            int minCol = -1;
            for (int col = 0; col < comparisonMatrix[row].length; col++) {
                double v = comparisonMatrix[row][col];
                if (v < minValue) {
                    minValue = v;
                    minCol = col;
                }
            }

            String linkedTrait = genotypeToTrait.get(gtSample);
            String linkedTraitScore = "N/A";
            if (linkedTrait == null) {
                linkedTrait = "N/A";
            } else {
                Integer linkedTraitID = traitToColIndex.get(linkedTrait);
                linkedTraitScore = "" + comparisonMatrix[row][linkedTraitID];
            }


            matchedGTOut.writeln(gtSample + "\t" + linkedTrait + "\t" + linkedTraitScore + "\t" + trColNames[minCol] + "\t" + comparisonMatrix[row][minCol] + "\t" + (!linkedTrait.equals("N/A") && !linkedTrait.equals(trColNames[minCol])));

        }
        matchedGTOut.close();

        // match trait samples
        TextFile matchedTrOut = new TextFile(outDir + "BestMatchPerTrait.txt", TextFile.W);
        matchedTrOut.writeln("Trait\tOriginalLinkedGenotype\tOriginalLinkedGenotypeScore\tBestMatchingGenotype\tBestMatchingGenotypeScore\tMixup");
        for (int col = 0; col < nrTraits; col++) {
            String trSample = trColNames[col];
            double minValue = Double.MAX_VALUE;
            int minRow = -1;
            for (int row = 0; row < comparisonMatrix.length; row++) {
                double v = comparisonMatrix[row][col];
                if (v < minValue) {
                    minValue = v;
                    minRow = row;
                }
            }

            String linkedGenotype = traitToGenotype.get(trSample);
            String linkedGenotypeScore = "N/A";
            if (linkedGenotype == null) {
                linkedGenotype = "N/A";
            } else {
                Integer linkedTraitID = genotypeToRowIndex.get(linkedGenotype);
                linkedGenotypeScore = "" + comparisonMatrix[linkedTraitID][col];
            }



            matchedTrOut.writeln(trSample + "\t" + linkedGenotype + "\t" + linkedGenotypeScore + "\t" + gtRowNames[minRow] + "\t" + comparisonMatrix[minRow][col] + "\t" + (!linkedGenotype.equals("N/A") && !linkedGenotype.equals(gtRowNames[minRow])));

        }
        matchedTrOut.close();


        // perform PCA


        // visualize results
        double[] expPC1EigenVector = new double[scoringMatrix.colObjects.size()];
        double[] genPC1EigenVector = new double[scoringMatrix.rowObjects.size()];
        double genvarPC1 = 0;
        double expvarPC1 = 0;

        MixupMapperVisualization v = new MixupMapperVisualization();
//        v.createROC(scoringMatrix, genotypeToTrait, traitToGenotype, outDir);


        String subtitle = null;
        if (expPC1EigenVector == null && genPC1EigenVector == null) {
        } else if (genPC1EigenVector != null && expPC1EigenVector == null) {
            subtitle = ("Genomic variance explained by PC1: " + (new java.text.DecimalFormat("#.###", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(genvarPC1));
        } else if (genPC1EigenVector == null && genPC1EigenVector != null) {
            subtitle = ("Gene expression variance explained by PC1: " + (new java.text.DecimalFormat("#.###", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(expvarPC1));
        } else {
            subtitle = ("Gene expression variance explained by PC1: " + (new java.text.DecimalFormat("#.###", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(expvarPC1) + " | Genomic variance explained by PC1: " + (new java.text.DecimalFormat("#.###", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(genvarPC1));
        }

        v.plotHeatMap(scoringMatrix, "HeatMap", subtitle, "Traits", "Genotypes", expPC1EigenVector, genPC1EigenVector, outDir);


        /*
         * String outfile;
         String plotTitle;

        
         */



    }
}
