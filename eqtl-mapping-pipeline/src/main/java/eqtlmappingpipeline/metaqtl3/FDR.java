/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import eqtlmappingpipeline.metaqtl3.graphics.QQPlot;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class FDR {

//    public static String permutationDir = null;
//    public static String outputDir = null;
    public enum FDRMethod {

        PROBELEVEL, GENELEVEL, FULL
    };

    public enum FileFormat {

        LARGE, REDUCED
    };

    /**
     * calculate the FalseDiscoveryRate for the discovered eQTLS
     *
     * @param eQTLTextFileLoc the location where the eQTL text files are stored
     * @param nrPermutationsFDR number of permutations performed
     * @param maxNrMostSignificantEQTLs maximum number of eQTLs to output
     * @param fdrcutoff the FDR cutoff
     * @param createQQPlot create a QQ plot after performing FDR calculations
     * @param outputDir set an alternate directory for output
     * @param permutationDir set an alternate directory for permutation files
     * @throws IOException
     */
    public static void calculateFDR(String eQTLTextFileLoc, int nrPermutationsFDR, int maxNrMostSignificantEQTLs, double fdrcutoff, boolean createQQPlot, String outputDir, String permutationDir) throws IOException {

        if (eQTLTextFileLoc == null || eQTLTextFileLoc.length() == 0) {
            throw new IllegalArgumentException("File containing real effects is not specified.");
        }
        if (nrPermutationsFDR < 1) {
            throw new IllegalArgumentException("Need at least one permutation to determine FDR");
        }
        if (maxNrMostSignificantEQTLs < 1) {
            throw new IllegalArgumentException("Need at least a single effect to perform FDR estimation");
        }
        if (fdrcutoff < 0 || fdrcutoff > 1) {
            throw new IllegalArgumentException("FDR threshold should be between 0.0 and 1.0! (Specified: " + fdrcutoff + ")");
        }

        //Load permuted data:
//        // load values for each permutation round:
        if (permutationDir == null) {
            permutationDir = eQTLTextFileLoc;
        }

        if (outputDir == null) {
            outputDir = eQTLTextFileLoc;
        }

        String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + 1 + ".txt.gz";
        TextFile tf = new TextFile(fileString, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int nrColsInPermutedFiles = elems.length;
        tf.close();

        System.out.println(nrColsInPermutedFiles + " columns in permuted QTL file.");
        // new permutationfile format requires different column layout...
        if (nrColsInPermutedFiles > 7) {
            System.out.println("Large permutation files detected.");
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.FULL, outputDir, permutationDir, createQQPlot);
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.PROBELEVEL, outputDir, permutationDir, createQQPlot);
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.GENELEVEL, outputDir, permutationDir, createQQPlot);
        } else {
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.FULL, outputDir, permutationDir, createQQPlot);
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.PROBELEVEL, outputDir, permutationDir, createQQPlot);
            if (nrColsInPermutedFiles >= 4) {
                runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.GENELEVEL, outputDir, permutationDir, createQQPlot);
            }
        }
    }

    private static void runFDR(String baseDir, int nrPermutationsFDR, int maxNrMostSignificantEQTLs,
            double fdrcutoff, FileFormat f, FDRMethod m, String outputDir, String permutationDir, boolean createQQPlot) throws IOException {
        //Load permuted data:
        // load values for each permutation round:
        System.out.println("");
        if (m == FDRMethod.GENELEVEL) {
            System.out.println("Performing gene level FDR");
        } else if (m == FDRMethod.PROBELEVEL) {
            System.out.println("Performing probe level FDR");
        } else if (m == FDRMethod.FULL) {
            System.out.println("Determining the FDR using all data");
        }

        double[][] permutedPValues = new double[nrPermutationsFDR][maxNrMostSignificantEQTLs];
//        ProgressBar pb = new ProgressBar(nrPermutationsFDR, "Reading permuted data:");
        int nrEQTLs = -1;
        System.out.println("Reading permuted files");
        for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
            String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + (permutationRound + 1) + ".txt.gz";
            System.out.println(fileString);
            // initialize the p-value matrix
            for (int s = 0; s < maxNrMostSignificantEQTLs; s++) {
                permutedPValues[permutationRound][s] = 1;
            }

            // read the permuted eqtl output
            TextFile gz = new TextFile(fileString, TextFile.R);

            String[] header = gz.readLineElems(TextFile.tab);
            int snpcol = -1;
            int pvalcol = -1;
            int probecol = -1;
            int genecol = -1;

            //PValue  SNP     Probe   Gene 
            for (int col = 0; col < header.length; col++) {
                if (header[col].equals("PValue")) {
                    pvalcol = col;
                }
                if (header[col].equals("SNP")) {
                    snpcol = col;
                }
                if (header[col].equals("Probe")) {
                    probecol = col;
                }
                if (header[col].equals("Gene")) {
                    genecol = col;
                }
            }
           //if (f == FileFormat.REDUCED) {
            //PValue  SNP     Probe   Gene
            if (snpcol == -1 || pvalcol == -1 || probecol == -1 && genecol == -1) {
                System.out.println("Column not found in permutation file: " + fileString);
                System.out.println("PValue: " + pvalcol);
                System.out.println("SNP: " + snpcol);
                System.out.println("Probe: " + probecol);
                System.out.println("Gene: " + genecol);
            }
            // }
            String[] data = gz.readLineElemsReturnReference(TextFile.tab);
            int itr = 0;

            HashSet<String> visitedEffects = new HashSet<String>();
            while (data != null) {

                if (data.length != 0) {
                    if (itr > maxNrMostSignificantEQTLs - 1) {
                        System.out.println("Breaking because: " + itr);
                        break;
                    } else {
                        int filteronColumn = -1;
                        String fdrId = "";
                        if (f == FileFormat.REDUCED) {
                            if (m == FDRMethod.FULL) {
                                fdrId = data[snpcol] + "-" + data[probecol];
                                filteronColumn = probecol;
                            } else if (m == FDRMethod.GENELEVEL && data.length > 3) {
                                fdrId = data[genecol];
                                filteronColumn = genecol;
                            } else {
                                fdrId = data[probecol];
                                filteronColumn = probecol;
                            }

                        } else {
                            if (m == FDRMethod.GENELEVEL) {
                                fdrId = data[eQTLTextFile.HUGO];
                                filteronColumn = eQTLTextFile.HUGO;
                            } else if (m == FDRMethod.PROBELEVEL) {
                                fdrId = data[4];
                                filteronColumn = 4;
                            } else {
                                fdrId = data[1] + "-" + data[4];
                                filteronColumn = 4;
                            }
                        }

                        // take top effect per gene / probe
                        if (data.length > filteronColumn) {

                            if (!fdrId.equals("-") && !visitedEffects.contains(fdrId)) {
                                permutedPValues[permutationRound][itr] = Double.parseDouble(data[0]);
                                visitedEffects.add(fdrId);
                                if (itr > 0 && permutedPValues[permutationRound][itr - 1] > permutedPValues[permutationRound][itr]) {
                                    System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
                                    System.exit(-1);
                                }
                                itr++;
                            }
                        } else {
                            System.out.println(Strings.concat(data, Strings.tab));
                        }
                        data = gz.readLineElemsReturnReference(TextFile.tab);
                    }
                }
            }
            gz.close();

            if (nrEQTLs == -1) {
                nrEQTLs = itr;
            }
        }
//        pb.close();
        maxNrMostSignificantEQTLs = nrEQTLs;
        double[][] actualPermutedPvals = new double[nrPermutationsFDR][nrEQTLs];
        for (int p = 0; p < nrPermutationsFDR; p++) {
            System.arraycopy(permutedPValues[p], 0, actualPermutedPvals[p], 0, nrEQTLs);
        }
        permutedPValues = actualPermutedPvals;

        //Load real data:
        double[] pValues = new double[maxNrMostSignificantEQTLs];
        for (int s = 0; s < maxNrMostSignificantEQTLs; s++) {
            pValues[s] = 1;
        }
        int nrRealDataEQTLs = 0;

        ArrayList<String> alEQTLData = new ArrayList<String>();

        String fileString = baseDir + "/eQTLs.txt.gz";
        if (!Gpio.exists(fileString)) {
            System.out.println("Could not find file: " + fileString + " trying un-GZipped file....");
            fileString = baseDir + "/eQTLs.txt";
        }
        if (!Gpio.exists(fileString)) {
            System.out.println("Could not find file: " + fileString);
            System.exit(0);
        }

        TextFile realEQTLs = new TextFile(fileString, TextFile.R);

        System.out.println(realEQTLs.countLines() + " lines in " + realEQTLs.getFileName());
        String header = realEQTLs.readLine();
        String str = realEQTLs.readLine();

// REAL DATA PROCESSING
        ProgressBar pb2 = new ProgressBar(maxNrMostSignificantEQTLs, "Now reading real data: " + fileString);
        int itr = 0;
        HashSet<String> visitedEffects = new HashSet<String>();
        while (str != null) {
            if (itr > maxNrMostSignificantEQTLs - 1) {
                break;
            } else {
                int filteronColumn = -1;
                String fdrId = "";
                String[] data = Strings.tab.split(str);

                if (m == FDRMethod.GENELEVEL) {
                    fdrId = data[eQTLTextFile.HUGO];
                    filteronColumn = eQTLTextFile.HUGO;
                } else if (m == FDRMethod.PROBELEVEL) {
                    fdrId = data[4];
                    filteronColumn = 4;
                } else {
                    fdrId = data[1] + "-" + data[4];
                    filteronColumn = 4;
                }

                if (data.length > filteronColumn) {
                    if (!fdrId.equals("-") && !visitedEffects.contains(fdrId)) {
                        alEQTLData.add(str);
                        pValues[itr] = Double.parseDouble(data[0]);
                        if (itr > 0 && pValues[itr - 1] > pValues[itr]) {
                            System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
                            System.exit(-1);
                        }
                        visitedEffects.add(fdrId);
                        itr++;
                    }
                }

                pb2.iterate();
                str = realEQTLs.readLine();
            }

        }

        pb2.close();
        nrRealDataEQTLs = itr;
        realEQTLs.close();

        //Process a certain P-Value, determine how many are significant:
        //System.out.println("\n");
        //System.out.println("Significant detected eQTLs:");
        int nrSignificantEQTLs = 0;
        boolean[] pValueSignificant = new boolean[nrRealDataEQTLs];

        String outFileName = "";
        String outFileNameAll = "";

        if (outputDir == null) {
            outputDir = baseDir;
        }

        if (m == FDRMethod.GENELEVEL) {
            outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + "-GeneLevel.txt";
            outFileNameAll = outputDir + "/eQTLsFDR-GeneLevel.txt.gz";
        } else if (m == FDRMethod.PROBELEVEL) {
            outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + "-ProbeLevel.txt";
            outFileNameAll = outputDir + "/eQTLsFDR-ProbeLevel.txt.gz";
        } else {
            outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + ".txt";
            outFileNameAll = outputDir + "/eQTLsFDR.txt.gz";
        }

        double pValueThresholdCorrespondingToRequestedFDR = -1;

        // for all p-values, determine how many eQTL we find below its p-value, in comparison to random data
        double previousPValueThreshold = -1;

        HashMap<Double, Integer> hashUniquePValues = new HashMap<Double, Integer>();
        ArrayList<Double> vecUniquePValues = new ArrayList<Double>();
        double previousPValue = -1;
        for (int p = 0; p < nrRealDataEQTLs; p++) {
            double pValue = pValues[p];
            if (previousPValue != pValue) {
                if (!hashUniquePValues.containsKey(pValue)) {
                    hashUniquePValues.put(pValue, null);
                    vecUniquePValues.add(pValue);
                }
                previousPValue = pValue;
            }
        }
        for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
            previousPValue = -1;
            for (int pPerm = 0; pPerm < maxNrMostSignificantEQTLs; pPerm++) {
                double pValue = permutedPValues[permutationRound][pPerm];
                if (previousPValue != pValue) {
                    if (!hashUniquePValues.containsKey(pValue)) {
                        hashUniquePValues.put(pValue, null);
                        vecUniquePValues.add(pValue);
                    }
                    previousPValue = pValue;
                }
            }
        }
        double[] uniquePValues = new double[hashUniquePValues.size()];
        Collections.sort(vecUniquePValues);
        hashUniquePValues.clear();
        for (int u = 0; u < vecUniquePValues.size(); u++) {
            uniquePValues[u] = vecUniquePValues.get(u);
            hashUniquePValues.put(uniquePValues[u], u);
        }

        System.out.println("Number of unique P Values:\t" + hashUniquePValues.size());

        long[] uniquePValuesNrEQTLsWithThisPValue = new long[hashUniquePValues.size()];
        long[] uniquePValuesNrEQTLsWithThisPValueCumulative = new long[hashUniquePValues.size()];
        previousPValue = -1;
        int pValueIndex = -1;
        for (int p = 0; p < nrRealDataEQTLs; p++) {
            double pValue = pValues[p];
            if (previousPValue != pValue || pValueIndex == -1) {
                pValueIndex = hashUniquePValues.get(pValue);
                previousPValue = pValue;
            }
            uniquePValuesNrEQTLsWithThisPValue[pValueIndex]++;
        }
        for (int p = 0; p < hashUniquePValues.size(); p++) {
            uniquePValuesNrEQTLsWithThisPValueCumulative[p] = uniquePValuesNrEQTLsWithThisPValue[p];
            if (p > 0) {
                uniquePValuesNrEQTLsWithThisPValueCumulative[p] += uniquePValuesNrEQTLsWithThisPValueCumulative[p - 1];
            }
        }

        int[][] permUniquePValuesNrEQTLsWithThisPValue = new int[hashUniquePValues.size()][nrPermutationsFDR];
        int[][] permUniquePValuesNrEQTLsWithThisPValueCumulative = new int[hashUniquePValues.size()][nrPermutationsFDR];

        for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
            previousPValue = -1;
            pValueIndex = -1;
            for (int p = 0; p < maxNrMostSignificantEQTLs; p++) {
                double pValue = permutedPValues[permutationRound][p];
                if (previousPValue != pValue || pValueIndex == -1) {
                    pValueIndex = hashUniquePValues.get(pValue);
                    previousPValue = pValue;
                }
                permUniquePValuesNrEQTLsWithThisPValue[pValueIndex][permutationRound]++;
            }
            for (int p = 0; p < hashUniquePValues.size(); p++) {
                permUniquePValuesNrEQTLsWithThisPValueCumulative[p][permutationRound] = permUniquePValuesNrEQTLsWithThisPValue[p][permutationRound];
                if (p > 0) {
                    permUniquePValuesNrEQTLsWithThisPValueCumulative[p][permutationRound] += permUniquePValuesNrEQTLsWithThisPValueCumulative[p - 1][permutationRound];
                }
            }
        }

        double[] fdrUniquePValues = new double[hashUniquePValues.size()];
        for (int p = 0; p < hashUniquePValues.size(); p++) {
            double meanNrEQTLsPerm = JSci.maths.ArrayMath.mean(permUniquePValuesNrEQTLsWithThisPValueCumulative[p]);
            double fdrVal = meanNrEQTLsPerm / (double) uniquePValuesNrEQTLsWithThisPValueCumulative[p];
            if (fdrVal > 1) {
                fdrVal = 1;
            }
            fdrUniquePValues[p] = fdrVal;

            // FDR = FP / TP;
            //System.out.println(p + "\t" + uniquePValues[p] + "\t" + uniquePValuesNrEQTLsWithThisPValueCumulative[p] + "\t" + meanNrEQTLsPerm + "\t" + fdrUniquePValues[p]);
        }

        // Write results
        if (m == FDRMethod.FULL) {
            TextFile outSignificant = new TextFile(outFileName, TextFile.W);
            outSignificant.write(header + "\tFDR\n");
            for (int p = 0; p < nrRealDataEQTLs; p++) {
                pValueIndex = hashUniquePValues.get(pValues[p]);
                String output = alEQTLData.get(p) + "\t" + fdrUniquePValues[pValueIndex];
                if (fdrUniquePValues[pValueIndex] <= fdrcutoff) {
                    pValueSignificant[p] = true;
                    outSignificant.writeln(output);
                    nrSignificantEQTLs++;
                }
            }
            outSignificant.close();
        } else {
            TextFile outSignificant = new TextFile(outFileName, TextFile.W);
            outSignificant.write(header + "\tFDR\n");
            realEQTLs.open();
            realEQTLs.readLine();
            String[] data = realEQTLs.readLineElems(TextFile.tab);
            Integer previousIndex = 0;
            while (data != null) {
                double p = Double.parseDouble(data[0]);
                Integer pIndex = hashUniquePValues.get(p);
                String output = "";

                if (pIndex == null) {
                    double fdrEstimate = 1;
                    if (previousIndex + 1 < fdrUniquePValues.length) {
                        fdrEstimate = (fdrUniquePValues[previousIndex] + fdrUniquePValues[previousIndex + 1]) / 2;
                    }
                    output = Strings.concat(data, Strings.tab) + "\t" + fdrEstimate;
                    if (fdrEstimate <= fdrcutoff) {
                        outSignificant.writeln(output);
                        nrSignificantEQTLs++;
                    }
                } else {
                    double fdr = fdrUniquePValues[pIndex];
                    output = Strings.concat(data, Strings.tab) + "\t" + fdr;
                    if (fdr <= fdrcutoff) {
                        outSignificant.writeln(output);
                        nrSignificantEQTLs++;
                    }
                    previousIndex = pIndex;
                }
                data = realEQTLs.readLineElems(TextFile.tab);
            }
            realEQTLs.close();
            outSignificant.close();
        }

        TextFile outAll = new TextFile(outFileNameAll, TextFile.W);
        outAll.writeln(header + "\tFDR");

        realEQTLs.open();
        realEQTLs.readLine();
        String[] data = realEQTLs.readLineElems(TextFile.tab);
        Integer previousIndex = 0;

        while (data != null) {
            double p = Double.parseDouble(data[0]);
            Integer pIndex = hashUniquePValues.get(p);
            String output = "";
            if (pIndex == null) {
                double fdrEstimate = 1;
                if (previousIndex + 1 < fdrUniquePValues.length) {
                    fdrEstimate = (fdrUniquePValues[previousIndex] + fdrUniquePValues[previousIndex + 1]) / 2;
                }
                output = Strings.concat(data, Strings.tab) + "\t" + fdrEstimate;
            } else {
                output = Strings.concat(data, Strings.tab) + "\t" + fdrUniquePValues[pIndex];
                previousIndex = pIndex;
            }

            outAll.writeln(output);
            data = realEQTLs.readLineElems(TextFile.tab);
        }
        realEQTLs.close();
        outAll.close();

        //System.out.println("");
        String output = "Number of significant eQTLs:\t" + nrSignificantEQTLs;
        System.out.println(output);

        String fileSuffix = "";
        if (m == FDRMethod.GENELEVEL) {
            fileSuffix = "-GeneLevel";
        } else if (m == FDRMethod.PROBELEVEL) {
            fileSuffix = "-ProbeLevel";
        }

        if (createQQPlot) {
            QQPlot qq = new QQPlot();
            String fileName = baseDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + "-QQPlot.pdf";
            qq.draw(fileName, fdrcutoff, nrPermutationsFDR,
                    maxNrMostSignificantEQTLs, permutedPValues, nrRealDataEQTLs, pValues,
                    pValueSignificant, nrSignificantEQTLs);
        }

        generateESNPsFile(outputDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + ".txt", outputDir + "/eQTLSNPsFDR" + fdrcutoff + fileSuffix + ".txt");
        generateEProbesFile(outputDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + ".txt", outputDir + "/eQTLProbesFDR" + fdrcutoff + fileSuffix + ".txt");

    }

    /**
     * Generates an eQTL SNPs file
     *
     * @param inputFile
     * @param outputFile
     */
    private static void generateESNPsFile(String inputFile, String outputFile) throws IOException {

        //Determine the number of unique SNPs and what the most significant eQTL for this SNP is:
        HashMap<String, Double> hashTopSNPsAbsZScore = new HashMap<String, Double>();
        HashMap<String, String> hashTopSNPsAnnotation = new HashMap<String, String>();
        ArrayList<String> vecTopSNPs = new ArrayList<String>();
        String header = "";

        TextFile in = new TextFile(inputFile, TextFile.R);
        String str = in.readLine();
        header = str;
        while ((str = in.readLine()) != null) {
            String[] data = str.split("\t");
            double absZScore = Math.abs(Double.parseDouble(data[10]));
            String rsName = data[1];
            if (hashTopSNPsAbsZScore.get(rsName) != null) {
                double absZScorePrevious = hashTopSNPsAbsZScore.get(rsName);
                if (absZScore > absZScorePrevious) {
                    hashTopSNPsAbsZScore.put(rsName, absZScore);
                    hashTopSNPsAnnotation.put(rsName, str);
                }
            } else {
                vecTopSNPs.add(rsName);
                hashTopSNPsAbsZScore.put(rsName, absZScore);
                hashTopSNPsAnnotation.put(rsName, str);
            }
        }
        in.close();

        //Per eSNP write the most significant eQTL that has been recorded to file:
        TextFile out = new TextFile(outputFile, TextFile.W);
        out.writeln(header);
        for (int s = 0; s < vecTopSNPs.size(); s++) {
            String rsName = vecTopSNPs.get(s);
            String annotation = hashTopSNPsAnnotation.get(rsName);
            out.writeln(annotation);
        }
        out.close();

        String output = " - Number of unique SNPs, constituting an eQTL:\t" + vecTopSNPs.size();
        System.out.println(output);

    }

    /**
     * Generates a file containing combinations of the most significant eQTL
     * with probes
     *
     * @param inputFile location of text file containing data on eQTLs
     * @param outputFile the location where the output will be written
     */
    private static void generateEProbesFile(String inputFile, String outputFile) throws IOException {

        //Determine the number of unique probes and what the most significant eQTL for each probe is:
        HashMap<String, Double> hashTopProbesAbsZScore = new HashMap<String, Double>();
        HashMap<String, String> hashTopProbesAnnotation = new HashMap<String, String>();
        ArrayList<String> vecTopProbes = new ArrayList<String>();
        String header = "";

        TextFile in = new TextFile(inputFile, TextFile.R);
        String str = in.readLine();
        header = str;
        while ((str = in.readLine()) != null) {
            String[] data = str.split("\t");
            double absZScore = Math.abs(Double.parseDouble(data[10]));
            String probe = data[4];
            if (hashTopProbesAbsZScore.get(probe) != null) {
                double absZScorePrevious = hashTopProbesAbsZScore.get(probe);
                if (absZScore > absZScorePrevious) {
                    hashTopProbesAbsZScore.put(probe, absZScore);
                    hashTopProbesAnnotation.put(probe, str);
                }
            } else {
                vecTopProbes.add(probe);
                hashTopProbesAbsZScore.put(probe, absZScore);
                hashTopProbesAnnotation.put(probe, str);
            }
        }
        in.close();

        //Per eSNP write the most significant eQTL that has been recorded to file:
        TextFile out = new TextFile(outputFile, TextFile.W);
        out.writeln(header);
        for (int s = 0; s < vecTopProbes.size(); s++) {
            String probe = vecTopProbes.get(s);
            String annotation = hashTopProbesAnnotation.get(probe);
            out.writeln(annotation);

        }

        out.close();

        String output = " - Number of unique probes, constituting an eQTL:\t" + vecTopProbes.size();
        System.out.println(output);
    }
}
