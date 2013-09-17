/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class eQTLFileCompare {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\t");
    private int nrShared = 0;
    private int nrOpposite = 0;

    public int getNrShared() {
        return nrShared;
    }

    public int getNrOpposite() {
        return nrOpposite;
    }

    public eQTLFileCompare() {
    }

    public eQTLFileCompare(String[] args) {


        String out = null;
        String file1 = null;
        String file2 = null;
        boolean matchOnGeneName = false;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--file1")) {
                file1 = val;
            } else if (arg.equals("--file2")) {
                file2 = val;
            } else if (arg.equals("--genebased")) {
                matchOnGeneName = true;
            }
        }

        if (out != null && file1 != null && file2 != null) {
            try {
                compareOverlapAndZScoreDirectionTwoEQTLFiles(file1, file2, out, matchOnGeneName);
            } catch (IOException ex) {
                Logger.getLogger(eQTLFileCompare.class.getName()).log(Level.SEVERE, null, ex);
            } catch (Exception ex) {
                Logger.getLogger(eQTLFileCompare.class.getName()).log(Level.SEVERE, null, ex);
            }

        } else {
            printUsage();
        }
    }

    private void printUsage() {
        System.out.print("QTL File comparison\n" + ConsoleGUIElems.LINE);
        System.out.println("Compares two eQTL files with each other.");

        System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);

        System.out.println("--out\t\tstring\t\tOutput file name\n"
                + "--file1\t\tstring\t\tLocation of file 1\n"
                + "--file2\t\tstring\t\tLocation of file 2\n"
                + "--genebased\t\t\tPerform comparison on the basis of gene names (optional, defaults to probe based comparison)\n");
    }

    public void compareOverlapAndZScoreDirectionTwoEQTLFiles(String file1, String file2, String outputFile, boolean matchOnGeneName) throws IOException, Exception {


        double filterOnFDR = -1; //Do we want to use another FDR measure? When set to -1 this is not used at all.

        HashMap<String, String> hashConvertProbeNames = new HashMap<String, String>(); //When comparing two eQTL files, run on different platforms, we can convert the probe names from one platform to the other, accommodating this comparison, example: hashConvertProbeNames.put(probeNameInFile1, equivalentProbeNameInFile2);
        HashSet<String> hashExcludeEQTLs = new HashSet<String>();   //We can exclude some eQTLs from the analysis. If requested, put the entire eQTL string in this HashMap for each eQTL
        HashSet<String> hashConfineAnalysisToSubsetOfProbes = new HashSet<String>(); //We can confine the analysis to only a subset of probes. If requested put the probe name in this HapMap
        HashSet<String> hashTestedSNPsThatPassedQC = null; //We can confine the analysis to only those eQTLs for which the SNP has been successfully passed QC, otherwise sometimes unfair comparisons are made. If requested, put the SNP name in this HashMap

        //Now load the eQTLs for file 1:
        HashMap<String, String[]> hashEQTLs = new HashMap<String, String[]>();
        HashSet<String> hashUniqueProbes = new HashSet<String>();
        HashSet<String> hashUniqueGenes = new HashSet<String>();


        TextFile in = new TextFile(file1, TextFile.R);
        in.readLine();
        String[] data = null;
        
        while ((data = SPLIT_ON_TAB.split(in.readLine())) != null) {
            if (hashConvertProbeNames.size() > 0) {
                if (hashConvertProbeNames.containsKey(data[4].trim())) {
                    data[4] = hashConvertProbeNames.get(data[4].trim());
                }
            }
            if (filterOnFDR == -1 || Double.parseDouble(data[18]) <= filterOnFDR) {
                if (hashConfineAnalysisToSubsetOfProbes.isEmpty() || hashConfineAnalysisToSubsetOfProbes.contains(data[4])) {
                    if (matchOnGeneName) {
                        if (data[16].length() > 1) {
                            if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
                                hashEQTLs.put(data[1] + "\t" + data[16], data);
                                hashUniqueProbes.add(data[4]);
                                hashUniqueGenes.add(data[16]);
                            }
                        }
                    } else {
                        if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
                            hashEQTLs.put(data[1] + "\t" + data[4], data);
                            hashUniqueProbes.add(data[4]);
                            hashUniqueGenes.add(data[16]);
                        }
                    }
                }
            }
        }
        in.close();

        //Initialize Graphics2D for the Z-Score allelic direction comparison:
        int width = 1000;
        int height = 1000;
        int margin = 100;
        int x0 = margin;
        int x1 = width - margin;
        int y0 = margin;
        int y1 = height - margin;

        ZScorePlot zs = new ZScorePlot();
        String zsOutFileName = outputFile + "-ZScoreComparison.png";
        zs.init(2, new String[]{"Dataset1", "Dataset2"}, false, zsOutFileName);

        //Variables holding variousStatistics:
        int nreQTLsIdenticalDirection = 0;
        int nreQTLsOppositeDirection = 0;
        HashMap<String, Integer> hashEQTLNrTimesAssessed = new HashMap<String, Integer>();
        ArrayList<String> vecEQTLNrTimesAssessed = new ArrayList<String>();

        HashMap<String, String[]> hashEQTLs2 = new HashMap<String, String[]>();
        HashSet<String> hashUniqueProbes2 = new HashSet<String>();
        HashSet<String> hashUniqueGenes2 = new HashSet<String>();
        HashSet<String> hashUniqueProbesOverlap = new HashSet<String>();
        HashSet<String> hashUniqueGenesOverlap = new HashSet<String>();

        int counterFile2 = 0;
        int overlap = 0;
        ArrayList<Double> vecX = new ArrayList<Double>();
        ArrayList<Double> vecY = new ArrayList<Double>();

        //Vector holding all opposite allelic effects:
        ArrayList<String> vecOppositeEQTLs = new ArrayList<String>();

        //Now process file 2:

        in = new TextFile(file2, TextFile.R);
        in.readLine();
        TextFile log = new TextFile(outputFile + "-eQTLComparisonLog.txt", TextFile.W);
        int lineno = 1;
        data = null;
        TextFile identicalOut = new TextFile(outputFile + "-eQTLsWithIdenticalDirecton.txt.gz", TextFile.W);
        while ((data = SPLIT_ON_TAB.split(in.readLine())) != null) {

            if (filterOnFDR == -1 || Double.parseDouble(data[18]) <= filterOnFDR) {

                if (hashConvertProbeNames.size() > 0) {
                    if (hashConvertProbeNames.containsKey(data[4].trim())) {
                        data[4] = hashConvertProbeNames.get(data[4].trim());
                    }
                }
                if (hashConfineAnalysisToSubsetOfProbes.isEmpty() || hashConfineAnalysisToSubsetOfProbes.contains(data[4])) {
                    if (matchOnGeneName) {
                        if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
                            if (data[16].length() > 1) {
                                hashUniqueProbes2.add(data[4]);
                                hashUniqueGenes2.add(data[16]);
                                if (!hashEQTLs2.containsKey(data[1] + "\t" + data[16])) {
                                    hashEQTLs2.put(data[1] + "\t" + data[16], data);
                                    counterFile2++;
                                }
                            }
                        }
                    } else {
                        if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
                            //hashEQTLs2.put(data[1] + "\t" + data[4], str);
                            hashUniqueProbes2.add(data[4]);
                            hashUniqueGenes2.add(data[16]);
                            counterFile2++;
                        }
                    }
                    String[] eQTL = null;
                    String identifier = null;
                    if (matchOnGeneName) {
                        //System.out.println(data.length + "\t" + str);
                        if (data.length > 16 && data[16].length() > 1) {
                            if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[16])) {
                                identifier = data[1] + "\t" + data[16];
                                if (hashEQTLs.containsKey(identifier)) {
                                    eQTL = hashEQTLs.get(identifier);
                                }
                            }
                        }
                    } else {
                        if (!hashExcludeEQTLs.contains(data[1] + "\t" + data[4])) {
                            identifier = data[1] + "\t" + data[4];
                            if (hashEQTLs.containsKey(identifier)) {
                                eQTL = hashEQTLs.get(identifier);
                            }
                        }
                    }

                    if (eQTL == null) {

                        //The eQTL, present in file 2 is not present in file 1:
                        double pValue = Double.parseDouble(data[0]);
                        //if (pValue < 1E-4) {
                        if (hashTestedSNPsThatPassedQC == null || hashTestedSNPsThatPassedQC.contains(data[1])) {
                            log.write("eQTL Present In New file But Not In Original File:\t" + identifier + "\t" + data[0] + "\t" + data[2] + "\t" + data[3] + "\t" + data[16] + "\n");
                        }
                        //}
                        double zScore2 = Double.parseDouble(data[10]);
                        int posX = 500 + (int) 0;
                        int posY = 500 - (int) Math.round(zScore2 * 10);
                        zs.draw(null, zScore2, 0, 1);

                    } else {
                        String[] eQtlData = eQTL;
                        boolean identicalProbe = true;
                        String probe = data[4];
                        String probeFound = eQtlData[4];
                        if (!probe.equals(probeFound)) {
                            identicalProbe = false;
                        }

                        overlap++;
                        hashUniqueProbesOverlap.add(data[4]);
                        hashUniqueGenesOverlap.add(data[16]);
                        if (!hashEQTLNrTimesAssessed.containsKey(identifier)) {
                            hashEQTLNrTimesAssessed.put(identifier, 1);
                            vecEQTLNrTimesAssessed.add(identifier);
                        } else {
                            hashEQTLNrTimesAssessed.put(identifier, 1 + hashEQTLNrTimesAssessed.get(identifier));
                        }
                        String alleles = eQtlData[8];
                        String alleleAssessed = eQtlData[9];

                        String correlations[] = (eQtlData[17]).split(";");
                        double correlation = 0;
                        int numCorr1 = 0;
                        for (int c = 0; c < correlations.length; c++) {
                            try {
                                if (!correlations[c].equals("-")) {
                                    correlation += Double.parseDouble(correlations[c]);
                                    numCorr1++;
                                }
                            } catch (Exception e) {
                            }
                        }

                        correlation /= (double) numCorr1;
//                       if(numCorr1 == 0){
//                           System.out.println("Warning: no correlations defined for eqtl file 1");
//                       }
                        double zScore = Double.parseDouble(eQtlData[10]);
//                        double pValue = Double.parseDouble(eQtlData[0]);
                        String alleles2 = data[8];
                        String alleleAssessed2 = data[9];
                        double zScore2 = Double.parseDouble(data[10]);
//                        double pValue2 = Double.parseDouble(data[0]);
                        String correlations2[] = data[17].split(";");
                        double correlation2 = 0;

                        boolean alleleflipped = false;
                        if (!alleleAssessed.equals(data[9])) {
                            if (data[9].equals(eQtlData[8].split("/")[0])) {
                                alleleflipped = true;
                            } else {
//                               System.out.println("WTF BBQ!");
                            }
                        }

                        int numCorr2 = 0;
                        for (int c = 0; c < correlations2.length; c++) {
                            try {
                                if (!correlations2[c].equals("-")) {

                                    correlation2 += (Double.parseDouble(correlations2[c]));


                                    numCorr2++;
                                }
                            } catch (NumberFormatException e) {
                            }
                        }
//                       if(numCorr2 == 0){
//                           System.out.println("Warning: no correlations defined for eqtl file 2");
//                       }
                        correlation2 /= (double) numCorr2;
                        if (alleleflipped) {
                            correlation2 = -correlation2;
                        }
                        boolean sameDirection = false;
                        int nrIdenticalAlleles = 0;
                        if (alleles.length() > 2 && alleles2.length() > 2) {
                            for (int a = 0; a < 3; a++) {
                                for (int b = 0; b < 3; b++) {
                                    if (a != 1 && b != 1) {
                                        if (alleles.getBytes()[a] == alleles2.getBytes()[b]) {
                                            nrIdenticalAlleles++;
                                        }
                                    }
                                }
                            }
                        }
                        if (nrIdenticalAlleles == 1) {
                            log.write("Error! SNPs have incompatible alleles!!:\t" + alleles + "\t" + alleles2 + "\t" + identifier + "\n");
                        }
                        if (nrIdenticalAlleles == 0) {
                            alleleAssessed2 = BaseAnnot.getComplement(alleleAssessed2);
                        }
                        if (!alleleAssessed.equals(alleleAssessed2)) {
                            zScore2 = -zScore2;
//                           correlation2 = -correlation2;
                            alleleAssessed2 = alleleAssessed;
                        }

                        //Recode alleles:
                        // if contains T, but no A, take complement
                        if (alleles.contains("T") && !alleles.contains("A")) {
                            alleles = BaseAnnot.getComplement(alleles);
                            alleleAssessed = BaseAnnot.getComplement(alleleAssessed);
                            alleleAssessed2 = BaseAnnot.getComplement(alleleAssessed2);
                        }

                        if (zScore2 * zScore > 0) {
                            sameDirection = true;
                        }

//                       if(correlation != correlation2 && (numCorr1 > 0 && numCorr2 > 0)){
//                           if(Math.abs(correlation - correlation2) > 0.00001){
//                               System.out.println("Correlations are different: "+lineno+"\t"+correlation +"\t"+correlation2+"\t"+str);
//                           }
//                           
//                       }
                        zs.draw(zScore, zScore2, 0, 1);
                        if (!sameDirection) {
                            nreQTLsOppositeDirection++;
                            if (matchOnGeneName) {
                                identifier = data[1] + "\t" + data[16];
                                if (!vecOppositeEQTLs.contains(identifier)) {
                                    vecOppositeEQTLs.add(identifier);
                                }
                            } else {
                                identifier = data[1] + "\t" + data[4];
                                if (!vecOppositeEQTLs.contains(identifier)) {
                                    vecOppositeEQTLs.add(identifier);
                                }
                            }
//                            int posX = 500 + (int) Math.round(zScore * 10);
//                            int posY = 500 - (int) Math.round(zScore2 * 10);
                            vecX.add(zScore);
                            vecY.add(zScore2);

                        } else {
                            // write to output
                            identicalOut.writeln(identifier);
                            nreQTLsIdenticalDirection++;
                            if (alleles.length() > 2 && !alleles.equals("A/T") && !alleles.equals("T/A") && !alleles.equals("C/G") && !alleles.equals("G/C")) {
//                                int posX = 500 + (int) Math.round(zScore * 10);
//                                int posY = 500 - (int) Math.round(zScore2 * 10);
                                vecX.add(zScore);
                                vecY.add(zScore2);
                            }
                        }
                    }
                }
            }
            lineno++;
        }
        identicalOut.close();
        log.close();
        in.close();

        double[] valsX = new double[vecX.size()];
        double[] valsY = new double[vecX.size()];
        for (int v = 0; v < valsX.length; v++) {
            valsX[v] = ((Double) vecX.get(v)).doubleValue();
            valsY[v] = ((Double) vecY.get(v)).doubleValue();
        }
        if (valsX.length > 2) {
            double correlation = JSci.maths.ArrayMath.correlation(valsX, valsY);
            double r2 = correlation * correlation;

            /*
             * randomEngine = new cern.jet.random.engine.DRand();
             tDistColt = new cern.jet.random.StudentT(olsY.length - 4, randomEngine);
             */
            cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = new cern.jet.random.tdouble.engine.DRand();
            cern.jet.random.tdouble.StudentT tDistColt = new cern.jet.random.tdouble.StudentT(valsX.length - 2, randomEngine);
            double pValuePearson = 1;
            double tValue = correlation / (Math.sqrt((1 - r2) / (double) (valsX.length - 2)));
            if (tValue < 0) {
                pValuePearson = tDistColt.cdf(tValue);
            } else {
                pValuePearson = tDistColt.cdf(-tValue);
            }
            pValuePearson *= 2;
            System.out.println("\nCorrelation between the Z-Scores of the overlapping set of eQTLs:\t" + correlation + "\tP-Value:\t" + pValuePearson);
        }



        TextFile out = new TextFile(outputFile + "-OppositeEQTLs.txt", TextFile.W);
        for (int o = 0; o < vecOppositeEQTLs.size(); o++) {
            out.write(vecOppositeEQTLs.get(o) + "\n");
        }
        out.close();

        in.close();

        zs.write(zsOutFileName);

        TextFile outSummary = new TextFile(outputFile + "-Summary.txt", TextFile.W);


        System.out.println("");
        System.out.println("Nr of eQTLs:\t" + hashEQTLs.size() + "\tin file:\t" + file1 + "\tNrUniqueProbes:\t" + hashUniqueProbes.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes.size());
        outSummary.writeln("Nr of eQTLs:\t" + hashEQTLs.size() + "\tin file:\t" + file1 + "\tNrUniqueProbes:\t" + hashUniqueProbes.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes.size());




        System.out.println("Nr of eQTLs:\t" + counterFile2 + "\tin file:\t" + file2 + "\tNrUniqueProbes:\t" + hashUniqueProbes2.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes2.size());
        outSummary.writeln("Nr of eQTLs:\t" + counterFile2 + "\tin file:\t" + file2 + "\tNrUniqueProbes:\t" + hashUniqueProbes2.size() + "\tNrUniqueGenes:\t" + hashUniqueGenes2.size());


        System.out.println("Overlap:\t" + overlap + "\tNrUniqueProbesOverlap:\t" + hashUniqueProbesOverlap.size() + "\tNrUniqueGenesOverlap:\t" + hashUniqueGenesOverlap.size());
        outSummary.writeln("Overlap:\t" + overlap + "\tNrUniqueProbesOverlap:\t" + hashUniqueProbesOverlap.size() + "\tNrUniqueGenesOverlap:\t" + hashUniqueGenesOverlap.size());


        System.out.println("");
        outSummary.writeln();

        System.out.println("Nr eQTLs with identical direction:\t" + nreQTLsIdenticalDirection);
        outSummary.writeln("Nr eQTLs with identical direction:\t" + nreQTLsIdenticalDirection);

        double proportionOppositeDirection = 100d * (double) nreQTLsOppositeDirection / (double) (nreQTLsOppositeDirection + nreQTLsIdenticalDirection);
        String proportionOppositeDirectionString = (new java.text.DecimalFormat("0.00;-0.00", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(proportionOppositeDirection);

        System.out.println("Nr eQTLs with opposite direction:\t" + nreQTLsOppositeDirection + "\t(" + proportionOppositeDirectionString + "%)");
        outSummary.writeln("Nr eQTLs with opposite direction:\t" + nreQTLsOppositeDirection + "\t(" + proportionOppositeDirectionString + "%)");


        outSummary.close();

        nrShared = hashUniqueProbesOverlap.size();
        nrOpposite = nreQTLsOppositeDirection;

    }
}
