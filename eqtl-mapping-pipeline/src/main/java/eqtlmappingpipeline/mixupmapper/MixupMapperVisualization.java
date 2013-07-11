/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.mixupmapper;

import eqtlmappingpipeline.graphics.curve.ROCCurve;
import eqtlmappingpipeline.graphics.histogram.MultiFrequencyDistributionHistogram;
import eqtlmappingpipeline.graphics.map.MixupMapperHeatMap;
import eqtlmappingpipeline.mixupmapper.containers.Trio;
import eqtlmappingpipeline.mixupmapper.stat.dist.Bin;
import eqtlmappingpipeline.mixupmapper.stat.dist.DiscreteDist;
import java.io.IOException;
import java.util.HashMap;
import java.util.ArrayList;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.TTest;

/**
 *
 * @author harmjan
 */
public class MixupMapperVisualization {

    private double[][] TPR = null;
    private double[][] FPR = null;
    private double[] AUC;
    private boolean loadedFamilyData;
    // For use in calculating log base 10. A log times this is a log base 10.
    private static final double LOG10SCALE = 1 / Math.log(10);
    private HashMap<String, Trio> sampleToTrio;
    private HashMap<String, String> etg;
    private HashMap<String, String> gte;
    private TextFile textOutput;

    // handy static methods
    public void createROC(DoubleMatrixDataset<String, String> result, HashMap<String, String> genotypeToTrait, HashMap<String, String> traitToGenotype, String outputLoc) throws IOException {
        textOutput = new TextFile(outputLoc + "ROC.txt", TextFile.W);

        etg = traitToGenotype;
        gte = genotypeToTrait;

        String[] rowNames = result.rowObjects.toArray(new String[0]);
        String[] colNames = result.colObjects.toArray(new String[0]);

        double[][] matrix = result.rawData;

        int numPositives = traitToGenotype.size();

        int matrixSize = rowNames.length * colNames.length;
        int numTP = numPositives;
        int numTN = matrixSize - numTP;

        int tp = 0;
        int tn = 0;

        double[] truePosValues = new double[numTP];
        double[] trueNegValues = new double[numTN];

        if (loadedFamilyData) {

            int numParents = 0;
            for (int r = 0; r < matrix.length; r++) {

                String rowName = rowNames[r];
                String coupledGenotype = etg.get(rowName);


                Trio relatives = sampleToTrio.get(coupledGenotype);

                boolean isChild = false;

                if (relatives != null) {

                    if (relatives.child.sampleName.equals(coupledGenotype)) {
                        isChild = true;
                    }

                    // check over all columns for values of relatives
                    for (int c = 0; c < matrix[matrix.length - 1].length; c++) {

                        String colName = colNames[c];
                        boolean valueIsFromRelative = false;
                        if (relatives != null) {
                            if (isChild) {
                                if (colName.equals(relatives.parent1.sampleName) || colName.equals(relatives.parent2.sampleName)) {
                                    valueIsFromRelative = true;
                                }
                            } else {
                                if (colName.equals(relatives.child.sampleName)) {
                                    valueIsFromRelative = true;
                                }
                            }
                        }

                        if (valueIsFromRelative) {
                            numParents++;
                        }
                    }
                }
            }

            System.out.println("Num relatives: " + numParents);
            // System.exit(0);


            if (numParents > 0) {
                double[] truePosMinFamilyValues = new double[numTP];
                double[] trueNegMinFamilyValues = new double[numTN - (numParents)];
                double[] parentValues = new double[numParents];

                numParents = 0;
                int tpc = 0;
                int tnc = 0;

//                outLog.log(matrixSize + " - " + truePosMinFamilyValues.length + " - "+ trueNegMinFamilyValues.length + " - "+ parentValues.length);

                for (int r = 0; r < matrix.length; r++) {

                    String rowName = rowNames[r];
                    String coupledGenotype = etg.get(rowName);


                    Trio relatives = sampleToTrio.get(coupledGenotype);


                    boolean isChild = false;

                    if (relatives.child.sampleName.equals(coupledGenotype)) {
                        isChild = true;
                    }

                    for (int c = 0; c < matrix[matrix.length - 1].length; c++) {
                        String colName = colNames[c];
                        boolean valueIsFromRelative = false;

                        if (relatives != null) {
                            if (isChild) {
                                if (colName.equals(relatives.parent1.sampleName) || colName.equals(relatives.parent2.sampleName)) {
                                    valueIsFromRelative = true;
                                }
                            } else {
                                if (colName.equals(relatives.child.sampleName)) {
                                    valueIsFromRelative = true;
                                }
                            }
                        }

                        if (valueIsFromRelative) {
                            parentValues[numParents] = matrix[r][c];
                            numParents++;
                        } else if (etg.get(rowName) != null && etg.get(rowName).equals(colName)) {
                            truePosMinFamilyValues[tpc] = matrix[r][c];
                            tpc++;
                        } else {
                            trueNegMinFamilyValues[tnc] = matrix[r][c];
                            tnc++;
                        }
                    }
                }
                drawDistributionsWithFamilyData(truePosMinFamilyValues, trueNegMinFamilyValues, parentValues, outputLoc);
            }
        }

        for (int row = 0; row < rowNames.length; row++) {
            String rowName = rowNames[row];

            for (int col = 0; col < colNames.length; col++) {
                String colName = colNames[col];

                if (etg.get(rowName) != null && etg.get(rowName).equals(colName)) {
                    truePosValues[tp] = matrix[row][col];
                    tp++;
                } else {
                    trueNegValues[tn] = matrix[row][col];
                    tn++;
                }
            }

        }

        // create distribution for the best matches of this matrix
        double[] lowestVals = new double[numPositives];
        double[] restVals = new double[numTN];
        int tnProcessed = 0;
        for (int col = 0; col < colNames.length; col++) {
            double lowestForThisCol = Double.MAX_VALUE;
            int lowestRow = -1;
            for (int row = 0; row < rowNames.length; row++) {
                if (matrix[row][col] < lowestForThisCol) {
                    lowestRow = row;
                    lowestForThisCol = matrix[row][col];
                }
            }
            for (int row = 0; row < rowNames.length; row++) {
                if (row == lowestRow) {
                    lowestVals[col] = matrix[row][col];
                } else {
                    restVals[tnProcessed] = matrix[row][col];
                    tnProcessed++;
                }
            }
        }

        drawROCCurve(lowestVals, restVals, true, outputLoc);
        drawROCCurve(truePosValues, trueNegValues, false, outputLoc);

        // now plot the ROC
        textOutput.close();
    }

    public void plotHeatMap(DoubleMatrixDataset<String, String> result, String plotTitle, String subTitle, String xAxisLabel, String yAxisLabel, double[] expPC1EigenVector, double[] genPC1EigenVector, String outdir) {
        int numRows = result.rowObjects.size();
        int numCols = result.colObjects.size();
        MixupMapperHeatMap hm = new MixupMapperHeatMap((numCols * 10) + 50, (numRows * 10) + 50, true, outdir+"HeatMap.pdf");

        hm.setLabels(result.rowObjects.toArray(new String[0]), result.colObjects.toArray(new String[0]));
        hm.setAxisLabels(xAxisLabel, yAxisLabel);
        // hm.sortByLabels(true, true);
        hm.setPC1EigenVector(genPC1EigenVector, expPC1EigenVector);

        hm.plot(result.rawData);
        hm.setTitle(plotTitle);
        if(subTitle!=null){
            hm.setSubTitle(subTitle);
        }
        

        hm.draw(outdir+"HeatMap.pdf");
        hm.close();
    }

    private void drawDistributionsWithFamilyData(double[] truePosValues, double[] trueNegValues, double[] parentValues, String outputLoc) {
        int numBins = 100;

        double[] tmpVals = new double[3];
        tmpVals[0] = JSci.maths.ArrayMath.max(truePosValues);
        tmpVals[1] = JSci.maths.ArrayMath.max(trueNegValues);
        tmpVals[2] = JSci.maths.ArrayMath.max(parentValues);

        // get the max for the distributions
        double distmax = JSci.maths.ArrayMath.max(tmpVals);

        tmpVals[0] = JSci.maths.ArrayMath.min(truePosValues);
        tmpVals[1] = JSci.maths.ArrayMath.min(trueNegValues);
        tmpVals[2] = JSci.maths.ArrayMath.min(parentValues);
        // get the min for the distributions

        double distmin = JSci.maths.ArrayMath.min(tmpVals);

        distmin = Math.floor(distmin);
        distmax = Math.ceil(distmax);

        double mu0 = JSci.maths.ArrayMath.mean(trueNegValues);
        double mu1 = JSci.maths.ArrayMath.mean(truePosValues);
        double mu2 = JSci.maths.ArrayMath.mean(parentValues);

        double sd0 = JSci.maths.ArrayMath.standardDeviation(trueNegValues);

        double var0 = JSci.maths.ArrayMath.variance(trueNegValues);
        double sd1 = JSci.maths.ArrayMath.standardDeviation(truePosValues);
        double var1 = JSci.maths.ArrayMath.variance(truePosValues);

        double sd2 = JSci.maths.ArrayMath.standardDeviation(parentValues);
        double var2 = JSci.maths.ArrayMath.variance(parentValues);


        double absmu = Math.abs(mu1 - mu0);
        double absmutpvsp = Math.abs(mu1 - mu2);

        double abssd = ((sd0 + sd1) / 2);
        double SNR = absmu / abssd;


        double absvar = (var0 + var1) / 2;
        double SNR2 = absmu / absvar;

        double absvartpvsp = (var1 + var2) / 2;
        double SNR2tpvsp = absmutpvsp / absvartpvsp;

        double SNRTPVsParents = absmutpvsp / ((sd0 + sd2) / 2);


        double test = TTest.test(trueNegValues, truePosValues);
        System.out.println("The SNR/distance between the original sample assignment frequency distribution and other sample pair frequency distribution:\t" + SNR2);
        System.out.println("The SNR/distance between the true positive and 1st degree family distributions:\t" + SNR2tpvsp);

        // make a true positive distribution
        DiscreteDist tp = new DiscreteDist();
        tp.createDist(distmin, distmax, numBins);

        // make a true negative distribution
        DiscreteDist tn = new DiscreteDist();
        tn.createDist(distmin, distmax, numBins);

        DiscreteDist pv = new DiscreteDist();
        pv.createDist(distmin, distmax, numBins);

        java.util.Arrays.sort(truePosValues);
        java.util.Arrays.sort(trueNegValues);
        java.util.Arrays.sort(parentValues);

        double range = 0.0;
        if (distmin < 0) {
            range = distmax + Math.abs(distmin);
        } else {
            range = distmax - distmin;
        }

        double binInterval = range / numBins;

        // put the true positive values in to a bin --> convert to frequency distribution
        for (int i = 0; i < truePosValues.length; i++) {
            double val = 0.0;
            if (distmin < 0.0) {
                val = truePosValues[i] + Math.abs(distmin);
            } else {
                val = truePosValues[i] - distmin;
            }

            int binNum = (int) Math.floor(val / binInterval);

            if (binNum == numBins) {
                binNum = numBins - 1;
            }
            tp.addToBin(binNum);
            // comp.addToBin(binNum);
        }

        // put the true negative values in to a bin --> convert to frequency distribution
        for (int i = 0; i < trueNegValues.length; i++) {
            double val = 0.0;
            if (distmin < 0.0) {
                val = trueNegValues[i] + Math.abs(distmin);
            } else {
                val = trueNegValues[i] - distmin;
            }
            int binNum = (int) Math.floor(val / binInterval);
            if (binNum == numBins) {
                binNum = numBins - 1;
            }
            tn.addToBin(binNum);
            // comp.addToBin(binNum);

        }

        for (int i = 0; i < parentValues.length; i++) {
            double val = 0.0;
            if (distmin < 0.0) {
                val = parentValues[i] + Math.abs(distmin);
            } else {
                val = parentValues[i] - distmin;
            }

            int binNum = (int) Math.floor(val / binInterval);

            if (binNum == numBins) {
                binNum = numBins - 1;
            }
            pv.addToBin(binNum);
            // comp.addToBin(binNum);
        }

        // create the distribution objects
        DiscreteDist[] dists = new DiscreteDist[3];
        tp.convertToFrequency(truePosValues.length);
        tn.convertToFrequency(trueNegValues.length);
        pv.convertToFrequency(parentValues.length);

        dists[0] = tp;
        dists[2] = pv;
        dists[1] = tn;

        // plot the distribution histogram based on frequency
        MultiFrequencyDistributionHistogram mfh = new MultiFrequencyDistributionHistogram(1100, 1100, true, outputLoc + "Distribution-Frequency-WithRelatives.pdf");
        mfh.plot(dists);
        mfh.draw(outputLoc + "Distribution-Frequency-WithRelatives.png");
        mfh.close();
    }

    private void drawROCCurve(double[] truePosValues, double[] trueNegValues, boolean bestmatches, String outputLoc) throws IOException {
        if (trueNegValues.length == 0 || truePosValues.length == 0) {
            System.out.println("WARNING: cannot determine ROCs since there are no values for true positives or true negatives: (" + truePosValues.length + " - " + trueNegValues.length + ")");
        } else {
            int numBins = 100;

            double[] tmpVals = new double[2];
            tmpVals[0] = JSci.maths.ArrayMath.max(truePosValues);
            tmpVals[1] = JSci.maths.ArrayMath.max(trueNegValues);

            // get the max for the distributions
            double distmax = JSci.maths.ArrayMath.max(tmpVals);

            tmpVals[0] = JSci.maths.ArrayMath.min(truePosValues);
            tmpVals[1] = JSci.maths.ArrayMath.min(trueNegValues);

            // get the min for the distributions
            double distmin = JSci.maths.ArrayMath.min(tmpVals);

            distmin = Math.floor(distmin);
            distmax = Math.ceil(distmax);

            double mu0 = JSci.maths.ArrayMath.mean(trueNegValues);
            double mu1 = JSci.maths.ArrayMath.mean(truePosValues);

            double sd0 = JSci.maths.ArrayMath.standardDeviation(trueNegValues);
            double var0 = JSci.maths.ArrayMath.variance(trueNegValues);

            double sd1 = JSci.maths.ArrayMath.standardDeviation(truePosValues);
            double var1 = JSci.maths.ArrayMath.variance(truePosValues);

            double absmu = Math.abs(mu1 - mu0);
            double abssd = ((sd0 + sd1) / 2);

            double SNR = absmu / abssd;

            double absvar = (var0 + var1) / 2;
            double SNR2 = absmu / absvar;


            double test = TTest.test(trueNegValues, truePosValues);
            System.out.println("P-value of T-test between original sample assignments and other sample pairs: " + test);
            if (bestmatches) {
                System.out.println("The SNR/distance between the best match frequency distribution and other sample pair frequency distribution:\t" + SNR2);
            } else {
                System.out.println("The SNR/distance between the original sample assignment frequency distribution and other sample pair frequency distribution:\t" + SNR2);
            }


            // make a true positive distribution
            DiscreteDist tp = new DiscreteDist();
            tp.createDist(distmin, distmax, numBins);

            // make a true negative distribution
            DiscreteDist tn = new DiscreteDist();
            tn.createDist(distmin, distmax, numBins);

//	    System.out.println("Sorting TP");
            java.util.Arrays.sort(truePosValues);

//	    System.out.println("Sorting FP");
            java.util.Arrays.sort(trueNegValues);

            double range = 0.0;
            if (distmin < 0) {
                range = distmax + Math.abs(distmin);
            } else {
                range = distmax - distmin;
            }

            double binInterval = range / numBins;

            // put the true positive values in to a bin --> convert to frequency distribution
//	    System.out.println("Binning TP");
            for (int i = 0; i < truePosValues.length; i++) {
                double val = 0.0;
                if (distmin < 0.0) {
                    val = truePosValues[i] + Math.abs(distmin);
                } else {
                    val = truePosValues[i] - distmin;
                }

                int binNum = (int) Math.floor(val / binInterval);

                if (binNum == numBins) {
                    binNum = numBins - 1;
                }
                tp.addToBin(binNum);
                // comp.addToBin(binNum);
            }

            // put the true negative values in to a bin --> convert to frequency distribution
//	    System.out.println("Binning FP");
            for (int i = 0; i < trueNegValues.length; i++) {
                double val = 0.0;
                if (distmin < 0.0) {
                    val = trueNegValues[i] + Math.abs(distmin);
                } else {
                    val = trueNegValues[i] - distmin;
                }
                int binNum = (int) Math.floor(val / binInterval);
                if (binNum == numBins) {
                    binNum = numBins - 1;
                }
                tn.addToBin(binNum);
                // comp.addToBin(binNum);

            }

            // create the distribution objects
            DiscreteDist[] dists = new DiscreteDist[2];
            tp.convertToFrequency(truePosValues.length);
            tn.convertToFrequency(trueNegValues.length);

            dists[0] = tp;
            dists[1] = tn;

            // plot the distribution histogram based on frequency
            MultiFrequencyDistributionHistogram mfh = null;
            if (bestmatches) {
                mfh = new MultiFrequencyDistributionHistogram(1100, 1100, true, outputLoc + "Distribution-Frequency-BestMatches.pdf");
            } else {
                mfh = new MultiFrequencyDistributionHistogram(1100, 1100, true, outputLoc + "Distribution-Frequency.pdf");
            }

            mfh.plot(dists);

            if (bestmatches) {
                mfh.draw(outputLoc + "Distribution-Frequency-BestMatches.pdf");

            } else {
                mfh.draw(outputLoc + "Distribution-Frequency.pdf");

            }
            mfh.close();

            // calculate the cumulative value for each bin
            tp.calcCumulative();
            tn.calcCumulative();

            // convert to frequencies
            tp.convertToFrequency(truePosValues.length);
            tn.convertToFrequency(trueNegValues.length);

            // reset the bin iterator
            tp.resetIterator();
            tn.resetIterator();

            int numPos = truePosValues.length;
            int numNeg = trueNegValues.length;

            // create new true postive rate and true negative rate bins
            TPR = new double[6][numBins];
            FPR = new double[6][numBins];
            AUC = new double[6];

            // determine the ROC for a priori chance between 0 and 1.0 with increments of 0.05
            for (int i = 1; i < 6; i++) {
                tp.resetIterator();
                tn.resetIterator();

                double apriori = i * 0.20;

                double[] tmpTPR = new double[numBins];
                double[] tmpFPR = new double[numBins];


                double tmpAUC = 0;
                int currentBinNum = 0;

                double prevTPR = 0.0; // y-axis
                double prevFPR = 0.0; // x-axis
                double curTPR = 0.0;
                double curFPR = 0.0;

                while (tp.hasNext()) {
                    Bin tpBin = tp.getNext();
                    Bin tnBin = tn.getNext();

                    double TP = tpBin.getCumulativeFrequency();
                    double FN = 1 - TP;

                    double FP = tnBin.getCumulativeFrequency();
                    double TN = 1 - FP;

                    curTPR = (TP * apriori) / (TP * apriori + FN);
                    // curTPR = (TP * apriori) / (TP * apriori + FN);
                    curFPR = FP / (FP + TN);

                    tmpTPR[currentBinNum] = curTPR;
                    tmpFPR[currentBinNum] = curFPR;



                    if (currentBinNum != 0) {
                        prevTPR = tmpTPR[currentBinNum - 1]; // y-axis
                        prevFPR = tmpFPR[currentBinNum - 1]; // x-axis
                        double squareArea = (curFPR - prevFPR) * (prevTPR); // square area
                        double triangleArea = (curFPR - prevFPR) * (curTPR - prevTPR) * 0.5;
                        tmpAUC += (squareArea + triangleArea);
                    }



                    textOutput.write(i + "\t" + apriori + "\t" + tmpTPR[currentBinNum] + "\t" + tmpFPR[currentBinNum]);
                    currentBinNum++;
                }

                // put the values for this a priori chance in the array
                FPR[i] = tmpFPR;
                TPR[i] = tmpTPR;
                AUC[i] = tmpAUC;
            }

            // plot curve based on the a priori chances
            if (!bestmatches) {
                ROCCurve roc = new ROCCurve(1024, 1024, true, outputLoc + "ROC.pdf");
                roc.plot(TPR, FPR, AUC, 0.20);
                roc.draw(outputLoc + "ROC.pdf");
                roc.close();
            }
        }
    }

    public void loadTruePositives(HashMap<String, String> genotypeToExpression, HashMap<String, String> expressionToGenotype) {
        etg = expressionToGenotype;
        gte = genotypeToExpression;
    }

    public void loadTruePositives(String inFile) throws IOException {
        gte = new HashMap<String, String>();
        etg = new HashMap<String, String>();

        TextFile in = new TextFile(inFile, TextFile.R);
        String line = "";
        while ((line = in.readLine()) != null) {
            String[] lineElems = line.split("\t");
            if (lineElems.length > 1) {
                gte.put(lineElems[0], lineElems[1]);
                etg.put(lineElems[1], lineElems[0]);
            }
        }
        in.close();

    }

    public void loadFamilyData(HashMap<String, Trio> sampleToTrio) {
        this.sampleToTrio = sampleToTrio;

        this.loadedFamilyData = true;
    }

    private boolean isRelative(String colName, String coupledGenotype, ArrayList<String> relatives) {
        boolean isRelative = false;
        for (int i = 0; i < relatives.size(); i++) {
            String relative = relatives.get(i);
            if (relative.equals(colName) && !relative.equals(coupledGenotype)) {
                isRelative = true;
                break;
            }
        }

        return isRelative;
    }
}
