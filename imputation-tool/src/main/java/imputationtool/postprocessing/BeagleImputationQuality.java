/*
 * This program checks for batch effects within beagle imputed dataset batches
 */
package imputationtool.postprocessing;

import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

/**
 *
 * @author harmjan
 */
public class BeagleImputationQuality {

    private static Pattern tab = Pattern.compile("\t");

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

    public void determineImputationQualityDistribution(String inputDir, String template, int numBatches, String outputLocation) throws IOException {

        String[] batchNames = getBatches(numBatches);
        int nrBatches = batchNames.length;
        int chrStart = 1;
        int chrEnd = 22;
        TextFile logGlobal = new TextFile(outputLocation + "/r2-distribution-global.txt", TextFile.W);

        int[] r2globalfreqdistribution = new int[11];

        boolean allFilesAvailable = true;
        int filesnotfound = 0;
        for (int batch = 0; batch < nrBatches; batch++) {
            int nrSamplesThisBatch = 0;
            for (int chr = chrStart; chr <= chrEnd; chr++) {
                // for(int chr: chromosomes){
                String templatecopy = new String(template);
                templatecopy = templatecopy.replace("BATCH", batchNames[batch]);
                templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);

                // Finn-CelGWAS2_Chr2-HM2-4.Finn-CelGWAS2_Chr2-4-BEAGLE.gprobs
                String fileName = inputDir + "/" + templatecopy + ".r2";

                File file = new File(fileName);
                if (!file.canRead()) {
                    System.out.println("Cannot open file:\t" + fileName);
                    allFilesAvailable = false;
                    filesnotfound++;
                }
            }
        }

        if (!allFilesAvailable) {
            System.out.println("Not all imputed dosage files are available!!! (" + filesnotfound + " out of " + (nrBatches * 22) + "). Exiting...");
            System.exit(-1);
        }

        int numsnpsglobal = 0;
        for (int chr = 1; chr < 23; chr++) {
            TextFile logChr = new TextFile(outputLocation + "/r2-distribution-" + chr + ".txt", TextFile.W);
            TextFile logWilChr = new TextFile(outputLocation + "/wilcoxon-" + chr + ".txt", TextFile.W);
            int numsnpschr = 0;
            int[] r2chrfreqdistribution = new int[11];

            for (int batch = 0; batch < batchNames.length; batch++) {
                TextFile logBatch = new TextFile(outputLocation + "/r2-distribution-" + chr + "-" + batchNames[batch] + ".txt", TextFile.W);
                int numsnpsbatch = 0;
                int[] r2batchfreqdistribution = new int[11];
                String templatecopy = new String(template);
                templatecopy = templatecopy.replace("BATCH", batchNames[batch]);
                templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);
                String fileName = inputDir + "/" + templatecopy + ".r2";




                TextFile reader = new TextFile(fileName, TextFile.R);
                String ln = "";

                while ((ln = reader.readLine()) != null) {
                    String[] elems = tab.split(ln);
                    if (elems.length == 2) {

                        double correlation = Double.parseDouble(elems[1]);
                        double absCor = Math.abs(correlation);

                        int binNumber = (int) (absCor * 10d);

                        if (Double.isNaN(absCor)) {
                            binNumber = 0;
                        }



                        r2globalfreqdistribution[binNumber]++;
                        r2chrfreqdistribution[binNumber]++;
                        r2batchfreqdistribution[binNumber]++;

                        numsnpsglobal++;
                        numsnpschr++;
                        numsnpsbatch++;

                    }
                }

                reader.close();


                System.out.println("Chr" + chr + "-" + batchNames[batch]);
                for (int i = 0; i < r2batchfreqdistribution.length; i++) {
                    // System.out.println( (i/10d) +"\t"+r2batchfreqdistribution[i] / numsnpsbatch);
                    logBatch.writeln((i / 10d) + "\t" + (double) r2batchfreqdistribution[i] / numsnpsbatch);
                }
                System.out.println("");

                logBatch.close();
            }


            // repeat batches, to compare against chromosome distribution...

            double[] r2chrfreqdistribution_d = new double[11];
            for (int i = 0; i < 11; i++) {
                r2chrfreqdistribution_d[i] = (double) r2chrfreqdistribution[i] / (double) numsnpschr;
            }

            for (int batch = 0; batch < batchNames.length; batch++) {

                int numsnpsbatch = 0;
                int[] r2batchfreqdistribution = new int[11];
                String templatecopy = new String(template);
                templatecopy = templatecopy.replace("BATCH", batchNames[batch]);
                templatecopy = templatecopy.replace("CHROMOSOME", "" + chr);
                String fileName = inputDir + "/" + templatecopy + ".r2";
                TextFile reader = new TextFile(fileName, TextFile.R);
                String ln = "";

                while ((ln = reader.readLine()) != null) {

                    String[] elems = tab.split(ln);
                    if (elems.length == 2) {


                        double correlation = Double.parseDouble(elems[1]);

                        double absCor = Math.abs(correlation);

                        int binNumber = (int) (absCor * 10d);
                        r2batchfreqdistribution[binNumber]++;
                        numsnpsbatch++;
                    }
                }

                reader.close();


                double[] r2batchfreqdistribution_d = new double[11];
                for (int i = 0; i < 11; i++) {
                    r2batchfreqdistribution_d[i] = (double) r2batchfreqdistribution[i] / (double) numsnpsbatch;
                }

                double pvalue = wilCoxon(r2batchfreqdistribution_d, r2chrfreqdistribution_d);

                //System.out.println("Chr"+chr+"-"+batchNames[batch]+"\t"+pvalue);

                logWilChr.writeln("Chr" + chr + "-" + batchNames[batch] + "\t" + pvalue);
            }

            System.out.println("Chr" + chr);
            for (int i = 0; i < r2chrfreqdistribution.length; i++) {
                // System.out.println( (i/10d) +"\t"+r2chrfreqdistribution[i] / numsnpschr);
                logChr.writeln((i / 10d) + "\t" + (double) r2chrfreqdistribution[i] / numsnpschr);
            }
            System.out.println("");
            logWilChr.close();
            logChr.close();
        }

        System.out.println("Global:");
        for (int i = 0; i < r2globalfreqdistribution.length; i++) {
            // System.out.println( (i/10d) +"\t"+r2globalfreqdistribution[i] / numsnpsglobal);
            logGlobal.writeln((i / 10d) + "\t" + (double) r2globalfreqdistribution[i] / numsnpsglobal);
        }
        System.out.println("");
        logGlobal.close();
    }

    private double wilCoxon(double[] vals1, double[] vals2) {
        WilcoxonMannWhitney wmw = new WilcoxonMannWhitney();
        double pValue = wmw.returnWilcoxonMannWhitneyPValue(vals1, vals2);
        return pValue;
    }
}
