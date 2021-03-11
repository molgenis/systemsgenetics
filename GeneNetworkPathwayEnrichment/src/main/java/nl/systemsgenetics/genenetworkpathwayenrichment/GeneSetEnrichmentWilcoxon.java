package nl.systemsgenetics.genenetworkpathwayenrichment;


import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.IntStream;

public class GeneSetEnrichmentWilcoxon {


    boolean includeForegroundGenesInBackgroundSet = true;

    public void runProperBackgroundCorrection(String matrix,
                                              String genesetfile,
                                              String annotationFile,
                                              String ensgToHugoFile,
                                              String limitGenes,
                                              String meanExpFile,
                                              int permutations,
                                              boolean squareValues,
                                              String outputfile) throws IOException {

        System.out.println("Gene set enrichment with proper bg correction.");
        // limit analysis to these genes
        HashSet<String> potentialBackgroundSet = null;
        if (limitGenes != null) {
            System.out.println("Limiting genes to those specified in " + limitGenes);
            potentialBackgroundSet = readSet(limitGenes);
            System.out.println("Limiting genes to " + potentialBackgroundSet.size());
        }

        // read geneset to test
        HashSet<String> potentialForegroundGeneset = readSet(genesetfile);
        System.out.println(potentialForegroundGeneset.size() + " genes in foreground geneset");
        potentialBackgroundSet.addAll(potentialForegroundGeneset);
        System.out.println(potentialBackgroundSet.size() + " genes in background geneset");

        // read average gene expression values from file
        ArrayList<ArrayList<String>> geneBins = null;
        HashMap<String, Integer> geneToBin = null;
        if (meanExpFile != null) {
            System.out.println("Reading expression data: " + meanExpFile);
            Pair<ArrayList<ArrayList<String>>, HashMap<String, Integer>> data = readAndBin(meanExpFile, potentialBackgroundSet);
            geneBins = data.getLeft();
            geneToBin = data.getRight();

            // limit to genes for which we have expression data
            HashSet<String> newGeneLimit = new HashSet<>();
            newGeneLimit.addAll(geneToBin.keySet());
            potentialBackgroundSet = newGeneLimit;
            System.out.println("Expression data read for: " + geneToBin.size() + " genes. Limiting analysis to this set of genes.");
        }

        // map ensembl ids to Hugo IDs
        HashMap<String, String> ensToHugo = null;
        if (ensgToHugoFile != null) {
            System.out.println("Reading Ensembl to HUGO annotation: " + ensgToHugoFile);
            ensToHugo = new HashMap<>();
            TextFile tf = new TextFile(ensgToHugoFile, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length > 1) {
                    String ensg = elems[0].split("\\.")[0];
                    ensToHugo.put(ensg, elems[1]);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            System.out.println(ensToHugo.size() + " gene annotations read.");
        }

        // Read pathway annotations
        HashMap<String, String> idToAnnotation = null;
        if (annotationFile != null) {
            System.out.println("Reading pathway annotation: " + annotationFile);
            idToAnnotation = new HashMap<String, String>();
            TextFile tf2 = new TextFile(annotationFile, TextFile.R);
            String[] elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length > 1) {
                    idToAnnotation.put(elems[0], elems[1]);
                }
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            System.out.println(idToAnnotation.size() + " pathway annotations read");
        }


        WilcoxonMannWhitney w = new WilcoxonMannWhitney();

        ArrayList<Pair<Double, Result>> results = new ArrayList<>();
        System.out.println("Using: " + matrix);
        // use a text based matrix
        System.out.println("Will square values: " + squareValues);
        System.out.println("Merging foreground genes into background: " + includeForegroundGenesInBackgroundSet);
        if (matrix.endsWith(".txt.gz") || matrix.endsWith(".txt")) {
            TextFile tf3 = new TextFile(matrix, TextFile.R);
            String[] header = tf3.readLineElems(TextFile.tab);
            ArrayList<String> allGenesInMatrix = new ArrayList<>();
            for (int i = 1; i < header.length; i++) {
                // remove ENSG version id
                String[] elems = header[i].split("\\.");
                if (elems.length > 1) {
                    allGenesInMatrix.add(elems[0]);
                } else {
                    allGenesInMatrix.add(header[i]);
                }
            }

            Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> genesets = getGeneIDs(allGenesInMatrix, potentialForegroundGeneset, potentialBackgroundSet, true);
            ArrayList<Integer> genesForeground = genesets.getLeft();
            ArrayList<Integer> genesBackground = genesets.getMiddle();
            ArrayList<String> presentForegroundGenes = genesets.getRight();
            ArrayList<Integer> allAllowedGeneIds = new ArrayList<>();
            allAllowedGeneIds.addAll(genesForeground);
            allAllowedGeneIds.addAll(genesBackground);

            String[] data = tf3.readLineElems(TextFile.tab);
            double[] foregroundGeneVals = new double[genesForeground.size()];
            double[] backgroundGeneVals = new double[genesBackground.size()];
            double[] pvalsperm = new double[permutations];
            double[] fcsPerm = new double[permutations];

            int lnctr = 0;
            while (data != null) {
                double[] rowOrig = convertToDouble(data);
                double[] row = null;
                if (squareValues) {
                    rowOrig = square(rowOrig);
                    row = rowOrig; // square(rowOrig);
                } else {
                    row = addMin(rowOrig);
                }


                copy(row, foregroundGeneVals, backgroundGeneVals, genesForeground, genesBackground);
                Pair<ArrayList<String>, ArrayList<String>> genesplit = split2(rowOrig, genesForeground, allGenesInMatrix, ensToHugo);

                double[] fgNoNaN = stripNaN(foregroundGeneVals);
                double[] bgNoNaN = stripNaN(backgroundGeneVals);

                double pvalWilcoxon = w.returnWilcoxonMannWhitneyPValue(fgNoNaN, bgNoNaN);
                double auc = w.getAUC();
                double meanFg = Descriptives.mean(fgNoNaN);
                double meanBg = Descriptives.mean(bgNoNaN);
                double fc = meanFg / meanBg;
                Result resultObj = new Result();
                String rowname = data[0];
                resultObj.rowname = rowname;
                resultObj.auc = auc;
                resultObj.meanFg = meanFg;
                resultObj.meanBg = meanBg;
                resultObj.fc = fc;
                resultObj.posAndNegGenes = genesplit;


                // permute
                ArrayList<ArrayList<String>> finalGeneBins1 = geneBins;
                HashMap<String, Integer> finalGeneToBin1 = geneToBin;
                HashSet<String> finalGeneLimit1 = potentialBackgroundSet;
                double[] finalRow = row;
                IntStream.range(0, permutations).parallel().forEach(p -> {

                    // match genes based on expression
                    ArrayList<Integer> permutedForegroundGenes = selectBackgroundFromBins(allGenesInMatrix, genesForeground, allAllowedGeneIds, finalGeneBins1, finalGeneToBin1);
                    HashSet<String> permFgGenes = new HashSet<>();
                    for (Integer i : permutedForegroundGenes) {
                        permFgGenes.add(allGenesInMatrix.get(i));
                    }
                    Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> genesetPerm = getGeneIDs(allGenesInMatrix, permFgGenes, finalGeneLimit1, false);

                    // select a new foreground

                    WilcoxonMannWhitney wp = new WilcoxonMannWhitney();
                    // select a new background
                    ArrayList<Integer> genesBackgroundPerm = genesetPerm.getMiddle();
                    double[] foregroundGeneValsPerm = new double[permutedForegroundGenes.size()];
                    double[] backgroundGeneValsPerm = new double[genesBackgroundPerm.size()];
                    copy(finalRow, foregroundGeneValsPerm, backgroundGeneValsPerm, permutedForegroundGenes, genesBackgroundPerm);

                    foregroundGeneValsPerm = stripNaN(foregroundGeneValsPerm);
                    backgroundGeneValsPerm = stripNaN(backgroundGeneValsPerm);

                    double pvalPerm = wp.returnWilcoxonMannWhitneyPValue(foregroundGeneValsPerm, backgroundGeneValsPerm);
                    double meanFgPerm = Descriptives.mean(foregroundGeneValsPerm);
                    double meanBgPerm = Descriptives.mean(backgroundGeneValsPerm);
                    double fcPerm = meanFgPerm / meanBgPerm;
                    pvalsperm[p] = pvalPerm;
                    fcsPerm[p] = fcPerm;
                });

                // determine permutation p-value
                Arrays.sort(pvalsperm);
                int nrLower = 0;
                for (int p = 0; p < pvalsperm.length; p++) {
                    if (pvalsperm[p] < pvalWilcoxon) {
                        nrLower++;
                    }
                }
                double meanFc = Descriptives.mean(fcsPerm);
                double sdFc = Math.sqrt(Descriptives.variance(fcsPerm));
                double fcZ = (resultObj.fc - meanFc) / sdFc;
                double empFDR = (double) nrLower / pvalsperm.length;

                resultObj.fcZ = fcZ;
                resultObj.fcZPval = ZScores.zToP(fcZ);
                resultObj.empFDR = empFDR;
                resultObj.pvalWilcoxon = pvalWilcoxon;

                results.add(new Pair<Double, Result>(pvalWilcoxon, resultObj, Pair.SORTBY.LEFT));
                data = tf3.readLineElems(TextFile.tab);
                lnctr++;
                if (lnctr % 10 == 0) {
                    System.out.print(lnctr + " terms processed\r");
                }
            }
            System.out.println("");
            System.out.println("Done");
            tf3.close();
            TextFile outputbg = new TextFile(outputfile + "-backgroundgenes.txt", TextFile.W);
            for (Integer i : genesBackground) {
                outputbg.writeln(allGenesInMatrix.get(i));
            }
            outputbg.close();
        }


        // sort results by p-value
        Collections.sort(results);
        double[] qvals = fdr(results, false);
        ArrayList<Pair<Double, Result>> results2 = new ArrayList<>();
        for (Pair<Double, Result> result : results) {
            results2.add(new Pair<Double, Result>(result.getRight().fcZPval, result.getRight()));
        }
        Collections.sort(results2);
        double[] qvalsFcP = fdr(results2, true);
        int q = 0;
        for (Pair<Double, Result> result : results2) {
            result.getRight().fcZPvalFDR = qvalsFcP[q];
            q++;
        }

        System.out.println("Output: " + outputfile);


        TextFile output = new TextFile(outputfile, TextFile.W);
        output.writeln("Pval\tFDR\tEmpiricalFDR\tTerm\tAnnotation\tAUC\tMeanFg\tMeanBg\tFoldChange\tFoldchangeZscore\tFoldChangePval\tFoldChangePvalFDR\t#PosGenes\t#NegGenes\tPosGenes\tNegGenes");
        int rank = 0;
        int sig = 0;

        // FDR/qval adjustment


        rank = 0;
        // write output
        for (Pair<Double, Result> result : results) {
            String id = result.getRight().rowname;
            double pval = result.getRight().pvalWilcoxon;
            double auc = result.getRight().auc;
            double meanfg = result.getRight().meanFg;
            double meanbg = result.getRight().meanBg;
            double foldchange = result.getRight().fc;
            double foldchangepval = result.getRight().fcZPval;
            double foldchangepvalfdr = result.getRight().fcZPvalFDR;
            String annotation = id;

            Pair<ArrayList<String>, ArrayList<String>> genes = result.getRight().posAndNegGenes;

            if (idToAnnotation != null) {
                annotation = idToAnnotation.get(id);
                if (annotation == null) {
                    annotation = id;
                }
            }
            double qvalue = qvals[rank];
            String posgenes = "";
            String neggenes = "";
            if (genes.getLeft().size() > 0) {
                posgenes = Strings.concat(genes.getLeft(), Strings.semicolon);
            }
            if (genes.getRight().size() > 0) {
                neggenes = Strings.concat(genes.getRight(), Strings.semicolon);
            }
            String outln = result.getLeft()
                    + "\t" + qvalue
                    + "\t" + result.getRight().empFDR
                    + "\t" + result.getRight().rowname
                    + "\t" + annotation
                    + "\t" + auc
                    + "\t" + meanfg
                    + "\t" + meanbg
                    + "\t" + foldchange
                    + "\t" + result.getRight().fcZ
                    + "\t" + foldchangepval
                    + "\t" + foldchangepvalfdr
                    + "\t" + genes.getLeft().size()
                    + "\t" + genes.getRight().size()
                    + "\t" + posgenes + "\t" + neggenes;

            if (qvalue < 0.05) {
//                System.out.println(outln);
                sig++;
            }
            rank++;
            output.writeln(outln);
        }

        output.close();
        System.out.println(sig + "/" + results.size() + " significant.");
        System.out.println();

    }

    private double[] fdr(ArrayList<Pair<Double, Result>> results, boolean fcPVal) {
        int rank = 0;
        double[] qvals = new double[results.size()];
        for (Pair<Double, Result> result : results) {
            if (fcPVal) {
                qvals[rank] = ((double) results.size() / (rank + 1)) * result.getRight().fcZPval;
            } else {
                qvals[rank] = ((double) results.size() / (rank + 1)) * result.getRight().pvalWilcoxon;
            }
            rank++;
        }
        for (int i = results.size() - 2; i > -1; i--) {
            qvals[i] = Math.min(qvals[i], qvals[i + 1]);
        }


        return qvals;
    }

    private Pair<ArrayList<String>, ArrayList<String>> split2(double[] row, ArrayList<Integer> genesForeground, ArrayList<String> allGenesInMatrix, HashMap<String, String> ensToHugo) {
        ArrayList<Pair<String, Double>> negP = new ArrayList<>();
        ArrayList<Pair<String, Double>> posP = new ArrayList<>();


        DecimalFormat f = new DecimalFormat("#.##");
        for (int d = 0; d < genesForeground.size(); d++) {
            Integer id = genesForeground.get(d);
            double v = row[id];
            String gene = allGenesInMatrix.get(id);
            String hugo = ensToHugo.get(gene);
            if (hugo == null) {
                hugo = gene;
            }
            if (v >= 0) {
                posP.add(new Pair<>(hugo, v, Pair.SORTBY.RIGHT));
//                pos.add(hugo + ":" + f.format(v));
            } else {
                negP.add(new Pair<>(hugo, v, Pair.SORTBY.RIGHT));
            }
        }

        Collections.sort(negP);
        Collections.sort(posP, Collections.reverseOrder());
        ArrayList<String> pos = new ArrayList<>();
        ArrayList<String> neg = new ArrayList<>();

        for (Pair<String, Double> p : negP) {
            neg.add(p.getLeft() + ":" + f.format(p.getRight()));
        }
        for (Pair<String, Double> p : posP) {
            pos.add(p.getLeft() + ":" + f.format(p.getRight()));
        }


        return new Pair<>(pos, neg);


    }

    class Result {
        public Pair<ArrayList<String>, ArrayList<String>> posAndNegGenes;
        public double fcZ = 0;
        public double empFDR = 1;
        public double fcZPval;
        public double pvalWilcoxon;
        public double fcZPvalFDR;
        double meanFg;
        double meanBg;
        double fc;
        String rowname = null;
        double auc = 0;
    }

//    public void run(String matrix,
//                    String genesetfile,
//                    String annotationFile,
//                    String ensgToHugoFile,
//                    String limitGenes,
//                    String meanExpFile,
//                    boolean randomlyMatchNumberOfGenes,
//                    boolean permute,
//                    int permutations,
//                    String outputfile) throws IOException {
//
//        // limit analysis to these genes
//        HashSet<String> allPotentialBackgroundGenesHash = null;
//        if (limitGenes != null) {
//            System.out.println("Limiting genes to those specified in " + limitGenes);
//            allPotentialBackgroundGenesHash = readSet(limitGenes);
//            System.out.println("Limiting genes to " + allPotentialBackgroundGenesHash.size());
//        }
//
//        // read geneset to test
//        HashSet<String> potentialForegroundGeneset = readSet(genesetfile);
//        System.out.println(potentialForegroundGeneset.size() + " genes in foreground geneset");
//        System.out.println(allPotentialBackgroundGenesHash.size() + " genes in background geneset");
//        allPotentialBackgroundGenesHash.addAll(potentialForegroundGeneset);
//        System.out.println(allPotentialBackgroundGenesHash.size() + " total genes");
//
//
//        // read average gene expression values from file
//        ArrayList<ArrayList<String>> geneBins = null;
//        HashMap<String, Integer> geneToBin = null;
//        if (meanExpFile != null) {
//            System.out.println("Reading expression data: " + meanExpFile);
//            Pair<ArrayList<ArrayList<String>>, HashMap<String, Integer>> data = readAndBin(meanExpFile, allPotentialBackgroundGenesHash);
//            geneBins = data.getLeft();
//            geneToBin = data.getRight();
//
//            // limit to genes for which we have expression data
//            HashSet<String> newGeneLimit = new HashSet<>();
//            newGeneLimit.addAll(geneToBin.keySet());
//            allPotentialBackgroundGenesHash = newGeneLimit;
//            System.out.println("Expression data read for: " + geneToBin.size() + " genes. Limiting analysis to this set of genes.");
//        }
//
//        // map ensembl ids to Hugo IDs
//        HashMap<String, String> ensToHugo = null;
//        if (ensgToHugoFile != null) {
//            System.out.println("Reading Ensembl to HUGO annotation: " + ensgToHugoFile);
//            ensToHugo = new HashMap<>();
//            TextFile tf = new TextFile(ensgToHugoFile, TextFile.R);
//            String[] elems = tf.readLineElems(TextFile.tab);
//            while (elems != null) {
//                if (elems.length > 1) {
//                    String ensg = elems[0].split("\\.")[0];
//                    ensToHugo.put(ensg, elems[1]);
//                }
//                elems = tf.readLineElems(TextFile.tab);
//            }
//            tf.close();
//            System.out.println(ensToHugo.size() + " gene annotations read.");
//        }
//
//        // Read pathway annotations
//        HashMap<String, String> idToAnnotation = null;
//        if (annotationFile != null) {
//            System.out.println("Reading pathway annotation: " + annotationFile);
//            idToAnnotation = new HashMap<String, String>();
//            TextFile tf2 = new TextFile(annotationFile, TextFile.R);
//            String[] elems = tf2.readLineElems(TextFile.tab);
//            while (elems != null) {
//                if (elems.length > 1) {
//                    idToAnnotation.put(elems[0], elems[1]);
//                }
//                elems = tf2.readLineElems(TextFile.tab);
//            }
//            tf2.close();
//            System.out.println(idToAnnotation.size() + " pathway annotations read");
//        }
//
//
//        WilcoxonMannWhitney w = new WilcoxonMannWhitney();
//
//        ArrayList<Pair<Double, Result>> results = new ArrayList<>();
//        System.out.println("Using: " + matrix);
//        // use a text based matrix
//        if (matrix.endsWith(".txt.gz") || matrix.endsWith(".txt")) {
//            TextFile tf3 = new TextFile(matrix, TextFile.R);
//            String[] header = tf3.readLineElems(TextFile.tab);
//            ArrayList<String> allGenesInMatrix = new ArrayList<>();
//            for (int i = 1; i < header.length; i++) {
//                allGenesInMatrix.add(header[i]);
//            }
//
//            Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> genesets = selectForeGroundAndBackground(allGenesInMatrix, potentialForegroundGeneset, allPotentialBackgroundGenesHash, geneBins, geneToBin, randomlyMatchNumberOfGenes, true);
//            ArrayList<Integer> genesForeground = genesets.getLeft();
//            ArrayList<Integer> genesBackground = genesets.getMiddle();
//            ArrayList<String> presentForegroundGenes = genesets.getRight();
//            String[] data = tf3.readLineElems(TextFile.tab);
//            double[] foregroundGeneVals = new double[genesForeground.size()];
//            double[] backgroundGeneVals = new double[genesBackground.size()];
//
//            System.out.println(genesBackground.size() + " genes in background geneset");
//
//            int lnctr = 0;
//            while (data != null) {
//                double[] row = convertToDouble(data);
//                copy(row, foregroundGeneVals, backgroundGeneVals, genesForeground, genesBackground);
//                Pair<ArrayList<String>, ArrayList<String>> genesplit = split(foregroundGeneVals, presentForegroundGenes, ensToHugo);
//
//                double[] noNaNforegroundGeneVals = stripNaN(foregroundGeneVals);
//                double[] noNaNbackgroundGeneVals = stripNaN(backgroundGeneVals);
//
//                double pval = w.returnWilcoxonMannWhitneyPValue(noNaNforegroundGeneVals, noNaNbackgroundGeneVals);
//                double auc = w.getAUC();
//                double meanFg = Descriptives.mean(noNaNforegroundGeneVals);
//                double meanBg = Descriptives.mean(noNaNbackgroundGeneVals);
//                double fc = meanFg / meanBg;
//                Result r = new Result();
//                String rowname = data[0];
//                r.rowname = rowname;
//                r.auc = auc;
//                r.meanFg = meanFg;
//                r.meanBg = meanBg;
//                r.fc = fc;
//                r.posAndNegGenes = genesplit;
//
//                if (permute) {
//                    double[] pvalsperm = new double[permutations];
//                    double[] fcsPerm = new double[permutations];
//                    HashMap<String, Integer> finalGeneToBin = geneToBin;
//                    ArrayList<ArrayList<String>> finalGeneBins = geneBins;
//                    HashSet<String> finalGeneLimit = allPotentialBackgroundGenesHash;
//                    IntStream.range(0, permutations).parallel().forEach(p -> {
//
//                        // select a new foreground
//                        Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> genesetPerm = selectForeGroundAndBackground(allGenesInMatrix, potentialForegroundGeneset, finalGeneLimit, finalGeneBins, finalGeneToBin, randomlyMatchNumberOfGenes, false);
//                        ArrayList<Integer> genesForegroundPerm = genesetPerm.getMiddle();
//                        HashSet<String> genesForegroundPermNames = new HashSet<>();
//                        for (Integer i : genesForegroundPerm) {
//                            genesForegroundPermNames.add(allGenesInMatrix.get(i));
//                        }
//
//                        WilcoxonMannWhitney wp = new WilcoxonMannWhitney();
//                        // select a new background
//                        genesetPerm = selectForeGroundAndBackground(allGenesInMatrix, genesForegroundPermNames, finalGeneLimit, finalGeneBins, finalGeneToBin, randomlyMatchNumberOfGenes, false);
//                        ArrayList<Integer> genesBackgroundPerm = genesetPerm.getMiddle();
//                        double[] foregroundGeneValsPerm = new double[genesForegroundPerm.size()];
//                        double[] backgroundGeneValsPerm = new double[genesBackgroundPerm.size()];
//                        copy(row, foregroundGeneValsPerm, backgroundGeneValsPerm, genesForegroundPerm, genesBackgroundPerm);
//
//                        foregroundGeneValsPerm = stripNaN(foregroundGeneValsPerm);
//                        backgroundGeneValsPerm = stripNaN(backgroundGeneValsPerm);
//
//                        double pvalPerm = wp.returnWilcoxonMannWhitneyPValue(foregroundGeneValsPerm, backgroundGeneValsPerm);
//                        double meanFgPerm = Descriptives.mean(foregroundGeneValsPerm);
//                        double meanBgPerm = Descriptives.mean(backgroundGeneValsPerm);
//                        double fcPerm = meanFgPerm / meanBgPerm;
//                        pvalsperm[p] = pvalPerm;
//                        fcsPerm[p] = fcPerm;
//                    });
//
//                    // determine permutation p-value
//                    Arrays.sort(pvalsperm);
//                    int nrLower = 0;
//                    for (int p = 0; p < pvalsperm.length; p++) {
//                        if (pvalsperm[p] < pval) {
//                            nrLower++;
//                        }
//                    }
//                    double meanFc = Descriptives.mean(fcsPerm);
//                    double sdFc = Math.sqrt(Descriptives.variance(fcsPerm));
//                    double fcZ = (r.fc - meanFc) / sdFc;
//                    double empFDR = (double) nrLower / pvalsperm.length;
//
////                    System.out.println(meanFc);
////                    System.out.println(sdFc);
////                    System.out.println(nrLower);
////                    System.out.println(empFDR);
////                    System.out.println(fcZ);
////                    System.exit(0);
//
//
//                    r.fcZ = fcZ;
//                    r.empFDR = empFDR;
//
//                }
//
//                results.add(new Pair<Double, Result>(pval, r, Pair.SORTBY.LEFT));
//                data = tf3.readLineElems(TextFile.tab);
//                lnctr++;
//                if (lnctr % 10 == 0) {
//                    System.out.print(lnctr + " terms processed\r");
//                }
//            }
//            System.out.println("");
//            System.out.println("Done");
//            tf3.close();
//            TextFile outputbg = new TextFile(outputfile + "-backgroundgenes.txt", TextFile.W);
//            for (Integer i : genesBackground) {
//                outputbg.writeln(allGenesInMatrix.get(i));
//            }
//            outputbg.close();
//        } else { // or a binary matrix
//            DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(matrix);
//            Set<String> cols = it.getCols();
//            ArrayList<String> allGenesInMatrix = new ArrayList<String>();
//            allGenesInMatrix.addAll(cols);
//            Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> genesets = selectForeGroundAndBackground(allGenesInMatrix, potentialForegroundGeneset, allPotentialBackgroundGenesHash, geneBins, geneToBin, randomlyMatchNumberOfGenes, true);
//            ArrayList<Integer> genesForeground = genesets.getLeft();
//            ArrayList<Integer> genesBackground = genesets.getMiddle();
//            ArrayList<String> presentForegroundGenes = genesets.getRight();
//            double[] foregroundGeneVals = new double[genesForeground.size()];
//            double[] backgroundGeneVals = new double[genesBackground.size()];
//
//            ArrayList<String> rowIds = new ArrayList<String>();
//            rowIds.addAll(it.getRows());
//            int rctr = 0;
//            for (double[] row : it) {
//                copy(row, foregroundGeneVals, backgroundGeneVals, genesForeground, genesBackground);
//
//                Pair<ArrayList<String>, ArrayList<String>> genesplit = split(foregroundGeneVals, presentForegroundGenes, ensToHugo);
//                double pval = w.returnWilcoxonMannWhitneyPValue(foregroundGeneVals, backgroundGeneVals);
//                double auc = w.getAUC();
//                double meanFg = Descriptives.mean(foregroundGeneVals);
//                double meanBg = Descriptives.mean(backgroundGeneVals);
//                double fc = meanFg / meanBg;
//                Result r = new Result();
//                String rowname = rowIds.get(rctr);
//                r.rowname = rowname;
//                r.auc = auc;
//                r.meanFg = meanFg;
//                r.meanBg = meanBg;
//                r.fc = fc;
//                r.posAndNegGenes = genesplit;
//
//                rctr++;
//
//                results.add(new Pair<Double, Result>(pval, r, Pair.SORTBY.LEFT));
//            }
//            it.close();
//
//            TextFile outputbg = new TextFile(outputfile + "-backgroundgenes.txt", TextFile.W);
//            for (Integer i : genesBackground) {
//                outputbg.writeln(allGenesInMatrix.get(i));
//            }
//            outputbg.close();
//        }
//
//        // sort results by p-value
//        Collections.sort(results);
//
//        System.out.println("Output: " + outputfile);
//
//
//        TextFile output = new TextFile(outputfile, TextFile.W);
//        output.writeln("Pval\tFDR\tEmpiricalFDR\tTerm\tAnnotation\tAUC\tMeanFg\tMeanBg\tFoldChange\tFoldchangeZscore\t#PosGenes\t#NegGenes\tPosGenes\tNegGenes");
//        int rank = 0;
//        int sig = 0;
//
//        // FDR/qval adjustment
//        double[] qvals = new double[results.size()];
//        for (Pair<Double, Result> result : results) {
//            qvals[rank] = ((double) results.size() / (rank + 1)) * result.getLeft();
//            rank++;
//        }
//        for (int i = results.size() - 2; i > -1; i--) {
//            qvals[i] = Math.min(qvals[i], qvals[i + 1]);
//        }
//        rank = 0;
//        // write output
//        for (Pair<Double, Result> result : results) {
//            String id = result.getRight().rowname;
//            double auc = result.getRight().auc;
//            double meanfg = result.getRight().meanFg;
//            double meanbg = result.getRight().meanBg;
//            double foldchange = result.getRight().fc;
//            String annotation = id;
//
//            Pair<ArrayList<String>, ArrayList<String>> genes = result.getRight().posAndNegGenes;
//
//            if (idToAnnotation != null) {
//                annotation = idToAnnotation.get(id);
//                if (annotation == null) {
//                    annotation = id;
//                }
//            }
//            double qvalue = qvals[rank];
//            String posgenes = "";
//            String neggenes = "";
//            if (genes.getLeft().size() > 0) {
//                posgenes = Strings.concat(genes.getLeft(), Strings.semicolon);
//            }
//            if (genes.getRight().size() > 0) {
//                neggenes = Strings.concat(genes.getRight(), Strings.semicolon);
//            }
//            String outln = result.getLeft()
//                    + "\t" + qvalue
//                    + "\t" + result.getRight().empFDR
//                    + "\t" + result.getRight().rowname
//                    + "\t" + annotation
//                    + "\t" + auc
//                    + "\t" + meanfg
//                    + "\t" + meanbg
//                    + "\t" + foldchange
//                    + "\t" + result.getRight().fcZ
//                    + "\t" + genes.getLeft().size()
//                    + "\t" + genes.getRight().size()
//                    + "\t" + posgenes + "\t" + neggenes;
//
//            if (qvalue < 0.05) {
////                System.out.println(outln);
//                sig++;
//            }
//            rank++;
//            output.writeln(outln);
//        }
//
//        output.close();
//        System.out.println(sig + "/" + results.size() + " significant.");
//        System.out.println();
//    }

    private double[] stripNaN(double[] vals) {
        ArrayList<Double> nnvs = new ArrayList<>();
        for (double d : vals) {
            if (!Double.isNaN(d)) {
                nnvs.add(d);
            }
        }
        return Primitives.toPrimitiveArr(nnvs);
    }

    private Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> selectForeGroundAndBackground(ArrayList<String> geneListInMatrix,
                                                                                                            HashSet<String> geneset,
                                                                                                            HashSet<String> geneLimit,
                                                                                                            ArrayList<ArrayList<String>> geneBins,
                                                                                                            HashMap<String, Integer> geneToBin,
                                                                                                            boolean randomlyMatchNumberOfGenes,
                                                                                                            boolean verbose) {
        ArrayList<Integer> potentialGenesBackground = new ArrayList<>();
        ArrayList<Integer> genesForeground = new ArrayList<>();

        ArrayList<String> presentForegroundGenes = new ArrayList<>();
        if (verbose) {
            System.out.println("Input matrix has " + geneListInMatrix.size() + " gene IDs.");
        }
        for (int c = 0; c < geneListInMatrix.size(); c++) {
            String gene = geneListInMatrix.get(c);
            if (geneLimit == null || geneLimit.contains(gene)) {
                if (geneset.contains(gene)) {
                    genesForeground.add(c);
                    presentForegroundGenes.add(gene);
                } else {
                    potentialGenesBackground.add(c);
                }
            }
        }

        if (genesForeground.size() == 0) {
            System.out.println("Error: no genes left in foreground");
            System.exit(0);
        }
        if (verbose) {
            System.out.println(genesForeground.size() + " out of " + geneset.size() + " foreground genes matched against " + (geneListInMatrix.size()) + " gene columns.");
            System.out.println("Potential background genes: " + potentialGenesBackground.size());
        }
        // select background genes
        ArrayList<Integer> genesBackground = new ArrayList<>();
        if (geneBins == null) {
            if (randomlyMatchNumberOfGenes) {
                if (verbose) {
                    System.out.println("Randomly selecting " + genesForeground.size() + " genes from " + potentialGenesBackground.size() + " background genes.");
                }
                genesBackground = randomlySelectBackground(potentialGenesBackground, genesForeground.size());
            } else {
                if (verbose) {
                    System.out.println("Including all " + potentialGenesBackground.size() + " potential background genes.");
                }
                genesBackground = potentialGenesBackground;
            }
        } else {
            // match on the basis of gene bins
            if (verbose) {
                System.out.println("Selecting  " + genesForeground.size() + " background genes using gene expression bins.");
            }
            genesBackground = selectBackgroundFromBins(geneListInMatrix, genesForeground, potentialGenesBackground, geneBins, geneToBin);
            if (verbose) {
                System.out.println(genesBackground.size() + " background genes selected.");
            }
        }

        if (genesBackground.size() == 0) {
            System.out.println("Error: no genes left in background");
            System.exit(0);
        }


        return new Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>>(genesForeground, genesBackground, presentForegroundGenes);
    }

    private Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> getGeneIDs(ArrayList<String> geneListInMatrix,
                                                                                         HashSet<String> geneset,
                                                                                         HashSet<String> geneLimit,
                                                                                         boolean verbose) {
        ArrayList<Integer> potentialGenesBackground = new ArrayList<>();
        ArrayList<Integer> genesForeground = new ArrayList<>();

        ArrayList<String> presentForegroundGenes = new ArrayList<>();
        if (verbose) {
            System.out.println("Input matrix has " + geneListInMatrix.size() + " gene IDs.");
        }
        for (int c = 0; c < geneListInMatrix.size(); c++) {
            String gene = geneListInMatrix.get(c);
            if (geneLimit == null || geneLimit.contains(gene)) {
                if (geneset.contains(gene)) {
                    genesForeground.add(c);
                    presentForegroundGenes.add(gene);
                    if (includeForegroundGenesInBackgroundSet) {
                        potentialGenesBackground.add(c);
                    }
                } else {
                    potentialGenesBackground.add(c);
                }
            }
        }

        return new Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>>(genesForeground, potentialGenesBackground, presentForegroundGenes);
    }

    private void copy(double[] row, double[] foregroundGeneVals, double[] backgroundGeneVals, ArrayList<Integer> genesForeground, ArrayList<Integer> genesBackground) {
        int ctr = 0;
        for (Integer i : genesForeground) {
            foregroundGeneVals[ctr] = row[i];
            ctr++;
        }
        ctr = 0;
        for (Integer i : genesBackground) {
            backgroundGeneVals[ctr] = row[i];
            ctr++;
        }
    }

    private ArrayList<Integer> selectBackgroundFromBins(ArrayList<String> allGenesInMatrix,
                                                        ArrayList<Integer> genesForeground,
                                                        ArrayList<Integer> genesBackground,
                                                        ArrayList<ArrayList<String>> geneBins,
                                                        HashMap<String, Integer> geneToBin) {
        int[] nrsPerBin = new int[geneBins.size()];
        HashSet<String> geneSet = new HashSet<>();
        HashSet<String> geneSetBg = new HashSet<>();
        for (Integer i : genesForeground) {
            String gene = allGenesInMatrix.get(i);
            geneSet.add(gene);
            Integer binno = geneToBin.get(gene);
            if (binno != null) {
                nrsPerBin[binno]++;
            } else {
                System.out.println(gene + " is not in any gene expression bin");
            }
        }

        for (Integer i : genesBackground) {
            String gene = allGenesInMatrix.get(i);
            geneSetBg.add(gene);
        }

        ArrayList<ArrayList<String>> backgroundBins = new ArrayList<>();
        for (ArrayList<String> bin : geneBins) {
            ArrayList<String> bgBin = new ArrayList<>();
            for (String s : bin) {
                if (geneSetBg.contains(s)) {
                    bgBin.add(s);
                }
            }
            backgroundBins.add(bgBin);
        }

        HashSet<String> selectedIds = new HashSet<>();
        int remainder = 0;
        for (int b = 0; b < nrsPerBin.length; b++) {
            ArrayList<String> bgBin = backgroundBins.get(b);
            int toSelect = nrsPerBin[b] + remainder;
            if (bgBin.size() <= toSelect) {
                // add all
                selectedIds.addAll(bgBin);
                if (remainder > 0) {
                    remainder -= bgBin.size();
                    if (remainder < 0) {
                        remainder = 0;
                    }
                }
            } else {
                HashSet<String> selected = new HashSet<>();
                while (selected.size() < toSelect) {
                    selected.add(bgBin.remove((int) Math.floor((Math.random() * bgBin.size()))));
                }
                selectedIds.addAll(selected);
                remainder = 0;
            }
        }

        ArrayList<Integer> bgOutput = new ArrayList<>();
        for (Integer i : genesBackground) {
            String gene = allGenesInMatrix.get(i);
            if (selectedIds.contains(gene)) {
                bgOutput.add(i);
            }
        }
        return bgOutput;
    }


    private ArrayList<Integer> randomlySelectBackground(ArrayList<Integer> genesBackground, int n) {
        ArrayList<Integer> output = new ArrayList<>();

        if (genesBackground.size() <= n) {
            return genesBackground;
        } else {
            while (output.size() < n) {
                output.add(genesBackground.remove((int) Math.floor((Math.random() * genesBackground.size()))));
            }
        }
        return output;
    }

    private Pair<ArrayList<ArrayList<String>>, HashMap<String, Integer>> readAndBin(String meanExpFile, HashSet<String> limitGenes) throws IOException {

        TextFile tf = new TextFile(meanExpFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<Pair<String, Double>> genes = new ArrayList<>();
        while (elems != null) {
            try {
                String gene = elems[0].split("\\.")[0];
                if (limitGenes == null || limitGenes.contains(gene)) {
                    Double val = Double.parseDouble(elems[1]);
                    Pair<String, Double> p = new Pair<String, Double>(gene, val, Pair.SORTBY.RIGHT);
                    genes.add(p);
                }
            } catch (NumberFormatException e) {

            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        Collections.sort(genes);

        ArrayList<ArrayList<String>> geneBins = new ArrayList<>();
        double nrGenesPerBin = Math.ceil((double) genes.size() / 10);
        HashMap<String, Integer> geneToBin = new HashMap<>();
        int bin = 0;
        int ctr = 0;
        ArrayList<String> currentBin = new ArrayList<>();

        for (Pair<String, Double> g : genes) {
            currentBin.add(g.getLeft());
            geneToBin.put(g.getLeft(), bin);
            ctr++;
            if (ctr == nrGenesPerBin) {
                geneBins.add(currentBin);
                currentBin = new ArrayList<>();
                ctr = 0;
                bin++;
            }
        }
        geneBins.add(currentBin);

        return new Pair<ArrayList<ArrayList<String>>, HashMap<String, Integer>>(geneBins, geneToBin);

    }

    private HashSet<String> readSet(String limitGenes) throws IOException {
        TextFile tf = new TextFile(limitGenes, TextFile.R);
        HashSet<String> output = new HashSet<String>();
        ArrayList<String> set = tf.readAsArrayList();
        for (String s : set) {
            String[] selems = s.split("\\.");
            output.add(selems[0]);
        }
        tf.close();
        return output;
    }

    private Pair<ArrayList<String>, ArrayList<String>> split(double[] valsForGeneSet, ArrayList<String> presentGenes, HashMap<String, String> ensToHugo) {


        ArrayList<Pair<String, Double>> negP = new ArrayList<>();
        ArrayList<Pair<String, Double>> posP = new ArrayList<>();


        DecimalFormat f = new DecimalFormat("#.##");
        for (int d = 0; d < valsForGeneSet.length; d++) {
            double v = valsForGeneSet[d];
            String gene = presentGenes.get(d);
            String hugo = ensToHugo.get(gene);
            if (hugo == null) {
                hugo = gene;
            }
            if (v >= 0) {
                posP.add(new Pair<>(hugo, v, Pair.SORTBY.RIGHT));
//                pos.add(hugo + ":" + f.format(v));
            } else {
                negP.add(new Pair<>(hugo, v, Pair.SORTBY.RIGHT));
            }
        }

        Collections.sort(negP);
        Collections.sort(posP, Collections.reverseOrder());
        ArrayList<String> pos = new ArrayList<>();
        ArrayList<String> neg = new ArrayList<>();

        for (Pair<String, Double> p : negP) {
            neg.add(p.getLeft() + ":" + f.format(p.getRight()));
        }
        for (Pair<String, Double> p : posP) {
            pos.add(p.getLeft() + ":" + f.format(p.getRight()));
        }


        return new Pair<>(pos, neg);
    }

    private double[] convertToDouble(String[] data) {
        double[] output = new double[data.length - 1];
        for (int i = 1; i < data.length; i++) {
            output[i - 1] = Double.parseDouble(data[i]);

        }


        return output;
    }

    public double[] addMin(double[] data) {
        double[] output = new double[data.length];
        double min = JSci.maths.ArrayMath.min(data);
        if (min < 0) {
            min = Math.abs(min);
        }
        for (int i = 0; i < output.length; i++) {
            output[i] = data[i] + min;
        }
        return output;
    }

    public double[] square(double[] data) {
        double[] output = new double[data.length];

        for (int i = 0; i < output.length; i++) {
            output[i] = data[i] * data[i];
        }
        return output;
    }

}
