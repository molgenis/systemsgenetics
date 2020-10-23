package nl.systemsgenetics.genenetworkpathwayenrichment;


import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class GeneSetEnrichmentWilcoxon {

    public void run(String matrix,
                    String genesetfile,
                    String annotationFile,
                    String ensgToHugoFile,
                    String limitGenes,
                    String meanExpFile,
                    boolean randomlyMatchNumberOfGenes,
                    String outputfile) throws IOException {


        // limit analysis to these genes
        HashSet<String> geneLimit = null;
        if (limitGenes != null) {
            geneLimit = readSet(limitGenes);
        }

        // read average gene expression values from file
        ArrayList<ArrayList<String>> geneBins = null;
        HashMap<String, Integer> geneToBin = null;
        if (meanExpFile != null) {
            Pair<ArrayList<ArrayList<String>>, HashMap<String, Integer>> data = readAndBin(meanExpFile, geneLimit);
            geneBins = data.getLeft();
            geneToBin = data.getRight();

            // limit to genes for which we have expression data
            HashSet<String> newGeneLimit = new HashSet<>();
            newGeneLimit.addAll(geneToBin.keySet());
            geneLimit = newGeneLimit;
        }

        // map ensembl ids to Hugo IDs
        HashMap<String, String> ensToHugo = null;
        if (ensgToHugoFile != null) {
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
        }

        // Read pathway annotations
        HashMap<String, String> idToAnnotation = null;
        if (annotationFile != null) {
            System.out.println("reading annotation: " + annotationFile);
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
            System.out.println(idToAnnotation.size() + " annotations read");
        }

        // read geneset to test
        HashSet<String> geneset = new HashSet<String>();
        System.out.println("Reading geneset: " + genesetfile);
        TextFile tf = new TextFile(genesetfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            geneset.add(elems[0].split("\\.")[0]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(geneset.size() + " genes in foreground geneset");

        HashMap<String, Pair<ArrayList<String>, ArrayList<String>>> posAndNegGenes = new HashMap<>();
        WilcoxonMannWhitney w = new WilcoxonMannWhitney();
        ArrayList<Pair<Double, String>> results = new ArrayList<>();
        System.out.println("Using: " + matrix);
        // use a text based matrix
        if (matrix.endsWith(".txt.gz") || matrix.endsWith(".txt")) {
            TextFile tf3 = new TextFile(matrix, TextFile.R);
            String[] header = tf3.readLineElems(TextFile.tab);
            ArrayList<String> allGenesInMatrix = new ArrayList<>();
            for (String s : header) {
                allGenesInMatrix.add(s);
            }
            Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> genesets = selectForeGroundAndBackground(allGenesInMatrix, geneset, geneLimit, geneBins, geneToBin, randomlyMatchNumberOfGenes);
            ArrayList<Integer> genesForeground = genesets.getLeft();
            ArrayList<Integer> genesBackground = genesets.getMiddle();
            ArrayList<String> presentForegroundGenes = genesets.getRight();
            String[] data = tf3.readLineElems(TextFile.tab);
            double[] foregroundGeneVals = new double[genesForeground.size()];
            double[] backgroundGeneVals = new double[genesBackground.size()];

            while (data != null) {
                double[] row = convertToDouble(data);
                copy(row, foregroundGeneVals, backgroundGeneVals, genesForeground, genesBackground);
                Pair<ArrayList<String>, ArrayList<String>> genesplit = split(foregroundGeneVals, presentForegroundGenes, ensToHugo);
                double pval = w.returnWilcoxonMannWhitneyPValue(foregroundGeneVals, backgroundGeneVals);
                String rowname = data[0];
                posAndNegGenes.put(rowname, genesplit);
                results.add(new Pair<Double, String>(pval, rowname, Pair.SORTBY.LEFT));
                data = tf3.readLineElems(TextFile.tab);
            }
            tf3.close();
        } else { // or a binary matrix
            DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(matrix);
            Set<String> cols = it.getCols();
            ArrayList<String> allGenesInMatrix = new ArrayList<String>();
            allGenesInMatrix.addAll(cols);
            Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> genesets = selectForeGroundAndBackground(allGenesInMatrix, geneset, geneLimit, geneBins, geneToBin, randomlyMatchNumberOfGenes);
            ArrayList<Integer> genesForeground = genesets.getLeft();
            ArrayList<Integer> genesBackground = genesets.getMiddle();
            ArrayList<String> presentForegroundGenes = genesets.getRight();
            double[] foregroundGeneVals = new double[genesForeground.size()];
            double[] backgroundGeneVals = new double[genesBackground.size()];

            ArrayList<String> rowIds = new ArrayList<String>();
            rowIds.addAll(it.getRows());
            int rctr = 0;
            for (double[] row : it) {
                copy(row, foregroundGeneVals, backgroundGeneVals, genesForeground, genesBackground);
                Pair<ArrayList<String>, ArrayList<String>> genesplit = split(foregroundGeneVals, presentForegroundGenes, ensToHugo);
                double pval = w.returnWilcoxonMannWhitneyPValue(foregroundGeneVals, backgroundGeneVals);
                String rowname = rowIds.get(rctr);
                posAndNegGenes.put(rowname, genesplit);
                rctr++;
                results.add(new Pair<Double, String>(pval, rowname, Pair.SORTBY.LEFT));
            }
            it.close();
        }

        // sort results by p-value
        Collections.sort(results);
        System.out.println("Output: " + outputfile);
        TextFile output = new TextFile(outputfile, TextFile.W);
        output.writeln("Pval\tFDR\tTerm\tAnnotation\tPosGenes\tNegGenes");
        int rank = 0;
        int sig = 0;

        // FDR/qval adjustment
        double[] qvals = new double[results.size()];
        for (Pair<Double, String> result : results) {
            qvals[rank] = ((double) results.size() / (rank + 1)) * result.getLeft();
            rank++;
        }
        for (int i = results.size() - 2; i > -1; i--) {
            qvals[i] = Math.min(qvals[i], qvals[i + 1]);
        }
        rank = 0;
        // write output
        for (Pair<Double, String> result : results) {
            String id = result.getRight();
            String annotation = id;

            Pair<ArrayList<String>, ArrayList<String>> genes = posAndNegGenes.get(id);

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
            String outln = result.getLeft() + "\t" + qvalue + "\t" + result.getRight() + "\t" + annotation + "\t" + posgenes + "\t" + neggenes;

            if (qvalue < 0.05) {
                System.out.println(outln);
                sig++;
            }
            rank++;
            output.writeln(outln);
        }

        output.close();
        System.out.println(sig + "/" + results.size() + " significant.");
        System.out.println();
    }

    private Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>> selectForeGroundAndBackground(ArrayList<String> geneListInMatrix, HashSet<String> geneset, HashSet<String> geneLimit, ArrayList<ArrayList<String>> geneBins,
                                                                                                            HashMap<String, Integer> geneToBin, boolean randomlyMatchNumberOfGenes) {
        ArrayList<Integer> potentialGenesBackground = new ArrayList<>();
        ArrayList<Integer> genesForeground = new ArrayList<>();

        ArrayList<String> presentForegroundGenes = new ArrayList<>();
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

        System.out.println("Potential background genes: " + potentialGenesBackground.size());

        // select background genes
        ArrayList<Integer> genesBackground = new ArrayList<>();
        if (geneBins == null) {
            if (randomlyMatchNumberOfGenes) {
                System.out.println("Randomly selecting " + genesForeground.size() + " genes from " + potentialGenesBackground.size() + " background genes.");
                genesBackground = randomlySelectBackground(potentialGenesBackground, genesForeground.size());
            } else {
                System.out.println("Including all " + potentialGenesBackground.size() + " potential background genes.");
                genesBackground = potentialGenesBackground;
            }
        } else {
            // match on the basis of gene bins
            System.out.println("Selecting  " + genesForeground.size() + " background genes using gene expression bins.");
            genesBackground = selectBackgroundFromBins(geneListInMatrix, genesForeground, potentialGenesBackground, geneBins, geneToBin);
            System.out.println(genesBackground.size() + " background genes selected.");
        }

        if (genesBackground.size() == 0) {
            System.out.println("Error: no genes left in background");
            System.exit(0);
        }

        System.out.println(genesForeground.size() + " out of " + geneset.size() + " matched against " + (geneListInMatrix.size() - 1) + " gene columns.");

        return new Triple<ArrayList<Integer>, ArrayList<Integer>, ArrayList<String>>(genesForeground, genesBackground, presentForegroundGenes);
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


}
