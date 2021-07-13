package nl.systemsgenetics.genenetworkpathwayenrichment;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

public class Utils {

    public static HashSet<String> readSet(String limitGenes) throws IOException {
        System.out.println("Reading set: " + limitGenes);
        TextFile tf = new TextFile(limitGenes, TextFile.R);
        HashSet<String> output = new HashSet<String>();
        ArrayList<String> set = tf.readAsArrayList();
        for (String s : set) {
            String[] selems = s.split("\\.");
            output.add(selems[0]);
        }
        tf.close();
        System.out.println(output.size() + " genes in set.");
        return output;
    }

    public static HashMap<String, String> readEnsToHugo(String ensgToHugoFile) throws IOException {
        System.out.println("Reading Ensembl to HUGO annotation: " + ensgToHugoFile);
        HashMap<String, String> ensgToHugo = new HashMap<>();
        TextFile tf = new TextFile(ensgToHugoFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > 1) {
                String ensg = elems[0].split("\\.")[0];
                ensgToHugo.put(ensg, elems[1]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(ensgToHugo.size() + " gene annotations read.");
        return ensgToHugo;
    }

    public static HashMap<String, String> readPathwayIDToAnnotation(String annotationFile) throws IOException {
        System.out.println("Reading pathway annotation: " + annotationFile);
        HashMap<String, String> pathwayIDToAnnotation = new HashMap<String, String>();
        TextFile tf2 = new TextFile(annotationFile, TextFile.R);
        String[] elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > 1) {
                pathwayIDToAnnotation.put(elems[0], elems[1]);
            }
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();
        System.out.println(pathwayIDToAnnotation.size() + " pathway annotations read");
        return pathwayIDToAnnotation;
    }

    private Random random;

    public static double[][] makeBackgroundEqualSize(double[] xOrig, double[] yOrig, ArrayList<ArrayList<Integer>> geneBinsAverageExpression, HashMap<Integer, Integer> geneToAverageExpressionBin) {

        ArrayList<Double> newX = new ArrayList<>();
        ArrayList<Double> newY = new ArrayList<>();

        for (int i = 0; i < yOrig.length; i++) {
            if (yOrig[i] == 1) {
                newX.add(xOrig[i]);
                newY.add(1d);

                // match a gene
                if (geneBinsAverageExpression != null) {
                    Integer binId = geneToAverageExpressionBin.get(i);
                    ArrayList<Integer> bin = geneBinsAverageExpression.get(binId);
                    int randomId = (int) Math.floor(Math.random() * bin.size());
                    newX.add(xOrig[randomId]);
                    newY.add(0d);
                } else {
                    // randomly pick a gene
                    int randomId = (int) Math.floor(Math.random() * xOrig.length);
                    newX.add(xOrig[randomId]);
                    newY.add(0d);
                }
            }
        }

        double[][] output = new double[2][];
        output[0] = Primitives.toPrimitiveArr(newX);
        output[1] = Primitives.toPrimitiveArr(newY);
        return output;
    }

    public static double[][] randomMatchForeground(double[] xOrig, double[] yOrig, ArrayList<ArrayList<Integer>> geneBinsAverageExpression, HashMap<Integer, Integer> geneToAverageExpressionBin) {

        ArrayList<Double> newX = new ArrayList<>();
        ArrayList<Double> newY = new ArrayList<>();

        for (int i = 0; i < yOrig.length; i++) {
            if (yOrig[i] == 1) {
                // match a gene
                if (geneBinsAverageExpression != null) {
                    Integer binId = geneToAverageExpressionBin.get(i);
                    ArrayList<Integer> bin = geneBinsAverageExpression.get(binId);
                    int randomId = (int) Math.floor(Math.random() * bin.size());
                    newX.add(xOrig[randomId]);
                    newY.add(1d);
                } else {
                    // randomly pick a gene
                    int randomId = (int) Math.floor(Math.random() * xOrig.length);
                    newX.add(xOrig[randomId]);
                    newY.add(1d);
                }
            } else {
                newX.add(xOrig[i]);
                newY.add(0d);
            }
        }

        double[][] output = new double[2][];
        output[0] = Primitives.toPrimitiveArr(newX);
        output[1] = Primitives.toPrimitiveArr(newY);
        return output;


    }

    /**
     * Code from method java.util.Collections.shuffle();
     */
    public void shuffle(double[] array) {
        if (random == null) random = new Random();
        int count = array.length;
        for (int i = count; i > 1; i--) {
            swap(array, i - 1, random.nextInt(i));
        }
    }

    private void swap(double[] array, int i, int j) {
        double temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }


    public static Pair<ArrayList<ArrayList<Integer>>, HashMap<Integer, Integer>> readAndBin2(String meanExpFile, ArrayList<String> intendedGeneOrder) throws IOException {
        HashMap<String, Integer> genemap = new HashMap<>();
        for (int i = 0; i < intendedGeneOrder.size(); i++) {
            genemap.put(intendedGeneOrder.get(i), i);
        }

        TextFile tf = new TextFile(meanExpFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<Pair<Integer, Double>> genes = new ArrayList<>();
        HashSet<Integer> found = new HashSet<>();
        while (elems != null) {
            try {
                String gene = elems[0].split("\\.")[0];
                if (intendedGeneOrder.contains(gene)) {
                    Double val = Double.parseDouble(elems[1]);
                    Integer geneId = genemap.get(gene);
                    found.add(geneId);
                    Pair<Integer, Double> p = new Pair<Integer, Double>(geneId, val, Pair.SORTBY.RIGHT);
                    genes.add(p);
                }
            } catch (NumberFormatException e) {

            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        if (found.size() != intendedGeneOrder.size()) {
            if (found.isEmpty()) {
                System.out.println("None of the IDs in " + meanExpFile + " match the ones in the matrix. Please check.");
                System.exit(0);
            }
            System.out.println("Warning not all genes present in average expression file: " + meanExpFile);
            System.out.println(found.size() + " out of " + intendedGeneOrder.size() + " gene ids matched.");
            System.out.println("Will randomly assign a bin for those.");
        }

        Collections.sort(genes);

        ArrayList<ArrayList<Integer>> geneBins = new ArrayList<>();
        int nrbins = 10;
        double nrGenesPerBin = Math.ceil((double) genes.size() / nrbins);
        HashMap<Integer, Integer> geneToBin = new HashMap<>();
        int bin = 0;
        int ctr = 0;
        ArrayList<Integer> currentBin = new ArrayList<>();

        for (Pair<Integer, Double> g : genes) {
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

        for (int i = 0; i < intendedGeneOrder.size(); i++) {
            if (!found.contains(i)) {
                int randombin = (int) Math.floor(Math.random() * nrbins);
                geneBins.get(randombin).add(i);
            }
        }

        return new Pair<ArrayList<ArrayList<Integer>>, HashMap<Integer, Integer>>(geneBins, geneToBin);
    }

    public static Pair<ArrayList<ArrayList<String>>, HashMap<String, Integer>> readAndBin(String meanExpFile, ArrayList<String> limitGeneList) throws IOException {
        HashSet<String> limitGenes = null;
        if (limitGeneList != null) {
            limitGenes = new HashSet<>();
            limitGenes.addAll(limitGeneList);
        }


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
}
