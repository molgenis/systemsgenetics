package nl.systemsgenetics.genenetworkpathwayenrichment;

import org.apache.commons.math3.distribution.TDistribution;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class GeneSetEnrichment {

    private TDistribution tdist;
    private ArrayList<ArrayList<Integer>> geneBinsAverageExpression;
    private HashMap<Integer, Integer> geneToAverageExpressionBin;

    public enum TESTTYPE {
        LOGISTIC,
        FISHEREXACT,
        WILCOXON
    }

    public void run(String matrixFile,
                    String foregroundFile,
                    String backgroundFile,
                    String annotationFile,
                    String ensgToHugoFile,
                    String meanExpFile,
                    int permutations,
                    TESTTYPE testtype,
                    String outputFile
    ) throws IOException {
        System.out.println("Gene set enrichment using logistic regression.");
        System.out.println("Using " + permutations + " permutations.");
        System.out.println("Using: " + testtype);
        System.out.println("Randomly selecting background genes to match foreground size: " + makeBackgroundEqualSize);

        //
        // initialization
        //

        HashSet<String> inputForegroundGeneset = Utils.readSet(foregroundFile);
        HashSet<String> inputBackroundGeneset = Utils.readSet(backgroundFile);

        HashSet<String> inputAllGenes = new HashSet<>();
        inputAllGenes.addAll(inputBackroundGeneset);
        inputAllGenes.addAll(inputForegroundGeneset);
        System.out.println(inputAllGenes.size() + " total genes.");


        // map ensembl ids to Hugo IDs
        HashMap<String, String> ensgToHugo = null;
        if (ensgToHugoFile != null) {
            ensgToHugo = Utils.readEnsToHugo(ensgToHugoFile);
        }

        // Read pathway annotations
        HashMap<String, String> pathwayIDToAnnotation = null;
        if (annotationFile != null) {
            pathwayIDToAnnotation = Utils.readPathwayIDToAnnotation(annotationFile);
        }

        //
        // testing
        //
        System.out.println("Processing: " + matrixFile);
        TextFile tf = new TextFile(matrixFile, TextFile.R);

        // get gene IDs from header
        String[] header = tf.readLineElems(TextFile.tab);
        ArrayList<Integer> colsToInclude = new ArrayList<>();
        ArrayList<Double> foreGroundOrBackground = new ArrayList<>();
        ArrayList<String> geneOrderInMatrix = new ArrayList<>();
        int nrOverlappingForeground = 0;
        for (int i = 1; i < header.length; i++) {
            // remove ENSG version id
            String[] elems = header[i].split("\\.");
            String gene = elems[0];
            if (inputAllGenes.contains(gene)) {
                colsToInclude.add(i);
                geneOrderInMatrix.add(gene);
                if (inputForegroundGeneset.contains(gene)) {
                    foreGroundOrBackground.add(1d);
                    nrOverlappingForeground++;
                } else {
                    foreGroundOrBackground.add(0d);
                }
            }
        }

        // read average gene expression values from file
        geneBinsAverageExpression = null;
        geneToAverageExpressionBin = null;
        if (meanExpFile != null) {
            Pair<ArrayList<ArrayList<Integer>>, HashMap<Integer, Integer>> genexp = Utils.readAndBin2(meanExpFile, geneOrderInMatrix);
            geneBinsAverageExpression = genexp.getLeft();
            geneToAverageExpressionBin = genexp.getRight();
        }

        System.out.println(colsToInclude.size() + " genes overlap between input genesets and pathway matrix.");
        if (colsToInclude.isEmpty()) {
            System.err.println("No overlapping genes.");
            System.exit(-1);
        }
        System.out.println(nrOverlappingForeground + " genes in foreground set present in pathway matrix.");

        TextFile outf = new TextFile(outputFile, TextFile.W);
        outf.writeln("P\tPermutationPvalue\tFDR\tFoldChange\tPathway\tPathwayAnnotation\tTestSpecificOutput");

        System.out.println("Output will be here: " + outputFile);
        int pwctr = 0;
        double[] y = Primitives.toPrimitiveArr(foreGroundOrBackground);
        tdist = new TDistribution(colsToInclude.size() - 1);
        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<Pair<Double, Result>> results = new ArrayList<>();
        while (elems != null) {
            String pw = elems[0];


            double[] xtmp = new double[y.length];
            int xctr = 0;
            for (Integer i : colsToInclude) {
                xtmp[xctr] = Double.parseDouble(elems[i]);
                xctr++;
            }

            /*
            y: 1 - gene is in foreground, 0 - gene is in background
            x: either binary or continuous (depending on test)
             */

            Result r = null;
            if (testtype.equals(TESTTYPE.FISHEREXACT)) {
                r = fetTest(xtmp, y, permutations);
            } else if (testtype.equals(TESTTYPE.LOGISTIC)) {
                // lrTest(xtmp, y, permutations);
            } else if (testtype.equals(TESTTYPE.WILCOXON)) {
                r = wilcoxonTest(xtmp, y, permutations);
            }

            if (r != null) {

                r.pw = pw;
                results.add(new Pair<Double, Result>(r.p, r));

//            outf.writeln(p + "\t" + fdr + "\t" + elems[0] + "\t" + annot + "\t" + betas[0][1] + "\t" + stder[0][1] + "\t" + t);

            }
            pwctr++;
            if (pwctr % 10 == 0) {
                System.out.print(pwctr + " pathways processed, " + results.size() + " have a result\r");
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        // calculate BH FDR
        Collections.sort(results);
        double[] qvals = fdr(results, false);
        for (int i = 0; i < results.size(); i++) {
            Result r = results.get(i).getRight();
            String annot = pathwayIDToAnnotation.get(r.pw);
            String outln = r.p + "\t" + r.termFDR + "\t" + qvals[i] + "\t" + r.pw + "\t" + r.fc + "\t" + annot + "\t" + r.testSpecific;
            outf.writeln(outln);
        }
        outf.close();

        System.out.println(pwctr + " total pathways processed.");


    }

    private double[] fdr(ArrayList<Pair<Double, Result>> results, boolean fcPVal) {
        int rank = 0;
        double[] qvals = new double[results.size()];
        for (Pair<Double, Result> result : results) {
            qvals[rank] = ((double) results.size() / (rank + 1)) * result.getRight().p;
            rank++;
        }
        for (int i = results.size() - 2; i > -1; i--) {
            qvals[i] = Math.min(qvals[i], qvals[i + 1]);
        }


        return qvals;
    }

    public void setMakeBackgroundEqualSize(boolean makeBackgroundEqualSize) {
        this.makeBackgroundEqualSize = makeBackgroundEqualSize;
    }

    private Result wilcoxonTest(double[] xOrig, double[] yOrig, int permutations) {

        double[] x = xOrig;
        double[] y = yOrig;
        if (makeBackgroundEqualSize) {
            double[][] xAndY = Utils.makeBackgroundEqualSize(xOrig, yOrig, geneBinsAverageExpression, geneToAverageExpressionBin);
            x = xAndY[0];
            y = xAndY[1];
        }

        ArrayList<Double> foreground = new ArrayList<>();
        ArrayList<Double> background = new ArrayList<>();
        for (int d = 0; d < x.length; d++) {
            if (y[d] == 1d) {
                foreground.add(x[d]);
                background.add(x[d]);
            } else {
                background.add(x[d]);
            }
        }

        WilcoxonMannWhitney wmw = new WilcoxonMannWhitney();
        double p = wmw.returnWilcoxonMannWhitneyPValue(Primitives.toPrimitiveArr(foreground), Primitives.toPrimitiveArr(background));
        double auc = wmw.getAUC();
        double fc = Descriptives.mean(Primitives.toPrimitiveArr(foreground)) / Descriptives.mean(Primitives.toPrimitiveArr(background));
        AtomicInteger c = new AtomicInteger();
        IntStream.range(0, permutations).parallel().forEach(perm -> {
            double[] xperm = null;
            double[] yperm = null;
            Utils u = new Utils();
            if (geneBinsAverageExpression != null) {
                double[][] xAndYPerm = Utils.randomMatchForeground(xOrig, yOrig, geneBinsAverageExpression, geneToAverageExpressionBin);
                xperm = xAndYPerm[0];
                yperm = xAndYPerm[1];
            } else {
                // randomly pick genes
                yperm = new double[yOrig.length];
                System.arraycopy(yOrig, 0, yperm, 0, yOrig.length);
                u.shuffle(yperm);
                xperm = xOrig;
            }

            if (makeBackgroundEqualSize) {
                double[][] xAndY = Utils.makeBackgroundEqualSize(xperm, yperm, geneBinsAverageExpression, geneToAverageExpressionBin);
                xperm = xAndY[0];
                yperm = xAndY[1];
            }

            ArrayList<Double> foregroundPerm = new ArrayList<>();
            ArrayList<Double> backgroundPerm = new ArrayList<>();
            for (int d = 0; d < xperm.length; d++) {
                if (yperm[d] == 1d) {
                    foregroundPerm.add(xperm[d]);
                    backgroundPerm.add(xperm[d]);
                } else {
                    backgroundPerm.add(xperm[d]);
                }
            }
            WilcoxonMannWhitney wmwp = new WilcoxonMannWhitney();
            double pperm = wmw.returnWilcoxonMannWhitneyPValue(Primitives.toPrimitiveArr(foregroundPerm), Primitives.toPrimitiveArr(backgroundPerm));

            if (pperm <= p) {
                c.getAndIncrement();
            }
        });

        double fdrForTerm = ((double) c.get()) / permutations;

        Result r = new Result();
        r.p = p;
        r.fc = fc;
        r.termFDR = fdrForTerm;
        return r;
    }

    private boolean makeBackgroundEqualSize;

    private class Result {

        public String pw;
        public String testSpecific;
        double p = 1;
        double fc = 0;
        double termFDR = 1;

    }

    private Result fetTest(double[] xOrig, double[] yOrig, int permutations) {

        double[] x = xOrig;
        double[] y = yOrig;
        if (makeBackgroundEqualSize) {
            double[][] xAndY = Utils.makeBackgroundEqualSize(xOrig, yOrig, geneBinsAverageExpression, geneToAverageExpressionBin);
            x = xAndY[0];
            y = xAndY[1];
        }

        FisherExactTest fet = new FisherExactTest();
        int n11 = 0;
        int n12 = 0;
        int n21 = 0;
        int n22 = 0;

        for (int d = 0; d < x.length; d++) {
            if (y[d] == 1) {
                if (x[d] == 0) {
                    n12++;
                } else {
                    n11++;
                }
            } else {
                if (x[d] == 0) {
                    n22++;
                } else {
                    n21++;
                }
            }
        }

        if (n11 + n21 < 10) { // n11 and n21 contain the counts of genes in each pathway..
            return null;
        } else {

            double p = fet.getFisherPValue(n11, n12, n21, n22);
            double fc = ((double) n11 / (n11 + n12)) / ((double) n21 / (n21 + n22));
            AtomicInteger c = new AtomicInteger();

            IntStream.range(0, permutations).parallel().forEach(perm -> {
                Utils u = new Utils();
                double[] yperm = null;
                double[] xperm = null;

                // match on the basis of average 'expression'
                if (geneBinsAverageExpression != null) {
                    double[][] xAndYPerm = Utils.randomMatchForeground(xOrig, yOrig, geneBinsAverageExpression, geneToAverageExpressionBin);
                    xperm = xAndYPerm[0];
                    yperm = xAndYPerm[1];
                } else {
                    // randomly pick genes
                    yperm = new double[yOrig.length];
                    System.arraycopy(yOrig, 0, yperm, 0, yOrig.length);
                    u.shuffle(yperm);
                    xperm = xOrig;
                }

                if (makeBackgroundEqualSize) {
                    double[][] xAndY = Utils.makeBackgroundEqualSize(xperm, yperm, geneBinsAverageExpression, geneToAverageExpressionBin);
                    xperm = xAndY[0];
                    yperm = xAndY[1];
                }

                int n11p = 0;
                int n12p = 0;
                int n21p = 0;
                int n22p = 0;

                for (int d = 0; d < xperm.length; d++) {
                    if (yperm[d] == 1) {
                        if (xperm[d] == 0) {
                            n12p++;
                        } else {
                            n11p++;
                        }
                    } else {
                        if (xperm[d] == 0) {
                            n22p++;
                        } else {
                            n21p++;
                        }
                    }
                }
                FisherExactTest fetp = new FisherExactTest();
                double pperm = fetp.getFisherPValue(n11p, n12p, n21p, n22p);
                if (pperm <= p) {
                    c.getAndIncrement();
                }
            });

            double fdrForTerm = ((double) c.get()) / permutations;

            Result r = new Result();
            r.p = p;
            r.fc = fc;
            r.termFDR = fdrForTerm;
            r.testSpecific = "n11: " + n11 + ", n12: " + n12 + ", n21: " + n21 + ", n22: " + n22 + ", group: " + (n11 + n21);

            return r;
        }
    }

//    private void lrTest(double[] xtmp, double[] y, int permutations) {
//        LogisticRegressionOptimized reg = new LogisticRegressionOptimized();
//
//        double[][] x = new double[y.length][2];
//        for (int i = 0; i < x.length; i++) {
//            x[i][0] = 1;
//            x[i][1] = xtmp[i];
//        }
//
//        LogisticRegressionResult result = reg.binomial(y, x);
//
//        double[][] betas = result.getBeta(); // format [0][0] == intercept, [0][1] == pathway z-scores
//        double[][] stder = result.getStderrs();
//
//        double t = betas[0][1] / stder[0][1];
//        double p = 2 * tdist.cumulativeProbability(-Math.abs(t));
//
//
//        // double[] pvalsPerm = new double[permutations];
//        AtomicInteger c = new AtomicInteger();
//        IntStream.range(0, permutations).parallel().forEach(perm -> {
//            LogisticRegressionOptimized regPerm = new LogisticRegressionOptimized();
//            // permute y
//            Utils u = new Utils();
//            double[] yperm = new double[y.length];
//            System.arraycopy(y, 0, yperm, 0, y.length);
//            u.shuffle(yperm);
//
//            LogisticRegressionResult resultPerm = regPerm.binomial(yperm, x);
//
//            double[][] betasPerm = resultPerm.getBeta(); // format [0][0] == intercept, [0][1] == pathway z-scores
//            double[][] stderPerm = resultPerm.getStderrs();
//
//            double tperm = betasPerm[0][1] / stderPerm[0][1];
//            double pperm = 2 * tdist.cumulativeProbability(-Math.abs(tperm));
////                pvalsPerm[perm] = pperm;
//            if (pperm <= p) {
//                c.getAndIncrement();
//            }
//        });
//
//        double fdr = ((double) c.get()) / permutations;
//
//    }


}
