package mbqtl.junk;

import mbqtl.stat.BetaDistributionMLE;
import umcg.genetica.io.text.TextFile;
import umontreal.iro.lecuyer.probdist.BetaDist;

import java.io.IOException;
import java.util.Arrays;

public class ShuangTest {

    public static void main(String[] args) throws IOException {

        String permfile = "D:\\Sync\\TMP\\shuang\\betadist\\rs2350629_FAM118A.perms.tsv";
        String topfile = "D:\\Sync\\TMP\\shuang\\betadist\\rs2350629.tsv";

        permfile = "D:\\Sync\\TMP\\shuang\\betadist\\rs67483168perm.txt";
        permfile = "D:\\Sync\\TMP\\shuang\\run3\\rs67483168perm.txt";
        topfile = "D:\\Sync\\TMP\\shuang\\betadist\\rs67483168_CST3.tsv";


        TextFile tf = new TextFile(topfile, TextFile.R);
        String header = tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        String lowest = null;
        double lowestP = 1;
        while (elems != null) {
            double metap = Double.parseDouble(elems[12]);
            if (metap < lowestP) {
                lowest = elems[0];
                lowestP = metap;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        System.out.println(lowest + "\t" + lowestP);

        TextFile tf2 = new TextFile(permfile, TextFile.R);
        String[] headerperm = tf2.readLineElems(TextFile.tab);
        double[] pvals = new double[headerperm.length - 2];
        Arrays.fill(pvals, 1);
        elems = tf2.readLineElems(TextFile.tab);
        int rctr = 0;
        while (elems != null) {
            for (int q = 2; q < elems.length; q++) {
                double v = Double.parseDouble(elems[q]);
                if (v == 0) {
                    System.out.println("Permutation: " + q + " gene " + elems[0] + "\t" + elems[1]);
                    v = 2.0E-323D;
                }
                if (v < pvals[q - 2]) {
                    pvals[q - 2] = v;
                }
            }
            rctr++;
            elems = tf2.readLineElems(TextFile.tab);
        }

        double propBetter = 0;
        for (int q = 0; q < pvals.length; q++) {
            if (pvals[q] <= lowestP) {
                propBetter++;
            }
            System.out.println(pvals[q]);
        }
        propBetter /= pvals.length;
        System.out.println("Proportion better pvals: " + propBetter);
        BetaDistributionMLE mle = new BetaDistributionMLE();
        double[] shape = mle.fit(pvals);
        BetaDist betaDistribution = new BetaDist(shape[0], shape[1]);
        double betaAdjPval = betaDistribution.cdf(lowestP);
        System.out.println("Alpha: " + shape[0]);
        System.out.println("Beta: " + shape[1]);
        System.out.println("BetaAdjP: " + betaAdjPval);

    }
}
