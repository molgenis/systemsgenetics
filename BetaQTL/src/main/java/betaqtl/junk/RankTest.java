package betaqtl.junk;

import betaqtl.stat.RankArray;

public class RankTest {

    public static void main(String[] args) {
        double[] data = new double[]{Double.NaN, Double.NaN, Double.NaN};
        RankArray ranker = new RankArray();
        double[] d = ranker.rank(data, true); // does this work with NaNs?
        for (int q = 0; q < data.length; q++) {
            System.out.println(q + "\t" + data[q] + "\t" + d[q]);
        }

        String gene = "ENSG00.1|ENSGE00.1";
        gene = gene.replaceAll("\\|", "_");
        gene = gene.replaceAll("\\.", "-");
        System.out.println(gene);
    }
}
