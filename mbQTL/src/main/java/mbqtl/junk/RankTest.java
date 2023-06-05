package mbqtl.junk;

import mbqtl.stat.RankArray;

public class RankTest {

	public static void main(String[] args) {
//        double[] data = new double[]{Double.NaN, Double.NaN, Double.NaN};
//        RankArray ranker = new RankArray();
//        double[] d = ranker.rank(data, true); // does this work with NaNs?
//        for (int q = 0; q < data.length; q++) {
//            System.out.println(q + "\t" + data[q] + "\t" + d[q]);
//        }
//
//        String gene = "ENSG00.1|ENSGE00.1";
//        gene = gene.replaceAll("\\|", "_");
//        gene = gene.replaceAll("\\.", "-");
//        System.out.println(gene);
//
//

		RankArray ranker = new RankArray();
		double[] data = new double[]{0.1, 0.3, 0.3, 0.4, 0.2};
		double[] datanan = new double[]{0.1, 0.3, 0.3, 0.4, Double.NaN};
		double[] dr = ranker.rank(data, true); // does this work with NaNs? answer: no
		double[] drnan = ranker.rank(datanan, true); // does this work with NaNs? answer: no
		for (int q = 0; q < data.length; q++) {
			System.out.println(q + "\t" + data[q] + "\t" + dr[q] + "\t" + drnan[q]);
		}
	}
}
