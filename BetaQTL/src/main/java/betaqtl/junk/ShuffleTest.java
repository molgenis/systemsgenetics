package betaqtl.junk;

import betaqtl.Util;

public class ShuffleTest {


	public static void main(String[] args) {

		double[] v = new double[100];
		for (int q = 0; q < v.length; q++) {
			v[q] = q;
		}

		double[][] z = new double[2][100];
		for (int p = 0; p < 2; p++) {
			double[] vcopy = new double[100];
			System.arraycopy(v, 0, vcopy, 0, vcopy.length);
			Util.shuffleArray(vcopy, p);
			z[p] = vcopy;
		}

		for (int p = 0; p < 100; p++) {
			System.out.println(p + "\t" + z[0][p] + "\t" + z[1][p]);
		}

	}

}
