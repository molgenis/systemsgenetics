package nl.systemsgenetics.genenetworkbackend.legacy;

/**
 *
 * @author ludefranke
 */
public class WilcoxonMannWhitney {

    private double k, n, kn, auc;

    public double getAUC() {
        return auc;
    }

    private double abs(double value) {
        return Math.abs(value);
    }

    private double k_out_n (double k, double n) {
        this.k = k;
        this.n = n;
        kn = 1.0d;
        while (k>0) {
            n--;
            k--;
            kn*=n/k;
        }
        return kn;
    }

    private double countSmallerRanks(double w, double sum, int m, int start, double[] rankList) {

        int i = 0;
        double temp = 0;
        double smaller = 0;
        int end = 0;
        int mminus1 = 0;

        if (sum > w) return 0;

        end = rankList.length - 1;
        if (m>0) {
            mminus1 = m - 1;
            for (i = start; i <= end - m; i++) {
                temp = sum + rankList[i];
                if (temp > w) return smaller;
                smaller += countSmallerRanks(w, temp, mminus1, i+1, rankList);
            }
        } else {
            if (sum + end + 1 <= w) return end - start + 1;
            for (i = start; i <= end; i++) {
                temp = sum + rankList[i];
                if (temp <= w) {
                    smaller++;
                } else {
                    return smaller;
                }
            }
        }
        return smaller;
    }

    private double normalZ (double z) {
        double x = z;
	double[] b = {0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
	double p = 0.2316419;
	double t = 1/(1+p*x);
	double fact = t;
	double sum = 0;
	for(int i=0; i <= b.length - 1; i++) {
  		sum += b[i]*fact;
  		fact *= t;
	}
	p = 2d * sum * Math.exp(-x*x/2.0d) / (Math.sqrt(2d*Math.PI));
        return p;
    }

  /** Creates a new instance of WilcoxonMannWhitney */
    public WilcoxonMannWhitney() {
    }


    public double returnWilcoxonMannWhitneyPValue(double[] A, double[] B) {
        java.util.Arrays.sort(A);
        java.util.Arrays.sort(B);
        double[] totalList = new double[A.length + B.length];
        for (int x=0; x<A.length; x++) totalList[x] = A[x];
        for (int y=0; y<B.length; y++) totalList[y + A.length] = B[y];
        java.util.Arrays.sort(totalList);

        double nA = A.length;
        double nB = B.length;
        double n = nA + nB;
        double maxSum = n * (n + 1d) / 2d;
        double h0 = maxSum / 2d;

        double previous = Double.MIN_VALUE;
        int start = 0;
        double[] totalRank = new double[totalList.length];
        for (int i = 0; i<totalList.length; i++) {
            if (totalList[i]==previous) {
                double meanRank = ((double) start + (double) i + 2d) / 2d;
                for (int j = start; j<=i; j++) {
                    totalRank[j] = meanRank;
                    //System.out.println(j + "\t" + totalList[j] + "\t" + totalRank[j]);
                }
            } else {
                totalRank[i] = i + 1;
                previous = totalList[i];
                start = i;
                //System.out.println(i + "\t" + totalList[i] + "\t" + totalRank[i]);
            }
        }

        double[] shortest = A;
        if (B.length < A.length) shortest = B;
        double nShortest = shortest.length;
        double w = 0;
        for (int a=0; a<shortest.length; a++) {
            int i = 0;
            while (1 == 1) {
                if (i >= totalList.length - 1 || shortest[a] == totalList[i]) break;
                i++;
            }
            w+=totalRank[i];
        }

        double nZ = nShortest; if (w>h0) nZ = n - nShortest;
        if (w>h0) w = maxSum - w;

        double p = 0;

        //Calculate AUC:
        double r1 = 0; int place = 0;
        for (int i=0; i<A.length; i++) {
            for (int j=place; j<totalList.length; j++) {
                if (A[i]==totalList[j]) {
                    r1+=(j + 1);
                    place = j + 1;
                    break;
                }
            }
        }
        double uA = r1 - nA * (nA + 1.0d) / 2.0d;
        auc = uA / (nA * nB);

        double permutations = k_out_n (nA, n);
        if (permutations >= 25000 || shortest.length >= 10) {
            double continuity = 0.5; if (w>=h0) continuity = -0.5;
            double z = Math.abs((w + continuity - nZ * (n + 1d) / 2d) / Math.sqrt(nA * nB * (n + 1d) / 12d));
            //System.out.println(maxSum + "\t" + w + "\t" + continuity + "\t" + z);
            p = normalZ(z);
            return p;
        }

        /*
        if (shortest.length < 10 && p < 0.25  && permutations < 60000) {
            double less = countSmallerRanks(w, 0, shortest.length, 0, totalRank);
            if (2 * less > permutations) {
                less = countSmallerRanks (w - 1, 0, shortest.length, 0, totalRank);
                less = permutations - less;
            }
            double sumFrequencies = permutations;
            p = 2.0d * less / sumFrequencies;
        }

        System.out.println("accurate p-value:\t" + p);

        if (shortest.length < 10 && p >= 0.25 || permutations >= 60000) {
            System.out.println("Cannot accurately determine P-Value!");
        }
         */

        return -1;

    }

}

