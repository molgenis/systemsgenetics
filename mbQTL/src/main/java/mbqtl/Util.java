package mbqtl;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class Util {

    public static HashMap<String, Integer> hash(ArrayList<String> list) {
        HashMap<String, Integer> hash = new HashMap<>();
        for (int i = 0; i < list.size(); i++) {
            hash.put(list.get(i), i);
        }
        return hash;
    }

    public static void shuffleArray(double[] ar, long seed) {
        Random rnd = new Random(seed);
        for (int i = ar.length - 1; i >= 0; i--) {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            double a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }

    public static double[] centerScale(double[] data) {
        double mean = JSci.maths.ArrayMath.mean(data);
        double sd = JSci.maths.ArrayMath.standardDeviation(data);
        double[] output = new double[data.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = (data[i] - mean) / sd;
        }
        return output;
    }

    public static double meanGenotype(double[] dosagesForDataset) {
        double sum = 0;
        int ctr = 0;
        for (int v = 0; v < dosagesForDataset.length; v++) {
            double q = dosagesForDataset[v];
            if (q != -1 && !Double.isNaN(q)) {
                sum += dosagesForDataset[v];
                ctr++;
            }
        }
        return sum / ctr;
    }


}
