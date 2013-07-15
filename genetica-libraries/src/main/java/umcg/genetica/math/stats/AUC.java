/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import java.util.Arrays;

/**
 *
 * @author juha
 */
public class AUC {

    public static double getAUC(double[] a1, double[] a2) {
        double[] totalList = new double[a1.length + a2.length];
        System.arraycopy(a1, 0, totalList, 0, a1.length);
        System.arraycopy(a2, 0, totalList, a1.length, a2.length);
        Arrays.sort(totalList);
        int r1 = 0;
        int place = 0;
        for (int i = 0; i < a1.length; i++) {
            for (int j = place; j < totalList.length; j++) {
                if (a1[i] == totalList[j]) {
                    r1 += (j + 1);
                    place = j + 1;
                    break;
                }
            }
        }
        double uA = r1 - a1.length * (a1.length + 1.0d) / 2.0d;
        double auc = uA / (a1.length * a2.length);
        return auc;
    }
}
