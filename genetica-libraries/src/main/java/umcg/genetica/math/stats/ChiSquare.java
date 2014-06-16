/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

import cern.jet.stat.tdouble.Probability;


/**
 *
 * @author harmjan
 */
public class ChiSquare {

    public static double getX(int a, int b, int c, int d) {
        return getX((double) a, (double) b, (double) c, (double) d);
    }

    public static double getX(double a, double b, double c, double d) {
        //     real  perm  total
        //in   a     c     ac
        //out  b     d     bd
        //     ab    cd    n
        double ac = a + c;
        double bd = b + d;
        double ab = a + b;
        double cd = c + d;
        double n = a + b + c + d;

        double expa = (ab * ac) / n;
        double expb = (bd * ab) / n;
        double expc = (ac * cd) / n;
        double expd = (bd * cd) / n;

        double aminexp = ((a - expa) * (a - expa)) / expa;
        double bminexp = ((b - expb) * (b - expb)) / expb;
        double cminexp = ((c - expc) * (c - expc)) / expc;
        double dminexp = ((d - expd) * (d - expd)) / expd;

        double xsq = aminexp + bminexp + cminexp + dminexp;
//        System.out.println(a + "\t" + b + "\t" + c + "\t" + d + "\t" + aminexp / expa + "\t" + bminexp / expb + "\t" + cminexp / expc + "\t" + dminexp / expd + "\t" + n + "\t" + xsq);

//        System.out.println(a + "\t" + b + "\t" + c + "\t" + d + "\t" + aminexp / expa + "\t" + bminexp / expb + "\t" + cminexp / expc + "\t" + dminexp / expd + "\t" + n + "\t" + xsq);


        return xsq;
    }

    public static double getP(int df, double x) {
        double p = Probability.chiSquareComplemented(df, x);
        return p;
    }
}
