package mbqtl.junk;

import umcg.genetica.math.stats.ZScores;

public class ZTest {

    public static void main(String[] args) {
//        double p = 3.679676003373988E-33;
//        double z = ZScores.pToZ(p);
//        double p2 = ZScores.zToP(z);
//
//        NormalDist d = new NormalDist();
//        double p3 = d.cdf(z);
//        JSci.maths.statistics.NormalDistribution distribution = new JSci.maths.statistics.NormalDistribution(0, 1);
//        double p4 = distribution.cumulative(z);
//        System.out.println(p + "\t" + z + "\t" + p2 + "\t" + p3 + "\t" + p4);
//
//        System.out.println();
//        System.out.println();
//        System.out.println();
//
//        double zq = -11.9893839;
//        double zqp = ZScores.zToP(zq);
//        double zqpz = ZScores.pToZTwoTailed(zqp);
//        System.out.println(zq + "\t" + zqp + "\t" + zqpz + "\t" + (zq - zqpz) + "\t" + ZScores.zToP(zqpz));
//
//
//        System.out.println(Math.round(0.4));
//        System.out.println(Math.round(0.5));
//        System.out.println(Math.round(1.5));

        System.out.println(ZScores.zToP(-2.858));

    }
}
