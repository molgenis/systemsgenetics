package mbqtl.junk;

import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

import java.text.DecimalFormat;

public class PTest {
    public static void main(String[] args) {
//        double spearman = -0.1396061331821292;
//        int nrSamples = 552;
//        double t = spearman / Math.sqrt((1.0D - spearman * spearman) / (double) (nrSamples - 2));
//        StudentT tDistColt = new StudentT((double) (nrSamples - 2), new DRand());
//        double p = tDistColt.cdf(t) * 2;
//        double z = ZScores.pToZ(p);

        DecimalFormat decimalFormat = new DecimalFormat("#.######");

        Correlation.correlationToZScore(554);

        for (double d = 0; d < 1; d += 0.01) {
            double corr = d;
            double df = 554;
            double t = corr / Math.sqrt((1.0D - corr * corr) / df);
            StudentT tDistColt = new StudentT((double) (df), new DRand());
            double pt = tDistColt.cdf(-t) * 2;
            JSci.maths.statistics.FDistribution dist = new JSci.maths.statistics.FDistribution(1, df);
            double fstat = df * corr * corr / (1 - corr * corr);
            double pf = dist.cumulative(fstat);
            double z = Correlation.convertCorrelationToZScore((int) df, corr);
            double pz = ZScores.zToP(z);
            System.out.println(d + "\t\t" + pf + "\t\t" + (1 - pf) + "\t\t" + pt + "\t\t" + pz);
        }


    }
}
