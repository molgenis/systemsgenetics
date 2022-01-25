package betaqtl.junk;

import org.apache.commons.math3.distribution.FDistribution;
import umontreal.iro.lecuyer.probdist.BetaDist;

public class BetaDistTest {

    public static void main(String[] args) {
        // ENSG00000096746.17
        // 5650
        // 1.05381
        // 330.113
        // 464.888 10:68246902:rs61854799:T_C -84273 0.0060772 0.223284 0.973735 0.977862
        double alpha = 1.05888;
        double beta = 1055.37;
        double pval = 0.0489563;
        BetaDist bdist = new BetaDist(alpha, beta);
        double pperm = bdist.cdf(pval);
        double pperm2 = bdist.inverseF(pval);
        System.out.println(pperm);
        System.out.println(pperm2);

        ;
        // 0.218019
        // 0.11515712908223462

        System.out.println(Double.parseDouble(""));


    }
}
