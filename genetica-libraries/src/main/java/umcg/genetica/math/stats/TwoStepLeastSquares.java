/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math.stats;

/**
 *
 * @author harmjan
 */
public class TwoStepLeastSquares extends Regression {

    // for IV analysis, this would be trans, cis, snp
    public static double[] tsls(double[] y, double[] x, double[] instruments) {
        // first do linear regression using instruments
        double[] coefficients = getLinearRegressionCoefficients(instruments, x);
        double[] predictedX = new double[y.length];
        for (int i = 0; i < predictedX.length; i++) {
            predictedX[i] = (coefficients[0] * instruments[i]) + coefficients[1];
        }

        double[] tslsCoefficients = getLinearRegressionCoefficients(predictedX, y);
        double[] output = new double[3];
        output[0] = tslsCoefficients[0];
        output[1] = tslsCoefficients[1];
        output[2] = tslsCoefficients[2];
        return output;
    }

    public static double[] tsls(double[] y, double[] x, double[] instruments, boolean b) {
        // first do linear regression using instruments
        double[] coefficients = getLinearRegressionCoefficients(instruments, x);
        double[] predictedX = new double[y.length];
        for (int i = 0; i < predictedX.length; i++) {
            predictedX[i] = (coefficients[0] * instruments[i]) + coefficients[1];
            System.out.println(y[i] + "\t" + x[i] + "\t" + instruments[i] + "\t" + predictedX[i]);
        }

        double[] tslsCoefficients = getLinearRegressionCoefficients(predictedX, y);
        double[] output = new double[3];
        output[0] = tslsCoefficients[0];
        output[1] = tslsCoefficients[1];
        output[2] = tslsCoefficients[2];
        return output;
    }
}
