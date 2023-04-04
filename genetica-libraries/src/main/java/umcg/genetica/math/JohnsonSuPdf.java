/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

/**
 *
 * @author patri
 */
public class JohnsonSuPdf {

	private static final double SQRT2 = FastMath.sqrt(2.0);
	
	public static double cumulativeProbability(final double x, final double mu, final double sigma, final double nu, final double tau) {
		if (sigma < 0) {
			throw new IllegalArgumentException("sigma must be positive");
		}
		if (tau < 0) {
			throw new IllegalArgumentException("tau must be positive");
		}
		
		final double rTau = 1 / tau;
		final double w = rTau < 1e-07 ? 1 : Math.exp(Math.pow(rTau, 2));
		final double omega = -nu * rTau;
		final double c = Math.pow(0.5 * (w - 1) * (w * Math.cosh(2 * omega) + 1), -0.5);
		final double z = (x - (mu + c * sigma * Math.sqrt(w) * Math.sinh(omega))) / (c * sigma);
		final double r = -nu + Math.log(z + Math.sqrt(1 + Math.pow(z, 2))) / rTau;
		
        if (FastMath.abs(r) > 40) {
            return r < 0 ? 0.0d : 1.0d;
        }
        return 0.5 * Erf.erfc(-r / SQRT2);
		
	}

}
