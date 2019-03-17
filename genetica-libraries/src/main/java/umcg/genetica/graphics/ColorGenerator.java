package umcg.genetica.graphics;

import java.awt.*;

/**
 * Created by hwestra on 2/10/15.
 */
public class ColorGenerator {

	public static Color generate(double s, double v) {
		double golden_ratio_conjugate = 0.618033988749895;
		double h = Math.random();
		h += golden_ratio_conjugate;
		h %= 1;


		return Color.getHSBColor((float) h, (float) s, (float) v);
	}

	public static Color generate() {

		return generate(0.5, 0.95);

	}

}
