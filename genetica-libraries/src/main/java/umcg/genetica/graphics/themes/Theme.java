package umcg.genetica.graphics.themes;

import java.awt.*;

/**
 * Created by hwestra on 7/16/15.
 */
public interface Theme {


	abstract Color getColor(int i);
	abstract Color getLightGrey();
	abstract Color getDarkGrey();
	abstract Font getLargeFont();
	abstract Font getMediumFont();
	abstract Font getLargeFontBold();
	abstract Font getMediumFontBold();
	abstract Font getSmallFont();
	abstract Font getSmallFontBold();
	abstract Stroke getStroke();
	abstract Stroke getStrokeDashed();
	abstract Stroke getThickStroke();
	abstract Stroke getThickStrokeDashed();

	abstract Color getColorSetOpacity(int i, float v);

	abstract Color getDarkerColor(Color color, double perc);


}
