package umcg.genetica.graphics.themes;

import java.awt.*;

/**
 * Created by hwestra on 8/26/15.
 */
public class ComplementaryColor implements Theme {

	public final Font LARGE_FONT = new Font("Helvetica", Font.PLAIN, 14);
	public final Font LARGE_FONT_BOLD = new Font("Helvetica", Font.BOLD, 14);
	public final Font SMALL_FONT = new Font("Helvetica", Font.PLAIN, 10);
	public final Font SMALL_FONT_BOLD = new Font("Helvetica", Font.BOLD, 10);

	public final Stroke strokeDashed = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
	public final Stroke stroke2pt = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	public final Stroke stroke2ptDashed = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
	public final Stroke stroke = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

	private final Color darkgrey = new Color(70, 67, 58);
	private final Color lightgrey = new Color(174, 164, 140);

	private final Color[] colors = new Color[]{
			new Color(255, 51, 51),
			new Color(209, 42, 148),
			new Color(122, 36, 180),
			new Color(37, 52, 185),
			new Color(39, 121, 193),
			new Color(42, 177, 208),
			new Color(50, 189, 38),
			new Color(160, 211, 42),
			new Color(255, 255, 51),
			new Color(255, 211, 51),
			new Color(255, 167, 51),
			new Color(255, 129, 51),
			new Color(255, 51, 51),

	};


	@Override
	public Color getColor(int i) {
		return colors[i % colors.length];
	}

	@Override
	public Color getLightGrey() {
		return lightgrey;
	}

	@Override
	public Color getDarkGrey() {
		return darkgrey;
	}

	@Override
	public Font getLargeFont() {
		return LARGE_FONT;
	}

	@Override
	public Font getMediumFont() {
		return null;
	}

	@Override
	public Font getLargeFontBold() {
		return LARGE_FONT_BOLD;
	}

	@Override
	public Font getMediumFontBold() {
		return null;
	}

	@Override
	public Font getSmallFont() {
		return SMALL_FONT;
	}

	@Override
	public Font getSmallFontBold() {
		return SMALL_FONT_BOLD;
	}

	@Override
	public Stroke getStroke() {
		return stroke;
	}

	@Override
	public Stroke getStrokeDashed() {
		return strokeDashed;
	}

	@Override
	public Stroke getThickStroke() {
		return stroke2pt;
	}

	@Override
	public Stroke getThickStrokeDashed() {
		return stroke2ptDashed;
	}

	@Override
	public Color getColorSetOpacity(int i, float v) {
		Color c = colors[i];
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
		return new Color(r, g, b, (v * 255));
	}

	@Override
	public Color getDarkerColor(Color color, double perc) {
		double delta = (1 - perc);
		int r = (int) Math.ceil(color.getRed() * delta);
		int g = (int) Math.ceil(color.getGreen() * delta);
		int b = (int) Math.ceil(color.getBlue() * delta);
		return new Color(r, g, b);

	}
}
