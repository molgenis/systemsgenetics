package umcg.genetica.graphics.themes;

import java.awt.*;

/**
 * Created by hwestra on 7/16/15.
 */
public class DefaultTheme implements Theme {

    public final Font LARGE_FONT = new Font("Helvetica", Font.PLAIN, 14);
    public final Font LARGE_FONT_BOLD = new Font("Helvetica", Font.BOLD, 14);
    public final Font MEDIUM_FONT = new Font("Helvetica", Font.PLAIN, 12);
    public final Font MEDIUM_FONT_BOLD = new Font("Helvetica", Font.BOLD, 12);
    public final Font SMALL_FONT = new Font("Helvetica", Font.PLAIN, 10);
    public final Font SMALL_FONT_BOLD = new Font("Helvetica", Font.BOLD, 10);

    public final Stroke strokeDashed = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
    public final Stroke stroke2pt = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
    public final Stroke stroke2ptDashed = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
    public final Stroke stroke = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

    private final Color darkgrey = new Color(70, 67, 58);
    private final Color lightgrey = new Color(174, 164, 140);

    private Color[] colors = new Color[]{
            new Color(70, 67, 58),
            new Color(208, 83, 77),
            new Color(98, 182, 177),
            new Color(116, 156, 80),
            new Color(124, 87, 147),
            new Color(174, 164, 140),
            new Color(185, 113, 65),
            new Color(63, 93, 126),
            new Color(79, 96, 64),
            new Color(109, 54, 96),
            new Color(110, 36, 30),
            new Color(32, 79, 74),
            new Color(81, 94, 28),
            new Color(67, 42, 82),
    };

    public void setColors(Color[] colors) {
        this.colors = colors;
    }

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
        return MEDIUM_FONT;
    }

    @Override
    public Font getLargeFontBold() {
        return LARGE_FONT_BOLD;
    }

    @Override
    public Font getMediumFontBold() {
        return MEDIUM_FONT_BOLD;
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
        Color c = getColor(i);
        int r = c.getRed();
        int g = c.getGreen();
        int b = c.getBlue();
        int a = (int) Math.floor((v * 255));
        return new Color(r, g, b, a);
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
