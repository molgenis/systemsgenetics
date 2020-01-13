package umcg.genetica.graphics.themes;


import java.awt.*;

public class ColorBlindTheme implements Theme {
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

    private final Color[] colors = new Color[]{
            new Color(70, 67, 58),
            new Color(109, 136, 196),
            new Color(220, 37, 127),
            new Color(242, 99, 34),
            new Color(108, 99, 172),
            new Color(252, 175, 23),

            new Color(152, 149, 133),
            new Color(159, 188, 255),
            new Color(242, 150, 100),
            new Color(108, 150, 172),
            new Color(252, 175, 23),
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
        if (i > colors.length - 1) {
            i = i % colors.length;
        }
        Color c = colors[i];
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
