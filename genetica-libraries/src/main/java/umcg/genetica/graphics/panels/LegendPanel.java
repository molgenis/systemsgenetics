package umcg.genetica.graphics.panels;

import umcg.genetica.graphics.DefaultGraphics;

import java.awt.*;
import java.text.DecimalFormat;

public class LegendPanel extends Panel {

    private double min;
    private double max;
    private double avg;
    private int n;


    public LegendPanel(int nrRows, int nrCols) {
        super(nrRows, nrCols);
    }


    private String title;

    public void setTitle(String title) {
        this.title = title;
    }

    public void setVals(double min, double max) {
        this.min = min;
        this.max = max;
    }

    public void setAvg(double avg, int n) {
        this.avg = avg;
        this.n = n;
    }

    @Override
    public void draw(DefaultGraphics g) {

        Graphics2D g2d = g.getG2d();
        int nrsteps = 10;
        int boxsize = 20;
        if (title != null) {
            g2d.setColor(theme.getDarkGrey());
            g2d.setFont(theme.getMediumFont());
            g2d.drawString(title, x0 + marginX, y0 + marginY - 10);
        }
        for (int i = 0; i < nrsteps + 1; i++) {
            int xPosYAxis = x0 + marginX;
            int yPosYAxis = y0 + marginY + ((nrsteps - i) * boxsize);
            double perc = (double) i / nrsteps;
            Color c = g.interpolateColor(theme.getColor(0), theme.getColor(1), perc);
            g2d.setColor(c);
            g2d.fillRect(xPosYAxis, yPosYAxis, boxsize, boxsize);

            g2d.setColor(theme.getDarkGrey());
            g2d.drawString("" + perc, xPosYAxis + boxsize + 10, yPosYAxis + boxsize);
        }


        int xPosYAxis = x0 + marginX;
        int yPosYAxis = y0 + marginY + ((nrsteps) * boxsize) + boxsize;
        DecimalFormat format = new DecimalFormat("#.###");
        g2d.drawString("Min: " + min, xPosYAxis + boxsize, yPosYAxis + boxsize + 10);
        g2d.drawString("Max: " + max, xPosYAxis + boxsize, yPosYAxis + boxsize + 20);
        g2d.drawString("Avg(%): " + format.format(avg) + " (n=" + n + ")", xPosYAxis + boxsize, yPosYAxis + boxsize + 30);


    }
}
