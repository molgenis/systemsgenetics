package umcg.genetica.graphics.panels;

import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class DensityGraphPanel extends Panel {

    public DensityGraphPanel(int nrRows, int nrCols) {
        super(nrRows, nrCols);
    }

    private Range dataRange;
    private String xAxisLabel;
    private String yAxisLabel;
    private String[] datasetLabels;
    private boolean plotLinearRegression = false;
    private boolean plothistogram = false;
    private boolean plotAxisTickLabels;
    private boolean plotLegend;
    private float alpha = 1.0f;

    double[] x;
    double[] y;

    @Override
    public void draw(DefaultGraphics defaultGraphics) {


        Graphics2D g2d = defaultGraphics.getG2d();

        // determine range
        if (dataRange == null) {
            // determine max and min using all data

            dataRange = new Range(Double.MAX_VALUE,
                    Double.MAX_VALUE,
                    -Double.MAX_VALUE,
                    -Double.MAX_VALUE);

            Range r = new Range(x, y);

            if (r.getMaxX() > dataRange.getMaxX()) {
                dataRange = new Range(dataRange.getMinX(),
                        dataRange.getMinY(),
                        r.getMaxX(),
                        dataRange.getMaxY());
            }
            if (r.getMaxY() > dataRange.getMaxY()) {
                dataRange = new Range(dataRange.getMinX(),
                        dataRange.getMinY(),
                        dataRange.getMaxX(),
                        r.getMaxY());
            }
            if (r.getMinX() < dataRange.getMinX()) {
                dataRange = new Range(r.getMinX(),
                        dataRange.getMinY(),
                        dataRange.getMaxX(),
                        dataRange.getMaxY());
            }
            if (r.getMinY() < dataRange.getMinY()) {
                dataRange = new Range(dataRange.getMinX(),
                        r.getMinY(),
                        dataRange.getMaxX(),
                        dataRange.getMaxY());
            }

            System.out.println(dataRange.toString());

            dataRange.round();
            System.out.println(dataRange.toString());
        }

        Range plotRange = dataRange;


        int nrPixelsMaxX = width - (2 * marginX);
        int nrPixelsMaxY = height - (2 * marginY);


        // axis labels
        // plot y-axis
        g2d.setColor(theme.getDarkGrey());
        double tickUnitY = plotRange.getRangeY() / 10;
        String pattern = "###,###,###.##";
        DecimalFormat decimalFormat = new DecimalFormat(pattern);

        g2d.setFont(theme.getMediumFont());
        FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

        int xPosYAxis = x0 + marginX - 10;
        int yPosYAxis = y0 + marginY;
        g2d.drawLine(xPosYAxis, yPosYAxis, xPosYAxis, yPosYAxis + nrPixelsMaxY);

        int maxlen = 0;
        for (double y = plotRange.getMinY(); y < plotRange.getMaxY() + (tickUnitY / 2); y += tickUnitY) {
            double yPerc = plotRange.getRelativePositionY(y);

            int ypos = y0 + marginY + (int) Math.ceil((1 - yPerc) * nrPixelsMaxY);
            int startx = xPosYAxis - 5;
            int stopx = startx + 5;
            g2d.drawLine(startx, ypos, stopx, ypos);
            String formattedStr = decimalFormat.format(y);
            int adv = metrics.stringWidth(formattedStr);
            if (adv > maxlen) {
                maxlen = adv;
            }
            g2d.setFont(theme.getMediumFont());
            if (plotAxisTickLabels) {
                if (formattedStr.equals("-0")) {
                    formattedStr = "0";
                }
                g2d.drawString(formattedStr, startx - adv - 5, ypos);
            }

        }

        if (yAxisLabel != null) {
            g2d.setFont(theme.getLargeFont());
            // determine middle of axis
            int middle = yPosYAxis + (nrPixelsMaxY / 2);
            int lengthOfAxisStr = metrics.stringWidth(yAxisLabel);
            int halfLength = lengthOfAxisStr / 2;
            int drawy = middle + halfLength;
            int drawx = xPosYAxis - maxlen - 20;

            drawRotate(g2d, drawx, drawy, -90, yAxisLabel);
        }

        // X axis
        // plot x-axis
        int yPosXAxis = y0 + marginY + nrPixelsMaxY + 10;

        int xPosXAxis = x0 + marginX;
        g2d.drawLine(xPosXAxis, yPosXAxis, xPosXAxis + nrPixelsMaxX, yPosXAxis);
        double tickUnitX = plotRange.getRangeX() / 10;

        for (double x = plotRange.getMinX(); x < plotRange.getMaxX() + (tickUnitX / 2); x += tickUnitX) {
            double xPerc = plotRange.getRelativePositionX(x);
            int xpos = xPosXAxis + (int) Math.ceil(xPerc * nrPixelsMaxX);
            int starty = yPosXAxis;
            int stopy = starty + 5;

            String formattedStr = decimalFormat.format(x);
            if (formattedStr.equals("-0")) {
                formattedStr = "0";
            }
            g2d.drawLine(xpos, starty, xpos, stopy);
            int adv = metrics.stringWidth(formattedStr);
            g2d.setFont(theme.getMediumFont());
            if (plotAxisTickLabels) {
                g2d.drawString(formattedStr, xpos - (adv / 2), stopy + 15);
            }

        }

        if (xAxisLabel != null) {
            g2d.setFont(theme.getLargeFont());
            // determine middle of axis
            int middle = xPosXAxis + (nrPixelsMaxX / 2);
            int lengthOfAxisStr = metrics.stringWidth(xAxisLabel);
            int halfLength = lengthOfAxisStr / 2;
            int drawx = middle - halfLength;
            int drawy = y0 + marginY + nrPixelsMaxY + 10 + (metrics.getHeight() * 2) + 10;
            g2d.drawString(xAxisLabel, drawx, drawy);
        }

        // draw title
        if (title != null) {
            int titlePosX = x0 + marginX;
            int titlePosY = y0 + marginY - theme.getLargeFont().getSize() - 5;
            g2d.setColor(theme.getDarkGrey());
            g2d.drawString(title, titlePosX, titlePosY);
        }

        // plot cross X
        g2d.setStroke(theme.getStrokeDashed());
        g2d.setColor(theme.getLightGrey());
        if (dataRange.getMinX() < 0 && dataRange.getMaxX() > 0) {
            double perc = dataRange.getRelativePositionX(0);
            int pixelX = x0 + marginX + (int) Math.ceil(nrPixelsMaxX * perc);
            int ystart = y0 + marginY + nrPixelsMaxY;
            int ystop = y0 + marginY;
            g2d.drawLine(pixelX, ystart, pixelX, ystop);
        }

        // plot cross Y
        if (dataRange.getMinY() < 0 && dataRange.getMaxY() > 0) {
            double perc = dataRange.getRelativePositionY(0);
            int pixelY = y0 + marginY + (int) Math.ceil(nrPixelsMaxY * perc);
            int xstart = x0 + marginX + nrPixelsMaxX;
            int xstop = x0 + marginX;
            g2d.drawLine(xstart, pixelY, xstop, pixelY);
        }

        // calculate the density profile
        int[][] grid = new int[nrPixelsMaxX][nrPixelsMaxY];
        g2d.setStroke(theme.getStroke());

        Color c = theme.getColorSetOpacity(0, alpha);

        g2d.setColor(c);
        long maxpix = 0;

        for (int j = 0; j < x.length; j++) {
            double xval = x[j];
            double yval = y[j];

            if (!Double.isNaN(xval) && !Double.isNaN(yval)) {
                double xperc = dataRange.getRelativePositionX(xval);
                double yperc = dataRange.getRelativePositionY(yval);

                // if yperc==0; this is actually the bottom of the grid (cg logic)


                int pixelX = (int) Math.ceil(nrPixelsMaxX * xperc);
                int pixelY = nrPixelsMaxY - (int) Math.ceil(nrPixelsMaxY * yperc);


                if (pixelX < 0) {
                    pixelX = 0;
                }
                if (pixelY < 0) {
                    pixelY = 0;
                }

                if (pixelY > nrPixelsMaxY - 1) {
                    pixelY = nrPixelsMaxY - 1;
                }

                if (pixelX > nrPixelsMaxX - 1) {
                    pixelX = nrPixelsMaxX - 1;
                }

                grid[pixelX][pixelY]++;
                int density = grid[pixelX][pixelY];
                if (density > maxpix) {
                    maxpix = density;

                }
            } else {
                // plot some markers somewhere?
            }
        }

        System.out.println("MAX: " + maxpix);


        ArrayList<Color> colortmp = new ArrayList<>();
        for (int i = 200; i > 0; i -= 25) {
            System.out.println(i);
            colortmp.add(new Color(i, i, i));
        }

        Color[] colors = colortmp.toArray(new Color[0]);
//                new Color[]{
//                new Color(37, 52, 185, alpha),
//                new Color(39, 121, 193, alpha),
//                new Color(42, 177, 208, alpha),
//                new Color(50, 189, 38, alpha),
//                new Color(160, 211, 42, alpha),
//                new Color(255, 255, 51, alpha),
//                new Color(255, 211, 51, alpha),
//                new Color(255, 167, 51, alpha),
//                new Color(255, 129, 51, alpha),
//                new Color(255, 51, 51, alpha)
//        };

        long sum = 0;
        for (int i = 0; i < grid.length; i++) {
            int xpos = x0 + i + marginX;
            for (int j = 0; j < grid[i].length; j++) {
                double val = grid[i][j];
                sum += val;
            }
        }
        double maxperc = (double) maxpix / sum;


        for (int i = 0; i < grid.length; i++) {
            int xpos = x0 + i + marginX;
            for (int j = 0; j < grid[i].length; j++) {
                double val = grid[i][j];
                int ypos = y0 + j + marginY;
                if (val > 0) {
//					double intensity = val / maxpix;

                    double perc = val / sum;
                    double percadj = perc / maxperc;

                    int colorbin = (int) Math.floor(percadj * colors.length);

                    if (colorbin < 0) {
                        colorbin = 0;
                    }
                    if (colorbin > colors.length - 1) {
                        colorbin = colors.length - 1;
                    }
                    Color col = colors[colorbin];
//					System.out.println(perc + "\t" + percadj + "\t" + colorbin);


//					Color col = defaultGraphics.interpolateColor(c1, c2, intensity);
                    g2d.setColor(col);
                    g2d.drawRect(xpos, ypos, 1, 1);


                }
            }
        }


        if (plotLinearRegression) {
            g2d.setColor(theme.getColor(0));
            double[] xval = x;
            double[] yval = y;

            Pair<double[], double[]> filtered = removeNulls(xval, yval);
            xval = filtered.getLeft();
            yval = filtered.getRight();

            // y = a + bx
            double[] params = Regression.getLinearRegressionCoefficients(xval, yval);
            double beta = params[0];
            double alpha = params[1];

            double minX = dataRange.getMinX();
            double maxX = dataRange.getMaxX();

            double y1 = alpha + (minX * beta);
            double y2 = alpha + (maxX * beta);


            double x1perc = dataRange.getRelativePositionX(minX);
            double y1perc = dataRange.getRelativePositionY(y1);
            double x2perc = dataRange.getRelativePositionX(maxX);
            double y2perc = dataRange.getRelativePositionY(y2);


            int pixelX1 = x0 + marginX + (int) Math.ceil(nrPixelsMaxX * x1perc);
            int pixelY1 = y0 + marginY + (int) Math.ceil(nrPixelsMaxY - (nrPixelsMaxY * y1perc));
            int pixelX2 = x0 + marginX + (int) Math.ceil(nrPixelsMaxX * x2perc);
            int pixelY2 = y0 + marginY + (int) Math.ceil(nrPixelsMaxY - (nrPixelsMaxY * y2perc));

            g2d.drawLine(pixelX1, pixelY1, pixelX2, pixelY2);

            // for shits and giggles, calculate the R-squared using the correlation coefficient
            double corr = Correlation.correlate(xval, yval);
            double rsq = corr * corr;

            DecimalFormat format = new DecimalFormat("#.##");

            g2d.drawString("" + format.format(rsq), pixelX2, pixelY2);

        }

        if (plotLegend && datasetLabels != null) {
            // plot a legend

            int maxStrWidth = 0;
            for (int i = 0; i < datasetLabels.length; i++) {
                String label = datasetLabels[i];
                int strwdth = getStringWidth(g2d, label);
                if (strwdth > maxStrWidth) {
                    maxStrWidth = strwdth;
                }
            }

            for (int i = 0; i < datasetLabels.length; i++) {
                g2d.setColor(theme.getColor(i));
                int pixelX = x0 + marginX + nrPixelsMaxX - maxStrWidth - 15;
                int pixelY = y0 + marginY;

                g2d.fillRect(pixelX, pixelY + (i * 15), 10, 10);
                g2d.drawString(datasetLabels[i], pixelX + 15, pixelY + (i * 15) + 10);

            }

        }

    }

    private Pair<double[], double[]> removeNulls(double[] xval, double[] yval) {
        boolean[] isnull = new boolean[xval.length];
        int ctrnotnull = 0;
        for (int i = 0; i < xval.length; i++) {
            if (Double.isNaN(xval[i]) || Double.isNaN(yval[i])) {
                isnull[i] = true;
            } else {
                ctrnotnull++;
            }
        }

        double[] outx = new double[ctrnotnull];
        double[] outy = new double[ctrnotnull];
        ctrnotnull = 0;
        for (int i = 0; i < xval.length; i++) {
            if (Double.isNaN(xval[i]) || Double.isNaN(yval[i])) {

            } else {
                outx[ctrnotnull] = xval[i];
                outy[ctrnotnull] = yval[i];
                ctrnotnull++;
            }
        }
        return new Pair<double[], double[]>(outx, outy);
    }

    public void setData(ArrayList<Double> x, ArrayList<Double> y) {
        this.x = Primitives.toPrimitiveArr(x);
        this.y = Primitives.toPrimitiveArr(y);
    }

    public void setPlotElems(boolean tickmarks, boolean legend) {
        this.plotAxisTickLabels = tickmarks;
        this.plotLegend = legend;
    }

    public void setLabels(String xlabel, String ylabel) {
        this.xAxisLabel = xlabel;
        this.yAxisLabel = ylabel;
    }

    public void setDataRange(Range range) {
        this.dataRange = range;
    }

    public void setData(double[] xarr, double[] yarr) {
        this.x = xarr;
        this.y = yarr;
    }

    public void setAlpha(float v) {

        this.alpha = v;
    }
}
