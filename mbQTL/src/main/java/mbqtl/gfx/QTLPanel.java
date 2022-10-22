package mbqtl.gfx;

import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.Panel;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.graphics.themes.Theme;
import umcg.genetica.math.stats.Regression;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class QTLPanel extends Panel {


    double[] x = null;
    double[] y = null;
    private Range dataRange;

    private int axisMargin = 10;
    private int intMargin = 50;
    private String[] alleles;
    private String datasetname;
    private String genename;
    private String rsid;
    private double z = Double.NaN;
    private double r = Double.NaN;
    private double p = Double.NaN;
    private boolean notTested = false;

    public QTLPanel(int nrRows, int nrCols) {
        super(nrRows, nrCols);
    }

    @Override
    public void draw(DefaultGraphics g) {
        Graphics2D g2d = g.getG2d();

        Color BLACK = new Color(0, 0, 0);
        Color LIGHTGRAY = new Color(200, 200, 200);
        Color DARKGRAY = new Color(100, 100, 100);

        // set some default colors, stroke, etc
        Theme t = new DefaultTheme();
        g2d.setStroke(t.getStroke());
        g2d.setColor(DARKGRAY);

        // draw dataset name
        if (datasetname != null) {
            g2d.setFont(t.getLargeFont());
            int x0p = x0 + axisMargin + 5;
            int y0p = y0 + axisMargin;
            g2d.drawString(datasetname, x0p, y0p);
        }

        if (notTested || x == null || y == null) {
            g2d.setColor(LIGHTGRAY);
            g2d.setFont(t.getSmallFont());
            int x0p = x0 + axisMargin + 5;
            int y0p = y0 + axisMargin + 10;
            g2d.drawString("Effect not tested in dataset", x0p, y0p);
        } else {

            int x0p = x0 + axisMargin;
            int y0p = y0 + axisMargin;

            int axisWidth = width - (2 * axisMargin);
            int axisHeight = height - (2 * axisMargin);

            int pWidth = width - (2 * intMargin);
            int pHeight = height - (2 * intMargin);


            double[][] filteredxy = stripMissing(x, y);
            x = filteredxy[0];
            y = filteredxy[1];

            int sumStatOffset = 10;
            if (rsid != null) {
                g2d.setFont(t.getSmallFont());
                // draw statistics, if defined
                g2d.setColor(LIGHTGRAY);
                int xp = x0 + axisMargin + 5;
                int yp = y0 + axisMargin + sumStatOffset;
                g2d.drawString(rsid, xp, yp);
                sumStatOffset += 10;
            }

            if (!Double.isNaN(r) || !Double.isNaN(z) || !Double.isNaN(p)) {
                g2d.setFont(t.getSmallFont());
                // draw statistics, if defined
                g2d.setColor(LIGHTGRAY);
                DecimalFormat f = new DecimalFormat("#.###");
                DecimalFormat pf = new DecimalFormat("0.0E0");
                String sumstatStr = "";
                if (!Double.isNaN(r)) {
                    sumstatStr += "r=" + f.format(r);
                }
                if (!Double.isNaN(z)) {
                    if (sumstatStr.length() > 0) {
                        sumstatStr += "; z=" + f.format(z);
                    } else {
                        sumstatStr += "z=" + f.format(z);
                    }
                }
                if (!Double.isNaN(p)) {
                    String pstr = "";
                    if (p > 0.01) {
                        pstr = f.format(p);
                    } else {
                        pstr = pf.format(p);
                    }
                    if (sumstatStr.length() > 0) {
                        sumstatStr += "; p=" + pstr;
                    } else {
                        sumstatStr += "p=" + pstr;
                    }
                }
                int xp = x0 + axisMargin + 5;
                int yp = y0 + axisMargin + sumStatOffset;
                g2d.drawString(sumstatStr, xp, yp);
            }

            g2d.setColor(DARKGRAY);
            dataRange = new Range(0, Primitives.min(y), 2, Primitives.max(y));


            // draw axes


            g2d.drawLine(x0p, y0p, x0p, y0p + axisHeight);
            g2d.drawLine(x0p, y0p + axisHeight, x0p + axisWidth, y0p + axisHeight);

            // draw y-axis labels
            int nrTicks = 2;
            double rangeY = dataRange.getRangeY();
            double unitY = dataRange.determineUnit(rangeY);
            double minYRemainder = dataRange.getMinY() % unitY;
            double maxYRemainder = dataRange.getMaxY() % unitY;
            double minYRound = Math.floor(dataRange.getMinY() - minYRemainder);
            double maxYRound = Math.floor(unitY + dataRange.getMaxY() - maxYRemainder);

            // this may need a fix later on
            if (rangeY / unitY < nrTicks) {
                unitY /= 4;
            }

            DecimalFormat format = new DecimalFormat("#.##");
            g2d.setColor(LIGHTGRAY);
            int maxTickLabelStrWidth = 0;
            for (double i = minYRound; i < maxYRound; i += unitY) {
                if (i >= dataRange.getMinY() && i <= dataRange.getMaxY()) {
                    double percY = dataRange.getRelativePositionY(i);
                    int intY = y0 + intMargin + pHeight - (int) Math.floor(percY * pHeight);
                    int intX1 = x0p;
                    int intX2 = x0p - 5;
                    g2d.setColor(DARKGRAY);
                    g2d.drawLine(intX1, intY, intX2, intY);
                    g2d.setFont(theme.getSmallFont());
                    String valStr = format.format(i);
                    int strw = this.getStringWidth(g2d, valStr);
                    if (strw > maxTickLabelStrWidth) {
                        maxTickLabelStrWidth = strw;
                    }
                    g2d.drawString(valStr, x0p - 5 - strw, intY + 4);
                }
            }

            // plot gene name, if any
            if (genename != null) {
                g2d.setColor(DARKGRAY);
                int xp = x0p - 5 - maxTickLabelStrWidth - 5;

                g2d.setFont(t.getSmallFont());
                int angle = -90;

                int strw = getStringWidth(g2d, genename);
                int yp = y0 + axisMargin + (axisHeight / 2) + (strw / 2);
                drawRotate(g2d, xp, yp, angle, genename);
            }


            // get regression line
            double[] coeff = Regression.getLinearRegressionCoefficients(x, y); // beta, alpha, se, t
            double xlmin = -0.1;
            double xlmax = 2.1;
            double yl1 = coeff[1] + (coeff[0] * xlmin);
            double yl2 = coeff[1] + (coeff[0] * xlmax);
            double yl1perc = dataRange.getRelativePositionY(yl1);
            double yl2perc = dataRange.getRelativePositionY(yl2);
            double xl1perc = dataRange.getRelativePositionX(xlmin);
            double xl2perc = dataRange.getRelativePositionX(xlmax);
            int yl1posY = y0 + intMargin + pHeight - (int) Math.floor(yl1perc * pHeight);
            int yl2posY = y0 + intMargin + pHeight - (int) Math.floor(yl2perc * pHeight);
            int xl1posX = x0 + intMargin + (int) Math.floor(xl1perc * pWidth);
            int xl2posX = x0 + intMargin + (int) Math.floor(xl2perc * pWidth);
            g2d.setColor(BLACK);
            g2d.drawLine(xl1posX, yl1posY, xl2posX, yl2posY);

            g2d.setColor(LIGHTGRAY);

            // group data by genotype labels
            ArrayList<Double> xAA = new ArrayList<>();
            ArrayList<Double> xAB = new ArrayList<>();
            ArrayList<Double> xBB = new ArrayList<>();
            ArrayList<Double> yAA = new ArrayList<>();
            ArrayList<Double> yAB = new ArrayList<>();
            ArrayList<Double> yBB = new ArrayList<>();
            for (int i = 0; i < x.length; i++) {
                if (x[i] < 0.5) {
                    xAA.add(x[i]);
                    yAA.add(y[i]);
                } else if (x[i] > 1.5) {
                    xBB.add(x[i]);
                    yBB.add(y[i]);
                } else {
                    xAB.add(x[i]);
                    yAB.add(y[i]);
                }
            }

            // draw x-axis labels
            g2d.setColor(LIGHTGRAY);
            Composite originalComposite = g2d.getComposite();
            for (int i = 0; i < 3; i++) {
                double percXTick = dataRange.getRelativePositionX(i);
                int intXTick = x0 + intMargin + (int) Math.floor(percXTick * pWidth);
                g2d.setColor(DARKGRAY);
                g2d.drawLine(intXTick, y0p + axisHeight, intXTick, y0p + axisHeight + 5);
                g2d.setFont(theme.getSmallFont());
                if (alleles == null) {
                    alleles = new String[2];
                    alleles[0] = "A";
                    alleles[1] = "B";
                }
                String alleleStr = alleles[0] + "/" + alleles[0];
                String sizeStr = "(n=" + xAA.size() + ")";
                if (i == 1) {
                    alleleStr = alleles[0] + "/" + alleles[1];
                    sizeStr = "(n=" + xAB.size() + ")";
                } else if (i == 2) {
                    alleleStr = alleles[1] + "/" + alleles[1];
                    sizeStr = "(n=" + xBB.size() + ")";
                }
                int strw = getStringWidth(g2d, alleleStr) / 2;
                // draw allleles
                g2d.setColor(DARKGRAY);
                g2d.drawString(alleleStr, intXTick - strw, y0p + axisHeight + 15);
                // draw sample size
                g2d.setColor(LIGHTGRAY);
                strw = getStringWidth(g2d, sizeStr) / 2;
                g2d.drawString(sizeStr, intXTick - strw, y0p + axisHeight + 25);

                // draw violin
                double[] vals = null;

                if (i == 0) {
                    vals = Primitives.toPrimitiveArr(yAA);

                } else if (i == 1) {
                    vals = Primitives.toPrimitiveArr(yAB);
                } else {
                    vals = Primitives.toPrimitiveArr(yBB);
                }

                int yViolin = y0 + intMargin;
                int widthViolin = 50 ; //(pWidth - (4 * intMargin)) / 3;
                int xViolin = intXTick - (widthViolin / 2);
                int heightViolin = pHeight;

                Color c = g2d.getColor();
                g2d.setColor(LIGHTGRAY);
                if (vals.length > 5) {
                    AlphaComposite alphaComposite10 = AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.2f);
                    g2d.setComposite(alphaComposite10);
                    drawViolinPlot(g2d, xViolin, yViolin, widthViolin, heightViolin, vals, dataRange.getMinY(), dataRange.getMaxY());
                    int widthBox = 24;
                    g2d.setColor(DARKGRAY);
                    drawBoxPlot(g2d, intXTick - 12, yViolin, widthBox, heightViolin, vals, dataRange.getMinY(), dataRange.getMaxY(), true);
                    AlphaComposite alphaComposite1 = AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 1f);
                    g2d.setComposite(alphaComposite1);
                }
            }
            g2d.setComposite(originalComposite);

            // plot data points
            g2d.setColor(DARKGRAY);
            for (int i = 0; i < x.length; i++) {
                double percX = dataRange.getRelativePositionX(x[i]);
                int intX = x0 + intMargin + (int) Math.floor(percX * pWidth);
                double percY = dataRange.getRelativePositionY(y[i]);
                int intY = y0 + intMargin + pHeight - (int) Math.floor(percY * pHeight);
                // System.out.println(i + "\t" + x[i] + "\t" + percX + "\t" + intX + "\t" + y[i] + "\t" + percY + "\t" + intY);
                g2d.fillOval(intX - 1, intY - 1, 2, 2);
            }

            System.out.println(dataRange.getMinX() + "\t" + dataRange.getMaxX() + "\t" + dataRange.getMinY() + "\t" + dataRange.getMaxY());


        }
    }

    private double[][] stripMissing(double[] x, double[] y) {
        int nonmissing = 0;
        for (int i = 0; i < x.length; i++) {
            if (!Double.isNaN(x[i]) && x[i] >= 0) {
                nonmissing++;
            }
        }
        double[][] output = new double[2][nonmissing];
        int ctr = 0;
        for (int i = 0; i < x.length; i++) {
            if (!Double.isNaN(x[i]) && x[i] >= 0) {
                output[0][ctr] = x[i];
                output[1][ctr] = y[i];
                ctr++;
            }
        }
        return output;
    }


    public void drawViolinPlot(Graphics2D g2d, int x, int y, int width, int height, double[] vals, double minValue, double maxValue) {

        int nrVals = vals.length;

        //Determine range of values:
        double minVals = Primitives.min(vals);
        double maxVals = Primitives.max(vals);

        //Make frequency distribution:
        int nrBins = 1 + (int) Math.round(Math.sqrt(nrVals) / 4d);
        int[] binCount = new int[nrBins];
        for (int n = 0; n < nrBins; n++) {
            double lower = minVals + (maxVals - minVals) * (double) n / (double) nrBins;
            double upper = minVals + (maxVals - minVals) * (double) (n + 1) / (double) nrBins;
            for (int v = 0; v < nrVals; v++) {
                if (vals[v] >= lower && vals[v] < upper) {
                    binCount[n]++;
                }
            }
        }

        //Smooth the distribution:
        int posYMin = y + height - (int) Math.round((double) height * (maxVals - minValue) / (maxValue - minValue));
        int posYMax = y + height - (int) Math.round((double) height * (minVals - minValue) / (maxValue - minValue));
        double[] posVal = new double[posYMax - posYMin + 1];
        for (int pos = posYMin; pos <= posYMax; pos++) {
            double value = (((-pos + y + height) * (maxValue - minValue)) / (double) height) + minValue;
            for (int n = 0; n < nrBins; n++) {
                double lower = minVals + (maxVals - minVals) * (double) n / (double) nrBins;
                double upper = minVals + (maxVals - minVals) * (double) (n + 1) / (double) nrBins;
                if (value >= lower && value < upper) {
                    posVal[pos - posYMin] = binCount[n];
                }

            }
        }
        double kernelWidth = (double) height / (double) nrBins / 4d;
        double[] kernelWeights = new double[201];
        for (int d = -100; d <= 100; d++) {
            double weight = Math.pow(Math.E, -((double) d / (double) kernelWidth) * ((double) d / (double) kernelWidth) / 2);
            kernelWeights[d + 100] = weight;
        }
        double[] posValSmoothed = new double[posYMax - posYMin + 1];
        for (int pos = posYMin; pos <= posYMax; pos++) {
            double valSmoothed = 0;
            double sumWeights = 0;
            for (int q = pos - 100; q <= pos + 100; q++) {
                if (q >= posYMin && q <= posYMax) {
                    sumWeights += kernelWeights[100 + q - pos];
                    valSmoothed += kernelWeights[100 + q - pos] * posVal[q - posYMin];
                }
            }
            posValSmoothed[pos - posYMin] = valSmoothed / sumWeights;
        }
        double maxSmoothedVal = Primitives.max(posValSmoothed);
        for (int pos = posYMin; pos <= posYMax; pos++) {
            posValSmoothed[pos - posYMin] /= maxSmoothedVal;
        }

        //Draw shape:
        GeneralPath path = new GeneralPath();
        path.moveTo(x + width / 2, posYMin);
        for (int pos = posYMin; pos <= posYMax; pos++) {
            path.lineTo(x + width / 2 - posValSmoothed[pos - posYMin] * width / 2d - 1, pos);
        }
        for (int pos = posYMax; pos >= posYMin; pos--) {
            path.lineTo(x + width / 2 + posValSmoothed[pos - posYMin] * width / 2d + 1, pos);
        }
        path.closePath();
        g2d.draw(path);
        g2d.fill(path);

    }

    public void drawBoxPlot(Graphics2D g2d, int x, int y, int width, int height, double[] vals, double minValue, double maxValue, boolean drawOutliers) {

        double median = JSci.maths.ArrayMath.percentile(vals, 0.50d);
        double q1 = JSci.maths.ArrayMath.percentile(vals, 0.25d);
        double q3 = JSci.maths.ArrayMath.percentile(vals, 0.75d);
        double iqr = q3 - q1;

        //Draw median:
        int posY = y + height - (int) Math.round((double) height * (median - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND));
        g2d.drawLine(x, posY, x + width, posY);
        //Draw IQR:
        int posY1 = y + height - (int) Math.round((double) height * (q3 - minValue) / (maxValue - minValue));
        int posY2 = y + height - (int) Math.round((double) height * (q1 - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND));
        g2d.drawRect(x, posY1, width, posY2 - posY1);

        //Draw whiskers:
        double whiskerTop = q3 + 1.5d * iqr;
        double whiskerBottom = q1 - 1.5d * iqr;
        double max = Primitives.max(vals);
        double min = Primitives.min(vals);
        if (min > whiskerBottom) {
            whiskerBottom = min;
        }
        if (max < whiskerTop) {
            whiskerTop = max;
        }
        posY = y + height - (int) Math.round((double) height * (whiskerTop - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, new float[]{2.0f}, 0.0f));
        g2d.drawLine(x + width / 2, posY, x + width / 2, posY1);
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND));
        g2d.drawLine(x + width / 2 - 5, posY, x + width / 2 + 5, posY);
        posY = y + height - (int) Math.round((double) height * (whiskerBottom - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, new float[]{2.0f}, 0.0f));
        g2d.drawLine(x + width / 2, posY2, x + width / 2, posY);
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND));
        g2d.drawLine(x + width / 2 - 5, posY, x + width / 2 + 5, posY);

        //Draw outliers:
        if (drawOutliers) {
            AlphaComposite alphaComposite10 = AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.10f);
            g2d.setComposite(alphaComposite10);
            for (int v = 0; v < vals.length; v++) {
                if (vals[v] > whiskerTop || vals[v] < whiskerBottom) {
                    posY = y + height - (int) Math.round((double) height * (vals[v] - minValue) / (maxValue - minValue));
                    g2d.drawOval(x + width / 2 - 2, posY - 2, 5, 5);
                }
            }
        }
    }

    public void setAlleles(String[] alleles) {
        this.alleles = alleles;
    }

    public void setData(double[] x, double[] y) {
        this.x = x;
        this.y = y;
    }

    private void determineDataRange() {
        this.dataRange = new Range(Double.MAX_VALUE, Double.MAX_VALUE, -1.7976931348623157E308, -1.7976931348623157E308);

        for (int i = 0; i < this.x.length; ++i) {
            if (this.x.length > 0) {
                Range r = new Range(this.x, this.y);
                r.round();
                if (r.getMaxX() > this.dataRange.getMaxX()) {
                    this.dataRange = new Range(this.dataRange.getMinX(), this.dataRange.getMinY(), r.getMaxX(), this.dataRange.getMaxY());
                }

                if (r.getMaxY() > this.dataRange.getMaxY()) {
                    this.dataRange = new Range(this.dataRange.getMinX(), this.dataRange.getMinY(), this.dataRange.getMaxX(), r.getMaxY());
                }

                if (r.getMinX() < this.dataRange.getMinX()) {
                    this.dataRange = new Range(r.getMinX(), this.dataRange.getMinY(), this.dataRange.getMaxX(), this.dataRange.getMaxY());
                }

                if (r.getMinY() < this.dataRange.getMinY()) {
                    this.dataRange = new Range(this.dataRange.getMinX(), r.getMinY(), this.dataRange.getMaxX(), this.dataRange.getMaxY());
                }
            }
        }

        this.dataRange.round();
    }

    public void setDatasetDetails(String datasetName, String geneName, String rsId, double z, double p, double r) {
        this.datasetname = datasetName;
        this.genename = geneName;
        this.rsid = rsId;
        this.z = z;
        this.p = p;
        this.r = r;
    }

    public void setNotTested() {
        notTested = true;
    }
}
