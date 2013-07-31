/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.graphics;

import com.lowagie.text.DocumentException;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.GeneralPath;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Locale;
import javax.imageio.ImageIO;
import umcg.genetica.math.Fmath;

/**
 *
 * @author harmjan
 */
public class ViolinBoxPlot {

    public enum Output {

        PDF, PNG
    };

    private double getWidth(String text, java.awt.Font font) {
        java.awt.Graphics2D g2d = (new java.awt.image.BufferedImage(1, 1, java.awt.image.BufferedImage.TYPE_INT_ARGB)).createGraphics();
        java.awt.font.TextLayout tL = new java.awt.font.TextLayout(text, font, g2d.getFontRenderContext());
        return (double) tL.getBounds().getWidth();
    }

    // draw multiple violinplots (for example for multiple datasets) format: vals[dataset][category][value]
    public void draw(double[][][] vals, String[] datasetNames, String[][] xLabels2, int width, int height, int margin, int betweenPlotMargin, Output output, String outputFileName) throws IOException {

        Locale defaultLocale = Locale.getDefault();
        Locale.setDefault(Locale.US);
        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.lowagie.text.Document document = null;
        com.lowagie.text.pdf.PdfWriter writer = null;
        com.lowagie.text.pdf.PdfContentByte cb = null;
        
        BufferedImage bi = null;
        bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();


        int marginLeft = 50;

        if (output == Output.PDF) {
            com.lowagie.text.Rectangle rectangle = new com.lowagie.text.Rectangle(width + marginLeft, height);
            document = new com.lowagie.text.Document(rectangle);
            try {
                writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFileName));
            } catch (DocumentException e) {
                throw new IOException(e.fillInStackTrace());
            }
            document.open();
            cb = writer.getDirectContent();
            cb.saveState();
            //com.lowagie.text.pdf.DefaultFontMapper fontMap = new com.lowagie.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(width + marginLeft, height);
        } else {
            bi = new BufferedImage(width + marginLeft, height, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }


        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, width + marginLeft, height);

        int nrPlots = 0;
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int ds = 0; ds < vals.length; ds++) {
            if (vals[ds] != null) {
                for (int j = 0; j < vals[ds].length; j++) {
                    double min = Fmath.minimum(vals[ds][j]);
                    double max = Fmath.maximum(vals[ds][j]);
                    if (min < minValue) {
                        minValue = min;
                    }
                    if (max > maxValue) {
                        maxValue = max;
                    }
                }
                nrPlots += vals[ds].length;
            }

        }

        int plotwidth = ((width - (margin * 2) - ((nrPlots - 1) * betweenPlotMargin)) / nrPlots);
        System.out.println(width + "\t" + plotwidth + "\t" + nrPlots);
        int plotheight = height - (2 * margin);

        // draw individual box/violinplots
        java.awt.AlphaComposite alphaComposite25 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f);
        java.awt.AlphaComposite alphaComposite50 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f);
        java.awt.AlphaComposite alphaComposite100 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC, 1.00f);



        float fontSize = 12f;
        java.awt.Font font = new java.awt.Font("Gill Sans MT", java.awt.Font.PLAIN, (int) fontSize);
        java.awt.Font fontBold = new java.awt.Font("Gill Sans MT", java.awt.Font.BOLD, (int) fontSize);
        java.awt.Font fontSmall = new java.awt.Font("Gill Sans MT", java.awt.Font.PLAIN, 8);
        java.awt.Font fontBoldSmall = new java.awt.Font("Gill Sans MT", java.awt.Font.BOLD, 8);

        int ctr = 0;

        for (int d = 0; d < vals.length; d++) {
            // plot datasetname
            if (vals[d] != null) {
                String dsName = datasetNames[d];
                int xLab = marginLeft + margin + (ctr * betweenPlotMargin) + (plotwidth * ctr);

                int yLab = margin;
                g2d.setFont(fontSmall);
                g2d.setComposite(alphaComposite100);
                g2d.setColor(Color.black);
                g2d.drawString(dsName, xLab, yLab);

                for (int j = 0; j < vals[d].length; j++) {
                    double[] valsDs = vals[d][j];

                    String name = xLabels2[d][j];
                    int x = marginLeft + margin + (ctr * betweenPlotMargin) + (plotwidth * ctr);
                    int y = margin;



                    g2d.setComposite(alphaComposite25);
                    g2d.setColor(new Color(223, 36, 20));
                    drawViolinPlot(g2d, x, y, plotwidth, plotheight, valsDs, minValue, maxValue);

                    g2d.setComposite(alphaComposite100);
                    g2d.setColor(new Color(0, 0, 0));
                    drawBoxPlot(g2d, x, y, plotwidth, plotheight, valsDs, minValue, maxValue, true, false);

                    g2d.setComposite(alphaComposite100);
                    g2d.setColor(Color.black);
                    g2d.drawString(name, x, height - margin);
                    ctr++;
                }
            } else {
                // plot that the dataset is not there.
            }


        }

        int marginTop = margin;
        int innerHeight = height - (2 * margin);

        //Draw expression level indicator:
        g2d.setComposite(alphaComposite100);
        g2d.setColor(new Color(100, 100, 100));
        int startVal = (int) Math.ceil(minValue);
        int endVal = (int) Math.floor(maxValue);
        int posY1 = marginTop + innerHeight - (int) Math.round((double) innerHeight * (startVal - minValue) / (maxValue - minValue));
        int posY2 = marginTop + innerHeight - (int) Math.round((double) innerHeight * (endVal - minValue) / (maxValue - minValue));
        g2d.drawLine(marginLeft - 10, posY1, marginLeft - 10, posY2);
        g2d.setFont(fontBold);
        for (int v = startVal; v <= endVal; v++) {
            int posY = marginTop + innerHeight - (int) Math.round((double) innerHeight * (v - minValue) / (maxValue - minValue));
            g2d.drawLine(marginLeft - 10, posY, marginLeft - 20, posY);
            g2d.drawString(String.valueOf(v), marginLeft - 25 - (int) getWidth(String.valueOf(v), g2d.getFont()), posY + 3);
        }

        g2d.translate(marginLeft - 40, marginTop + innerHeight / 2);
        g2d.rotate(-0.5 * Math.PI);
        g2d.drawString("Relative gene expression", -(int) getWidth("Relative gene expression", g2d.getFont()) / 2, 0);
        g2d.rotate(+0.5 * Math.PI);
        g2d.translate(-(marginLeft - 40), -(marginTop + innerHeight / 2));

        g2d.dispose();
        if (output == Output.PDF) {
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(outputFileName));
        }


        Locale.setDefault(defaultLocale);
    }

    public void draw(double[][] vals, String[] xLabels, int width, int height, int margin, int betweenPlotMargin, Output output, String outputFileName) throws IOException, DocumentException {

        Locale defaultLocale = Locale.getDefault();
        Locale.setDefault(Locale.US);
        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.lowagie.text.Document document = null;
        com.lowagie.text.pdf.PdfWriter writer = null;
        BufferedImage bi = null;
        bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();

com.lowagie.text.pdf.PdfContentByte cb = null;
        if (output == Output.PDF) {
            com.lowagie.text.Rectangle rectangle = new com.lowagie.text.Rectangle(width, height);
            document = new com.lowagie.text.Document(rectangle);
            writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFileName));

            document.open();
            cb =  writer.getDirectContent();
            cb.saveState();
            //com.lowagie.text.pdf.DefaultFontMapper fontMap = new com.lowagie.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(width, height);
        } else {
            bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }


        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, width, height);


        int nrPlots = xLabels.length;
        int plotwidth = width - (margin * 2) - ((nrPlots - 1) * betweenPlotMargin) / nrPlots;
        int plotheight = height - (2 * margin);

        double minVal = Double.MAX_VALUE;
        double maxVal = -Double.MAX_VALUE;
        for (int i = 0; i < xLabels.length; i++) {

            double min = Fmath.minimum(vals[i]);
            double max = Fmath.maximum(vals[i]);
            if (min < minVal) {
                minVal = min;
            }
            if (max > maxVal) {
                maxVal = max;
            }

        }

        // draw individual box/violinplots
        for (int i = 0; i < xLabels.length; i++) {
            double[] valsDs = vals[i];
            String name = xLabels[i];
            int x = margin + (i * betweenPlotMargin) + (plotwidth * i);
            int y = margin;

            drawViolinPlot(g2d, x, y, plotwidth, plotheight, valsDs, minVal, maxVal);

            drawBoxPlot(g2d, x, y, plotwidth, plotheight, valsDs, minVal, maxVal, true, false);
        }

        g2d.dispose();
        if (output == Output.PDF) {
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(outputFileName));
        }


        Locale.setDefault(defaultLocale);
    }

    public void drawViolinPlot(Graphics2D g2d, int x, int y, int width, int height, double[] vals, double minValue, double maxValue) {

        int nrVals = vals.length;

        //Determine range of values:
        double minVals = JSci.maths.ArrayMath.min(vals);
        double maxVals = JSci.maths.ArrayMath.max(vals);

        //Make frequency distribution:
        int nrBins = 1 + (int) Math.round(Math.sqrt(nrVals) / 2d);
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
        double kernelWidth = 10;
        double[] kernelWeights = new double[201];
        for (int d = -100; d <= 100; d++) {
            double weight = java.lang.Math.pow(java.lang.Math.E, -((double) d / (double) kernelWidth) * ((double) d / (double) kernelWidth) / 2);
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
        double maxSmoothedVal = JSci.maths.ArrayMath.max(posValSmoothed);
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

    public void drawBoxPlot(Graphics2D g2d, int x, int y, int width, int height, double[] vals, double minValue, double maxValue, boolean drawOutliers, boolean fill) {

        double median = JSci.maths.ArrayMath.median(vals);
        double q1 = JSci.maths.ArrayMath.percentile(vals, 0.25d);
        double q3 = JSci.maths.ArrayMath.percentile(vals, 0.75d);
        double iqr = q3 - q1;


        //Draw median:
        int posY = y + height - (int) Math.round((double) height * (median - minValue) / (maxValue - minValue));
        g2d.setStroke(new java.awt.BasicStroke(2.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawLine(x, posY, x + width, posY);
        //Draw IQR:
        int posY1 = y + height - (int) Math.round((double) height * (q3 - minValue) / (maxValue - minValue));
        int posY2 = y + height - (int) Math.round((double) height * (q1 - minValue) / (maxValue - minValue));

        if (fill) {
            Color c = g2d.getColor();
            Color c2 = new Color(255, 255, 255, 128);
            g2d.setColor(c2);
            g2d.fillRect(x, posY1, width, posY2 - posY1);
            g2d.setColor(c);
        }
        g2d.setStroke(new java.awt.BasicStroke(1.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawRect(x, posY1, width, posY2 - posY1);

        //Draw whiskers:
        double whiskerTop = q3 + 1.5d * iqr;
        double whiskerBottom = q1 - 1.5d * iqr;
        double max = JSci.maths.ArrayMath.max(vals);
        double min = JSci.maths.ArrayMath.min(vals);
        if (min > whiskerBottom) {
            whiskerBottom = min;
        }
        if (max < whiskerTop) {
            whiskerTop = max;
        }
        posY = y + height - (int) Math.round((double) height * (whiskerTop - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, new float[]{2.0f}, 0.0f));
        g2d.drawLine(x + width / 2, posY, x + width / 2, posY1);
        g2d.setStroke(new java.awt.BasicStroke(1.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawLine(x + width / 2 - 5, posY, x + width / 2 + 5, posY);
        posY = y + height - (int) Math.round((double) height * (whiskerBottom - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, new float[]{2.0f}, 0.0f));
        g2d.drawLine(x + width / 2, posY2, x + width / 2, posY);
        g2d.setStroke(new java.awt.BasicStroke(1.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawLine(x + width / 2 - 5, posY, x + width / 2 + 5, posY);

        //Draw outliers:
        if (drawOutliers) {
            java.awt.AlphaComposite alphaComposite10 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.10f);
            g2d.setComposite(alphaComposite10);
            for (int v = 0; v < vals.length; v++) {
                if (vals[v] > whiskerTop || vals[v] < whiskerBottom) {
                    posY = y + height - (int) Math.round((double) height * (vals[v] - minValue) / (maxValue - minValue));
                    g2d.drawOval(x + width / 2 - 2, posY - 2, 5, 5);
                }
            }
        }


    }
}
