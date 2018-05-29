/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.graphics;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import umcg.genetica.math.stats.Regression;

/**
 *
 * @author harm-jan
 */
public class QQPlot {

    private final static int FILE_TYPE_PNG = 1;
    private final static int FILE_TYPE_PDF = 2;
    private int outputPlotsFileType = FILE_TYPE_PDF;

    public void draw(String fileName, double fdrCutOff, int nrPermutationsFDR, int maxNrMostSignificantEQTLs, double[][] permutedPValues, double[] pValues, boolean[] pValueSignificant, int nrSignificantEQTLs) {

        System.setProperty("java.awt.headless", "true");

        //Draw QQ plot, color the significant eQTLs with a different color, and calculate lambda inflation statistic
        int width = 600;
        int height = 600;
        int margin = 50;
        int x0 = margin;
        int x1 = width - margin;
        int innerWidth = x1 - x0;
        int y0 = margin + 40;
        int y1 = height - margin - 12 - 10;
        int innerHeight = y1 - y0;

        File fileQQPlot = null;
        if (outputPlotsFileType == FILE_TYPE_PNG) {
            fileQQPlot = new File(fileName);
        } else {
            fileQQPlot = new File(fileName);
        }
        Graphics2D g2d = null;
        BufferedImage bi = null;
        com.itextpdf.text.Document document = null;
        com.itextpdf.text.pdf.PdfContentByte cb = null;
        com.itextpdf.text.pdf.PdfWriter writer = null;
        if (outputPlotsFileType == FILE_TYPE_PNG) {
            bi = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        } else {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);

            document = new com.itextpdf.text.Document(rectangle);

            try {
                writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(fileQQPlot));
                document.open();
                cb = writer.getDirectContent();
                cb.saveState();
                g2d = cb.createGraphics(width, height);
            } catch (Exception e) {
                System.out.println("Cannot write to PDF file!:\t" + fileQQPlot.getAbsolutePath());
                System.exit(-1);
            }

        }
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height);

        double capLog = 16;
        g2d.setFont(new java.awt.Font("Arial", java.awt.Font.PLAIN, 11));
        for (int l = 0; l <= capLog; l++) {
            double log10 = l;
            int posX = x0 + (int) Math.round((double) innerWidth * (log10 / capLog));
            int posY = y1 - (int) Math.round((double) innerHeight * (log10 / capLog));
            g2d.setColor(new Color(220, 220, 220));
            g2d.drawLine(x0, posY, x1, posY);
            g2d.drawLine(posX, y0, posX, y1);
            g2d.setColor(new Color(0, 0, 0));
            g2d.drawLine(x0 - 3, posY, x0, posY);
            g2d.drawLine(posX, y1, posX, y1 + 3);
            //Draw axis:
            g2d.setColor(new Color(0, 0, 0));
            String logPValueString = (new java.text.DecimalFormat("##;-##", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(log10);
            g2d.drawString(logPValueString, posX - 1 - getWidth(logPValueString, g2d.getFont()) / 2, y1 + 13);
            g2d.drawString(logPValueString, x0 - 7 - getWidth(logPValueString, g2d.getFont()), posY + 4);
        }

        g2d.setColor(new Color(200, 200, 200));
        g2d.drawLine(x0, y1, x1, y0);
        g2d.setColor(new Color(0, 0, 0));
        g2d.drawLine(x0, y1, x0, y0);
        g2d.drawLine(x0, y1, x1, y1);

        double[] distLog10Observed = new double[maxNrMostSignificantEQTLs];
        double[] distLog10Null = new double[maxNrMostSignificantEQTLs];
        for (int p = 0; p < maxNrMostSignificantEQTLs; p++) {
            double log10Observed = -Math.log10(pValues[p]);
            distLog10Observed[p] = log10Observed;
            if (log10Observed > capLog) {
                log10Observed = capLog;
            }
            double log10Null = 0;
            for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
                log10Null += permutedPValues[permutationRound][p];
            }
            
            log10Null = -Math.log10(log10Null / ((double) nrPermutationsFDR));
            distLog10Null[p] = log10Null;
            if (log10Null > capLog) {
                log10Null = capLog;
            }
            if (pValueSignificant[p]) {
                g2d.setColor(new Color(255, 0, 0));
            } else {
                g2d.setColor(new Color(0, 0, 0));
            }
            int posX = x0 + (int) Math.round((double) innerWidth * (log10Null / capLog));
            int posY = y1 - (int) Math.round((double) innerHeight * (log10Observed / capLog));
            g2d.fillOval(posX - 1, posY - 1, 3, 3);
        }
        double[] rc = Regression.getLinearRegressionCoefficients(distLog10Null, distLog10Observed);
        g2d.setColor(new Color(0, 0, 0));
        String slopeString = "Lambda inflation: " + (new java.text.DecimalFormat("##.##;-##.##", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(rc[0]);
        g2d.drawString(slopeString, x0, y0 - 15);
        String fdrString = (new java.text.DecimalFormat("##.##;-##.##", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(fdrCutOff);
        String nrEQTLsString = "Number of significant eQTLs (FDR " + fdrString + "):\t" + nrSignificantEQTLs;
        g2d.drawString(nrEQTLsString, x0, y0 - 30);

        //double correlation = JSci.maths.ArrayMath.correlation(distLog10Null, distLog10Observed);
        //String correlationString = "Correlation: " + (new java.text.DecimalFormat("##.##;-##.##", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(correlation);
        //g2d.drawString(correlationString, x0, y0 - 30);
        //System.out.println(rc[0] + "\t" + rc[1] + "\t" + correlation);

        //Save image:
        if (outputPlotsFileType == FILE_TYPE_PNG) {
            try {
                javax.imageio.ImageIO.write(bi, "png", fileQQPlot);
            } catch (Exception e) {
                System.out.println(e.getMessage());
                System.out.println(e.getStackTrace());
            }
        } else {
            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        }
    }

    private int getWidth(String text, java.awt.Font font) {
        java.awt.Graphics2D g2d = (new java.awt.image.BufferedImage(1, 1, java.awt.image.BufferedImage.TYPE_INT_ARGB)).createGraphics();
        java.awt.font.TextLayout tL = new java.awt.font.TextLayout(text, font, g2d.getFontRenderContext());
        return (int) tL.getBounds().getWidth();
    }
}
