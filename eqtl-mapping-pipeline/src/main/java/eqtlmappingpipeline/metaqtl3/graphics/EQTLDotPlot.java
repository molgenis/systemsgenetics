/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.graphics;

import com.itextpdf.text.DocumentException;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harm-jan
 *
 *
 *
 */
public class EQTLDotPlot {

    private long[] cumChrPos;
    private int[] chromosomeLength = {0,
        247199719, 242751149, 199446827, 191263063, 180837866,
        170896992, 158821424, 146274826, 140273252, 135374737,
        134452384, 132289534, 114127980, 106360585, 100338915,
        88822254, 78654742, 76117153, 63806651, 62435964,
        46944323, 49591432, 154913754, 57772954, 154913754, 0};

    public enum Output {

        PDF, PNG
    };

    public void draw(String inputFile, String outputFile, Output output) throws IOException, DocumentException {

        System.setProperty("java.awt.headless", "true");

        //Init chromosomal relative location:
        cumChrPos = new long[25];
        for (int chr = 0; chr < 25; chr++) {
            if (chr > 0) {
                cumChrPos[chr] = cumChrPos[chr - 1];
            }
            cumChrPos[chr] += chromosomeLength[chr];
        }

        //Init image:
        int width = 1200;
        int height = 1200;
        int margin = 100;
        int x0 = margin;
        int x1 = width - margin;
        int innerWidth = x1 - x0;
        int y0 = margin;
        int y1 = height - margin;
        int innerHeight = y1 - y0;

        Graphics2D g2d = null;
        com.itextpdf.text.Document document = null;
        com.itextpdf.text.pdf.PdfWriter writer = null;
        com.itextpdf.text.pdf.PdfContentByte cb = null;
        BufferedImage bi = null;
        if (output == Output.PDF) {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);
            writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFile));

            document.open();
             cb = writer.getDirectContent();
             cb.saveState();
            //com.itextpdf.text.pdf.DefaultFontMapper fontMap = new com.itextpdf.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(width, height);
        } else {
            bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }

//        BufferedImage bi = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_RGB);
//        Graphics2D g2d = bi.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height);
        g2d.setFont(new java.awt.Font(g2d.getFont().getFontName(), java.awt.Font.PLAIN, 8));

        for (int chr = 1; chr <= 24; chr++) {
            int x = getPlotPosition(chr, 0, innerWidth);
            int x2 = getPlotPosition(chr, chromosomeLength[chr], innerWidth);
            if (chr % 2 == 1) {
                g2d.setColor(new Color(150, 150, 150));
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.02f));
                g2d.fillRect(x0 + x, y0, x2 - x, innerHeight);
                g2d.fillRect(x0, y1 - x2, innerWidth, x2 - x);
                g2d.setColor(new Color(150, 150, 150));
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f));
                g2d.drawLine(x0 + x, y0, x0 + x, x1);
                g2d.drawLine(x0 + x2, y0, x0 + x2, x1);
                g2d.drawLine(x0, y1 - x2, x1, y1 - x2);
                g2d.drawLine(x0, y1 - x, x1, y1 - x);
            }

            String chromosomeName = String.valueOf(chr);
            if (chr == 23) {
                chromosomeName = "X";
            }
            if (chr == 24) {
                chromosomeName = "Y";
            }
            g2d.setColor(new Color(0, 0, 0));
            g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
            g2d.drawString(chromosomeName, x0 + ((x2 + x) / 2) - getWidth(chromosomeName, g2d.getFont()) / 2, y1 + 15);
            g2d.drawString(chromosomeName, x0 - 5 - getWidth(chromosomeName, g2d.getFont()), y1 - (x2 + x) / 2 + 3);
        }

        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.0f));
        g2d.setColor(new Color(180, 180, 180));
        g2d.drawRect(x0, y0, innerWidth, innerHeight);

        double pValueMin = 1;
        double pValueMax = 0;
        double log10Min = -1;
        double log10Max = -1;

        TextFile in = new TextFile(inputFile, TextFile.R);
        
        String str = in.readLine();
        while ((str = in.readLine()) != null) {
            String[] data = str.split("\t");
            double pValue = Double.parseDouble(data[0]);
            if (pValue < pValueMin) {
                pValueMin = pValue;
            }
            if (pValue > pValueMax) {
                pValueMax = pValue;
            }
        }
        in.close();
        log10Min = -Math.log10(pValueMax);
        log10Max = -Math.log10(pValueMin);
        if (log10Max > 65.90) {
            log10Max = 65.90;
        }
        
        in.open();
        str = in.readLine();
        int[] distOverChr = new int[1000];
        while ((str = in.readLine()) != null) {
            String[] data = str.split("\t");
            double pValue = Double.parseDouble(data[0]);
            double log10PValue = -Math.log10(pValue);
            if (log10PValue > 65.90) {
                log10PValue = 65.90;
            }
            int chrSNP = Integer.parseInt(data[2]);
            int chrPosSNP = Integer.parseInt(data[3]);
            int chrProbe = Integer.parseInt(data[5]);
            int chrPosProbe = Integer.parseInt(data[6]);
            int x = getPlotPosition(chrSNP, chrPosSNP, innerWidth);
            distOverChr[x / 10]++;
            int y = getPlotPosition(chrProbe, chrPosProbe, innerWidth);
            double proportion = (log10PValue - log10Min) / (log10Max - log10Min);
            int ovalSize = (int) Math.round(4.0d + proportion * 10.0d);
            g2d.setColor(new Color(0, 0, 0));
            g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f));
            g2d.fillOval(margin + x - ovalSize / 2, y1 - y - ovalSize / 2, ovalSize, ovalSize);
            g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.75f));
            g2d.drawOval(margin + x - ovalSize / 2, y1 - y - ovalSize / 2, ovalSize, ovalSize);
        }

        //Draw distribution of eQTLs:
        int maxDist = 0;
        for (int d = 0; d < 1000; d++) {
            if (distOverChr[d] > maxDist) {
                maxDist = distOverChr[d];
            }
        }
        g2d.setColor(new Color(0, 0, 0));
        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f));
        g2d.drawLine(x0, y1 + 80, x1, y1 + 80);
        for (int d = 0; d < 1000; d++) {
            if (distOverChr[d] > 0) {
                int x = d * 10;
                int y = 80;
                g2d.setColor(new Color(0, 0, 0));
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f));
                int barHeight = (int) Math.round(70.0d * (double) distOverChr[d] / (double) maxDist);
                g2d.fillRect(margin + x - 3, y1 + y - barHeight, 6, barHeight);
                g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f));
                g2d.drawRect(margin + x - 3, y1 + y - barHeight, 6, barHeight);
            }
        }

        in.close();


        //Save image:
        if (output == Output.PDF) {
            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(outputFile));
        }
    }

    private int getPlotPosition(int chr, int chrPos, int graphWidth) {
        int locus = chr;
        if (locus != -1 && locus < 25) {
            long cumPos = chrPos;
            if (chr > 0) {
                cumPos += cumChrPos[chr - 1];
            }
            return (int) ((double) graphWidth * ((double) cumPos / (double) cumChrPos[cumChrPos.length - 1]));
        }
        return -1;
    }

    private int getWidth(String text, java.awt.Font font) {
        java.awt.Graphics2D g2d = (new java.awt.image.BufferedImage(1, 1, java.awt.image.BufferedImage.TYPE_INT_ARGB)).createGraphics();
        java.awt.font.TextLayout tL = new java.awt.font.TextLayout(text, font, g2d.getFontRenderContext());
        return (int) tL.getBounds().getWidth();
    }
}
