/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.graphics;

import com.lowagie.text.DocumentException;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfTemplate;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class ScatterPlot {

    private BufferedImage bi;
    private Graphics2D g2d;
    private int graphHeight;
    private int graphWidth;
    private int drawWidth;
    private int drawHeight;
    private int margin = 50;
    private Color color;
    private Font font;
    private double unitX = Double.NaN;
    private double unitY = Double.NaN;
    private double[] x;
    private double[] y;
    private double maxX;
    private double maxY;
    private double minY;
    private double minX;
    private double rangeX;
    private double rangeY;
    private int fontheight;
    private OUTPUTFORMAT format;
    private String outfilename;
    private com.lowagie.text.Document document = null;
    private com.lowagie.text.pdf.PdfWriter writer = null;
    private PdfContentByte cb;

    public enum OUTPUTFORMAT {

        PDF, PNG, JPG
    }

    public ScatterPlot(int sizeX, int sizeY, double[] x, double[] y, OUTPUTFORMAT format, String outfile) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Error initializing scatterplot: X and Y do not have the same length! " + x.length + "x" + y.length);
        }
        graphHeight = sizeY;
        graphWidth = sizeX;

        this.outfilename = outfile;
        this.format = format;

        this.x = x;
        this.y = y;
        init();
        plot();
        draw();
    }

    public ScatterPlot(int sizeX, int sizeY, double unitX, double unitY, double[] x, double[] y, OUTPUTFORMAT format, String outfile) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Error initializing scatterplot: X and Y do not have the same length! " + x.length + "x" + y.length);
        }
        graphHeight = sizeY;
        graphWidth = sizeX;
        this.format = format;
        this.outfilename = outfile;

        this.unitX = unitX;
        this.unitY = unitY;

        this.x = x;
        this.y = y;

        init();
        plot();
        draw();
    }

    private void init() {

        if (margin > graphHeight || margin > graphWidth) {
            throw new IllegalArgumentException("Size of graph should be > " + (margin * 2) + " pixels in both directions");
        }

        drawWidth = graphWidth - (2 * margin);
        drawHeight = graphHeight - (2 * margin);

        if (format == OUTPUTFORMAT.PDF) {
            Rectangle rectangle = new Rectangle(graphWidth, graphHeight);
            document = new com.lowagie.text.Document(rectangle);
            
            if (!outfilename.toLowerCase().endsWith(".pdf")) {
                outfilename += ".pdf";
            }
            try {
                writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outfilename));
                
            } catch (DocumentException e) {
                e.printStackTrace();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            document.open();
            cb = writer.getDirectContent();
            cb.saveState();
            
//            com.lowagie.text.pdf.DefaultFontMapper fontMap = new com.lowagie.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(graphWidth, graphHeight);
        } else {
            bi = new java.awt.image.BufferedImage(graphWidth, graphHeight, java.awt.image.BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }

        color = new Color(255, 255, 255);

        g2d.setColor(color);
        g2d.setFont(new Font("Verdana", Font.PLAIN, 10));
        FontMetrics fontmetrics = g2d.getFontMetrics();
        fontheight = fontmetrics.getHeight();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        g2d.fillRect(0, 0, graphWidth, graphHeight);

        // draw axis
        color = new Color(0, 0, 0);
        g2d.setColor(color);
//        g2d.drawRect(margin, margin, graphHeight - (2 * margin), graphWidth - (2 * margin));

        g2d.setColor(color);

        // 
        determineRange();

        // draw axes
        drawAxis();
    }

    private void plot() {
        g2d.setColor(new Color(0, 128, 255, 64));
        for (int i = 0; i < x.length; i++) {

            double xval = x[i];
            double yval = y[i];

            if (Double.isInfinite(xval)) {
                if (xval < 0) {
                    xval = -Double.MAX_VALUE;
                } else {
                    xval = -Double.MAX_VALUE;
                }
            }

            if (Double.isInfinite(yval)) {
                if (yval < 0) {
                    yval = -Double.MAX_VALUE;
                } else {
                    yval = -Double.MAX_VALUE;
                }
            }


            int posX = margin + (int) Math.ceil((Math.abs(minX - xval) / rangeX) * drawWidth);
            int posY = margin + drawHeight - (int) Math.ceil((Math.abs(minY - yval) / rangeY) * drawHeight);

            g2d.fillRect(posX - 1, posY, 3, 3);
        }
    }

    private void draw() {
        try {
            g2d.dispose();
            if (format == OUTPUTFORMAT.JPG) {
                if (!outfilename.toLowerCase().endsWith(".jpg")) {
                    outfilename += ".jpg";
                }
                javax.imageio.ImageIO.write(bi, "jpg", new File(outfilename));
            } else if (format == OUTPUTFORMAT.PNG) {
                if (!outfilename.toLowerCase().endsWith(".png")) {
                    outfilename += ".png";
                }
                javax.imageio.ImageIO.write(bi, "png", new File(outfilename));
            } else {
                cb.restoreState();
                document.close();
                writer.close();
            }
        } catch (Exception e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }

    }

    private void drawAxis() {
        Color originalColor = g2d.getColor();
        Font oriFont = g2d.getFont();

        FontMetrics fontmetrics = g2d.getFontMetrics();
        int tickFontHeight = fontmetrics.getHeight();

        g2d.setFont(new Font("SansSerif", Font.PLAIN, 9));
        g2d.setColor(Color.gray);
        // x-axis
        int xposYAxis = 0;
        int yposXAxis = 0;
        if (minY >= 0) {
            // all Y values above 0, X axis starts at bottom left: |_
            g2d.drawLine(margin, margin + drawHeight, margin + drawWidth, margin + drawHeight);
            yposXAxis = margin + drawHeight;
        } else if (maxY <= 0) {
            // all Y values below 0, X axis starts at top left 
            g2d.drawLine(margin, margin, margin + drawWidth, margin);
            yposXAxis = margin;
        } else {
            // X-axis crosses the y-axis somewhere, at DIFF
            // this method does show some rounding errors at the moment..
            int diff = (int) Math.ceil((Math.abs(minY) / rangeY) * drawHeight);
            yposXAxis = margin + drawHeight - diff;
            g2d.drawLine(margin, yposXAxis, margin + drawWidth, yposXAxis);

        }

        // y-axis
        if (minX > 0) {
            // all X values above 0, Y-axis starts at top left
            g2d.drawLine(margin, margin, margin, margin + drawHeight);
            xposYAxis = margin;
        } else if (maxX <= 0) {
            // all X values below 0, axis starts at top right 
            g2d.drawLine(margin + drawWidth, margin, margin + drawWidth, margin + drawHeight);
            xposYAxis = margin + drawWidth;
        } else {
            // Y-axis crosses the X-axis somewhere, at DIFF
            // this method does show some rounding errors at the moment..
            int diff = (int) Math.ceil((Math.abs(minX) / rangeX) * drawWidth);
            xposYAxis = margin + diff;
            g2d.drawLine(margin + diff, margin, margin + diff, margin + drawHeight);

        }

        // X-ticks
        if (!Double.isNaN(unitX)) {
            // start drawing from minX, add one unitX at a time..
            // first round minX using the unitX
            double tickX = minX - (minX % unitX);

            while (tickX <= maxX) {

                double nonroundedPerc = Math.abs(minX - tickX) / rangeX;
//                double roundedPerc = roundToDecimals((nonroundedPerc), 2);
                int diff = (int) Math.floor(nonroundedPerc * drawWidth); // for the position relative to min and maxX

                String tickLabelFormatted = null;
                if (unitX > 10000 || unitX < 0.001) {
                    tickLabelFormatted = new DecimalFormat("0.#E0").format(tickX);
                } else {
                    tickLabelFormatted = new DecimalFormat("#.###").format(tickX);
                }

                if (tickX == 0) {
                    g2d.drawLine(xposYAxis, yposXAxis - 3, xposYAxis, yposXAxis + 3);
                    int nrPixelsString = g2d.getFontMetrics().stringWidth(tickLabelFormatted) / 2;
                    g2d.drawString(tickLabelFormatted, margin + diff - nrPixelsString, yposXAxis + tickFontHeight + 3);

                } else {
                    g2d.drawLine(margin + diff, yposXAxis - 3, margin + diff, yposXAxis + 3);
                    int nrPixelsString = g2d.getFontMetrics().stringWidth(tickLabelFormatted) / 2;
                    g2d.drawString(tickLabelFormatted, margin + diff - nrPixelsString, yposXAxis + tickFontHeight + 3);

                }

                tickX += unitX;

            }
        }

        // Y-ticks
        if (!Double.isNaN(unitY)) {
            double tickY = minY - (minY % unitY);
            while (tickY <= maxY) {

                double nonroundedPerc = Math.abs(minY - tickY) / rangeY;

//                double perc = roundToDecimals((Math.abs(minY - tickY) / rangeY), 2);
                int diff = (int) Math.floor(nonroundedPerc * drawHeight);
                String tickLabelFormatted = null;
                if (unitY > 100000 || unitY < 0.0001) {
                    tickLabelFormatted = new DecimalFormat("0.#E0").format(tickY);
                } else {
                    tickLabelFormatted = new DecimalFormat("#.###").format(tickY);
                }

                if (tickY == 0) {
                    g2d.drawLine(xposYAxis - 3, yposXAxis, xposYAxis + 3, yposXAxis);
                    g2d.drawString(tickLabelFormatted, xposYAxis - g2d.getFontMetrics().stringWidth(tickLabelFormatted) - 5, margin + drawHeight - diff + (tickFontHeight / 2) - 3);
                } else {
                    g2d.drawLine(xposYAxis - 3, margin + drawHeight - diff, xposYAxis + 3, margin + drawHeight - diff);
                    g2d.drawString(tickLabelFormatted, xposYAxis - g2d.getFontMetrics().stringWidth(tickLabelFormatted) - 5, margin + drawHeight - diff + (tickFontHeight / 2) - 3);
                }


                tickY += unitY;
            }
        }
        g2d.setFont(oriFont);
        g2d.setColor(originalColor);
    }

    private double roundToDecimals(double d, int c) {
        int temp = (int) ((d * Math.pow(10, c)));
        return (((double) temp) / Math.pow(10, c));
    }

    private void determineRange() {
        maxX = Primitives.max(x);
        maxY = Primitives.max(y);
        minX = Primitives.min(x);
        minY = Primitives.min(y);

        if (minY > 0) {
            minY = 0;
        }
        if (minX > 0) {
            minX = 0;
        }

        if (maxY < 0) {
            maxY = 0;
        }
        if (maxX < 0) {
            maxX = 0;
        }

        if (Double.isInfinite(maxX)) {
            maxX = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minX)) {
            minX = -Double.MAX_VALUE;
        }

        if (Double.isInfinite(maxY)) {
            maxY = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minY)) {
            minY = -Double.MAX_VALUE;
        }

        rangeX = Math.abs(minX - maxX);
        rangeY = Math.abs(minY - maxY);

//        System.out.println("MinX: " + minX + "\nMaxX: " + maxX + "\nMinY: " + minY + "\nMaxY: " + maxY + "\nRangeX: " + rangeX + "\nRangeY: " + rangeY + "\nUnitX: " + unitX + "\nUnitY: " + unitY);

        if (Double.isNaN(unitX)) {
            unitX = determineUnit(rangeX);
            if (unitX >= Math.abs(minX) && unitX >= Math.abs(maxX)) {
                unitX /= 10;
            }
        }
        if (Double.isNaN(unitY)) {
            unitY = determineUnit(rangeY);
            // prevent excessive tickmarking..
            if (unitY >= Math.abs(minY) && unitY >= Math.abs(maxY)) {
                unitY /= 10;
            }
        }


        double remainder = Math.abs(maxX) % unitX;
        if (remainder > 0d && maxX != 0) {
            double diff = unitX - remainder;
            maxX += diff;
        }

        remainder = Math.abs(minX) % unitX;
//        System.out.println(unitX + " min rem: " + remainder);
        if (remainder > 0d && minX != 0) {
            double diff = unitX - remainder;
            minX -= diff;
        }

        rangeX = Math.abs(minX - maxX);
        if (rangeX == 0) {
            maxX = unitX;
        }

        // ensure the max Y is rounded off by to the next unitY
        remainder = Math.abs(maxY) % unitY;
        if (remainder > 0d && maxY != 0d) {
            double diff = unitY - remainder;
            maxY += diff;
        }


        remainder = Math.abs(minY) % unitY; // -8 , unit == 10, remainder = 8 // diff = 2
//        System.out.println(remainder);
        if (remainder > 0d && minY != 0) {
            double diff = unitY - remainder;
            minY -= diff;
//            if (Math.abs(minY) - remainder == 0) {
//                minY = 0;
//            } else if (minY < 0) {
//                minY -= (unitY - remainder);
//            } else {
//                minY += (unitY - remainder);
//            }
        }

        rangeY = Math.abs(minY - maxY);
        if (rangeY == 0) {
            maxY = unitY;
        }

        if (Double.isInfinite(maxY)) {
            maxY = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minY)) {
            minY = -Double.MAX_VALUE;
        }

        if (Double.isInfinite(maxX)) {
            maxX = Double.MAX_VALUE;
        }

        if (Double.isInfinite(minX)) {
            minX = -Double.MAX_VALUE;
        }

        rangeX = Math.abs(minX - maxX);
        rangeY = Math.abs(minY - maxY);

//        System.out.println("");
//        System.out.println("MinX: " + minX + "\nMaxX: " + maxX + "\nMinY: " + minY + "\nMaxY: " + maxY + "\nRangeX: " + rangeX + "\nRangeY: " + rangeY + "\nUnitX: " + unitX + "\nUnitY: " + unitY);

    }

    private double determineUnit(double range) {

        double divisor = Math.log10(range);
        divisor = Math.floor(divisor);
        divisor = Math.pow(10, divisor);
        return divisor;
    }
}
