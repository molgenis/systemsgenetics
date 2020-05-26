/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.graphics;

import JSci.maths.ArrayMath;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Rectangle;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 *
 * @author juha
 */
public class Heatmap {

    private static final Font ROW_FONT = new Font("Verdana", Font.PLAIN, 14);
    private static final Font COL_FONT = new Font("Verdana", Font.PLAIN, 10);
    private static final Font LEGEND_FONT = new Font("Verdana", Font.PLAIN, 14);
//    private static final int LEGEND_X = 10;
    private static final int LEGEND_Y = 10;
    private static final int LEGEND_WIDTH = 200;
    private static final int LEGEND_HEIGHT = 100;
    private static final Logger LOGGER = Logger.getLogger(Heatmap.class.getName());

    public enum Output {

        PDF, PNG
    };

    /**
     * Original heatmap implementation Juha
     *
     * @param values
     * @param rowHeaders
     * @param colHeaders
     * @param width
     * @param height
     * @param filename
     * @param output
     * @throws IOException
     * @throws DocumentException
     */
    public static void drawHeatmap(double[][] values, String[] rowHeaders, String[] colHeaders, int width, int height, String filename, Output output) throws IOException, DocumentException {

        if (values.length != rowHeaders.length) {
            throw new IllegalArgumentException("Data length and number of row headers differ!");
        }
        if (values[0].length != colHeaders.length) {
            throw new IllegalArgumentException("Data length and number of column headers differ!");
        }
        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.itextpdf.text.Document document = null;
        com.itextpdf.text.pdf.PdfWriter writer = null;
        com.itextpdf.text.pdf.PdfContentByte cb = null;
        BufferedImage bi = null;
        if (output == Output.PDF) {
            Rectangle rectangle = new Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);
            writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(filename));

            document.open();
            cb = writer.getDirectContent();
            //com.itextpdf.text.pdf.DefaultFontMapper fontMap = new com.itextpdf.text.pdf.DefaultFontMapper();
            cb.saveState();
            g2d = cb.createGraphics(width, height);
        } else {
            bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, width, height);

        g2d.setFont(ROW_FONT);
        FontMetrics fontmetrics = g2d.getFontMetrics();
        int leftMargin = 0;
        for (String s : rowHeaders) {
            leftMargin = Math.max(leftMargin, fontmetrics.stringWidth(s));
        }
        leftMargin += 10;
        g2d.setFont(COL_FONT);
        fontmetrics = g2d.getFontMetrics();
        int topMargin = 0;
        for (String s : colHeaders) {
            topMargin = Math.max(topMargin, fontmetrics.stringWidth(s));
        }
        topMargin += 10;
        topMargin += 200;

        int plotWidth = width - leftMargin;
        int plotHeight = height - topMargin;

        int tileW = (int) (double) plotWidth / values[0].length; // tile width in pixels
        int tileH = (int) (double) plotHeight / values.length; // tile height in pixels

        if (tileW < 3) {
            throw new IllegalArgumentException("Map tiles would become less than 3 pixels wide. Output a wider image or use less columns.");
        }
        if (tileH < 3) {
            throw new IllegalArgumentException("Map tiles would become less than 3 pixels tall. Output a wider image or use less rows.");
        }

        // draw row headers
        g2d.setFont(ROW_FONT);
        g2d.setColor(Color.black);
        for (int row = 0; row < values.length; row++) {
            g2d.drawString(rowHeaders[row], 0, topMargin + row * tileH + tileH / 2 + 5);
        }

        // draw col headers
        g2d.setFont(COL_FONT);
        g2d.setColor(Color.black);
        for (int col = 0; col < values[0].length; col++) {
            g2d.rotate(-Math.PI / 2.0);
//            g2d.drawString(colHeaders[col], leftMargin + col * tileW + tileW / 2, 20);
            g2d.drawString(colHeaders[col], -topMargin + 10, leftMargin + col * tileW + tileW / 2);
            g2d.rotate(Math.PI / 2.0);
        }

        // draw heatmap
        double min = ArrayMath.min(values);
        double max = ArrayMath.max(values);
        normalize(values);
        for (int row = 0; row < values.length; row++) {
            for (int col = 0; col < values[0].length; col++) {
//                g2d.setColor(getColor(values[row][col]));
                g2d.setColor(getRGB(values[row][col]));
                g2d.fillRect(leftMargin + col * tileW, topMargin + row * tileH, tileW, tileH);
            }
        }

        // draw scale legend
        int legendX = width - LEGEND_WIDTH - 20;
        GradientPaint gp = new GradientPaint(
                legendX, LEGEND_Y, Color.blue,
                legendX + LEGEND_WIDTH / 2, LEGEND_Y, Color.white);
        g2d.setPaint(gp);
        g2d.fillRect(legendX, LEGEND_Y, LEGEND_WIDTH / 2, LEGEND_HEIGHT);
        gp = new GradientPaint(
                legendX + LEGEND_WIDTH / 2, LEGEND_Y, Color.white,
                legendX + LEGEND_WIDTH, LEGEND_Y, Color.red);
        g2d.setPaint(gp);
        g2d.fillRect(legendX + LEGEND_WIDTH / 2, LEGEND_Y, LEGEND_WIDTH / 2, LEGEND_HEIGHT);
        g2d.setColor(Color.black);
        g2d.setFont(LEGEND_FONT);
        String minStr = Math.round(min) + "";
        int stringWidth = fontmetrics.stringWidth(minStr);
        g2d.drawString(minStr, legendX - stringWidth / 2, LEGEND_Y + LEGEND_HEIGHT + 20);
        String avgStr = Math.round((max - min) / 2) + "";
        stringWidth = fontmetrics.stringWidth(avgStr);
        g2d.drawString(avgStr, legendX + LEGEND_WIDTH / 2 - stringWidth / 2, LEGEND_Y + LEGEND_HEIGHT + 20);
        String maxStr = Math.round(max) + "";
        stringWidth = fontmetrics.stringWidth(maxStr);
        g2d.drawString(maxStr, legendX + LEGEND_WIDTH - stringWidth / 2, LEGEND_Y + LEGEND_HEIGHT + 20);

        g2d.dispose();
        if (output == Output.PDF) {
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(filename));
        }
    }

    /**
     * Correlation heatmap implementation Marc Jan Always use minimal value:-1
     * and maximum value: 1
     *
     * @param values
     * @param rowHeaders
     * @param colHeaders
     * @param width
     * @param height
     * @param filename
     * @param output
     * @throws IOException
     * @throws DocumentException
     */
    public static void drawCorrelationHeatmap(double[][] values, String[] rowHeaders, String[] colHeaders, int width, int height, String filename, Output output) throws IOException, DocumentException {

        if (values.length != rowHeaders.length) {
            throw new IllegalArgumentException("Data length and number of row headers differ!");
        }
        if (values[0].length != colHeaders.length) {
            throw new IllegalArgumentException("Data length and number of column headers differ!");
        }

        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.itextpdf.text.Document document = null;
        com.itextpdf.text.pdf.PdfWriter writer = null;
        com.itextpdf.text.pdf.PdfContentByte cb = null;
        BufferedImage bi = null;
        if (output == Output.PDF) {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);
            writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(filename));

            document.open();
            cb = writer.getDirectContent();
            cb.saveState();
            //com.itextpdf.text.pdf.DefaultFontMapper fontMap = new com.itextpdf.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(width, height);
        } else {
            bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, width, height);

        g2d.setFont(ROW_FONT);
        FontMetrics fontmetrics = g2d.getFontMetrics();
        int leftMargin = 0;
        for (String s : rowHeaders) {
            leftMargin = Math.max(leftMargin, fontmetrics.stringWidth(s));
        }
        leftMargin += 10;
        g2d.setFont(COL_FONT);
        fontmetrics = g2d.getFontMetrics();
        int topMargin = 0;
        for (String s : colHeaders) {
            topMargin = Math.max(topMargin, fontmetrics.stringWidth(s));
        }
        topMargin += 10;
        topMargin += 200;

        int plotWidth = width - leftMargin;
        int plotHeight = height - topMargin;

        int tileW = (int) (double) plotWidth / values[0].length; // tile width in pixels
        int tileH = (int) (double) plotHeight / values.length; // tile height in pixels

        if (tileW < 3) {
            throw new IllegalArgumentException("Map tiles would become less than 3 pixels wide (" + plotWidth + "/" + values[0].length + " = " + ((double) plotWidth / values[0].length) + "). Output a wider image or use less columns.");
        }
        if (tileH < 3) {
            throw new IllegalArgumentException("Map tiles would become less than 3 pixels tall (" + plotHeight + "/" + values.length + " = " + ((double) plotHeight / values.length) + "). Output a wider image or use less rows.");
        }

        // draw row headers
        g2d.setFont(ROW_FONT);
        g2d.setColor(Color.black);
        for (int row = 0; row < values.length; row++) {
            g2d.drawString(rowHeaders[row], 0, topMargin + row * tileH + tileH / 2 + 5);
        }

        // draw col headers
        g2d.setFont(COL_FONT);
        g2d.setColor(Color.black);
        for (int col = 0; col < values[0].length; col++) {
            g2d.rotate(-Math.PI / 2.0);
//            g2d.drawString(colHeaders[col], leftMargin + col * tileW + tileW / 2, 20);
            g2d.drawString(colHeaders[col], -topMargin + 10, leftMargin + col * tileW + tileW / 2);
            g2d.rotate(Math.PI / 2.0);
        }

        // draw heatmap
        double min = -1;
        double max = 1;
        normalizeCorrelations(values);
        for (int row = 0; row < values.length; row++) {
            for (int col = 0; col < values[0].length; col++) {
//                g2d.setColor(getColor(values[row][col]));
                g2d.setColor(getRGB(values[row][col]));
                g2d.fillRect(leftMargin + col * tileW, topMargin + row * tileH, tileW, tileH);
            }
        }

        // draw scale legend
        int legendX = width - LEGEND_WIDTH - 20;
        GradientPaint gp = new GradientPaint(
                legendX, LEGEND_Y, Color.blue,
                legendX + LEGEND_WIDTH / 2, LEGEND_Y, Color.white);
        g2d.setPaint(gp);
        g2d.fillRect(legendX, LEGEND_Y, LEGEND_WIDTH / 2, LEGEND_HEIGHT);
        gp = new GradientPaint(
                legendX + LEGEND_WIDTH / 2, LEGEND_Y, Color.white,
                legendX + LEGEND_WIDTH, LEGEND_Y, Color.red);
        g2d.setPaint(gp);
        g2d.fillRect(legendX + LEGEND_WIDTH / 2, LEGEND_Y, LEGEND_WIDTH / 2, LEGEND_HEIGHT);
        g2d.setColor(Color.black);
        g2d.setFont(LEGEND_FONT);
        String minStr = Math.round(min) + "";
        int stringWidth = fontmetrics.stringWidth(minStr);
        g2d.drawString(minStr, legendX - stringWidth / 2, LEGEND_Y + LEGEND_HEIGHT + 20);
        String avgStr = "0";
        stringWidth = fontmetrics.stringWidth(avgStr);
        g2d.drawString(avgStr, legendX + LEGEND_WIDTH / 2 - stringWidth / 2, LEGEND_Y + LEGEND_HEIGHT + 20);
        String maxStr = Math.round(max) + "";
        stringWidth = fontmetrics.stringWidth(maxStr);
        g2d.drawString(maxStr, legendX + LEGEND_WIDTH - stringWidth / 2, LEGEND_Y + LEGEND_HEIGHT + 20);

        g2d.dispose();
        if (output == Output.PDF) {
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(filename));
        }
    }

    /**
     * Converts the given array to 0...1 scale linearly.
     *
     * @param values
     */
    private static void normalize(double[][] values) {

        double min = ArrayMath.min(values);
        double max = ArrayMath.max(values) + min;
        System.out.println("Range: " + min + " -- " + max);
        for (int row = 0; row < values.length; row++) {
            for (int col = 0; col < values[0].length; col++) {
                values[row][col] += min;
                values[row][col] /= max;
                System.out.println(values[row][col]);
            }
        }
        System.out.println("min " + ArrayMath.min(values));
        System.out.println("max " + ArrayMath.max(values));

    }

    /**
     * Converts the given correlation (-1 to 1) array to 0...1 scale linearly.
     * Add 1 and divide by 2.
     *
     * @param values
     */
    private static void normalizeCorrelations(double[][] values) {

        for (int row = 0; row < values.length; row++) {
            for (int col = 0; col < values[0].length; col++) {
                values[row][col] = (values[row][col] + 1) / 2;
            }
        }

    }

    private Rectangle getScaleGradient(int width, int height) {
        com.itextpdf.text.Rectangle r = new com.itextpdf.text.Rectangle(width, height);
        return r;

    }

    private static Color getRGB(double power) {
        if (power < 0.5) {
            power *= 2;
            return new Color((float) power, (float) power, 1f); // blue to white, 0 -> 1: (0, 0, 1) -> (1, 1, 1)
        } else {
            power -= 0.5;
            power *= 2;
            return new Color(1f, (float) (1 - power), (float) (1 - power)); // white to red, 0 -> 1: (1, 1, 1) -> (1, 0, 0)
        }
//        return new Color(0f, 0f, (float) (1 - power)); // blue to black
//        return new Color(1f, 1f, (float) (1 - power)); // white to yellow
//        return new Color(1f, 1f, (float) power); // yellow to white
    }

    private static Color getColor(double power) {

        float h = (float) (power * 0.3); // Hue
//        float h = (float) (240d / 360d); // Blue Hue
        float s = (float) 0.9; // Saturation
        float b = (float) 0.9; // Brightness

        return Color.getHSBColor(h, s, b);
    }
}
