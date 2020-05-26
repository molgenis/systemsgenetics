/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.graphics;

import java.awt.*;

import java.io.File;


import java.awt.geom.*;
import java.awt.image.*;

/**
 *
 * @author harm-jan
 */
public class Graphics {

    private com.itextpdf.text.Document document;
    private boolean usePDF = false;
    protected BufferedImage bi;
    protected Graphics2D g2d;
    protected Color color;
    protected Font font;
    protected int graphHeight, graphWidth;
    protected int drawWidth, drawHeight;
    protected int marginTop, marginBottom, marginLeft, marginRight;
    protected double scalingX, scalingY;
    protected int FILE_TYPE;
    protected com.itextpdf.text.pdf.PdfContentByte cb;
    protected String outputLoc = "";
    protected com.itextpdf.text.pdf.PdfWriter writer;

    public Graphics() {
        bi = new java.awt.image.BufferedImage(100, 100, java.awt.image.BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();
    }

    public Graphics(int width, int height, boolean outputPDF, String outputLoc) {
        usePDF = outputPDF;
        this.outputLoc = outputLoc;
        init(width, height);
    }

    public Graphics(int width, int height) {
        init(width, height);
    }

    protected void init(int width, int height) {

        if (usePDF) {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);
            writer = null;
            try {
                writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputLoc));
                document.open();
                cb = writer.getDirectContent();
                cb.saveState();
            } catch (Exception e) {
                System.out.println("Cannot write to PDF file!:\t" + outputLoc);
                System.exit(-1);
            }


            this.usePDF = true;
        }

        graphHeight = height;
        graphWidth = width;

        drawWidth = width;
        drawHeight = height;

        if (usePDF) {
            g2d = cb.createGraphics(width, height);
        } else {
            bi = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }

        color = new Color(255, 255, 255);
        font = new Font("Georgia", Font.PLAIN, 10);

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2d.setColor(color);
        g2d.fillRect(0, 0, graphWidth, graphHeight);

        marginTop = 0;
        marginBottom = 0;
        marginLeft = 0;
        marginRight = 0;

        scalingX = 1;
        scalingY = 1;
    }

    // draw (device)
    public void draw(String outputFile) {

        if (usePDF) {
            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            try {
                javax.imageio.ImageIO.write(bi, "png", new File(outputFile));
            } catch (Exception e) {
                System.out.println(e.getMessage());
                System.out.println(e.getStackTrace());
            }
        }
    }

    public void close() {
        g2d = null;
        document = null;
    }

    // protected functions
    protected void setGraphBackgroundColor(int r, int g, int b) {
        g2d.setColor(new Color(r, g, b));
        g2d.fillRect(0, 0, graphWidth, graphHeight);
        g2d.setColor(color);
    }

    protected void setDrawBackgroundColor(int r, int g, int b) {
        g2d.setColor(new Color(r, g, b));
        g2d.fillRect(marginLeft, marginTop, drawWidth, drawHeight);
        g2d.setColor(color);
    }

    protected void setColor(Color newColor) {
        color = newColor;
        g2d.setColor(color);
    }

    protected void setColor(int r, int g, int b, int a) {
        color = new Color(r, g, b, a);
        g2d.setColor(color);
    }

    protected void setMargins(int margin) {
        marginTop = margin;
        marginBottom = margin;
        marginLeft = margin;
        marginRight = margin;
        calcDrawArea();
    }

    protected void setMargins(int top, int bottom, int left, int right) {
        marginTop = top;
        marginBottom = bottom;
        marginLeft = left;
        marginRight = right;
        calcDrawArea();
    }

    public void setTitle(String title) {
        if (marginTop > 0) {
            Font oldFont = font;
            Color oldColor = color;

            setColor(126, 126, 126, 255);
            setFont(20, Font.BOLD, "Georgia");

            int height = getStringHeight(title);
            int ymid = (int) Math.floor(marginTop / 2);

            if (height <= marginTop) {
                drawText(marginLeft, ymid, title);
            }

            setFont(oldFont);
            setColor(oldColor);
        }
    }

    public void setSubTitle(String title) {
        if (marginTop > 0) {
            Font oldFont = font;
            Color oldColor = color;

            setColor(126, 126, 126, 255);

            setFont(20, Font.PLAIN, "Georgia");
            int titleheight = getStringHeight("AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz0123456789[]./';.,{}");

            setFont(12, Font.PLAIN, "Georgia");
            int height = getStringHeight(title);
            int ymid = marginTop - (int) Math.floor(height / 2);

            if (height <= marginTop) {
                drawText(marginLeft, ymid, title);
            }

            setFont(oldFont);
            setColor(oldColor);
        }
    }

    public void setAxisLabels(String xlab, String ylab) {

        if (marginBottom > 0 && marginRight > 0) {
            Font oldFont = font;
            Color oldColor = color;

            setColor(126, 126, 126, 255);

            setFont(15, Font.BOLD, "Georgia");

            // x label
            int width = (int) Math.floor((double) getStringWidth(xlab) / 2);
            int xmid = (int) Math.floor(drawWidth / 2);
            int xpos = marginLeft + xmid - width;

            int height = getStringHeight(xlab);
            int ymid = (int) Math.floor(marginBottom / 2);
            int ypos = graphHeight - (ymid - height);

            if (height <= marginBottom) {
                drawText(xpos, ypos, xlab);
            }

            // y label
            width = (int) Math.floor((double) getStringWidth(ylab) / 2);
            ymid = (int) Math.floor(drawHeight / 2);
            ypos = marginTop + ymid + (width);

            height = getStringHeight(ylab);
            xmid = (int) Math.floor(marginLeft / 2);
            xpos = xmid - (xmid - height);

            if (height <= marginRight) {
                drawText(xpos, ypos, -90, ylab);
            }

            setFont(oldFont);
            setColor(oldColor);
        }
    }

    protected void setXAxisLabelPos(int x, int y, int rotation) {
    }

    protected void setYAxisLabelPos(int x, int y, int rotation) {
    }

    protected void setAxis(boolean x, boolean y) {
    }

    protected void setFont(int fontSize, int style, String fontFace) {
        font = new Font(fontFace, style, fontSize);
        g2d.setFont(font);
    }

    protected void setFont(int fontSize, String fontFace) {
        font = new Font(fontFace, Font.PLAIN, fontSize);
        g2d.setFont(font);
    }

    protected void setFont(Font newFont) {
        font = newFont;
        g2d.setFont(font);
    }

    protected void setStroke(int width) {
        g2d.setStroke(new BasicStroke(width, BasicStroke.CAP_BUTT, 1));
    }

    protected void setStroke(int width, int join) {
        g2d.setStroke(new BasicStroke(width, BasicStroke.CAP_BUTT, join));
    }

    protected int calcStringWidth(String input) {

        return 0;
    }

    protected int calcStringHeight(String input) {

        return 0;
    }

    protected void drawArea(int x, int y) {
    }

    protected void drawAxis() {
    }

    protected void drawRect(int xCoord1, int yCoord1, int width, int height, boolean filled) {
        if (filled) {
            g2d.fillRect(xCoord1, yCoord1, width, height);
        } else {
            g2d.draw(new Rectangle2D.Double(xCoord1, yCoord1, width, height));
        }

    }

    protected void drawLine(int xCoord1, int yCoord1, int xCoord2, int yCoord2) {
        g2d.draw(new Line2D.Double(xCoord1, yCoord1, xCoord2, yCoord2));
    }

    protected void drawPolyLine(int[] x, int[] y, int width) {
        g2d.draw(new Polygon(x, y, x.length));
    }

    protected void drawText(int xCoord, int yCoord, String text) {
        g2d.drawString(text, xCoord, yCoord);
    }

    protected void drawText(int xCoord, int yCoord, int angle, String text) {

        // Create a rotation transformation for the font.
        AffineTransform fontAT = new AffineTransform();

        // Derive a new font using a rotatation transform
        fontAT.rotate(Math.toRadians(angle));
        Font theDerivedFont = font.deriveFont(fontAT);

        // set the derived font in the Graphics2D context
        g2d.setFont(theDerivedFont);

        // Render a string using the derived font
        g2d.drawString(text, xCoord, yCoord);

        // put the original font back
        g2d.setFont(font);


    }

    protected void calcDrawArea() {
        drawHeight = graphHeight - marginTop - marginBottom;
        drawWidth = graphWidth - marginLeft - marginRight;
    }

    protected void calcXScaling(double maxX) {
        calcDrawArea();
        scalingX = (double) drawWidth / maxX;
    }

    protected void calcYScaling(double maxY) {
        calcDrawArea();
        scalingY = (double) drawHeight / maxY;
    }

    protected void calcScaling(double maxX, double maxY) {
        calcDrawArea();
        scalingY = (double) drawWidth / maxY;
        scalingX = (double) drawHeight / maxX;
    }

    protected void calcScalingDoNotRecalculateDrawArea(double maxX, double maxY) {
        scalingY = (double) drawWidth / maxY;
        scalingX = (double) drawHeight / maxX;
    }

    protected int getXCoord(double value) {
        return marginLeft + (int) Math.round(scalingX * value);
    }

    protected int getYCoord(double value) {
        return (drawHeight - (int) Math.round(scalingY * value)) + marginTop;

    }

    protected int getXCoord(int value) {
        return marginLeft + (int) Math.round(scalingX * value);
    }

    protected int getYCoord(int value) {
        return (drawHeight - (int) Math.round(scalingY * value)) + marginTop;
    }

    protected void setOpacity() {
    }

    protected int getStringWidth(String input) {
        // get metrics from the graphics
        FontMetrics metrics = g2d.getFontMetrics(font);
        // get the height of a line of text in this font and render context
        // int hgt = metrics.getHeight();
        // get the advance of my text in this font and render context
        int adv = metrics.stringWidth(input);
        return adv;
        // calculate the size of a box to hold the text with some padding.
        //Dimension size = new Dimension(adv+2, hgt+2);

    }

    protected int getStringHeight(String input) {
        // get metrics from the graphics
        FontMetrics metrics = g2d.getFontMetrics(font);
        // get the height of a line of text in this font and render context
        // int hgt = metrics.getHeight();
        // get the advance of my text in this font and render context
        int adv = metrics.getHeight();
        return adv;
        // calculate the size of a box to hold the text with some padding.
        //Dimension size = new Dimension(adv+2, hgt+2);

    }

    protected Dimension stringBoundingBox(String input) {
        // get metrics from the graphics
        FontMetrics metrics = g2d.getFontMetrics(font);
        // get the height of a line of text in this font and render context
        int hgt = metrics.getHeight();
        // get the advance of my text in this font and render context
        int adv = metrics.getHeight();

        // calculate the size of a box to hold the text with some padding.
        Dimension size = new Dimension(adv + 2, hgt + 2);
        return size;
    }

    protected int getMaxStringLength(String[] data) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < data.length; i++) {
            if (getStringWidth(data[i]) > max) {
                max = getStringWidth(data[i]);
            }
        }
        return max;
    }
    // private functions
}
