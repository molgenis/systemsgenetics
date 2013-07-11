/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.graphics;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import umcg.genetica.containers.Pair;

/**
 *
 * @author harm-jan
 */
public class MultiVenn {

    int width = 0;
    int height = 0;
    int spacer = 50;
    BufferedImage bimage;
    Graphics2D g2d;

    public static void main(String[] args) {
        String[] setnames = new String[3];
        setnames[0] = "SetA";
        setnames[1] = "SetB";
        setnames[2] = "SetC";

        double[] setsizes = new double[3];
        setsizes[0] = 100;
        setsizes[1] = 100;
        setsizes[2] = 25;

        double[][] overlaps = new double[3][3];
        overlaps[0][1] = 50;
        overlaps[0][2] = 10;
        overlaps[1][2] = 10;
        MultiVenn m = new MultiVenn();
        m.plot(setnames, setsizes, overlaps);
        System.out.println("drawing");
        m.draw("d:\\work\\multivenn.png");
    }
    private int requestedsize;
    private FontMetrics fontmetrics;
    private int margin;
    private int plotwidth;
    private int plotheight;

    public void plot(String[] setnames, double[] setsizes, double[][] overlaps) {
        int columns = setnames.length - 1;
        System.out.println(columns);

        this.requestedsize = 200;
        margin = 75;
        bimage = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        g2d = bimage.createGraphics();
        fontmetrics = g2d.getFontMetrics();


        plotwidth = 200;
        plotheight = 200;
        width = (columns * plotwidth) + ((columns + 1) * margin);
        height = width;

        bimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        g2d = bimage.createGraphics();

//        System.out.println(requestedsize + "\t" + width);
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height);
        g2d.setColor(new Color(0, 0, 0));


        for (int i = 0; i < columns; i++) {
            for (int j = i + 1; j < setnames.length; j++) {
                int row = i;
                int col = j - 1;
                int offsetY = margin + (row * plotwidth) + (row * margin);
                int offsetX = margin + (col * plotwidth) + (col * margin);

                System.out.println(offsetX + "\t" + offsetY + "\t" + setnames[i] + "\t" + setnames[j]);
                plot(setsizes[i], setsizes[j], overlaps[i][j], setnames[i], setnames[j], offsetX, offsetY);
            }
        }


        System.out.println(width + "\t" + height);


//        width = (margin * 2) + 10
//                + fontmetrics.stringWidth(setAName)
//                + fontmetrics.stringWidth(setBName)
//                + (requestedsize - (margin * 2));
//        height = width / 2;

    }

    private void plot(double sizeSetA, double sizeSetB, double overlap, String setAName, String setBName, int offsetX, int offsetY) {
        double distanceradians = getDistance(sizeSetA, sizeSetB, overlap);

        g2d.setStroke(new BasicStroke(2f));

        int widthofstringA = fontmetrics.stringWidth(setAName);
        int widthofstringB = fontmetrics.stringWidth(setBName);

        int maxCircleSize = plotwidth / 2;

        System.out.println("max circle\t" + maxCircleSize);
        double width1radians = (getR(sizeSetA) * 2);
        double width2radians = (getR(sizeSetB) * 2);
        double maxwidthradians = width1radians;
        if (width2radians > width1radians) {
            maxwidthradians = width2radians;
        }

        int circle1sizepixels = (int) Math.ceil(maxCircleSize * (width1radians / maxwidthradians));
        int circle2sizepixels = (int) Math.ceil(maxCircleSize * (width2radians / maxwidthradians));
        int distancepixels = (int) Math.ceil(maxCircleSize * (distanceradians / maxwidthradians));

        System.out.println(circle1sizepixels + "\t" + circle2sizepixels + "\t" + distancepixels);
        // label
        double radius1 = (double) circle1sizepixels / 2;
        double radius2 = (double) circle2sizepixels / 2;

        //int totaldistance = distance + (int)Math.ceil(height1/2) + (int)Math.ceil(height2/2);
        Color black = new Color(0, 0, 0);
        Color white = new Color(255, 255, 255);
        Color grey = new Color(255, 255, 255, 125);

        int totaldistance = distancepixels + (int) Math.ceil(radius1) + (int) Math.ceil(radius2);
        int remaining = plotwidth - totaldistance;
        int marginx = remaining / 2;

        System.out.println(totaldistance + "\t" + remaining + "\t" + marginx);
        // circle 1
        g2d.setColor(black);
        g2d.fillOval(offsetX + marginx, offsetY + (plotheight / 2) - (int) Math.floor(radius1), circle1sizepixels, circle1sizepixels);

        int circle1origin = offsetX + marginx + (int) Math.ceil(radius1);
        int circle2origin = circle1origin + distancepixels;
        int circle2start = circle2origin - (int) Math.ceil(radius2);
        g2d.setColor(grey);
        g2d.fillOval(circle2start, offsetY + (plotheight / 2) - (int) Math.floor(radius2), circle2sizepixels, circle2sizepixels);

        g2d.setColor(black);
        g2d.drawOval(circle2start, offsetY + (plotheight / 2) - (int) Math.floor(radius2), circle2sizepixels, circle2sizepixels);






        g2d.setColor(black);
        int originX1 = offsetX + (int) Math.ceil(marginx + radius1);
        int originY1 = offsetY + (int) Math.ceil(plotheight / 2);

//        System.out.println(margin);

        int originX2 = (int) Math.ceil(circle2start + radius2);
        int originY2 = offsetY + (int) Math.ceil((plotheight / 2));
//        System.out.println(originY2 + "\t" + (height / 2) + "");
//        for (int deg = 0; deg < 360; deg += 5) {
        int deg = 225;
        Pair<Integer, Integer> xy0 = calcPosOnCircle(radius1, 0, 0, deg);
        Pair<Integer, Integer> xy1 = calcPosOnCircle(radius1 + 10, 0, 0, deg);

        g2d.drawLine(originX1 + xy0.getLeft(),
                originY1 + xy0.getRight(),
                originX1 + xy1.getLeft(),
                originY1 + xy1.getRight());

        String setADesc = "(" + sizeSetA + ")";
        int setADescWidth = fontmetrics.stringWidth(setADesc);
        int descheight = fontmetrics.getHeight() * 2;
        g2d.drawString(setAName, originX1 + xy1.getLeft() - widthofstringA - 5, originY1 + xy1.getRight() - (descheight / 2));
        g2d.drawString(setADesc, originX1 + xy1.getLeft() - setADescWidth - 5, originY1 + xy1.getRight() + fontmetrics.getHeight() - (descheight / 2));

        deg = 270 + 45;
//        deg =0;
        xy0 = calcPosOnCircle(radius2, 0, 0, deg);
        xy1 = calcPosOnCircle(radius2 + 10, 0, 0, deg);

        g2d.drawLine(originX2 + xy0.getLeft(),
                originY2 + xy0.getRight(),
                originX2 + xy1.getLeft(),
                originY2 + xy1.getRight());

        String setBDesc = "(" + sizeSetB + ")";

        g2d.drawString(setBName, originX2 + xy1.getLeft() + 5, originY2 + xy1.getRight() - (descheight / 2));
        g2d.drawString(setBDesc, originX2 + xy1.getLeft() + 5, originY2 + xy1.getRight() + fontmetrics.getHeight() - (descheight / 2));

        String overlapstr = "Overlap";
        String overlapdesc = "(" + overlap + ")";

        int overlapstrwidth = fontmetrics.stringWidth(overlapstr);
        int overlapstrwidthhalway = overlapstrwidth / 2;
        int overlapdescstrwidth = fontmetrics.stringWidth(overlapdesc);
        int overlapdescstrwidthhalway = overlapdescstrwidth / 2;


        int dd2 = (distancepixels - (int) Math.floor(radius1));
        int dd1 = (distancepixels - (int) Math.floor(radius2));
        int ddx = distancepixels - dd1 - dd2;
        int overlapmid = offsetX + marginx + (int) Math.floor(radius1) + dd1 + (int) Math.floor(ddx / 2);


        System.out.println(overlapmid);


        int halfmaxcircleheight = maxCircleSize / 2;

        g2d.drawLine(overlapmid, offsetY + (plotheight / 2), overlapmid, offsetY + (plotheight / 2) - halfmaxcircleheight - 10);
        g2d.fillOval(overlapmid - 5, offsetY + (plotheight / 2) - 5, 10, 10);
        g2d.drawString(overlapstr, overlapmid - overlapstrwidthhalway, offsetY + (plotheight / 2) - halfmaxcircleheight - (fontmetrics.getHeight() * 2));
        g2d.drawString(overlapdesc, overlapmid - overlapdescstrwidthhalway, offsetY + (plotheight / 2) - halfmaxcircleheight - fontmetrics.getHeight());

    }

    public void draw(String loc) {
        try {
            javax.imageio.ImageIO.write(bimage, "png", new File(loc));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private double calc_overlap(double r1, double r2, double d) {
        double alpha = 2 * Math.acos((Math.pow(d, 2) + Math.pow(r1, 2) - Math.pow(r2, 2)) / (2 * r1 * d));
        double beta = 2 * Math.acos((Math.pow(d, 2) + Math.pow(r2, 2) - Math.pow(r1, 2)) / (2 * r2 * d));
        double s = (0.5 * Math.pow(r1, 2) * (alpha - Math.sin(alpha))) + (0.5 * Math.pow(r2, 2) * (beta - Math.sin(beta)));
        return s;
    }

    private Pair<Integer, Integer> calcPosOnCircle(double r, double originX, double originY, double angle) {
//        x = cx + r * cos(a)
//y = cy + r * sin(a)
        double radians = angle * Math.PI / 180;
        Integer x = (int) Math.floor(originX + (r * Math.cos(radians)));
        Integer y = (int) Math.floor(originY + (r * Math.sin(radians)));

        return new Pair<Integer, Integer>(x, y);
    }

    private double getR(double r) {
        return Math.sqrt((double) r / Math.PI);
    }

    private double getDistance(double sizeSetA, double sizeSetB, double overlap) {
        double r1 = getR(sizeSetA);
        double r2 = getR(sizeSetB);

        double minOverlap = Math.abs(r1 - r2);
        double maxOverlap = r1 + r2;

        double upper = maxOverlap;
        double lower = minOverlap;
        double testpos = 0;
        double circleOverlap = 0;
        double allowed_difference = overlap / 100;

        int i = 0;
        while (lower != upper) {
            testpos = (lower + upper) / 2;
            circleOverlap = calc_overlap(r1, r2, testpos);
            if (Math.abs(circleOverlap - overlap) <= allowed_difference) {
                // position found

                break;
            } else {
                if (circleOverlap > overlap) {
                    lower = testpos;
                } else {
                    upper = testpos;
                }
            }

//            System.out.println(r1 + "\t" + r2 + "\t" + testpos + "\t" + circleOverlap + "\t" + overlap + "\t" + lower + "\t" + upper + "\t" + Math.abs(circleOverlap - overlap) + "\t" + allowed_difference);
            if (i == 100) {
                break;
            }
            i++;

        }

        return testpos;
    }
}
