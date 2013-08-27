/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.graphics;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.Goniometry;

/**
 *
 * @author harm-jan
 */
public class VennDiagram {

    private int width = 0;
    private int height = 0;
    private int spacer = 50;
    private BufferedImage bimage;
    private Graphics2D g2d;
    private int requestedsize;

    public VennDiagram(int size) {
        this.requestedsize = size;
    }

    public void plot(double sizeSetA, double sizeSetB, double overlap, String setAName, String setBName) {


        // first calculate the actual minimal size of the image:
        int margin = 50;
        bimage = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        g2d = bimage.createGraphics();
        FontMetrics fontmetrics = g2d.getFontMetrics();

        width = (margin * 2) + 10 + fontmetrics.stringWidth(setAName) + fontmetrics.stringWidth(setBName) + (requestedsize - (margin * 2));
        height = width / 2;

        bimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        g2d = bimage.createGraphics();

//        System.out.println(requestedsize + "\t" + width);
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height);
        g2d.setColor(new Color(0, 0, 0));

        double d = getDistance(sizeSetA, sizeSetB, overlap);

        int widthofstringA = fontmetrics.stringWidth(setAName);
        int widthofstringB = fontmetrics.stringWidth(setBName);

        int meanwidthOfString = widthofstringA + widthofstringB;

        int maxCircleSize = (requestedsize - (margin * 2)) / 2;

        double width1radians = (getR(sizeSetA) * 2);
        double width2radians = (getR(sizeSetB) * 2);
        double maxwidth = width1radians;
        if (width2radians > width1radians) {
            maxwidth = width2radians;
        }

        int circle1sizepixels = (int) Math.ceil(maxCircleSize * (width1radians / maxwidth));
        int circle2sizepixels = (int) Math.ceil(maxCircleSize * (width2radians / maxwidth));
        int distancepixels = (int) Math.ceil(maxCircleSize * (d / maxwidth));

        // label
        double radius1 = (double) circle1sizepixels / 2;
        double radius2 = (double) circle2sizepixels / 2;

        //int totaldistance = distance + (int)Math.ceil(height1/2) + (int)Math.ceil(height2/2);
//        Color black = new Color(0, 0, 0);
//        Color white = new Color(255, 255, 255);
//        Color grey = new Color(255, 255, 255, 125);
        Color black = new Color(137, 142, 216);
        Color grey = new Color(73, 78, 134, 125);
//        Color white = new Color(220, 220, 220);//, 125);

        int totaldistance = distancepixels + (int) Math.ceil(radius1) + (int) Math.floor(radius2);
        int remaining = width - totaldistance;
        int marginx = remaining / 2;

        // circle 1
        g2d.setColor(black);
        g2d.fillOval(marginx, (height / 2) - (int) Math.floor(radius1), circle1sizepixels, circle1sizepixels);

        int circle1origin = marginx + (int) Math.ceil(radius1);
        int circle2origin = circle1origin + distancepixels;
        int circle2start = circle2origin - (int) Math.ceil(radius2);
        g2d.setColor(grey);
        g2d.fillOval(circle2start, (height / 2) - (int) Math.floor(radius2), circle2sizepixels, circle2sizepixels);

        g2d.setColor(black);
        g2d.drawOval(circle2start, (height / 2) - (int) Math.floor(radius2), circle2sizepixels, circle2sizepixels);

        g2d.setStroke(new BasicStroke(2f));




        g2d.setColor(black);
        int originX1 = (int) Math.ceil(marginx + radius1);
        int originY1 = (int) Math.ceil(height / 2.0d);

//        System.out.println(margin);

        int originX2 = (int) Math.ceil(circle2start + radius2);
        int originY2 = (int) Math.ceil((height / 2));
//        System.out.println(originY2 + "\t" + (height / 2) + "");
//        for (int deg = 0; deg < 360; deg += 5) {
        int deg = 225;
        Pair<Integer, Integer> xy0 = Goniometry.calcPosOnCircle(radius1, 0, 0, deg);
        Pair<Integer, Integer> xy1 = Goniometry.calcPosOnCircle(radius1 + 10, 0, 0, deg);

//        g2d.drawLine(originX1 + xy0.getLeft(),
//                originY1 + xy0.getRight(),
//                originX1 + xy1.getLeft(),
//                originY1 + xy1.getRight());

        String setADesc = "(" + sizeSetA + ")";
        int setADescWidth = fontmetrics.stringWidth(setADesc);
        int descheight = fontmetrics.getHeight() * 2;
        g2d.drawString(setAName, originX1 + xy1.getLeft() - widthofstringA - 5, originY1 + xy1.getRight() - (descheight / 2));
        g2d.drawString(setADesc, originX1 + xy1.getLeft() - setADescWidth - 5, originY1 + xy1.getRight() + fontmetrics.getHeight() - (descheight / 2));

        deg = 270 + 45;
//        deg =0;
        xy0 = Goniometry.calcPosOnCircle(radius2, 0, 0, deg);
        xy1 = Goniometry.calcPosOnCircle(radius2 + 10, 0, 0, deg);

//        g2d.drawLine(originX2 + xy0.getLeft(),
//                originY2 + xy0.getRight(),
//                originX2 + xy1.getLeft(),
//                originY2 + xy1.getRight());

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
        int overlapmid = marginx + (int) Math.floor(radius1) + dd1 + (int) Math.floor(ddx / 2);




        int halfmaxcircleheight = maxCircleSize / 2;

//        g2d.drawLine(overlapmid, (height / 2), overlapmid, (height / 2) - halfmaxcircleheight - 10);
//        g2d.fillOval(overlapmid - 5, (height / 2) - 5, 10, 10);
        g2d.drawString(overlapstr, overlapmid - overlapstrwidthhalway, (height / 2) - halfmaxcircleheight - (fontmetrics.getHeight() * 2));
        g2d.drawString(overlapdesc, overlapmid - overlapdescstrwidthhalway, (height / 2) - halfmaxcircleheight - fontmetrics.getHeight());

    }

    public void draw(String loc) throws IOException {

        javax.imageio.ImageIO.write(bimage, "png", new File(loc));

    }

    private double calc_overlap(double r1, double r2, double d) {
        double alpha = 2 * Math.acos((Math.pow(d, 2) + Math.pow(r1, 2) - Math.pow(r2, 2)) / (2 * r1 * d));
        double beta = 2 * Math.acos((Math.pow(d, 2) + Math.pow(r2, 2) - Math.pow(r1, 2)) / (2 * r2 * d));
        double s = (0.5 * Math.pow(r1, 2) * (alpha - Math.sin(alpha))) + (0.5 * Math.pow(r2, 2) * (beta - Math.sin(beta)));
        return s;
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
        while (Math.abs(lower - upper) < .0000001 ) {
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
