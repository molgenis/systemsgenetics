/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.graphics;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harm-jan
 */
public class ScatterPlot {
    private Color[] colors = new Color[10];

    public ScatterPlot(){
        colors[0] = new Color(122,217,76,200);
        colors[1] = new Color(72,199,207,200);
        colors[2] = new Color(194,161,77,200);
        colors[3] = new Color(190,190,190,130);
        colors[4] = new Color(250,38,207,200);
        colors[5] = new Color(0,0,0,200);
    }

//    private void drawPlot(String outDir, int i, int j, HashMap<String, Double> groupPC1T, HashMap<String, Double> groupPC1W, HashMap<String, Double> groupPC2T, HashMap<String, Double> groupPC2W, HashMap<String, Vector<Double>> groupValuesPC1, HashMap<String, Vector<Double>> groupValuesPC2, double maleTPC1, double maleTPC2, Vector<Double> maleValsPC1, Vector<Double> maleValsPC2, double maleWPC1, double maleWPC2, HashMap<String, Vector<String>> genders, double maxXVal, double minXVal, double maxYVal, double minYVal) {
//
//        int width = 1000;
//        int height = 1000;
//        int margin = 100;
//
//        int x0 = margin;
//        int x1 = width - margin;
//
//        int innerWidth = x1 - x0;
//
//        int y0 = margin;
//        int y1 = height - margin;
//
//        int innerHeight = y1 - y0;
//
//        BufferedImage bimage = new BufferedImage(width, height + 100, BufferedImage.TYPE_INT_RGB);
//        Graphics2D g2d = bimage.createGraphics();
//
//        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
//        g2d.setColor(new Color(255, 255, 255));
//        g2d.fillRect(0, 0, width, height+100);
//        g2d.setColor(new Color(128,128,128));
//
//        g2d.drawLine(x0, (int) Math.floor((double)height/2), innerWidth+x0, (int) Math.floor((double)height/2));       // PC1-axis
//
//        g2d.drawString("PC "+(i), x0, (int) Math.floor((double)height/2) + 15);
//
//        g2d.drawLine((int) Math.floor((double) width / 2), y0, (int) Math.floor((double) width / 2), height-y0);       // PC2-axis
//
//        int strWidth = getStringWidth("BMI", g2d, g2d.getFont());
//        drawText((int) Math.floor((double) width / 2) + 12, y0 + strWidth, -90, "BMI", g2d, g2d.getFont());
//
//        DecimalFormat df = new DecimalFormat("#.###");
//        DecimalFormat df2 = new DecimalFormat("0.##E0");
//
//        for(int g=0; g<phenotypeGroups.size(); g++){
//            String group = phenotypeGroups.get(g);
//            Vector<Double> val1 = groupValuesPC1.get(group);
//            Vector<Double> val2 = groupValuesPC2.get(group);
//            Vector<String> gen = genders.get(group);
//
//            g2d.setStroke(new BasicStroke(1.0f));
//            g2d.setPaint(new Color(190,190,190,190));
//
//            int x = 580;
//
//            for (int v=0; v<val1.size(); v++) {
//
//                int plotX = x0 +     (int) Math.round((val1.get(v) - minXVal) / (maxXVal - minXVal) * (double) innerWidth) ;
//                int plotY = margin + (int) Math.round((val2.get(v) - minYVal) / (maxYVal - minYVal) * (double) innerHeight);
//
//                String sex = gen.get(v);
//                if(sex == null || !(sex.equals("0") || sex.equals("1")) ){
//                    g2d.drawRect(plotX, plotY, 12, 12);
//                    g2d.setColor(colors[g]);
//                    g2d.fillRect(plotX, plotY, 12, 12);
//                } else if(sex.equals("0")){
//                    g2d.drawOval(plotX, plotY, 12, 12);
//                    g2d.setColor(colors[g]);
//                    g2d.fillOval(plotX, plotY, 12, 12);
//                } else if (sex.equals("1")) {
//                    Polygon poly = new Polygon();
//                    poly.addPoint(plotX, plotY);
//                    poly.addPoint(plotX, plotY+12);
//                    poly.addPoint(plotX + 12, plotY);
//                    g2d.drawPolygon(poly);
//                    g2d.setColor(colors[g]);
//                    g2d.fillPolygon(poly);
//                }
//
//            }
//
//            int incr = 30 + 15 * (g + 0);
//            x = 580;
//            g2d.setColor( colors[g].darker() );
//            String groupName = phenotypeGroups.get(g);
//            strWidth = getStringWidth(groupName, g2d, g2d.getFont());
//            g2d.drawString(groupName, x - 10 - strWidth, y1 + incr);
//
//            double pc1t = groupPC1T.get(group);
//            double pc2t = groupPC2T.get(group);
//            double pc1w = groupPC1W.get(group);
//            double pc2w = groupPC2W.get(group);
//
//            incr = 30 + 15 * (g + 0);
//
//            if(pc1t < 0.001){
//                g2d.drawString(df2.format(pc1t), x,  y1 + incr);
//            } else {
//                g2d.drawString(df.format(pc1t), x,  y1 + incr);
//            }
//
//            x += 10 + getStringWidth("T-test PC1", g2d, g2d.getFont());
//            if(pc2t < 0.001){
//                g2d.drawString(df2.format(pc2t), x,  y1 + incr);
//            } else {
//                g2d.drawString(df.format(pc2t), x,  y1 + incr);
//            }
//
//            x += 10 + getStringWidth("T-test BMI", g2d, g2d.getFont());
//            if(pc1w < 0.001){
//                g2d.drawString(df2.format(pc1w), x,  y1 + incr);
//            } else {
//                g2d.drawString(df.format(pc1w), x,  y1 + incr);
//            }
//
//            x += 10 + getStringWidth("Wilcoxon PC1", g2d, g2d.getFont());
//            if(pc2w < 0.001){
//                g2d.drawString(df2.format(pc2w), x,  y1 + incr);
//            } else {
//                g2d.drawString(df.format(pc2w), x,  y1 + incr);
//            }
//
//        }
//
//
//
//        int incr = 30 + (phenotypeGroups.size() * 15);
//        int x = 580;
//        g2d.setColor(Color.black);
//
//        strWidth = getStringWidth("sex", g2d, g2d.getFont());
//        g2d.drawString("sex",  x - 10 - strWidth, y1 + incr);
//
//
//        if(maleTPC1 < 0.001){
//            g2d.drawString(df2.format(maleTPC1), x,  y1 + incr);
//        } else {
//            g2d.drawString(df.format(maleTPC1), x,  y1 + incr);
//        }
//
//        x += 10 + getStringWidth("T-test PC1", g2d, g2d.getFont());
//        if(maleTPC2 < 0.001){
//            g2d.drawString(df2.format(maleTPC2), x,  y1 + incr);
//        } else {
//            g2d.drawString(df.format(maleTPC2), x,  y1 + incr);
//        }
//
//        x += 10 + getStringWidth("T-test BMI", g2d, g2d.getFont());
//        if(maleWPC1 < 0.001){
//            g2d.drawString(df2.format(maleWPC1), x,  y1 + incr);
//        } else {
//            g2d.drawString(df.format(maleWPC1), x,  y1 + incr);
//        }
//
//        x += 10 + getStringWidth("Wilcoxon PC1", g2d, g2d.getFont());
//        if(maleWPC2 < 0.001){
//            g2d.drawString(df2.format(maleWPC2), x,  y1 + incr);
//        } else {
//            g2d.drawString(df.format(maleWPC2), x,  y1 + incr);
//        }
//
//        incr = 15;
//        x = 580;
//        g2d.drawString("T-test PC1", x,    y1 + incr);
//        x += 10 + getStringWidth("T-test PC1", g2d, g2d.getFont());
//
//        g2d.drawString("T-test BMI", x,    y1 + incr);
//
//        x += 10 + getStringWidth("T-test BMI", g2d, g2d.getFont());
//        g2d.drawString("Wilcoxon PC1", x,  y1 + incr);
//
//        x += 10 + getStringWidth("Wilcoxon PC1", g2d, g2d.getFont());
//        g2d.drawString("Wilcoxon BMI", x,  y1 + incr);
//
//        int plotY = y1 + 15;
//        int plotX = x0;
//
//        g2d.setColor(new Color(0,0,0));
//        g2d.fillOval(plotX,  plotY - 8, 8, 8);
//        g2d.drawString("Sex 0 (Female)", plotX + 12, plotY);
//
//        plotY = y1 + 30;
//        Polygon poly = new Polygon();
//        poly.addPoint(plotX, plotY - 8);
//        poly.addPoint(plotX, plotY);
//        poly.addPoint(plotX + 8, plotY-8);
//        g2d.fillPolygon(poly);
//        g2d.drawString("Sex 1 (Male)", plotX + 12, plotY);
//
//        plotY = y1 + 45;
//        g2d.fillRect(plotX, plotY - 8, 8, 8);
//        g2d.drawString("Sex Undetermined", plotX + 12, plotY);
//
//        try {
//            javax.imageio.ImageIO.write(bimage, "png", new File(outDir+"Comp"+i+"-BMI.png"));
//        } catch (IOException e) {
//            System.out.println(e.getMessage());
//            e.printStackTrace();
//        }
//    }

     protected int getStringWidth(String input, Graphics2D g2d, Font font) {
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

    protected int getStringHeight(String input, Graphics2D g2d, Font font) {
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

    protected Dimension stringBoundingBox(String input, Graphics2D g2d, Font font) {
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

    protected void drawText(int xCoord, int yCoord, int angle, String text, Graphics2D g2d, Font font) {

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

    public void draw(double[] x, double[] y, String xDesc, String yDesc, String subtitle, String outLoc) {

        double maxX = Primitives.max(x);
        double minX = Primitives.min(x);
        
        double maxY = Primitives.max(y);
        double minY = Primitives.min(y);

        int width   = 1000;
        int height  = 1000;
        int margin  = 100;

        int x0      = margin;
        int x1      = width - margin;

        int innerWidth = x1 - x0;

        int y0 = margin;
        int y1 = height - margin;

        int innerHeight = y1 - y0;

        BufferedImage bimage = new BufferedImage(width, height + 100, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = bimage.createGraphics();

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height+100);
        g2d.setColor(new Color(128,128,128));

        DecimalFormat df = new DecimalFormat("#.###");

        double xAbs = Math.abs(minX) + Math.abs(maxX);
        double distOfMinToMid = Math.abs(minX) / xAbs;
        double distOfMaxToMid = Math.abs(maxX) / xAbs;
        int distToYAxis = margin + (int) Math.floor( distOfMinToMid * innerWidth );
        
        if(minY >= 0){
            g2d.drawLine(x0, y1, innerWidth+x0, y1);       // x-axis
            g2d.drawString(xDesc, x0, y1 - 12);
            g2d.drawString(df.format(minX), x0, y1 + 20);
            g2d.drawString(df.format(maxX), x1 - getStringWidth(df.format(maxX), g2d, g2d.getFont()), y1 + 20);
            g2d.drawString("0", distToYAxis + 10, y1 + 20); // plot 0;
            g2d.drawString(df.format(minY), distToYAxis + 10, y1 - 12); // plot 0;
            minY = 0;
        } else {
            g2d.drawLine(x0, (int) Math.floor((double)height/2), innerWidth+x0, (int) Math.floor((double)height/2));       // x-axis
            g2d.drawString(xDesc, x0, (int) Math.floor((double)height/2) - 15);
            g2d.drawString(df.format(minX), x0, (int) Math.floor((double)height/2) + 20);
            g2d.drawString("0", distToYAxis + 10, (int) Math.floor((double)height/2) + 20); // plot 0;
            g2d.drawString(df.format(minY), distToYAxis + 10, (int) Math.floor((double)height/2) - 20); // plot 0;
        }
    
        g2d.drawLine(distToYAxis, y0, distToYAxis, y1 );       // y-axis
        

        int strWidth = getStringWidth(yDesc, g2d, g2d.getFont());
        drawText(distToYAxis - 12, y0 + strWidth, -90, yDesc, g2d, g2d.getFont());

        Font f = g2d.getFont();

        Font font = new Font("Calibri", Font.BOLD, 30);
        g2d.setFont(font);
        String title = xDesc + " - " + yDesc;
        int theight = getStringHeight(title, g2d, font);
        g2d.drawString(title, x0, 0 + 10 + theight);

        font = new Font("Georgia", Font.BOLD, 20);
        g2d.setFont(font);
        g2d.drawString(subtitle, x0, 0 + 10 + theight + 30);


        g2d.setFont(f);
        //        DecimalFormat df = new DecimalFormat("#.###");
//        DecimalFormat df2 = new DecimalFormat("0.###E0");

        for(int i=0; i<x.length; i++){
            int plotX = x0 + (int) Math.round((x[i] - minX) / (maxX - minX) * (double) innerWidth);
            int plotY = y1 - (int) Math.round((y[i] - minY) / (maxY - minY) * (double) innerHeight);

            g2d.drawRect(plotX - 6, plotY - 6, 12, 12);
            g2d.setColor(colors[5]);
            g2d.fillRect(plotX - 6, plotY - 6, 12, 12);

        }


        try {
            javax.imageio.ImageIO.write(bimage, "png", new File(outLoc+xDesc+"-"+yDesc+".png"));
        } catch (IOException e) {
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }

    public void draw(double[] x, double[] y, String xDesc, String yDesc, int[] gender, int[] group, String subtitle, String outLoc) {

        double maxX = Primitives.max(x);
        double minX = Primitives.min(x);

        double maxY = Primitives.max(y);
        double minY = Primitives.min(y);

        int width   = 1000;
        int height  = 1000;
        int margin  = 100;

        int x0      = margin;
        int x1      = width - margin;

        int innerWidth = x1 - x0;

        int y0 = margin;
        int y1 = height - margin;

        int innerHeight = y1 - y0;

        BufferedImage bimage = new BufferedImage(width, height + 100, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = bimage.createGraphics();

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height+100);
        g2d.setColor(new Color(128,128,128));

        DecimalFormat df = new DecimalFormat("#.###");

        double xAbs = Math.abs(minX) + Math.abs(maxX);
        double distOfMinToMid = Math.abs(minX) / xAbs;
        double distOfMaxToMid = Math.abs(maxX) / xAbs;
        int distToYAxis = margin + (int) Math.floor( distOfMinToMid * innerWidth );

        if(minY >= 0){
            minY = 0;
            g2d.drawLine(x0, y1, innerWidth+x0, y1);       // x-axis
            g2d.drawString(xDesc, x0, y1 - 12);
            g2d.drawString(df.format(minX), x0, y1 + 20);
            g2d.drawString(df.format(maxX), x1 - getStringWidth(df.format(maxX), g2d, g2d.getFont()), y1 + 20);
            g2d.drawString("0", distToYAxis + 10, y1 + 20); // plot 0;
            g2d.drawString(df.format(minY), distToYAxis + 10, y1 - 12); // plot 0;
        } else {
            g2d.drawLine(x0, (int) Math.floor((double)height/2), innerWidth+x0, (int) Math.floor((double)height/2));       // x-axis
            g2d.drawString(xDesc, x0, (int) Math.floor((double)height/2) - 15);
            g2d.drawString(df.format(minX), x0, (int) Math.floor((double)height/2) + 20);
            g2d.drawString("0", distToYAxis + 10, (int) Math.floor((double)height/2) + 20); // plot 0;
            g2d.drawString(df.format(minY), distToYAxis + 10, (int) Math.floor((double)height/2) - 20); // plot 0;
        }

        g2d.drawLine(distToYAxis, y0, distToYAxis, y1 );       // y-axis

        int strWidth = getStringWidth(yDesc, g2d, g2d.getFont());
        drawText(distToYAxis - 12, y0 + strWidth, -90, yDesc, g2d, g2d.getFont());
        drawText(distToYAxis + 10, y0 + 10, 0, df.format(maxY), g2d, g2d.getFont());

        Font f = g2d.getFont();

        Font font = new Font("Calibri", Font.BOLD, 30);
        g2d.setFont(font);
        String title = xDesc + " - " + yDesc;
        int theight = getStringHeight(title, g2d, font);
        g2d.drawString(title, x0, 0 + 10 + theight);

        font = new Font("Georgia", Font.BOLD, 20);
        g2d.setFont(font);
        g2d.drawString(subtitle, x0, 0 + 10 + theight + 30);


        g2d.setFont(f);
        //        DecimalFormat df = new DecimalFormat("#.###");
//        DecimalFormat df2 = new DecimalFormat("0.###E0");

        for(int i=0; i<x.length; i++){
            int plotX = x0 + (int) Math.round((x[i] - minX) / (maxX - minX) * (double) innerWidth);
            int plotY = y1 - (int) Math.round((y[i] - minY) / (maxY - minY) * (double) innerHeight);

            Integer sex = gender[i];
            Integer g = group[i];
            
            if( sex == null || !(sex == 0 || sex == 1 ) ){
                g2d.drawRect(plotX -6, plotY-6, 12, 12);
                g2d.setColor(colors[g]);
                g2d.fillRect(plotX, plotY, 12, 12);
            } else if( sex == 0 ){
                g2d.drawOval(plotX-6, plotY-6, 12, 12);
                g2d.setColor(colors[g]);
                g2d.fillOval(plotX-6, plotY-6, 12, 12);
            } else if ( sex == 1 ) {
                Polygon poly = new Polygon();
                poly.addPoint(plotX-6, plotY-6);
                poly.addPoint(plotX-6, plotY-6+12);
                poly.addPoint(plotX-6 + 12, plotY-6);
                g2d.drawPolygon(poly);
                g2d.setColor(colors[g]);
                g2d.fillPolygon(poly);
            }
            

        }


        try {
            javax.imageio.ImageIO.write(bimage, "png", new File(outLoc+xDesc+"-"+yDesc+".png"));
        } catch (IOException e) {
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }
}
