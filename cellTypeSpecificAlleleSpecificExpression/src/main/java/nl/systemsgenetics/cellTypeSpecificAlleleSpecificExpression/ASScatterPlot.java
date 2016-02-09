/*
 * Copyright (C) 2015 adriaan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.xml.bind.helpers.AbstractMarshallerImpl;


/**
 *
 * @author adriaan
 * adapted from harm jan, : imputation-tool/src/main/java/imputationtool/postprocessing/ScatterPlot.java
 */
public class ASScatterPlot {
    
    private BufferedImage bi;
    private Graphics2D g2d;
    private int graphHeight;
    private int graphWidth;
    private int drawWidth;
    private int drawHeight;
    private Color color;
    private Font font;
    private FontMetrics fontmetrics ;
    
    private int xStart;
    private int xEnd;
    private int yStart;
    private int yEnd;
    
    private int xLength;
    private int yLength;
    
    private double[] margins = {0.2, 0.1, 0.15, 0.15};
    
    public ASScatterPlot(int size ) {
        init(size, size);
    }

    protected void init(int width, int height) {
        
        
        graphHeight = height;
        graphWidth = width;
        
        //margins in the following order: top, right, bottom, left
        
        
        int  xTickLength = (int) (0.015 * ((double) graphWidth)); // ratio of tick length compared to canvas size
        int  yTickLength = (int) (0.015 * ((double) graphHeight));
        
        bi = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();

        color = new Color(255, 255, 255);
        font = new Font("Georgia", Font.PLAIN, 10);
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2d.setColor(color);
        g2d.fillRect(0, 0, graphWidth, graphHeight);
        
        
        
        fontmetrics = g2d.getFontMetrics();
        
        color = new Color(0,0,0);
        g2d.setColor(color);
        
        //definePlotBoundary
        xStart = (int)Math.floor( ((double) graphWidth) *  margins[3] ) ;
        xEnd =  (int)Math.floor(((double) graphWidth) *  (1 - margins[1]) );
        yStart = (int)Math.floor( ((double) graphHeight) *  (1 - margins[2]));
        yEnd = (int)Math.floor( ((double) graphHeight) *  margins[0]) ;
        
        //drawAxis
        g2d.draw(new Line2D.Double(xStart, yStart, xStart,yEnd));
        g2d.draw(new Line2D.Double(xStart, yStart, xEnd, yStart));
        
        xLength = xEnd - xStart;
        yLength = yStart - yEnd;
        //draw ticks and a value at every 0.1 interval
        for(int i = 0; i< 11;i++){
            //
            int xPos = xStart + (int) (( (double) xLength )* ( ((double) i) / 10) );  
            g2d.draw(new Line2D.Double(xPos, yStart , xPos, (yStart + yTickLength)));
            drawStringInMiddlePosition(String.format("%1.1f",(double) i / 10), xPos  , yStart + 3 * yTickLength);
            
            int yPos = yStart - (int) (( (double) yLength )* ( ((double) i) / 10) );  
            g2d.draw(new Line2D.Double(xStart, yPos , (xStart - xTickLength), yPos));
            drawStringInMiddlePosition(String.format("%1.1f",(double) i / 10), xStart - 3 * xTickLength, yPos );
            
        }
        
        //draw grid lines and a value at every 0.2 interval
        color = new Color(0, 0, 0, 30);
        g2d.setColor(color);
        
        for(int i = 0; i < 11 ;i++){
            
            int xPos = xStart + (int) (( (double) xLength )* ( ((double) i) / 10) );  
            g2d.draw(new Line2D.Double(xPos, yStart , xPos, yEnd));
            
            int yPos = yStart - (int) (( (double) yLength )* ( ((double) i) / 10) );  
            g2d.draw(new Line2D.Double(xStart, yPos , xEnd, yPos));

            
        }
        
        color = new Color(0, 0, 0);
        g2d.setColor(color);
        
       
        //print the x- axis
       
        drawStringInMiddlePosition("ASratio", xStart + (int) (( (double) xLength ) / 2) , yStart + (int)Math.floor((double) graphHeight * ( margins[2] / 2.0) ) );
        
        //print the y-axis, rotation and all
        AffineTransform orig = g2d.getTransform();
        g2d.rotate(-Math.PI * 0.5);
        drawStringInMiddlePosition("Phenotype proportion", -xStart - (int) (( (double) xLength ) / 2), yEnd - (int)Math.floor((double) graphWidth * ( margins[3] ) ));
        
        g2d.setTransform(orig);
    }

    //take the middle of the string and plot it at this exact position
    private void drawStringInMiddlePosition(String text, int x, int y){
        
        int widthOfString = fontmetrics.stringWidth(text); 
        int heightOfString = fontmetrics.getAscent();
        
        int xNew = (int) Math.floor( ((double) x ) - ((double) widthOfString) / 2.0 );
        int yNew = (int) Math.floor( ((double) y ) + ((double) heightOfString) / 2.0 );
        g2d.drawString(text, xNew, yNew);
    
    }
    
    
    
    void plot(int x, int y) {
    }

    // x == correlation
    // y == r2 score
    public void plot(double x, double y) {
        color = new Color(0,0,0, 150);
        g2d.setColor(color);
        
        int xPos = xStart + (int) Math.floor(((double) x) * ((double) xLength) );
        int yPos = yStart - (int) Math.floor(((double) y) * ((double) yLength) );
        
        g2d.fillRect(xPos , yPos , 4, 4);
        
    }

    public void draw(String outputFile) {
        
        //print the title, currently just the file name:
        drawStringInMiddlePosition(outputFile, xStart + (int) (( (double) xLength ) / 2) , yEnd - (int)Math.floor((double) graphHeight * ( margins[0] / 2.0) ) );
        
        try {
            javax.imageio.ImageIO.write(bi, "png", new File(outputFile));
        } catch (Exception e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }

    }
    
}
