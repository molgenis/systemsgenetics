/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.pcaoptimum;

import eqtlmappingpipeline.graphics.Graphics;

/**
 *
 * @author harmjan
 */
public class PCAOptimumPlot extends Graphics {

    PCAOptimumPlot(int i, int i0, boolean b, String out) {
        super(i, i0, b, out);
    }
    public void plot(double[] cisx, double[] cisy, double[] transx, double[] transy, int cisOrig, int cisNew, int cisShared, int transOrig, int transNew, int transShared ){
        
        int width  = graphWidth;
        int height = graphHeight;
        int margin = 50;
        
        int plotwidth  = (width - (3 * margin))/ 2; 
        int plotheight = (height - (3 * margin))/ 2; 
        
        scatterplot(cisx, cisy, plotwidth, plotheight, margin, margin);
        scatterplot(transx, transy, plotwidth, plotheight, margin, (2*margin)+plotheight);
        
    }
    
    private void scatterplot(double[] x, double[] y, int plotwidth, int plotheight, int offsetx, int offsety){
        // cisplot
        double maxX = max(x);
        double minX = min(x);
        double maxY = max(y);
        double minY = min(y);
        
        double rangeX = getrange(minX, maxX);
        double rangeY = getrange(minY, maxY);
//        double xPixelsPerUnit = plotwidth / rangeX;
//        double yPixelsPerUnit = plotwidth / rangeX;
//        
        // draw axis:
        
        // x-axis + y-axis
        setColor(0, 0, 0, 255);
        int xCoord = (int) (Math.floor( (Math.abs(0-minX) / rangeX) * plotwidth)) + offsetx;
        int yCoord = (int) (Math.floor( (Math.abs(0-minY) / rangeY) * plotheight))+ offsety;
        drawLine(offsetx, yCoord, offsetx+plotwidth, yCoord);
        drawLine(xCoord, offsety, xCoord, plotheight+offsety);
        
        // plot dots
        setColor(0, 0, 0, 128);
        for(int i=0; i<x.length; i++){
            xCoord = (int)Math.floor( (Math.abs(x[i]-minX)/rangeX)*plotwidth );
            yCoord = (int)Math.floor( (Math.abs(y[i]-minY)/rangeY)*plotwidth );
            
            g2d.fillOval(xCoord, yCoord, 2, 2);
        }
        
    }
    
    public double max(double[] x){
        double max = Double.MIN_VALUE;
        for(double d: x){
            if(d>max){
                max = d;
            }
        }
        return max;
        
    }
    
    public double min(double[] x){
        double min = Double.MAX_VALUE;
        for(double d: x){
            if(d<min){
                min = d;
            }
        }
        return min;

    }
    
    public double getrange(double min, double max){
        double range = 0;
        if(min < 0 && max < 0){
           range = Math.abs(min) - Math.abs(max);
        } 
        if(min < 0 && max > 0){
           range = Math.abs(min) + max;
        }
        if(min > 0 && max > 0){
           range = max - min;
        }
        return range;
    }
}
