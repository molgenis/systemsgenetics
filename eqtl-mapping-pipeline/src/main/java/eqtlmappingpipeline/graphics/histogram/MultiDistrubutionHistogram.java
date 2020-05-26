/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.graphics.histogram;

import eqtlmappingpipeline.mixupmapper.stat.dist.Bin;
import eqtlmappingpipeline.mixupmapper.stat.dist.DiscreteDist;
import java.awt.Color;

/**
 *
 * @author harm-jan
 */
public class MultiDistrubutionHistogram extends Histogram {

    public MultiDistrubutionHistogram(int width, int height){
        super(width, height);
    }

    public MultiDistrubutionHistogram(int width, int height, boolean outputPDF, String outputLoc){
        super(width, height, outputPDF, outputLoc);
    }

    

    public void plot(DiscreteDist[] dists){
         if (graphWidth > 3 * 50 && graphHeight > 3 * 50) {
            setMargins(100);


            boolean sameXRange = true;
            double prevMax = 0.0;
            double prevMin = 0.0;
            double newRangeMax = Double.MIN_VALUE;
            double newRangeMin = Double.MAX_VALUE;
            for (int i = 0; i < dists.length; i++) {
                // check if all dists have the same range on the x-axis
                double minX = dists[i].getMin();
                double maxX = dists[i].getMax();

                if(minX < newRangeMin){
                    newRangeMin = minX;
                }

                if(maxX > newRangeMax){
                    newRangeMax = maxX;
                }

                if (i == 0) {
                    prevMax = maxX;
                    prevMin = minX;
                } else {
                    if (prevMax != maxX) {
                        sameXRange = sameXRange && false;
                    }
                    if (prevMin != minX) {
                        sameXRange = sameXRange && false;
                    }
                }
            }

            if (sameXRange) {
                 // get the bounds for this plot

                int yRange          = 0;

                double maxY = Double.MIN_VALUE;
                double minY = 0;


                for(int i=0; i<dists.length; i++){
                    dists[i].resetIterator();
                    Bin currentBin = dists[i].getNext();
                    int count = currentBin.getCount();
                    if(count > maxY){
                        maxY = count;
                    }
                    while(dists[i].hasNext()){
                        currentBin = dists[i].getNext();
                        count = currentBin.getCount();
                        if(count > maxY){
                            maxY = count;
                        }
                      
                    }

                    
                }

                yRange = (int) maxY;

                for (int i = 0; i < dists.length; i++) {
                    int numBins = dists[i].getNumBins();

                    double barsize = Math.floor(drawWidth / numBins);
                    drawWidth = (int) Math.floor(barsize * numBins);

                    calcScalingDoNotRecalculateDrawArea(numBins, yRange);
                    determineColor(i);
                    plotDist(dists[i]);
                }

                // draw the axis and axis labels
                // draw the axis and axis labels
                double yticksevery = 0.05;
                if(yRange < 1 ){
                    yticksevery = 0.1;
                } else if( yRange < 10){
                    yticksevery = 0.5;
                } else if( yRange < 100) {
                    yticksevery = 5;
                } else if( yRange < 1000){
                    yticksevery = 50;
                } else if( yRange < 10000){
                    yticksevery = 500;
                } else if( yRange < 1000000){
                    yticksevery = 50000;
                }


                plotYAxis(yticksevery,0.0,yRange);

                double xticksevery = 0;
                if(prevMin > 0.0 && prevMax < 1.0){
                    xticksevery = 0.1;
                } else {
                    xticksevery = 0.5;
                }

                plotXAxis(xticksevery,prevMin,prevMax);

                setAxisLabels("Score", "Count");

            } else {
                // create a new distribution for each original distribution, but now in the range of newRangeMin and newRangeMax
            }
            
         }
    }

    private void determineColor(int i) {
        Color[] colors = {
            new Color(0, 0, 0, 128),
            new Color(128, 128, 128, 128),
            new Color(192, 192, 192, 64),
            new Color(0, 0, 255, 128),
            new Color(255, 255, 0, 128),
            new Color(0, 255, 255, 128),
            new Color(255, 0, 255, 128),};

        while (i > colors.length - 1) {
            i -= colors.length - 1;
        }

        setColor(colors[i]);

    }
    
    private void plotDist(DiscreteDist dist){
        int barWidth = (int) Math.floor((double) drawWidth / dist.getNumBins());
        int iterator = 0;
        dist.resetIterator();
	System.out.println(dist.getNumBins());
        while(dist.hasNext()) {
            Bin currentBin =  dist.getNext();
            int count     = currentBin.getCount();
            double xCoord   = getXCoord(iterator);
            int iXCoord     = (int) Math.floor(xCoord);
            int height      = (int) Math.floor(scalingY * count);
            int yCoord      = (getYCoord(0) - height);

            // draw the bar
            drawRect(iXCoord, yCoord, barWidth, height, true);
            // draw the frequency
            // format the freqyency first
            
            iterator++;
        }
        dist.resetIterator();
    }


}
