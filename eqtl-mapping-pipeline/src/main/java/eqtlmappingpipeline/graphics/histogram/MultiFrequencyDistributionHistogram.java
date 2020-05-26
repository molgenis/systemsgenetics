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
public class MultiFrequencyDistributionHistogram extends FrequencyDistributionHistogram {
    private boolean plotCumulative;

    public MultiFrequencyDistributionHistogram(int width, int height) {
        super(width, height);
    }

    public MultiFrequencyDistributionHistogram(int width, int height, boolean outputPDF, String outputLoc){
        super(width, height, outputPDF, outputLoc);
    }

    public void plot(DiscreteDist[] a) {
        if (graphWidth > 3 * 50 && graphHeight > 3 * 50) {
            setMargins(100);

            boolean sameXRange = true;
            double prevMax = 0.0;
            double prevMin = 0.0;
            double newRangeMax = Double.MIN_VALUE;
            double newRangeMin = Double.MAX_VALUE;

            for (int i = 0; i < a.length; i++) {
                // check if all dists have the same range on the x-axis
                double minX = a[i].getMin();
                double maxX = a[i].getMax();

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
                 double yRange = Double.MIN_VALUE;
                for(int i=0; i<a.length; i++){
                    a[i].resetIterator();
                    Bin currentbin = a[i].getNext();
                    if(currentbin.getFrequency() > yRange){
                        yRange = currentbin.getFrequency();
                    }
                    while(a[i].hasNext()){
                        currentbin = a[i].getNext();
                        if(currentbin.getFrequency() > yRange){
                            yRange = currentbin.getFrequency();
                        }
                    }
                }
                

                for (int i = 0; i < a.length; i++) {
                    int numBins = a[i].getNumBins();
                    double barsize = Math.floor(drawWidth / numBins);
                    drawWidth = (int) Math.floor(barsize * numBins);

                    determineColor(i);
                    
                    if(plotCumulative){
                        calcScalingDoNotRecalculateDrawArea(numBins, 1.0);
                        plotCumulative(a[i]);
                    } else {
                        
                        calcScalingDoNotRecalculateDrawArea(numBins, yRange);
                    }

                    plotDist(a[i]);
                }

                // draw the axis and axis labels

                plotYAxis(0.01,0.0,yRange);

                double xticksevery = 0;
                if(prevMin > 0.0 && prevMax < 1.0){
                    xticksevery = 0.1;
                } else {
                    xticksevery = 0.5;
                }

                plotXAxis(xticksevery,prevMin,prevMax);
                
                setAxisLabels("Score", "Frequency");

            } else {
                // create a new distribution for each original distribution, but now in the range of newRangeMin and newRangeMax
            }

        }


        // create a histogram on the basis of the relative frequency of each bin


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
}
