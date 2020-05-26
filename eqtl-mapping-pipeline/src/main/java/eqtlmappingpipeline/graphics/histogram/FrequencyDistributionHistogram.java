/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.graphics.histogram;

import eqtlmappingpipeline.mixupmapper.stat.dist.Bin;
import eqtlmappingpipeline.mixupmapper.stat.dist.DiscreteDist;

/**
 *
 * @author harm-jan
 */
public class FrequencyDistributionHistogram extends Histogram {

    public FrequencyDistributionHistogram(int width, int height){
        super(width, height);
    }

    public FrequencyDistributionHistogram(int width, int height, boolean outputPDF, String outputLoc){
        super(width, height, outputPDF, outputLoc);
    }

    @Override
    public void plot(DiscreteDist a){
         if (graphWidth > 3 * 50 && graphHeight > 3 * 50) {
            setMargins(50);
            int yRange          = 1;

            int numBins = a.getNumBins();
            double barsize = Math.floor(drawWidth / numBins);
            drawWidth = (int) Math.floor(barsize * numBins);

            calcScalingDoNotRecalculateDrawArea(numBins, yRange);
            setColor(0,0,0,128);
            plotDist(a);
            plotCumulative(a);
            

            // draw the axis and axis labels
            plotYAxis(4,0.0,1.0);
            plotXAxis(4,a.getMin(),a.getMin());
         }
    }

    protected void plotDist(DiscreteDist dist){
        int barWidth = (int) Math.ceil((double) drawWidth / dist.getNumBins());
        int iterator = 0;
        dist.resetIterator();
        setStroke(0);
        while(dist.hasNext()) {
            Bin currentBin =  dist.getNext();
            double freq     = currentBin.getFrequency();
            double xCoord   = getXCoord(iterator);
            int iXCoord     = (int) Math.floor(xCoord);
            int height      = (int) Math.floor(scalingY * freq);
            int yCoord      = (getYCoord(0) - height);

            // draw the bar
            drawRect(iXCoord, yCoord, barWidth, height, true);
            // draw the frequency
            // format the freqyency first
                      
            iterator++;
        }
        dist.resetIterator();
    }

    protected void plotCumulative(DiscreteDist dist){
        int iterator = 0;
        dist.resetIterator();

        int[] x = new int[dist.getNumBins()];
        int[] y = new int[dist.getNumBins()];

        double freq = 0;

        while( dist.hasNext() ) {
            Bin currentBin  = dist.getNext();

            double xCoord   = getXCoord(iterator);
            int iXCoord     = (int) Math.floor(xCoord);

            freq            = currentBin.getCumulativeFrequency();
            int height      = (int) Math.floor(scalingY * freq);
            int yCoord      = (getYCoord(0) - height);
            if(iterator == 0){
                yCoord      = getYCoord(0);
            }

            x[iterator] = iXCoord;
            y[iterator] = yCoord;

            iterator++;
        }

        int xCoord1 = 0;
        int xCoord2 = 0;
        int yCoord1 = 0;
        int yCoord2 = 0;

        for(int i=0; i<x.length; i++){
            if(i == 0){
                xCoord1 = getXCoord(0);
                yCoord1 = getYCoord(0);
            } else {
                xCoord1 = x[i-1];
                yCoord1 = y[i-1];
            }

            xCoord2 = x[i];
            yCoord2 = y[i];

            drawLine(xCoord1, yCoord1, xCoord2, yCoord2);
        }

        dist.resetIterator();
    }
}
