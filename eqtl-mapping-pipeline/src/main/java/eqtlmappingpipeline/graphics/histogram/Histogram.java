/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.graphics.histogram;

import eqtlmappingpipeline.graphics.Graphics;
import eqtlmappingpipeline.mixupmapper.stat.dist.Bin;
import eqtlmappingpipeline.mixupmapper.stat.dist.DiscreteDist;

/**
 *
 * @author harm-jan
 */
public class Histogram extends Graphics {

    

    public Histogram(int width, int height) {
        super(width, height);
    }

    public Histogram(int width, int height, boolean outputPDF, String outputLoc){
        super(width, height, outputPDF, outputLoc);
        
    }

    public void plot(DiscreteDist dist) {
        setMargins(50);

        double maxX = dist.getMax();
        double minX = dist.getMin();
        int numBins = dist.getNumBins();
        double binInterval = dist.getBinInterval();
        int yRange = dist.getTotalNumValues();

        calcScaling(numBins, yRange);

        int xpos = 0;
        int ypos = 0;

        setStroke(3, 1);

        setColor(196, 196, 196, 64);
        plotDist(dist, xpos, ypos);
        double median = dist.getMedian();
        double medianpos = median / binInterval;
        setColor(196, 196, 196, 64);
        // median
        drawLine(getXCoord(medianpos), getYCoord((yRange / 2) - (yRange / 8)), getXCoord(medianpos), getYCoord((yRange / 2) + (yRange / 8)));

        // cumulative discordant
        setColor(196, 196, 196, 128);
        plotCumulative(dist, xpos, ypos);

        // axis
        setColor(0, 0, 0, 200);

        plotXAxis(binInterval, numBins, minX, xpos, ypos);
        plotYAxis(xpos, yRange);

    }

    protected void plotDist(DiscreteDist dist, int startX, int startY) {

        int barWidth = (int) Math.floor((double) drawWidth / dist.getNumBins());
        int iterator = 0;
        while (dist.hasNext()) {
            Bin currentBin = dist.getNext();
            int freq = currentBin.getCount();
            double xCoord = getXCoord(iterator);
            int iXCoord = (int) Math.floor(xCoord + startX);
            int height = (int) Math.floor(scalingY * freq);
            int yCoord = (getYCoord(startY) - height);

            // draw the bar
            drawRect(iXCoord, yCoord, barWidth, height, true);
            // draw the frequency
            drawText(iXCoord, yCoord - 20, String.valueOf(freq));
            iterator++;
        }
        dist.resetIterator();
    }

    protected void plotCumulative(DiscreteDist dist, int startX, int startY) {

        int iterator = 0;

        dist.calcCumulative();

        int[] x = new int[dist.getNumBins()];
        int[] y = new int[dist.getNumBins()];

        int freq = 0;

        while (dist.hasNext()) {
            Bin currentBin = (Bin) dist.getNext();

            double xCoord = getXCoord(iterator);
            int iXCoord = (int) Math.floor(xCoord + startX);

            freq = currentBin.getCumulative();
            int height = (int) Math.floor(scalingY * freq);
            int yCoord = (getYCoord(startY) - height);

            x[iterator] = iXCoord;
            y[iterator] = yCoord;

            iterator++;
        }

        int xCoord1 = 0;
        int xCoord2 = 0;
        int yCoord1 = 0;
        int yCoord2 = 0;

        for (int i = 0; i < x.length; i++) {
            if (i == 0) {
                xCoord1 = getXCoord(0);
                yCoord1 = getYCoord(0);
            } else {
                xCoord1 = x[i - 1];
                yCoord1 = y[i - 1];
            }

            xCoord2 = x[i];
            yCoord2 = y[i];

            drawLine(xCoord1, yCoord1, xCoord2, yCoord2);
        }

        dist.resetIterator();

    }

    protected void plotYAxis(int numTicks, double minY, double maxY) {
        setColor(0, 0, 0, 255);
        setStroke(1);

        int xCoord1 = marginLeft;
        int xCoord2 = marginLeft;
        int yCoord1 = marginTop;
        int yCoord2 = graphHeight - marginBottom;

        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

        java.text.DecimalFormat dformat = new java.text.DecimalFormat("0.00");

        int tickInterval = (int) Math.floor((double) (yCoord2 - yCoord1) / numTicks);
        double interval  = Math.abs(maxY - minY) / numTicks;
        for (int i = 0; i < numTicks; i++) {
            int yPos = yCoord2 - (tickInterval * i);
            int xPos = xCoord1 - 2;

            if (i != 0) {
                drawLine(xPos, yPos, xPos + 4, yPos);   
            }
            // draw some text labels
            double value = minY + (interval * i);
            int stringWidth = getStringWidth(dformat.format(value));
            drawText(xPos - stringWidth - 2, yPos, dformat.format(value));

        }

        int yPos = yCoord2 - (tickInterval * numTicks);
        int xPos = xCoord1 - 2;
        drawLine(xPos, yPos, xPos + 4, yPos);
        double value = minY + (interval * numTicks);
        int stringWidth = getStringWidth(dformat.format(value));
        drawText(xPos - stringWidth - 2, yPos, dformat.format(value));
    }

    protected void plotYAxis(double ticksEvery, double minY, double maxY) {
        setColor(0, 0, 0, 255);
        setStroke(1);
        setFont(12, "Verdana");
        int xCoord1 = marginLeft;
        int xCoord2 = marginLeft;
        int yCoord1 = marginTop;
        int yCoord2 = graphHeight - marginBottom;

        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);
        java.text.DecimalFormat dformat = new java.text.DecimalFormat("0.00");

        double yrange = Math.abs(maxY) - Math.abs(minY);
        if(minY < 0.0 && maxY > 0.0){
            yrange = Math.abs(maxY) + Math.abs(minY);
        }

        double yscaling = drawHeight / yrange;

        double prevTick = 0.0;
        while(prevTick <= maxY){

            int yPos = yCoord2 - (int) Math.floor( prevTick * yscaling );
            int xPos = xCoord1 - 2;
            drawLine(xPos, yPos, xPos + 4, yPos);

            // draw some text labels
            double value = prevTick;
            int stringWidth = getStringWidth(dformat.format(value));
            drawText(xPos - stringWidth - 12, yPos, dformat.format(value));

            prevTick += ticksEvery;
        }

        prevTick = 0.0 - ticksEvery;
        while(prevTick >= minY){

            int yPos = yCoord2 - (int) Math.floor( prevTick * yscaling );
            int xPos = xCoord1 - 2;
            drawLine(xPos, yPos, xPos + 4, yPos);

            // draw some text labels
            double value = prevTick;
            int stringWidth = getStringWidth(dformat.format(value));
            drawText(xPos - stringWidth - 12, yPos, dformat.format(value));

            prevTick -= ticksEvery;
        }

    }

    protected void plotXAxis(double ticksEvery, double minX, double maxX) {

        setColor(0, 0, 0, 255);
        setStroke(1);
        setFont(12, "Verdana");
        int xCoord1 = marginLeft;
        int xCoord2 = graphWidth - marginRight;
        int yCoord1 = graphHeight - marginBottom;
        int yCoord2 = graphHeight - marginBottom;

        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);
        java.text.DecimalFormat dformat = new java.text.DecimalFormat("0.00");

        double xrange = Math.abs(maxX) - Math.abs(minX);
        if(minX < 0.0 && maxX > 0.0){
            xrange = Math.abs(maxX) + Math.abs(minX);
        }

        double xscaling = drawWidth / xrange;

        double prevTick = 0.0;
        while(prevTick <= maxX){

            int xPos = (int) Math.floor(xCoord1 + (xscaling * (prevTick + Math.abs(minX) ) ) );
            int yPos = yCoord1 - 2;
            drawLine(xPos, yPos, xPos, yPos + 4);

            // draw some text labels
            double value = prevTick;
            int stringWidth = getStringWidth(dformat.format(value));
            drawText(xPos + 6 , yPos + stringWidth + 10 , -90, dformat.format(value));
            
            prevTick += ticksEvery;
        }

        prevTick = 0.0 - ticksEvery;
        while(prevTick >= minX){

            int xPos = (int) Math.floor(xCoord1 + (xscaling * (prevTick + Math.abs(minX) ) ) );
            int yPos = yCoord1 - 2;
            drawLine(xPos, yPos, xPos, yPos + 4);

            // draw some text labels
            double value = prevTick;
            int stringWidth = getStringWidth(dformat.format(value));
            drawText(xPos + 6, yPos + stringWidth + 10 , -90, dformat.format(value));

            prevTick -= ticksEvery;
        }

    }

    protected void plotXAxis(int numTicks, double minX, double maxX) {
        setColor(0, 0, 0, 255);
        setStroke(1);

        int xCoord1 = marginLeft;
        int xCoord2 = graphWidth - marginRight;
        int yCoord1 = graphHeight - marginBottom;
        int yCoord2 = graphHeight - marginBottom;

        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

        java.text.DecimalFormat dformat = new java.text.DecimalFormat("0.00");

        int tickInterval = (int) Math.floor((double) (xCoord2 - xCoord1) / numTicks);
        double interval  = Math.abs(maxX - minX) / numTicks;
        for (int i = 0; i < numTicks; i++) {
            if (i != 0) {
                int xPos = xCoord1 + (tickInterval * i);
                int yPos = yCoord1 - 2;
                drawLine(xPos, yPos, xPos, yPos + 4);

                // draw some text labels
                double value = minX + (interval * i);
                int stringWidth = getStringWidth(dformat.format(value));
                drawText(xPos, yPos + stringWidth + 5 , -90, dformat.format(value));
            }
        }

        int xPos = xCoord1 + (tickInterval * numTicks);
        int yPos = yCoord1 - 2;
        drawLine(xPos, yPos, xPos, yPos + 4);

        double value = minX + (interval * numTicks);
        int stringWidth = getStringWidth(dformat.format(value));
        drawText(xPos, yPos + stringWidth + 5 , -90, dformat.format(value));
    }

    protected void plotYAxis(int startX, int stopY) {

        drawLine(getXCoord(startX) - 5, getYCoord(0) + 5, getXCoord(startX) - 5, getYCoord(stopY) - 5);

        // 0
        drawText(getXCoord(startX) - 30, getYCoord(0), "0");
        drawText(getXCoord(startX) - 30, getYCoord(stopY / 2), String.valueOf((int) Math.floor(stopY / 2)));
        drawText(getXCoord(startX) - 30, getYCoord(stopY), String.valueOf(stopY));



    }

    

    protected void plotXAxis(double binInterval, int numBins, double minX, int startX, int startY) {
        // plot a Line

        drawLine(getXCoord(startX) - 5, getYCoord(startY) + 5, drawWidth + marginLeft + 5, getYCoord(startY) + 5);

        // axis value
        for (int i = 0; i < numBins; i++) {
            if (i == 0 || i == numBins - 1 || i == (int) Math.floor((double) numBins / 2)) {
                double axisValue = minX + (binInterval * i);
                double xCoord = getXCoord(i);
                int iXCoord = (int) Math.floor(xCoord + startX);
                drawText(iXCoord, getYCoord(startY) + 25, String.valueOf(axisValue));
            }
        }
    }
}
