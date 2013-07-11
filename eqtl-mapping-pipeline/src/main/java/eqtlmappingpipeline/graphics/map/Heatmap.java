/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.graphics.map;

import eqtlmappingpipeline.graphics.Graphics;

/**
 *
 * @author harm-jan
 */
public class Heatmap extends Graphics {

    protected int r, g, b, a;
    protected boolean useAlphaGradient, colSorted, rowSorted;
    protected String[] rowNames, colNames;
    protected int[] rowOrder, colOrder;
    protected int desiredWidth, desiredHeight;
    protected boolean invertGradient;

    public Heatmap(int width, int height) {

        super();

        desiredWidth    = width;
        desiredHeight   = height;
        r = 0;
        g = 0;
        b = 255;
        a = 0;
        useAlphaGradient = true;
        invertGradient = false;
        colSorted = false;
        rowSorted = false;
        rowNames = null;
        colNames = null;
        rowOrder = null;
        colOrder = null;

        setMargins(50);
    }

     public Heatmap(int width, int height, boolean outputPDF, String outputLoc){
        super(width, height, outputPDF, outputLoc);
        desiredWidth    = width;
        desiredHeight   = height;
        r = 0;
        g = 0;
        b = 255;
        a = 0;
        useAlphaGradient = true;
        invertGradient = false;
        colSorted = false;
        rowSorted = false;
        rowNames = null;
        colNames = null;
        rowOrder = null;
        colOrder = null;

        setMargins(50);
        
    }

    public void plot(int[][] matrix) {}

    public void plot(double[][] matrix) {
        
        // we have set the desired width, and desired height, although this might not fit...
        graphWidth  = desiredWidth;
        graphHeight = desiredHeight;
        calcDrawArea();

        int numRows = matrix.length;
        int numCols = matrix[numRows - 1].length;

        int boxWidth  = (int) Math.floor((double) drawWidth  / numCols);
        int boxHeight = (int) Math.floor((double) drawHeight / numRows);

        double min  = Double.MAX_VALUE;
        double max  = Double.MIN_VALUE;
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                double val = matrix[row][col];
                if (val > max) {
                    max = val;
                }
                if (val < min) {
                    min = val;
                }
            }
        }

        if(min < 0.0){
            min = Math.abs(min);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {
                    double val = matrix[row][col];
                    matrix[row][col] = val + min;
                }
            }
            max = max + min;
            min = 0.0;
        }

        double range = max;

        if (boxWidth == 0 || boxHeight == 0) {
            System.out.println("Cannot plot a heatmap on " + drawWidth + "x" + drawHeight + " for a matrix of " + numRows + "x" + numCols);
        } else {
            // i want symmetrical boxes, so take the lowest value...
            if (boxWidth > boxHeight) {
                boxWidth = boxHeight;
            } else {
                boxHeight = boxWidth;
            }

            setFont(boxWidth, "Monospaced");
            
            int maxRowNameLength = 0;
            int maxColNameLength = 0;
            
            if(rowNames != null){
                // find the maximal length of the string within the rownames
                maxRowNameLength = getMaxStringLength(rowNames);
            }

            if(colNames != null){
                // find the maximal length of the string within the rownames
                maxColNameLength = getMaxStringLength(colNames);
            }

            int calculatedWidth  = marginLeft + (boxWidth * numCols) + maxColNameLength + boxWidth + marginRight;
            int calculatedHeight = marginTop +  (boxWidth * numRows) + maxRowNameLength + boxWidth + marginBottom;

            // initialize the graphing objects
            init(calculatedWidth, calculatedHeight);
            setMargins(50);

            // plot the title
            setFont(15, "Georgia");

            // draw the grid
            drawGrid(numRows, numCols, boxWidth);

            // draw colored boxes in a certain order 
            if(colSorted && rowSorted){
                for(int row=0; row < numRows; row++){
                    int rowPosition = rowOrder[row];
                    for(int col=0; col < numCols; col++){
                        int colPosition = colOrder[col];
                        double val = matrix[row][col];
                        determineColor(val, min, range);
                        int y = (rowPosition * boxWidth) + marginTop;
                        int x = (colPosition * boxWidth) + marginLeft;
                        drawRect(x, y, boxWidth, boxHeight, true);
                    }
                }
            } else if (colSorted && !rowSorted){
                for(int row=0; row < numRows; row++){
                    int rowPosition = row;
                    for(int col=0; col < numCols; col++){
                        int colPosition = colOrder[col];
                        double val = matrix[row][col];
                        determineColor(val, min, range);
                        int y = (rowPosition * boxWidth) + marginTop;
                        int x = (colPosition * boxWidth) + marginLeft;
                        drawRect(x, y, boxWidth, boxHeight, true);
                    }
                }
            } else if (rowSorted && !colSorted){
                for(int row=0; row < numRows; row++){
                    int rowPosition = rowOrder[row];
                    for(int col=0; col < numCols; col++){
                        int colPosition = col;
                        double val = matrix[row][col];
                        determineColor(val, min, range);
                        int y = (rowPosition * boxWidth) + marginTop;
                        int x = (colPosition * boxWidth) + marginLeft;
                        drawRect(x, y, boxWidth, boxHeight, true);
                    }
                }
            } else {
                for (int row = 0; row < numRows; row++) {
                    for (int col = 0; col < numCols; col++) {
                        double val = matrix[row][col];
                        determineColor(val, min, range);
                        int y = (row * boxWidth) + marginTop;
                        int x = (col * boxWidth) + marginLeft;
                        drawRect(x, y, boxWidth, boxHeight, true);
                    }
                }
            }

            drawLabels(numRows, numCols, boxWidth, 0);
        }
    }

    protected void drawLabels(int numRows, int numCols, int boxWidth, int margin){
        setColor(128,128,128,255);
        setFont(boxWidth, "Verdana");
        // draw labels
        if(rowNames != null){
            int x = (numCols * boxWidth) + marginLeft + boxWidth + margin;
            if(rowSorted){
                for(int row=0; row<numRows; row++){
                    int yPos = rowOrder[row];
                    int y = (yPos * boxWidth) + marginTop + boxWidth - 2;
                    drawText(x, y,rowNames[row]);
                }
            } else {
                 for(int row=0; row<numRows; row++){
                    int y = (row * boxWidth) + marginTop + boxWidth - 2;
                    drawText(x, y,rowNames[row]);
                }
            }
        }

        if(colNames != null){
            if(colSorted){
                for (int col = 0; col < numCols; col++) {
                    int xPos = colOrder[col];
                    int x = (xPos * boxWidth) + marginTop + boxWidth - 1;
                    int stringWidth = getStringWidth(colNames[col]);
                    int y = (numRows * boxWidth) + boxWidth + stringWidth + marginTop + margin;
                    drawText(x, y,-90, colNames[col]);
                }
            } else {
                for (int col = 0; col < numCols; col++) {
                    int x = (col * boxWidth) + marginTop + boxWidth - 1;
                    int stringWidth = getStringWidth(colNames[col]);
                    int y = (numRows * boxWidth) + boxWidth + stringWidth + marginTop + margin;
                    drawText(x, y,-90, colNames[col]);
                }
            }
        }
    }

    public void drawLegend(int numRows, int numCols, int boxWidth, double min, double range){

        double increments = range / 20;

        int tickheight = 3;

        double val = min - increments;
        for(int i=0; i<20; i++){
            val += increments;
            determineColor(val, min, range);
            int x = (numCols * boxWidth) + (2*boxWidth) + marginLeft;
            int y = marginTop + (numRows*boxWidth) + (i*tickheight) + boxWidth;
            drawRect(x, y, boxWidth, tickheight, true);
        }

        java.text.DecimalFormat df = new java.text.DecimalFormat("#.##");
        setFont(boxWidth, "Verdana");
        setColor(128,128,128,255);

        int y = marginTop + (numRows*boxWidth) + (2*boxWidth);
        int x = (numCols * boxWidth) + boxWidth + marginLeft + (2*boxWidth);
        
        drawText(x,y, df.format(min));

        y = marginTop + (numRows*boxWidth) + (20*tickheight) + (2*boxWidth);
        drawText(x,y, df.format(min+range));
    }

    public void setLabels(String[] rowLabels, String[] colLabels){
        rowNames = rowLabels;
        colNames = colLabels;
    }

    public void sortByLabels(boolean sortRows, boolean sortCols){
        if(sortRows && rowNames != null){
            String[] sortedRows = java.util.Arrays.copyOf(rowNames, rowNames.length);

            java.util.Arrays.sort(sortedRows, String.CASE_INSENSITIVE_ORDER);

            rowOrder = new int[rowNames.length];
            for(int i=0; i<rowNames.length; i++){
                rowOrder[i] = getIndexOfValue(sortedRows, rowNames[i]);
            }

            rowSorted = true;

        }

        if(sortCols && colNames != null){
            String[] sortedCols = java.util.Arrays.copyOf(colNames, colNames.length);
            java.util.Arrays.sort(sortedCols, String.CASE_INSENSITIVE_ORDER);

            colOrder = new int[colNames.length];
            for(int i=0; i<colNames.length; i++){
                colOrder[i] = getIndexOfValue(sortedCols, colNames[i]);
            }

            colSorted = true;
        }
    }

    protected int getIndexOfValue(String[] data, String value){
        for(int i=0; i<data.length; i++){
            if(data[i].equals(value)){
                return i;
            }
        }
        return -1;
    }

    protected void drawGrid(int numRows, int numCols, int boxWidth){
        setColor(196,196,196,255);
        setStroke(1);
        // draw a grid...
        // horizontal lines
        for(int row=0; row < numRows; row = row + 5){
            if(row % 25 == 0){
                setColor(128,128,128,255);
                setStroke(2);
            }

            int y1 = (row * boxWidth) + marginTop;
            int x1 = marginLeft;
            int y2 = (row * boxWidth) + marginTop;
            int x2 = graphWidth - marginRight + boxWidth;

            drawLine(x1, y1, x2, y2);
            setColor(196,196,196,255);
            setStroke(1);
        }

        setColor(196,196,196,255);
        setStroke(1);
        // vertical lines
        for(int col=0; col < numCols; col = col + 5){
            if(col % 25 == 0){
                setColor(128,128,128,255);
                setStroke(2);
            }

            int y1 = marginTop;
            int x1 = (col * boxWidth) + marginLeft;
            int y2 = graphHeight - marginBottom + boxWidth;
            int x2 = (col * boxWidth) + marginLeft;

            drawLine(x1, y1, x2, y2);
            setStroke(1);
            setColor(196,196,196,255);
        }
    }

    public void useAlphaGradient() {
        useAlphaGradient = useAlphaGradient & false;
    }

    public void setBaseColor(int r, int g, int b, int a) {

    }

    public void setMaxColor(int r, int g, int b, int a) {
        
    }

    protected void determineColor(double value, double min, double range) {

        if (useAlphaGradient) {
            // determine alpha gradient
            double perc = 0;
            if( min < 0){
                perc = (value + Math.abs(min)) / range;
            } else {
                perc = Math.abs(value - min) / range;
            }
            
            int alpha = (int) Math.floor(perc * 255);
            if(invertGradient){
                alpha = 255 - alpha;
            }
            if(alpha > 255){
                alpha = 255;
		System.out.println("Warming: value for alpha higher than max!: "+ value + "\t" + range);
            }
	    if(alpha < 0){
		alpha = 0;
		System.out.println("Warming: value for alpha lower than min!: "+ value + "\t" + range);
	    }
            setColor(r, g, b, alpha);
        } else {
            // determine color gradient from baseColor
        }

    }
}
