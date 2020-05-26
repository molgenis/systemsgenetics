/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.graphics.map;

import java.awt.Color;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.ArrayList;

/**
 *
 * @author harm-jan
 */
public class MixupMapperHeatMap extends Heatmap {

    private double[] rowEigenVector, colEigenVector;
    private int[] bestGenotypeMatches;
    private boolean familyDataLoaded;
    private HashMap<String, String> expressionToGenotypeCouplingNames;
    private String xLab;
    private String yLab;

    public MixupMapperHeatMap(int width, int height) {
        super(width, height);
        rowEigenVector = null;
        colEigenVector = null;
        bestGenotypeMatches = null;
    }

    public MixupMapperHeatMap(int width, int height, boolean outputPDF, String outputLoc) {
        super(width, height, outputPDF, outputLoc);
        rowEigenVector = null;
        colEigenVector = null;
        bestGenotypeMatches = null;
    }

    @Override
    public void plot(double[][] matrix) {
        // we have set the desired width, and desired height, although this might not fit...
        graphWidth = desiredWidth;
        graphHeight = desiredHeight;
        calcDrawArea();

        int numRows = matrix.length;
        int numCols = matrix[numRows - 1].length;

        int boxWidth = (int) Math.floor((double) drawWidth / numCols);
        int boxHeight = (int) Math.floor((double) drawHeight / numRows);

        if (boxWidth == 0 || boxHeight == 0) {
            boxWidth = 1;
            boxHeight = 1;
        }

        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
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

        // determine range for gradient
        double range;
        if (min < 0 && max > 0) {
            range = max + Math.abs(min);
        } else {
            range = Math.abs(max - min);
        }

        if (boxWidth == 0 || boxHeight == 0) {
            System.out.println("Cannot plot a heatmap on " + drawWidth + "x" + drawHeight + " for a matrix of " + numRows + "x" + numCols);
        } else {
            // i want symmetrical boxes, so take the lowest value...
            if (boxWidth > boxHeight) {
                boxWidth = boxHeight;
            } else {
                boxHeight = boxWidth;
            }

            setFont(boxWidth, "Verdana");

            int maxRowNameLength = 0;
            int maxColNameLength = 0;

            if (rowNames != null) {
                // find the maximal length of the string within the rownames
                maxRowNameLength = getMaxStringLength(rowNames);
            }

            if (colNames != null) {
                // find the maximal length of the string within the rownames
                maxColNameLength = getMaxStringLength(colNames);
            }

            // correlation box

            int correlationBoxWidth = 100;
            if (rowEigenVector == null && colEigenVector == null) {
                correlationBoxWidth = 0;
            }

            int calculatedWidth;
            int calculatedHeight;

            double minRowEigenVectorValue = Double.MAX_VALUE;
            double maxRowEigenVectorValue = Double.MIN_VALUE;

            double minColEigenVectorValue = Double.MAX_VALUE;
            double maxColEigenVectorValue = Double.MIN_VALUE;

            if (rowEigenVector != null) {
                calculatedWidth = marginLeft + (boxWidth * numCols) + (2 * maxColNameLength) + boxWidth + correlationBoxWidth + boxWidth + marginRight;
                for (Double d : rowEigenVector) {
                    if (d > maxRowEigenVectorValue) {
                        maxRowEigenVectorValue = d;
                    }
                    if (d < minRowEigenVectorValue) {
                        minRowEigenVectorValue = d;
                    }
                }
            } else {
                calculatedWidth = marginLeft + (boxWidth * numCols) + (2 * maxColNameLength) + boxWidth + marginRight;
            }

            if (colEigenVector != null) {
                calculatedHeight = marginTop + (boxWidth * numRows) + maxColNameLength + boxWidth + correlationBoxWidth + boxWidth + marginBottom;

                for (Double d : colEigenVector) {
                    if (d > maxColEigenVectorValue) {
                        maxColEigenVectorValue = d;
                    }
                    if (d < minColEigenVectorValue) {
                        minColEigenVectorValue = d;
                    }
                }
            } else {
                calculatedHeight = marginTop + (boxWidth * numRows) + maxColNameLength + boxWidth + marginBottom;
            }

            // initialize the graphing objects
            init(calculatedWidth, calculatedHeight);
            setMargins(50);

            // plot the title
            setFont(15, "Verdana");

            setColor(128, 128, 128, 255);
            drawText(graphWidth - marginRight + (int) Math.floor(marginRight / 2) + 10, (int) Math.floor(graphWidth / 2), -90, yLab); //+ (int) Math.floor(getStringWidth("Gene expression sample")/2)
            drawText((int) Math.floor(graphWidth / 2) - (int) Math.floor(getStringWidth(xLab)), graphHeight - (int) Math.floor(marginBottom / 2), xLab);

            // draw the grid
            drawGrid(numRows, numCols, boxWidth);

            if (!colSorted) {
                colOrder = new int[numCols];
                for (int i = 0; i < numCols; i++) {
                    colOrder[i] = i;
                }
            }

            if (!rowSorted) {
                rowOrder = new int[numRows];
                for (int rownum = 0; rownum < numRows; rownum++) {
                    rowOrder[rownum] = rownum;
                }
            }

            for (int row = 0; row < numRows; row++) {

                int rowPosition = rowOrder[row];

//                String coupledGenotype = expressionToGenotypeCouplingNames.get(rowNames[row]);

                if (rowEigenVector != null) {
                    // draw the sample correlation boxes
                    r = 0;
                    g = 0;
                    b = 0;
                    a = 0;
                    setColor(128, 128, 128, 128);
                    // determineColor(eigenExpVector[row], minEigenExpVectorValue, maxEigenExpVectorValue - minEigenExpVectorValue);

                    int xCoord = marginLeft + (boxWidth * numCols);
                    int yCoord = (rowPosition * boxWidth) + marginTop;

                    int wdth = (int) Math.floor(correlationBoxWidth * rowEigenVector[row]);

                    g2d.fillRect(xCoord, yCoord, wdth, boxHeight);

                    Stroke oldStroke = g2d.getStroke();
                    Color oldColor = g2d.getColor();
                    g2d.setColor(new Color(200, 200, 200));
                    setStroke(1);
                    g2d.draw(new Rectangle2D.Double(xCoord, yCoord, wdth, boxHeight));
                    g2d.setColor(oldColor);
                    g2d.setStroke(oldStroke);

                    setColor(0, 0, 0, 255);
                    setFont(boxWidth, "Verdana");
                    java.text.DecimalFormat df = new java.text.DecimalFormat("0.00");
                    drawText(xCoord + boxWidth, yCoord + boxWidth, df.format(rowEigenVector[row]));
                }

                for (int col = 0; col < numCols; col++) {
                    int colPosition = colOrder[col];
                    double val = matrix[row][col];



                    invertGradient = false;
                    r = 0;
                    g = 0;
                    b = 0;
                    a = 0;

                    if (bestGenotypeMatches != null) {
                        if (bestGenotypeMatches[col] == row) {
                            invertGradient = true;
                            r = 0;
                            g = 0;
                            b = 0;
                            a = 0;
                        }
                    }

                    determineColor(val, min, range);
                    int y = (rowPosition * boxWidth) + marginTop;
                    int x = (colPosition * boxWidth) + marginLeft;
                    drawRect(x, y, boxWidth, boxHeight, true);
                    invertGradient = false;
                }
            }

            for (int col = 0; col < numCols; col++) {
                int colPosition = colOrder[col];
                if (colEigenVector != null) {
                    r = 0;
                    g = 0;
                    b = 0;
                    a = 0;
                    setColor(128, 128, 128, 128);
                    // determineColor(eigenGenVector[col], minEigenGenVectorValue, maxEigenGenVectorValue - minEigenGenVectorValue);

                    int xCoord = marginLeft + (colPosition * boxWidth);
                    int yCoord = marginTop + (boxWidth * numRows);

                    int wdth = (int) Math.floor(correlationBoxWidth * colEigenVector[col]);

                    g2d.fillRect(xCoord, yCoord, boxWidth, wdth);

                    Stroke oldStroke = g2d.getStroke();
                    Color oldColor = g2d.getColor();
                    g2d.setColor(new Color(200, 200, 200));
                    setStroke(1);
                    g2d.draw(new Rectangle2D.Double(xCoord, yCoord, boxWidth, wdth));
                    g2d.setColor(oldColor);
                    g2d.setStroke(oldStroke);

                    setColor(0, 0, 0, 255);
                    setFont(boxWidth, "Verdana");
                    java.text.DecimalFormat df = new java.text.DecimalFormat("0.00");
                    drawText(xCoord + boxWidth, yCoord + 3 * boxWidth, -90, df.format(colEigenVector[col]));
                }
            }


            r = 0;
            g = 0;
            b = 0;
            a = 255;
            setColor(0, 0, 0, 255);
            drawLabels(numRows, numCols, boxWidth, correlationBoxWidth);

            invertGradient = false;
            r = 0;
            g = 0;
            b = 0;
            a = 0;
            drawLegend(numRows, numCols, boxWidth, min, range);
        }


    }
    
    @Override
    public void setAxisLabels(String xlab, String ylab) {
        this.xLab = xlab;
        this.yLab = ylab;
    }

    public void setBestGenotypeMatches(int[] matches) {
        bestGenotypeMatches = matches;
    }

    public void setPC1EigenVector(double[] rowEig, double[] colEig) {
        rowEigenVector = rowEig;
        colEigenVector = colEig;
    }

    public void setFamilyData(HashMap<String, String> sampleToFamilyId, HashMap<String, ArrayList<String>> familyToSampleIds, HashMap<String, String> etgNames) {
        familyDataLoaded = true;
        this.expressionToGenotypeCouplingNames = etgNames;
    }
}
