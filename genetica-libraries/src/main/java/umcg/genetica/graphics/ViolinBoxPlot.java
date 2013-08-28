/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.graphics;

import com.lowagie.text.DocumentException;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.GeneralPath;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Locale;
import javax.imageio.ImageIO;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

/**
 *
 * @author harmjan
 */
public class ViolinBoxPlot {

    public enum Output {

        PDF, PNG
    };

    private double getWidth(String text, java.awt.Font font) {
        java.awt.Graphics2D g2d = (new java.awt.image.BufferedImage(1, 1, java.awt.image.BufferedImage.TYPE_INT_ARGB)).createGraphics();
        java.awt.font.TextLayout tL = new java.awt.font.TextLayout(text, font, g2d.getFontRenderContext());
        return (double) tL.getBounds().getWidth();
    }

    // draw multiple violinplots (for example for multiple datasets) format: vals[dataset][category][value]
    // xlabels format: xlabels2[dataset][category]
    public void draw(double[][][] vals, String[] datasetNames, String[][] xLabels2, Output output, String outputFileName) throws IOException {

        Locale defaultLocale = Locale.getDefault();
        Locale.setDefault(Locale.US);
        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.lowagie.text.Document document = null;
        com.lowagie.text.pdf.PdfWriter writer = null;
        com.lowagie.text.pdf.PdfContentByte cb = null;

        BufferedImage bi = null;

        // draw individual box/violinplots
        java.awt.AlphaComposite alphaComposite25 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.25f);
        java.awt.AlphaComposite alphaComposite50 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f);
        java.awt.AlphaComposite alphaComposite100 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC, 1.00f);

        float fontSize = 12f;
        java.awt.Font font = new java.awt.Font("Gill Sans MT", java.awt.Font.PLAIN, (int) fontSize);
        java.awt.Font fontBold = new java.awt.Font("Gill Sans MT", java.awt.Font.BOLD, (int) fontSize);
        java.awt.Font fontSmall = new java.awt.Font("Gill Sans MT", java.awt.Font.PLAIN, 8);
        java.awt.Font fontBoldSmall = new java.awt.Font("Gill Sans MT", java.awt.Font.BOLD, 8);

        int[] datasetWitdths = new int[datasetNames.length];
        int individualPlotWidth = 70;
        int individualPlotMarginLeft = 10;
        int individualPlotMarginRight = 10;
        int betweenDatasetMargin = 20;
        int marginLeft = 100;
        int marginRight = 100;
        int marginTop = 100;
        int marginBottom = 100;
        int innerHeight = 500;

        int totalDatasetWidth = 0;
        for (int x = 0; x < datasetNames.length; x++) {
            int nrCategoriesInDataset = vals[x].length;
            datasetWitdths[x] = nrCategoriesInDataset * (individualPlotWidth + individualPlotMarginLeft + individualPlotMarginRight);
            totalDatasetWidth += datasetWitdths[x];
        }

        totalDatasetWidth += (datasetNames.length - 1) * betweenDatasetMargin;

        int docWidth = marginLeft + marginRight + totalDatasetWidth;
        int docHeight = marginTop + marginBottom + innerHeight;

        if (output == Output.PDF) {
            com.lowagie.text.Rectangle rectangle = new com.lowagie.text.Rectangle(docWidth, docHeight);
            document = new com.lowagie.text.Document(rectangle);
            try {
                writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFileName));
            } catch (DocumentException e) {
                throw new IOException(e.fillInStackTrace());
            }
            document.open();
            cb = writer.getDirectContent();
            cb.saveState();
            //com.lowagie.text.pdf.DefaultFontMapper fontMap = new com.lowagie.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(docWidth, docHeight);
        } else {
            bi = new BufferedImage(docWidth, docHeight, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }


        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, docWidth, docHeight);



        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;

        cern.jet.random.engine.RandomEngine randomEngine = new cern.jet.random.engine.DRand();

        WilcoxonMannWhitney wmw = new WilcoxonMannWhitney();


        double[][][] aucs = new double[vals.length][vals[0].length][vals[0].length];
        double[][][] pvals = new double[vals.length][vals[0].length][vals[0].length];

        for (int dataset = 0; dataset < datasetNames.length; dataset++) {
            double[][] valsForDs = vals[dataset];
            // test between categories
            for (int category1 = 0; category1 < valsForDs.length; category1++) {
                double[] vals1 = valsForDs[category1];
                double min1 = JSci.maths.ArrayMath.min(vals1);
                double max1 = JSci.maths.ArrayMath.max(vals1);
                if (min1 < minValue) {
                    minValue = min1;
                }
                if (max1 > maxValue) {
                    maxValue = max1;
                }

                for (int category2 = category1 + 1; category2 < valsForDs.length; category2++) {

                    double[] vals2 = valsForDs[category2];
                    double pValueWilcoxon = wmw.returnWilcoxonMannWhitneyPValue(vals2, vals1);
                    double auc = wmw.getAUC();
                    double min2 = JSci.maths.ArrayMath.min(vals2);
                    double max2 = JSci.maths.ArrayMath.max(vals2);

                    aucs[dataset][category1][category2] = auc;
                    pvals[dataset][category1][category2] = pValueWilcoxon;


                    if (max2 > maxValue) {
                        maxValue = max2;
                    }
                    if (min2 < minValue) {
                        minValue = min2;
                    }
                }
            }

        }

        ArrayList<Pair<Double, Integer>> sortedPValuesPerDataset = new ArrayList<Pair<Double, Integer>>();
        ArrayList<Pair<Double, Triple<Integer, Integer, Integer>>> sortedPValuesPerCategory = new ArrayList<Pair<Double, Triple<Integer, Integer, Integer>>>();
        for (int dataset = 0; dataset < vals.length; dataset++) {

            // sort categories on the basis of their WMW results
            double[][] aucsfordataset = aucs[dataset];
            double[][] pvaluesfordataset = pvals[dataset];

            double minPvalueForDs = Double.MAX_VALUE;

            for (int category1 = 0; category1 < pvaluesfordataset.length; category1++) {
                for (int category2 = category1 + 1; category2 < pvaluesfordataset.length; category2++) {
                    double valueToSort = pvaluesfordataset[category1][category2];
//                    if(sortByAUC){
//                        valueToSort = aucsfordataset[category1][category2];
//                    }
                    sortedPValuesPerCategory.add(new Pair<Double, Triple<Integer, Integer, Integer>>(valueToSort, new Triple<Integer, Integer, Integer>(dataset, category1, category2), Pair.SORTBY.LEFT));
                    if (valueToSort < minPvalueForDs) {
                        minPvalueForDs = valueToSort;
                    }
                }
            }

            sortedPValuesPerDataset.add(new Pair<Double, Integer>(minPvalueForDs, dataset, Pair.SORTBY.LEFT));

        }

        Collections.sort(sortedPValuesPerDataset, Collections.reverseOrder());
        Collections.sort(sortedPValuesPerCategory, Collections.reverseOrder());



        int datasetCounter = 0;
        int deltaX = 0;
        for (Pair<Double, Integer> pvalueDatasetPair : sortedPValuesPerDataset) {

            Integer datasetNumber = pvalueDatasetPair.getRight();
            int datasetwidth = datasetWitdths[datasetNumber];

            // draw the gradient
            g2d.setComposite(alphaComposite100);
            int x = marginLeft + deltaX;
            int height = innerHeight;
            java.awt.GradientPaint gradient = new java.awt.GradientPaint(0, marginTop, new Color(230, 230, 230), 0, height, new Color(250, 250, 250));
            g2d.setPaint(gradient);
            g2d.fillRect(x - individualPlotMarginLeft / 2, marginTop - 90, datasetwidth, height + 90);

            // now get the sorted results per category within this dataset
            ArrayList<Triple<Integer, Integer, Integer>> sortedPValuesForDataset = new ArrayList<Triple<Integer, Integer, Integer>>();
            for (Pair<Double, Triple<Integer, Integer, Integer>> pvalueTriplePair : sortedPValuesPerCategory) {
                System.out.println(pvalueTriplePair.toString());
                if (pvalueTriplePair.getRight().getLeft().equals(datasetNumber)) {
                    // this result belongs to this dataset
                    sortedPValuesForDataset.add(pvalueTriplePair.getRight());
                }
            }

            // reorder categories..
            HashSet<Integer> visitedCategories = new HashSet<Integer>();
            ArrayList<Integer> categoryOrder = new ArrayList<Integer>();
            for (Triple<Integer, Integer, Integer> datasetCategoryCombo : sortedPValuesForDataset) {
                Integer category1 = datasetCategoryCombo.getMiddle();
                Integer category2 = datasetCategoryCombo.getRight();

                if (!visitedCategories.contains(category1)) {
                    categoryOrder.add(category1);
                    visitedCategories.add(category1);
                }

                if (!visitedCategories.contains(category2)) {
                    categoryOrder.add(category2);
                    visitedCategories.add(category2);
                }

            }

            String tissueName = datasetNames[datasetNumber];
            int y = marginTop;
            g2d.setColor(new Color(0, 0, 0));
            g2d.setFont(fontBold);

            //            Vector vecSentences = multi.getSentencesText(tissueName, tissueWidth, g2d);
//            for (int t = 0; t < vecSentences.size(); t++) {
//                String sentence = ((String) vecSentences.get(t)).trim();
//                g2d.drawString(sentence, x + tissueWidth / 2 - (int) Math.round(getWidth(sentence, g2d.getFont()) / 2), y + t * 13 + 13 - 90);
//            }
            g2d.setFont(font);

            // now we can plot the combinations between the categories (and their pvalues and aucs, yay!)
            int[] categoryIndex = new int[vals[datasetNumber].length];
            int categoryCounter = 0;
            int plotStart = x + individualPlotMarginLeft;
            for (Integer category : categoryOrder) {
                double[] vals1 = vals[datasetNumber][category];
                categoryIndex[category] = categoryCounter;

                // plot the individual box and violin plots
                g2d.setComposite(alphaComposite25);
                g2d.setColor(new Color(223, 36, 20));

                int xposViolin = plotStart + (individualPlotMarginLeft + individualPlotWidth + individualPlotMarginRight) * categoryCounter;
                int xposBoxPlot = xposViolin; // + (individualPlotWidth / 2 - 3) / 2;
                drawViolinPlot(g2d, xposViolin, y, individualPlotWidth, height, vals1, minValue, maxValue);
                g2d.setComposite(alphaComposite100);
                g2d.setColor(new Color(0, 0, 0));
                drawBoxPlot(g2d, xposBoxPlot, y, individualPlotWidth, height, vals1, minValue, maxValue, false);

                // Draw bottom discrimator between candidate genes and other genes:
                double minVal1 = JSci.maths.ArrayMath.min(vals1);

                int posY1 = y + height - (int) Math.round((double) height * (minVal1 - minValue) / (maxValue - minValue));


                int linePos = xposViolin + (individualPlotWidth / 2);
                g2d.setComposite(alphaComposite25);
                g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, new float[]{2.0f}, 0.0f));
                g2d.setColor(new Color(223, 36, 20));
                g2d.drawLine(linePos, posY1 + 5, linePos, docHeight - marginBottom + 15);

                g2d.setComposite(alphaComposite100);
                g2d.setFont(fontBoldSmall);
                g2d.setColor(new Color(223, 36, 20));

                g2d.drawString(xLabels2[datasetNumber][category],
                        linePos - (int) Math.round(getWidth(xLabels2[datasetNumber][category], g2d.getFont()) / 2),
                        docHeight - marginBottom + 25);

                categoryCounter++;
            }

            // plot the wilcoxon result lines..
            int comboCounter = 0;

            int previousHeight = 0;

            for (Triple<Integer, Integer, Integer> datasetCategoryCombo : sortedPValuesForDataset) {
                Integer category1 = datasetCategoryCombo.getMiddle();
                Integer category2 = datasetCategoryCombo.getRight();

                int category1Index1 = categoryIndex[category1];
                int category1Index2 = categoryIndex[category2];

                int xpos1 = plotStart + (individualPlotMarginLeft + individualPlotWidth + individualPlotMarginRight) * category1Index1 + (individualPlotWidth / 2);
                int xpos2 = plotStart + (individualPlotMarginLeft + individualPlotWidth + individualPlotMarginRight) * category1Index2 + (individualPlotWidth / 2);

                double maxVal1 = JSci.maths.ArrayMath.max(vals[datasetNumber][category1]);
                double maxVal2 = JSci.maths.ArrayMath.max(vals[datasetNumber][category2]);
                int posY1 = y + height - (int) Math.round((double) height * (maxVal1 - minValue) / (maxValue - minValue));
                int posY2 = y + height - (int) Math.round((double) height * (maxVal2 - minValue) / (maxValue - minValue));
                
                int horizontalLineYPos = 0;
                if(posY1 < posY2){
                    horizontalLineYPos = posY1 - 25;
                } else {
                    horizontalLineYPos = posY2 - 25;
                }
                
                g2d.setComposite(alphaComposite100);
                g2d.setColor(new Color(100, 100, 100));

                int comboMargin = comboCounter * 20;

                g2d.drawLine(xpos1, posY1, xpos1, horizontalLineYPos); // vertical line 1
                g2d.drawLine(xpos1, horizontalLineYPos, xpos2, horizontalLineYPos); // horizontal line
                g2d.drawLine(xpos2, posY2, xpos2, horizontalLineYPos); // vertical line 2

                int midpos = xpos1 + ((xpos2 - xpos1) / 2);



                double pValueWilcoxon = pvals[datasetNumber][category1][category2];
                String pValueWilcoxonString = (new java.text.DecimalFormat("0.#E0", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(pValueWilcoxon);
                if (pValueWilcoxon > 0.001) {
                    pValueWilcoxonString = (new java.text.DecimalFormat("##.###;-##.###", new java.text.DecimalFormatSymbols(java.util.Locale.US))).format(pValueWilcoxon);
                }
                //  + midpos 
                g2d.drawString(pValueWilcoxonString, midpos - (int) Math.round(getWidth(pValueWilcoxonString, g2d.getFont()) / 2), horizontalLineYPos - 5);

                comboCounter++;
            }
            deltaX += datasetwidth + betweenDatasetMargin;
            datasetCounter++;
        }

        //Draw expression level indicator:
        g2d.setComposite(alphaComposite100);

        g2d.setColor(new Color(100, 100, 100));
        int startVal = (int) Math.ceil(minValue);
        int endVal = (int) Math.floor(maxValue);
        int posY1 = marginTop + innerHeight - (int) Math.round((double) innerHeight * (startVal - minValue) / (maxValue - minValue));
        int posY2 = marginTop + innerHeight - (int) Math.round((double) innerHeight * (endVal - minValue) / (maxValue - minValue));

        g2d.drawLine(marginLeft
                - 10, posY1, marginLeft - 10, posY2);
        g2d.setFont(fontBold);
        for (int v = startVal;
                v <= endVal;
                v++) {
            int posY = marginTop + innerHeight - (int) Math.round((double) innerHeight * (v - minValue) / (maxValue - minValue));
            g2d.drawLine(marginLeft - 10, posY, marginLeft - 20, posY);
            g2d.drawString(String.valueOf(v), marginLeft - 25 - (int) getWidth(String.valueOf(v), g2d.getFont()), posY + 3);
        }

        g2d.translate(marginLeft
                - 60, marginTop + innerHeight / 2);
        g2d.rotate(
                -0.5 * Math.PI);
        g2d.drawString(
                "Relative gene expression", -(int) getWidth(
                "Relative gene expression", g2d.getFont()) / 2, 0);
        g2d.rotate(
                +0.5 * Math.PI);
        g2d.translate(
                -(marginLeft - 60), -(marginTop + innerHeight / 2));




        g2d.dispose();
        if (output == Output.PDF) {
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(outputFileName));
        }

        Locale.setDefault(defaultLocale);
    }

    public void drawViolinPlot(Graphics2D g2d, int x, int y, int width, int height, double[] vals, double minValue, double maxValue) {

        int nrVals = vals.length;

        //Determine range of values:
        double minVals = JSci.maths.ArrayMath.min(vals);
        double maxVals = JSci.maths.ArrayMath.max(vals);

        //Make frequency distribution:
        int nrBins = 1 + (int) Math.round(Math.sqrt(nrVals) / 2d);
        int[] binCount = new int[nrBins];
        for (int n = 0; n < nrBins; n++) {
            double lower = minVals + (maxVals - minVals) * (double) n / (double) nrBins;
            double upper = minVals + (maxVals - minVals) * (double) (n + 1) / (double) nrBins;
            for (int v = 0; v < nrVals; v++) {
                if (vals[v] >= lower && vals[v] < upper) {
                    binCount[n]++;
                }
            }
        }

        //Smooth the distribution:
        int posYMin = y + height - (int) Math.round((double) height * (maxVals - minValue) / (maxValue - minValue));
        int posYMax = y + height - (int) Math.round((double) height * (minVals - minValue) / (maxValue - minValue));
        double[] posVal = new double[posYMax - posYMin + 1];
        for (int pos = posYMin; pos <= posYMax; pos++) {
            double value = (((-pos + y + height) * (maxValue - minValue)) / (double) height) + minValue;
            for (int n = 0; n < nrBins; n++) {
                double lower = minVals + (maxVals - minVals) * (double) n / (double) nrBins;
                double upper = minVals + (maxVals - minVals) * (double) (n + 1) / (double) nrBins;
                if (value >= lower && value < upper) {
                    posVal[pos - posYMin] = binCount[n];
                }

            }
        }
        double kernelWidth = 10;
        double[] kernelWeights = new double[201];
        for (int d = -100; d <= 100; d++) {
            double weight = java.lang.Math.pow(java.lang.Math.E, -((double) d / (double) kernelWidth) * ((double) d / (double) kernelWidth) / 2);
            kernelWeights[d + 100] = weight;
        }
        double[] posValSmoothed = new double[posYMax - posYMin + 1];
        for (int pos = posYMin; pos <= posYMax; pos++) {
            double valSmoothed = 0;
            double sumWeights = 0;
            for (int q = pos - 100; q <= pos + 100; q++) {
                if (q >= posYMin && q <= posYMax) {
                    sumWeights += kernelWeights[100 + q - pos];
                    valSmoothed += kernelWeights[100 + q - pos] * posVal[q - posYMin];
                }
            }
            posValSmoothed[pos - posYMin] = valSmoothed / sumWeights;
        }
        double maxSmoothedVal = JSci.maths.ArrayMath.max(posValSmoothed);
        for (int pos = posYMin; pos <= posYMax; pos++) {
            posValSmoothed[pos - posYMin] /= maxSmoothedVal;
        }

        //Draw shape:
        GeneralPath path = new GeneralPath();
        path.moveTo(x + width / 2, posYMin);
        for (int pos = posYMin; pos <= posYMax; pos++) {
            path.lineTo(x + width / 2 - posValSmoothed[pos - posYMin] * width / 2d - 1, pos);
        }
        for (int pos = posYMax; pos >= posYMin; pos--) {
            path.lineTo(x + width / 2 + posValSmoothed[pos - posYMin] * width / 2d + 1, pos);
        }
        path.closePath();
        g2d.draw(path);
        g2d.fill(path);

    }

    public void drawBoxPlot(Graphics2D g2d, int x, int y, int width, int height, double[] vals, double minValue, double maxValue, boolean drawOutliers) {

        double median = JSci.maths.ArrayMath.percentile(vals, 0.50d);
        double q1 = JSci.maths.ArrayMath.percentile(vals, 0.25d);
        double q3 = JSci.maths.ArrayMath.percentile(vals, 0.75d);
        double iqr = q3 - q1;


        //Draw median:
        int posY = y + height - (int) Math.round((double) height * (median - minValue) / (maxValue - minValue));
        g2d.setStroke(new java.awt.BasicStroke(2.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawLine(x, posY, x + width, posY);
        //Draw IQR:
        int posY1 = y + height - (int) Math.round((double) height * (q3 - minValue) / (maxValue - minValue));
        int posY2 = y + height - (int) Math.round((double) height * (q1 - minValue) / (maxValue - minValue));
        g2d.setStroke(new java.awt.BasicStroke(1.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawRect(x, posY1, width, posY2 - posY1);

        //Draw whiskers:
        double whiskerTop = q3 + 1.5d * iqr;
        double whiskerBottom = q1 - 1.5d * iqr;
        double max = JSci.maths.ArrayMath.max(vals);
        double min = JSci.maths.ArrayMath.min(vals);
        if (min > whiskerBottom) {
            whiskerBottom = min;
        }
        if (max < whiskerTop) {
            whiskerTop = max;
        }
        posY = y + height - (int) Math.round((double) height * (whiskerTop - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, new float[]{2.0f}, 0.0f));
        g2d.drawLine(x + width / 2, posY, x + width / 2, posY1);
        g2d.setStroke(new java.awt.BasicStroke(1.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawLine(x + width / 2 - 5, posY, x + width / 2 + 5, posY);
        posY = y + height - (int) Math.round((double) height * (whiskerBottom - minValue) / (maxValue - minValue));
        g2d.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, new float[]{2.0f}, 0.0f));
        g2d.drawLine(x + width / 2, posY2, x + width / 2, posY);
        g2d.setStroke(new java.awt.BasicStroke(1.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_ROUND));
        g2d.drawLine(x + width / 2 - 5, posY, x + width / 2 + 5, posY);

        //Draw outliers:
        if (drawOutliers) {
            java.awt.AlphaComposite alphaComposite10 = java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.10f);
            g2d.setComposite(alphaComposite10);
            for (int v = 0; v < vals.length; v++) {
                if (vals[v] > whiskerTop || vals[v] < whiskerBottom) {
                    posY = y + height - (int) Math.round((double) height * (vals[v] - minValue) / (maxValue - minValue));
                    g2d.drawOval(x + width / 2 - 2, posY - 2, 5, 5);
                }
            }
        }
    }
}
