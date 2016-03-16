/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.graphics;

import com.itextpdf.text.DocumentException;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Locale;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harm-jan
 */
public class ForestPlot {

    private String[] geneNames;

    public enum Output {

        PDF, PNG
    };
    private static final Font LARGE_FONT = new Font("Verdana", Font.PLAIN, 14);
    private static final Font LARGE_FONT_BOLD = new Font("Verdana", Font.BOLD, 14);
    private static final Font SMALL_FONT = new Font("Verdana", Font.PLAIN, 10);
    private static final int upMaxR = 108;
    private static final int upMaxG = 189;
    private static final int upMaxB = 69;
    
    private static final int downMaxR = 0;
    private static final int downMaxG = 174;
    private static final int downMaxB = 239;
    
    private static final float downMaxH = 196f/360f;
    private static final float downMaxS = 1f;
    private static final float downMaxBr = 0.93f;
    
    private static final float upMaxH = 100f/360f;
    private static final float upMaxS = 0.63f;
    private static final float upMaxBr = 0.74f;
    
    private static final Color red1Color = new Color(242, 101, 94, 50);
    private static final Color red2Color = new Color(242, 101, 94, 25);
    private static final Color gray1Color = new Color(237, 237, 237);
    private static final Color gray2Color = new Color(247, 247, 247);
    private static final Logger LOGGER = Logger.getLogger(Heatmap.class.getName());
    private static final Stroke dashed = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
    private static final Stroke line2pt = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
    private static final Stroke line = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

    public void drawForrestPlot(String xAxisName, String[] yAxisNames, Double[] xValues, String filename, Output output, double[] significanceThresholds, Double minX, Double maxX, int[] weights, byte[] chr, int[] chrpos, int metaRow) throws IOException, DocumentException {
        Double[][] workvalues = new Double[1][xValues.length];
        System.arraycopy(xValues, 0, workvalues[0], 0, xValues.length);
        drawMultiForrestPlot(xAxisName, yAxisNames, workvalues, filename, output, significanceThresholds, minX, maxX, weights, chr, chrpos, metaRow);
    }

    public void setGeneNames(String[] geneNames) {
        this.geneNames = geneNames;
    }

    public void drawMultiForrestPlot(String xAxisName, String[] yAxisNames, Double[][] xValues, String filename, Output output, double[] significanceThresholds, Double minX, Double maxX, int[] weights, byte[] chr, int[] chrpos, int metaRow) throws IOException, DocumentException {
        if (yAxisNames.length != xValues[0].length) {
            throw new IllegalArgumentException("Data length and number of row headers differ!");
        }

        // set locale for digit grouping

        Locale defaultLocale = Locale.getDefault();
        Locale.setDefault(Locale.US);
        // set up Graphics2D depending on required format using iText in case PDF
        Graphics2D g2d = null;
        com.itextpdf.text.Document document = null;
        com.itextpdf.text.pdf.PdfWriter writer = null;
        BufferedImage bi = null;
        int width = 1;
        int height = 1;
        bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();


        g2d.setFont(LARGE_FONT);
        FontMetrics fontmetrics = g2d.getFontMetrics();


        // set left margin
        int leftMargin = 10;
        int maxStringSize = 0;
        for (String s : yAxisNames) {
            maxStringSize = Math.max(maxStringSize, fontmetrics.stringWidth(s));
        }

        int topMargin = 10;
        int plotwidth = 100;

        int plotspacer = 20;

        int nrPlots = xValues.length;


        // width = plotwidth + max string size + 3 margins
        width = (plotwidth * nrPlots) + (nrPlots * plotspacer) + maxStringSize + (leftMargin * 3);
        System.out.println(width);
        // height is topmargin * 2 + yAxisNames * fontsize
        int textpadding = 5;
        int maxBoxSize = 20;
        int minBoxSize = 5;
        int fontheight = fontmetrics.getHeight();
        int geneNameMargin = 0;
        int geneBoxHeight = 15;

        if (geneNames != null) {
            geneNameMargin = (leftMargin * 2) + geneBoxHeight + fontheight * 3;
        }
        height = (yAxisNames.length * textpadding) + (2 * textpadding) + (fontheight * yAxisNames.length) + (topMargin * 2) + geneNameMargin + fontheight + topMargin;
        System.out.println(height);
        // initialize plot
        com.itextpdf.text.pdf.PdfContentByte cb = null;
        if (output == ForestPlot.Output.PDF) {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);
            writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(filename));

            document.open();
            cb = writer.getDirectContent();
            cb.saveState();
            //com.itextpdf.text.pdf.DefaultFontMapper fontMap = new com.itextpdf.text.pdf.DefaultFontMapper();
            g2d = cb.createGraphics(width, height);
        } else {
            bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2d = bi.createGraphics();
        }


        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, width, height);
        g2d.setStroke(line);
        // draw datasetnames
        g2d.setColor(Color.gray);



        for (int row = 0; row < yAxisNames.length; row++) {
            int textStartY = topMargin + fontheight + ((row * fontheight) + (row * textpadding)) + geneNameMargin;
            g2d.setFont(LARGE_FONT);
            g2d.drawString(yAxisNames[row], leftMargin, textStartY);
        }
        g2d.setColor(Color.black);

        for (int plotNr = 0; plotNr < nrPlots; plotNr++) {
            // draw forrest plot
            Double[] workValues = xValues[plotNr];

            double min = Double.MAX_VALUE;
            double max = Double.MIN_VALUE;

            if (maxX != null && minX != null) {
                max = maxX;
                min = minX;
            } else {
                for (int d = 0; d < workValues.length; d++) {
                    if (workValues[d] != null) {
                        if (workValues[d] > max) {
                            max = workValues[d];
                        }
                        if (workValues[d] < min) {
                            min = workValues[d];
                        }
                    }
                }
            }

            double origMin = min;
            double origMax = max;

            double absmin = 0;
            boolean correctValues = false;
            if (min < 0) {
                absmin = Math.abs(min);
                correctValues = true;
            }
            System.out.println("");
            if (correctValues) {
                min += absmin;
                max += absmin;
            } else {
                min -= absmin;
                max -= absmin;
            }

            Double[] correctedValues = new Double[workValues.length];
            for (int i = 0; i < correctedValues.length; i++) {
                // normalize to 0
                if (workValues[i] != null) {
                    System.out.println(i + "\t" + workValues[i]);
                    if (correctValues) {
                        correctedValues[i] = workValues[i] + absmin;

                    } else {
                        correctedValues[i] = workValues[i] - absmin;

                    }
                    System.out.println(i + "\t" + correctedValues[i]);
                    correctedValues[i] /= max;
                    System.out.println(i + "\t" + correctedValues[i]);
                }
            }


            int plotStartX = leftMargin + maxStringSize + leftMargin + (plotNr * plotspacer) + (plotNr * plotwidth);
            int plotStopX = plotStartX + plotwidth;
            int plotStartY = topMargin + geneNameMargin;
            int plotStopY = height - 2 * topMargin - fontheight;
            DecimalFormat df = new DecimalFormat("#.##");

            g2d.setColor(Color.gray);

            g2d.setFont(SMALL_FONT);
            int strWidth = getWidth(df.format(origMin), LARGE_FONT);
            g2d.drawString(df.format(origMin), plotStartX, plotStopY + 2 * leftMargin);
            strWidth = getWidth(df.format(origMax), LARGE_FONT);
            g2d.drawString(df.format(origMax), plotStopX - strWidth + 2, plotStopY + 2 * leftMargin);

            strWidth = getWidth(df.format(0d), LARGE_FONT);
            g2d.drawString(df.format(0d), plotStopX - (plotwidth / 2) - strWidth + 3, plotStopY + 2 * leftMargin);

            g2d.setFont(LARGE_FONT);

            if (geneNames != null) {

                int geneBoxStartY = plotStartY - leftMargin * 2 - geneBoxHeight - fontheight * 3;
                g2d.setColor(Color.gray);
                double metaZ = workValues[metaRow];
                double upperSignificance = significanceThresholds[significanceThresholds.length - 1];
                int alpha = 255;
                if (Math.abs(metaZ) <= Math.abs(upperSignificance)) {
                    alpha = 126;
                }
                if (metaZ > 0d) {
                    g2d.setColor(new Color(upMaxR, upMaxG, upMaxB, alpha));
                } else {
                    g2d.setColor(new Color(downMaxR, downMaxG, downMaxB, alpha));
                }

                // get MetaZ





                Color c = null;





                g2d.fillRect(plotStartX, geneBoxStartY, plotwidth, geneBoxHeight);

                int textStartY = geneBoxStartY + geneBoxHeight + leftMargin + leftMargin;
                Rectangle2D pos = fontmetrics.getStringBounds(geneNames[plotNr], g2d);



                g2d.setColor(Color.gray);
                String chrStr = "Chr: " + ChrAnnotation.parseByte(chr[plotNr]);
                DecimalFormat df2 = new DecimalFormat("###,###,###,###,###");
                String chrPos = "Pos: " + df2.format(chrpos[plotNr]);

                g2d.setFont(LARGE_FONT_BOLD);
                g2d.drawString(geneNames[plotNr], plotStartX, textStartY);
                g2d.setFont(SMALL_FONT);
                g2d.drawString(chrStr, plotStartX, textStartY + fontheight);
                g2d.drawString(chrPos, plotStartX, textStartY + fontheight * 2);
                g2d.setFont(LARGE_FONT);
            }



            System.out.println("Min " + min + "\tmax " + max);
            // draw axes
//            g2d.drawLine(plotStartX, plotStartY, plotStartX, plotStopY); // y-axis 1
            g2d.setColor(Color.gray);
//            g2d.drawLine(plotStartX + plotwidth / 2, plotStartY, plotStartX + plotwidth / 2, plotStopY); // y-axis at 0
//            g2d.drawLine(plotStartX, plotStopY, plotStopX, plotStopY); // x-axis

            int plotheight = plotStopY - plotStartY;

            g2d.setColor(gray1Color);
            g2d.fillRect(plotStartX, plotStartY, plotwidth, plotheight);

            g2d.setColor(Color.black);
            System.out.println("plot height: " + plotheight);
            System.out.println("plot width: " + plotwidth);
            System.out.println("start " + plotStartX);
            System.out.println("");
            if (significanceThresholds != null) {

                Color[] colors = new Color[4];
                colors[0] = new Color(120, 120, 120, 128);
                colors[1] = new Color(138, 138, 138, 128);
                colors[2] = new Color(197, 197, 197, 128);
                colors[3] = new Color(249, 249, 249, 128);

                Arrays.sort(significanceThresholds);



                g2d.setStroke(dashed);
                for (int i = significanceThresholds.length - 1; i > -1; i--) {
                    double val = significanceThresholds[i];

                    System.out.println(i + "\t" + val);
                    // correct to new scale
                    if (correctValues) {
                        val += absmin;
                    } else {
                        val -= absmin;
                    }


                    System.out.println(i + "\tminmax " + val);
                    val /= max;

                    System.out.println("Perc: " + val);
                    int nrPixels1 = (int) Math.floor(val * plotwidth);
                    int nrPixels2 = plotwidth - nrPixels1;
                    System.out.println(i + "\tdivmax " + nrPixels1);
                    System.out.println(i + "\tdivmax " + nrPixels2);



                    // thresholds are on both sides in this plot
                    int startX = plotStartX + nrPixels1;
                    int stopX = plotStartX + nrPixels2;
                    int startY = plotStartY;
                    int stopY = plotStopY;

                    int boxwidth = startX - stopX;

                    Color selected = null;
                    if (i > colors.length) {
                        selected = colors[0];
                    } else {
                        selected = colors[i];
                    }
                    g2d.setColor(selected);
                    if (i > significanceThresholds.length - 2) {
                        g2d.setStroke(line2pt);
                        g2d.setColor(red1Color);
                    } else {
                        g2d.setStroke(dashed);
                        g2d.setColor(red2Color);
                    }

                    g2d.drawLine(startX, startY, startX, stopY);
                    g2d.drawLine(stopX, startY, stopX, stopY);




                    System.out.println(startX);
                    System.out.println(startY);
                    System.out.println(boxwidth);
                    System.out.println(plotheight);

                }
                g2d.setColor(Color.black);
            }
            g2d.setStroke(line);


            // draw values

            double[] sqrtwghts = new double[weights.length];
            double maxWght = 0;

            // determine the total sample
            int nrTotalSamples = 0;
            for (int i = 0; i < metaRow; i++) {
                if (workValues[i] != null) {
                    nrTotalSamples += weights[i];
                }
            }
            sqrtwghts[metaRow] = Math.sqrt(nrTotalSamples);
            System.out.println("total samples: " + nrTotalSamples);
            maxWght = sqrtwghts[metaRow];
            for (int i = 0; i < metaRow; i++) {
                sqrtwghts[i] = Math.sqrt(weights[i]);
            }
            
            // for replication studies
            for (int i = metaRow+1; i < weights.length; i++) {
                sqrtwghts[i] = Math.sqrt(weights[i]);
            }

            for (int i = 0; i < correctedValues.length; i++) {
                Double val = correctedValues[i];
                if (val != null) {
                    int posX = plotStartX + (int) Math.floor(val * plotwidth); // original threshold (right side of 0)


                    int posY = topMargin + fontheight + (fontheight * i) + (i * textpadding) - (textpadding / 2) + geneNameMargin;

                    System.out.println(i + "\t" + val + "\t" + posX + "\t" + posY);

//                    GeneralPath gp = new GeneralPath();
//                    gp.moveTo(-5,0);
//                    gp.lineTo(+5, 0);
//                    gp.





                    double relativeWeight = sqrtwghts[i] / maxWght;
                    double range = maxBoxSize - minBoxSize;
                    System.out.println("rel. weight: " + weights.length + "\t" + sqrtwghts[i] + "\t" + relativeWeight);
                    int boxSize = minBoxSize + (int) Math.ceil(range * relativeWeight);
                    int halfbox = boxSize / 2;
                    int x = posX - halfbox;
                    int y = posY - halfbox;

                    // determine fraction of max
                    double realZ = Math.abs(workValues[i]);

                    double percOfMax = realZ / origMax;
                    int minAlpha = 75;
                    int maxAlpha = 255;
                    int actualAlpha = minAlpha + (int) Math.ceil((maxAlpha - minAlpha) * percOfMax);

                    Color c = null;
                    float minBr = 0.5f;
                    if (workValues[i] >= 0) {
                        System.out.println(yAxisNames[i] + "\tUp color: " + realZ);
                        
                        float s = (float) (upMaxS * percOfMax);
                        float br = minBr +(float) (upMaxBr * percOfMax);
                        if(br > upMaxBr){
                            br = upMaxBr;
                        }
                        c = Color.getHSBColor(upMaxH, s, br);
                    } else {
                        float s = (float) (downMaxS * percOfMax);
                        float br = minBr + (float) (downMaxBr * percOfMax);
                        if(br > downMaxBr){
                            br = downMaxBr;
                        }
                        c = Color.getHSBColor(downMaxH, s,  br);
                    }

                    g2d.setColor(c);
                    if (i == metaRow) {
                        // draw diamond
                        g2d.fill(drawDiamond(x, y, boxSize, boxSize));
                    } else {
                        g2d.fillRect(x, y, boxSize, boxSize);
                    }
                    g2d.setColor(Color.black);
                }
            }

            g2d.setColor(Color.gray);
            g2d.setStroke(line2pt);
            g2d.drawLine(plotStartX + plotwidth / 2, plotStartY, plotStartX + plotwidth / 2, plotStopY); // y-axis at 0
            g2d.setStroke(line);
            g2d.setColor(Color.black);
        }




        g2d.dispose();
        if (output == ForestPlot.Output.PDF) {
            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            bi.flush();
            ImageIO.write(bi, output.toString().toLowerCase(), new File(filename));
        }

        Locale.setDefault(defaultLocale);
    }

    public int getWidth(String text, java.awt.Font font) {
        java.awt.Graphics2D g2d = (new java.awt.image.BufferedImage(1, 1, java.awt.image.BufferedImage.TYPE_INT_ARGB)).createGraphics();
        java.awt.font.TextLayout tL = new java.awt.font.TextLayout(text, font, g2d.getFontRenderContext());
        return (int) tL.getBounds().getWidth();
    }

    public GeneralPath drawDiamond(int x, int y, int width, int height) {
        GeneralPath diamond = new GeneralPath(GeneralPath.WIND_EVEN_ODD, 4);

        x += width / 2;
        diamond.moveTo(x, y);
        x += width / 2;
        y += height / 2;
        diamond.lineTo(x, y);
        x -= width / 2;
        y += height / 2;
        diamond.lineTo(x, y);
        x -= width / 2;
        y -= height / 2;
        diamond.lineTo(x, y);
        return diamond;
    }
}
