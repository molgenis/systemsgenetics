/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.chromosomeyexpressionplotter;

import java.awt.*;
import java.awt.image.*;

import java.io.File;
import java.util.ArrayList;

import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.util.RankDoubleArray;

/**
 *
 * @author harm-jan
 */
public class ChromosomeYExpressionPlot {
    public void draw(TriTyperGeneticalGenomicsDataset[] gg, String outputFile){

        System.setProperty("java.awt.headless", "true");

        //Load all genetical genomics datasets:
        System.out.println("Generating chromosome Y expression plot:");
        int numDatasets = gg.length;
        int[] cumSamples = new int[numDatasets];
        int totalNrSamples = 0;

        
        for (int d = 0; d < numDatasets; d++) {
            cumSamples[d] = totalNrSamples;
            totalNrSamples += gg[d].getTotalGGSamples();
        }

        //Init image:
        int height = 450;
        int margin = 20;
        int width = totalNrSamples * 10 + margin * 2;
        int x0 = margin;
        int x1 = width - margin;
        int innerWidth = x1 - x0;
        int y0 = margin;
        int y1 = height - margin - 100;
        int innerHeight = y1 - y0;
        BufferedImage bi = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = bi.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height);


        //Iterate through probes that map to the Y chromosome:
        for (int d = 0; d < gg.length; d++) {
            int nrYProbes = 0;
            ArrayList<Double> vecProbeMeans = new ArrayList<Double>();
            int probeCount = gg[d].getExpressionData().getProbes().length;
            for (int probe = 0; probe < probeCount; probe++) {

                if (gg[d].getExpressionData().getChr()[probe] == 24) {
                    vecProbeMeans.add(gg[d].getExpressionData().getProbeMean()[probe]);
                }
            }
            java.util.Collections.sort(vecProbeMeans);
            double meanChrYProbeExpressionThreshold = 0;
            if (vecProbeMeans.size() > 0) {
                meanChrYProbeExpressionThreshold = vecProbeMeans.get((int) (vecProbeMeans.size() * 0.75d));
            }
            for (int probe = 0; probe < probeCount; probe++) {
                if (gg[d].getExpressionData().getChr()[probe] == 24) {
                    if (gg[d].getExpressionData().getProbeMean()[probe] > meanChrYProbeExpressionThreshold) {
                        nrYProbes++;
                    }
                }
            }

            double[][] yProbeData = new double[gg[d].getTotalGGSamples()][nrYProbes];
            nrYProbes = 0;
            double[][] rawData = gg[d].getExpressionData().getMatrix();

            for (int probe = 0; probe < probeCount; probe++) {
                if (gg[d].getExpressionData().getChr()[probe] == 24) {
                    if (gg[d].getExpressionData().getProbeMean()[probe] > meanChrYProbeExpressionThreshold) {
                        //float[] chrYProbeRanked = gg[d].rankReturnFloats(gg[d].rawData[probe]);
                        for (int ind = 0; ind < gg[d].getTotalGGSamples(); ind++) {
                            yProbeData[ind][nrYProbes] = rawData[probe][ind];
                            //yProbeData[ind][nrYProbes]=chrYProbeRanked[ind];
                        }
                        nrYProbes++;
                    }
                }
                /*
                if (gg[d].probeChrName[probe]==23) {
                if (gg[d].probeHUGO[probe].equals("XIST")) {
                for (int ind=0; ind<gg[d].sampleCount; ind++) {
                int indWGA = gg[d].indWGA[ind];
                String sex = "Male"; if (gg[d].model.indIsFemale[indWGA]==1) sex = "Female";
                String sampleName = (String) gg[d].model.vectorInd.get(indWGA);
                System.out.println(gg[d].probeName[probe] + "\t" + gg[d].probeHUGO[probe] + "\t" + sampleName + "\t" + gg[d].rawData[probe][ind] + "\t" + sex);
                }
                }
                }
                 */
            }

            System.out.println(" - Dataset:\t" + gg[d].getSettings().name + "\tNumber of probes mapping to chromosome Y:\t" + nrYProbes);

            //No Y probes are present, quit method:
            if (nrYProbes > 0) {

                double[] meanArray = new double[gg[d].getTotalGGSamples()];
                for (int ind = 0; ind < gg[d].getTotalGGSamples(); ind++) {
                    meanArray[ind] = JSci.maths.ArrayMath.mean(yProbeData[ind]);
                }
                double min = JSci.maths.ArrayMath.min(meanArray);
                double max = JSci.maths.ArrayMath.max(meanArray);

                //Get the rank of each individual:
                RankDoubleArray rda = new RankDoubleArray();
                double[] meanArrayRanked = rda.rank(meanArray);
                int[] intIndex = new int[gg[d].getTotalGGSamples()];
                for (int ind = 0; ind < gg[d].getTotalGGSamples(); ind++) {
                    intIndex[ind] = (int) meanArrayRanked[ind];
                }
                for (int ind = 0; ind < gg[d].getTotalGGSamples(); ind++) {
                    int indIncrement = 0;
                    for (int ind2 = ind + 1; ind2 < gg[d].getTotalGGSamples(); ind2++) {
                        if (intIndex[ind] == intIndex[ind2]) {
                            indIncrement++;
                            intIndex[ind2] += indIncrement;
                        }
                    }
                }

                //Define dashed stroke:
                float dash[] = {3.0f};
                java.awt.BasicStroke strokeDashed = new java.awt.BasicStroke(3.0f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_MITER, 3.0f, dash, 0.0f);

                //Now plot the individual median Y probe intensities:
                g2d.setFont(new java.awt.Font(g2d.getFont().getFontName(), java.awt.Font.PLAIN, 8));
                int[] ggindWGA = gg[d].getExpressionToGenotypeIdArray();
                for (int ind = 0; ind < gg[d].getTotalGGSamples(); ind++) {
                    int indWGA = ggindWGA[ind];
                    g2d.setColor(new Color(0, 0, 255));
                    g2d.setStroke(strokeDashed);
                    String sex = "Male";
                    if (gg[d].getGenotypeData().getIsFemale()[indWGA]) {
                        sex = "Female";
                        g2d.setStroke(new java.awt.BasicStroke(3));
                        g2d.setColor(new Color(255, 0, 0));
                    }
                    int plotX = x0 + intIndex[ind] * 10 + cumSamples[d] * 10;
                    int plotY = y1 - (int) ((double) innerHeight * (meanArray[ind] - min) / (max - min));
                    g2d.drawLine(plotX, y1 + 1, plotX, plotY);

                    g2d.setColor(new Color(90, 90, 90));
                    g2d.setStroke(new java.awt.BasicStroke(2));
                    g2d.drawLine(plotX, y1 + 2, plotX, y1 + 5);

                    g2d.setColor(new Color(0, 0, 0));
                    g2d.setFont(new java.awt.Font(g2d.getFont().getFontName(), java.awt.Font.PLAIN, 8));
                    String sampleName = gg[d].getGenotypeData().getIndividuals()[indWGA];

                    g2d.translate(plotX + 3, height - 5);
                    g2d.rotate(-Math.PI / 2.0);
                    int sampleNamePixelWidth = getPixelWidthOfString(sampleName, g2d.getFont());
                    g2d.drawString(sampleName, 100 - sampleNamePixelWidth, 0);
                    g2d.rotate(+Math.PI / 2.0);
                    g2d.translate(-(plotX + 3), -(height - 5));
                    //System.out.println(sampleName + "\t" + meanArray[ind] + "\t" + sex);

                }

                int plotX = x0 + cumSamples[d] * 10 + 10;
                g2d.setColor(new Color(0, 0, 0));
                g2d.setFont(new java.awt.Font(g2d.getFont().getFontName(), java.awt.Font.BOLD, 12));
                g2d.drawString(gg[d].getSettings().name, plotX, y0 + 4 + 3);
                g2d.setFont(new java.awt.Font(g2d.getFont().getFontName(), java.awt.Font.BOLD, 8));
                g2d.setStroke(new java.awt.BasicStroke(3));
                g2d.setColor(new Color(255, 0, 0));
                g2d.drawLine(plotX, y0 + 20, plotX + 20, y0 + 20);
                g2d.setColor(new Color(0, 0, 0));
                g2d.drawString("Female average chromosome Y probe expression", plotX + 25, y0 + 20 + 3);
                g2d.setColor(new Color(0, 0, 255));
                g2d.setStroke(strokeDashed);
                g2d.drawLine(plotX, y0 + 35, plotX + 20, y0 + 35);
                g2d.setColor(new Color(0, 0, 0));
                g2d.drawString("Male average chromosome Y probe expression", plotX + 25, y0 + 35 + 3);

                g2d.setColor(new Color(90, 90, 90));
                g2d.setStroke(new java.awt.BasicStroke(2));
                g2d.drawLine(x0 + cumSamples[d] * 10 - 3, y1 + 2, x0 + (cumSamples[d] + gg[d].getTotalGGSamples()) * 10 - 3, y1 + 2);
                g2d.drawLine(x0 + cumSamples[d] * 10 - 3, y1 + 2, x0 + cumSamples[d] * 10 - 3, y0);
            }

        }


        //Save image:
        try {
            javax.imageio.ImageIO.write(bi, "png", new File(outputFile));
	    
        } catch (Exception e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }

        System.out.println("");
    }


     /**
     * gets the width in pixels of a string
     * @param text the text of which the width requires to be measured
     * @param font the font in which the text will be printed
     * @return width in pixels
     */
    public int getPixelWidthOfString(String text, Font font) {
        java.awt.Graphics2D g2d = (new java.awt.image.BufferedImage(1, 1, java.awt.image.BufferedImage.TYPE_INT_ARGB)).createGraphics();
        java.awt.font.TextLayout tL = new java.awt.font.TextLayout(text, font, g2d.getFontRenderContext());
        return (int) tL.getBounds().getWidth();
    }
}
