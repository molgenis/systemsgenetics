/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.graphics;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.text.DecimalFormat;
import java.util.concurrent.ArrayBlockingQueue;

/**
 *
 * @author harm-jan
 */
public class ZScorePlotThread extends Thread {

    int width = 0;
    int height = 0;
    int spacer = 50;
    BufferedImage bimage;
    Graphics2D g2d;
    int numDatasets = 0;
    private int unitsize = 0;
    private int halfsize;
    private int plotsize;
    ArrayBlockingQueue<PlotPackage> queue;
    private static final Color defaultColor = new Color(0, 0, 0, 128);
    private static final Color absentColor = new Color(255, 0, 0, 225);
    private final int[] individualSize;
    private final int[][] shared;
    private final int[][] identicaldirection;

    public ZScorePlotThread(ArrayBlockingQueue<PlotPackage> wqueue, BufferedImage image, Graphics2D g2dd,
            int width, int height, int spacer, int numdatasets, int unitsize, int halfsize, int plotsize) {
        this.queue = wqueue;
        this.bimage = image;
        this.g2d = g2dd;
        this.width = width;
        this.height = height;
        this.spacer = spacer;
        this.numDatasets = numdatasets;
        this.unitsize = unitsize;
        this.halfsize = halfsize;
        this.plotsize = plotsize;
        g2d.setColor(defaultColor);
        individualSize = new int[numdatasets];
        shared = new int[numdatasets][numdatasets];
        identicaldirection = new int[numdatasets][numdatasets];
//        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.15f));
    }

    @Override
    public void run() {
        boolean poison = false;
        while (!poison) {
            try {
                PlotPackage pack = queue.take();
                if (!pack.poison) {
                    plot(pack.zscore1, pack.zscore2, pack.dataset1, pack.dataset2);
                } else {
                    poison = pack.poison;
                }

                pack = null;
            } catch (InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        // draw final numbers
//        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0f));
        g2d.setColor(defaultColor);
        DecimalFormat df = new DecimalFormat("##");
        DecimalFormat df2 = new DecimalFormat("###,###,###,###,###");
        for (int dataset1 = 0; dataset1 < numDatasets; dataset1++) {
            for (int dataset2 = dataset1 + 1; dataset2 < numDatasets; dataset2++) {
                int leftmargin = (dataset2 + 1) * spacer + (dataset2 * plotsize) - (plotsize + spacer);
                int topmargin = (dataset1 + 1) * spacer + (dataset1 * plotsize);
                int posX = leftmargin;
                int posY = topmargin;
//                double percD1 = (shared[dataset1][dataset2]/)*100;
//                double percD2 = (shared[dataset1][dataset2]/)*100;
                double perc = ((double) identicaldirection[dataset1][dataset2] / shared[dataset1][dataset2]) * 100;
//                g2d.drawString("Nr QTLs in Dataset 1: " + + " (" +df.format(percD1)+")", posX, posY + 10);
//                g2d.drawString("Nr QTLs in Dataset 2: " + + " (" +df.format(percD2)+")", posX, posY + 22);
                g2d.drawString("Shared: " + df2.format(shared[dataset1][dataset2]), posX, posY + 10);
                g2d.drawString("Identical direction: " + df2.format(identicaldirection[dataset1][dataset2]) + " (" + df.format(perc) + "%)", posX, posY + 22);
            }
        }

    }

    private void plot(Double zScore, Double zScore2, int dataset1, int dataset2) {
        boolean absent1 = false;
        boolean absent2 = false;
        if (zScore == null) {
            zScore = 0d;
            absent1 = true;
        }
        if (zScore2 == null) {
            zScore2 = 0d;
            absent2 = true;
        }




        zScore = maxOut(zScore);
        zScore2 = maxOut(zScore2);
        int leftmargin = (dataset2 + 1) * spacer + (dataset2 * plotsize) - (plotsize + spacer);
        int topmargin = (dataset1 + 1) * spacer + (dataset1 * plotsize);
        int posX = leftmargin + halfsize + (int) Math.round(zScore * unitsize);
        int posY = topmargin + halfsize - (int) Math.round(zScore2 * unitsize);

        if (absent1 && absent2) {
            // don't do anything
        } else if (absent1 || absent2) {
            // Zscores absent in one or the other datasets get a different color :)
            g2d.setColor(absentColor);
            g2d.fillOval(posX - 3, posY - 3, 6, 6);
            g2d.setColor(defaultColor);
        } else {
            // shared
            shared[dataset1][dataset2]++;
            if (zScore >= 0 && zScore2 >= 0 || zScore < 0 && zScore2 < 0) {
                identicaldirection[dataset1][dataset2]++;
            }
            g2d.fillOval(posX - 3, posY - 3, 6, 6);
        }


    }

    private double maxOut(double zScore) {
        if (zScore > 40) {
            return 40;
        }
        if (zScore < -40) {
            return -40;
        }
        return zScore;
    }
}
