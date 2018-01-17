/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.binarymeta.meta.graphics;

import com.itextpdf.text.Document;
import com.itextpdf.text.pdf.PdfContentByte;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.concurrent.ArrayBlockingQueue;

/**
 *
 * @author harmjan
 */
public class ZScorePlot {

    int width = 0;
    int height = 0;
    int spacer = 50;
    BufferedImage bimage;
    Graphics2D g2d;
    int numDatasets = 0;
    private int unitsize = 0;
    private int halfsize;
    private int plotsize;
    private ZScorePlotThread thread;
    private ArrayBlockingQueue<PlotPackage> queue;
    private String[] datasets;
    private boolean pdfOutput = false;
    private String outfilename = "";
    private Document document;
    private PdfContentByte cb;
    private com.itextpdf.text.pdf.PdfWriter writer;

    public void init(int numdatasets, String[] datasets, boolean pdf, String filename) {

        if (pdf) {
            pdfOutput = true;
            outfilename = filename + ".pdf";
        }


        this.datasets = datasets;
        plotsize = 400;
        numDatasets = numdatasets;
        halfsize = 200;
        unitsize = 5;
        spacer = 5 * unitsize * 3;

        width = (plotsize * numDatasets) + ((numdatasets + 1) * spacer) - (plotsize + spacer);
        height = (plotsize * numDatasets) + ((numdatasets + 1) * spacer) - (plotsize + spacer);

        if (pdfOutput) {
            com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
            document = new com.itextpdf.text.Document(rectangle);
            writer = null;
            try {
                writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(filename));
                document.open();
                cb = writer.getDirectContent();
                cb.saveState();
                g2d = cb.createGraphics(width, height);
            } catch (Exception e) {
                System.out.println("Cannot write to PDF file!:\t" + filename);
                System.exit(-1);
            }

        } else {
            bimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2d = bimage.createGraphics();
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2d.setColor(new Color(255, 255, 255));
            g2d.fillRect(0, 0, width, height);
            g2d.setColor(new Color(0, 0, 0));
        }

        queue = new ArrayBlockingQueue<PlotPackage>(8);

        thread = new ZScorePlotThread(queue, bimage, g2d, width, height, spacer, numdatasets, unitsize, halfsize, plotsize);

        thread.setName("Plotter");
        thread.start();
//        
    }

    public synchronized void draw(Double zScore, Double zScore2, int dataset1, int dataset2) {
//        if(dataset1 > dataset2){
//            int tmp = dataset2;
//            dataset2 = dataset1;
//            dataset1 = tmp;
//        }

        try {
            PlotPackage p = new PlotPackage();
            p.zscore1 = zScore;
            p.zscore2 = zScore2;
            p.dataset1 = dataset1;
            p.dataset2 = dataset2;
            queue.put(p);
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }
    }

    public void write(String outputFile) throws Exception {
        try {
            PlotPackage p = new PlotPackage();
            p.poison = true;
            queue.put(p);
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }

        thread.join();
        g2d.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1.00f));
        Color black = new Color(0, 0, 0);
        Color grey = new Color(150, 150, 150);
        g2d.setColor(grey);
        g2d.setStroke(new java.awt.BasicStroke(2.0f, java.awt.BasicStroke.CAP_SQUARE, java.awt.BasicStroke.JOIN_ROUND));
        AffineTransform atNormal = new AffineTransform();
        AffineTransform at90deg = new AffineTransform();
        atNormal.setToRotation(0);
        at90deg.setToRotation(-Math.PI / 2.0);

        Font f = new Font("Sans-Serif", Font.BOLD, 15);
        g2d.setFont(f);
        AffineTransform fontAT = new AffineTransform();
        fontAT.rotate(Math.toRadians(-90));
        Font font = g2d.getFont();
        Font theDerivedFont = font.deriveFont(fontAT);


        g2d.setFont(font);



        for (int dataset1 = 0; dataset1 < numDatasets; dataset1++) {
//           int dataset2 = 1;
            for (int dataset2 = dataset1 + 1; dataset2 < numDatasets; dataset2++) {

                g2d.setColor(grey);
                int leftmargin = (dataset2 + 1) * spacer + (dataset2 * plotsize) - (plotsize + spacer);
                int topmargin = (dataset1 + 1) * spacer + (dataset1 * plotsize);

                int startX = leftmargin + 0;
                int stopX = leftmargin + plotsize;
                int midX = leftmargin + halfsize;
                int startY = topmargin + 0;
                int stopY = topmargin + plotsize;
                int midY = topmargin + halfsize;

                g2d.drawLine(startX, midY, stopX, midY); // x-axis
                g2d.drawLine(midX, startY, midX, stopY); // x-axis

                for (int z = -40; z <= 40; z += 5) {
                    g2d.drawLine(midX + z * unitsize, midY, midX + z * unitsize, midY + unitsize); // x-ticks
                    g2d.drawLine(midX, midY - z * unitsize, midX + unitsize, midY - z * unitsize);
                }

                g2d.setColor(black);
                // axis labels
                g2d.setFont(font);

                g2d.drawString(datasets[dataset1], startX, midY - unitsize);
                g2d.setFont(theDerivedFont);

                g2d.drawString(datasets[dataset2], midX - unitsize, stopY);

            }
        }

        g2d.setFont(font);
        g2d.setColor(grey);
        for (int dataset1 = 0; dataset1 < numDatasets; dataset1++) {
//           int dataset2 = 1;
            for (int dataset2 = dataset1 + 1; dataset2 < numDatasets; dataset2++) {

                int leftmargin = (dataset2 + 1) * spacer + (dataset2 * plotsize) - (plotsize + spacer);
                int topmargin = (dataset1 + 1) * spacer + (dataset1 * plotsize);

                int boxdist = 5 * unitsize;

                g2d.setStroke(new java.awt.BasicStroke(1.0f, java.awt.BasicStroke.CAP_SQUARE, java.awt.BasicStroke.JOIN_ROUND));
                g2d.setColor(new Color(0, 0, 0));
                g2d.drawRect(leftmargin - boxdist, topmargin - boxdist, plotsize + boxdist * 2, plotsize + boxdist * 2);
            }
        }



        /*
         * drawLine(int x1, int y1, int x2, int y2)
         */


        if (pdfOutput) {
            g2d.dispose();
            cb.restoreState();
            document.close();
            writer.close();
        } else {
            g2d.dispose();
            javax.imageio.ImageIO.write(bimage, "png", new File(outputFile + "-ZScoreComparison.png"));
        }
    }
}