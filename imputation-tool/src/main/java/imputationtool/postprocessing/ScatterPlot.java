/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputationtool.postprocessing;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.io.File;

/**
 *
 * @author harmjan
 */
public class ScatterPlot {

    private BufferedImage bi;
    private Graphics2D g2d;
    private int graphHeight;
    private int graphWidth;
    private int drawWidth;
    private int drawHeight;
    private Color color;
    private Font font;

    public ScatterPlot(int size) {
        init(size, size);
    }

    protected void init(int width, int height) {


        graphHeight = height;
        graphWidth = width;

        bi = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_RGB);
        g2d = bi.createGraphics();

        color = new Color(255, 255, 255);
        font = new Font("Georgia", Font.PLAIN, 10);

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2d.setColor(color);
        g2d.fillRect(0, 0, graphWidth, graphHeight);
        g2d.draw(new Line2D.Double(((double) graphWidth / 2), 0, ((double) graphWidth / 2), graphHeight));
        color = new Color(0, 0, 0);
        g2d.setColor(color);
    }

    void plot(int x, int y) {
    }

    // x == correlation
    // y == r2 score
    public void plot(double x, double y) {
        double xpos = ((double) graphWidth / 2) + (((double) graphWidth / 2) * x);
        double ypos = graphHeight - (graphHeight * y);
        int ixpos = (int) Math.floor(xpos);
        int iypos = (int) Math.floor(ypos);
        g2d.fillRect(ixpos - 1, iypos + 1, 2, 2);


    }

    public void draw(String outputFile) {

        try {
            javax.imageio.ImageIO.write(bi, "png", new File(outputFile));
        } catch (Exception e) {
            System.out.println(e.getMessage());
            System.out.println(e.getStackTrace());
        }

    }
}
