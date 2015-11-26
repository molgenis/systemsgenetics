/*
 * Copyright (C) 2015 adriaan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.io.File;


/**
 *
 * @author adriaan
 * adapted from harm jan, : imputation-tool/src/main/java/imputationtool/postprocessing/ScatterPlot.java
 */
public class ASScatterPlot {
    
    private BufferedImage bi;
    private Graphics2D g2d;
    private int graphHeight;
    private int graphWidth;
    private int drawWidth;
    private int drawHeight;
    private Color color;
    private Font font;

    public ASScatterPlot( ) {
        init(1000, 1000);
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

        double xpos = ((double) graphWidth ) + (((double) graphWidth ) * x);
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
