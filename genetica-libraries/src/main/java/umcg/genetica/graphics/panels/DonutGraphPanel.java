package umcg.genetica.graphics.panels;


import umcg.genetica.graphics.DefaultGraphics;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;

public class DonutGraphPanel extends Panel {
	
	
	private double[] data;
	private String text;
	
	public DonutGraphPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}
	
	class DonutGraphPanelData {
	
	}
	
	public void setData(double[] vals) {
		this.data = vals;
	}
	
	public void setText(String text) {
		this.text = text;
	}
	
	@Override
	public void draw(DefaultGraphics defaultGraphics) {
		
		Graphics2D g2d = defaultGraphics.getG2d();
		g2d.setStroke(new BasicStroke(0));
		
		Ellipse2D.Double circle = new Ellipse2D.Double();
		circle.setFrame(x0, y0, width, width);
		
		Area circlearea = new Area(circle);
		
		Ellipse2D.Double innercircle = new Ellipse2D.Double();
		int innerciclewidh = 60;
		int halfinnercirclewidth = innerciclewidh / 2;
		innercircle.setFrame(x0 + halfinnercirclewidth, y0 + halfinnercirclewidth, width - innerciclewidh, height - innerciclewidh);
		
		Area innercirclearea = new Area(innercircle);
		circlearea.subtract(innercirclearea);
		
		
		g2d.setColor(theme.getLightGrey());
		g2d.fill(circlearea);
		
		
		// plot the percentages
		
		for (int p = 0; p < data.length; p++) {
			double d = data[p];
			int radius = width;
			int start = 90;
			int stop = (int) Math.ceil(-360 * d); // start + (int) Math.ceil(d * 360);
			
			Arc2D.Double arc = new Arc2D.Double(x0, y0, radius, radius, start, stop, Arc2D.PIE);
			Area arcarea = new Area(arc);
			arcarea.subtract(innercirclearea);
			
			g2d.setColor(theme.getColor(p + 1));
			
			g2d.fill(arcarea);
			
		}
		
		
		Arc2D arc = new Arc2D.Double();
		
		g2d.setColor(theme.getDarkGrey());
		g2d.setFont(theme.getLargeFont());
		FontMetrics metrics = g2d.getFontMetrics();
		if (text != null) {
			int x = x0 + (width / 2);
			int y = y0 + (width / 2);
			int lnctr = 0;
			for (String line : text.split("\n")) {
				
				int tx = x - (metrics.stringWidth(line) / 2);
				
				if (lnctr > 0) {
					y += g2d.getFontMetrics().getHeight();
				}
				g2d.drawString(line, tx, y);
				lnctr++;
			}
		}
//		g2d.draw(circle);
		
		
	}
}
