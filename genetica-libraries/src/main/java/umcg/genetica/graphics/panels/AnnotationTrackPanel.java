package umcg.genetica.graphics.panels;

import umcg.genetica.features.Feature;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;

import java.awt.*;

/**
 * Created by hwestra on 11/27/16.
 */
public class AnnotationTrackPanel extends Panel {

	private double[][][] data; // [dataset][group][bp]
	private Feature region;
	private int[] highlight;
	private String[][] groupnames;

	public AnnotationTrackPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(double[][][] data, // [dataset][group][bp]
						int[] highlight, // [bp]
						Feature region,
						String[][] groupnames) {
		this.data = data;
		this.highlight = highlight;
		this.region = region;
		this.groupnames = groupnames; // [dataset][group]
	}

	@Override
	public void draw(DefaultGraphics g) {


		Graphics2D g2d = g.getG2d();

		int nrbins = data[0][0].length;
		int regionsize = region.getStop() - region.getStart();

		double bpPerBin = regionsize / nrbins;


		int figureWidth = width;
		int regionSize = nrbins;
		int nrPixelsX = figureWidth - (2 * marginX);
		double pixelPerBin = (double) nrPixelsX / nrbins;
		int marginBetween = 5;

		Color defaultLightGrey = new Color(175, 175, 175);
		Color defaultColor = new Color(90, 90, 90);
		Color highlightColor = new Color(208, 83, 77);


		int pixelsPerBp = (int) Math.floor((double) nrPixelsX / regionSize);
		if (pixelsPerBp < 1) {
			pixelsPerBp = 1;
		}


		int nrgroups = 0;
		for (int dataset = 0; dataset < data.length; dataset++) {
			for (int group = 0; group < data[dataset].length; group++) {
				nrgroups++;
			}
		}

		int groupctr = 0;
		int trackheight = 25;
		int trackmargin = 5;

		for (int dataset = 0; dataset < data.length; dataset++) {
			for (int group = 0; group < data[dataset].length; group++) {
				// get max
				double max = 0;
				for (int bp = 0; bp < data[dataset][group].length; bp++) {
					if (data[dataset][group][bp] > max) {
						max = data[dataset][group][bp];
					}
				}


				int y0 = this.y0 + marginY + (groupctr * trackheight) + (groupctr * trackmargin);
				g2d.setColor(defaultLightGrey);
				for (int bin = 0; bin < data[dataset][group].length; bin++) {
					int x1 = x0 + (int) Math.floor(marginX + (bin * pixelPerBin));
					double v = data[dataset][group][bin];
					if (v > 0) {
						v /= max;
						int h = (int) Math.floor(trackheight * v);
						int y1 = y0 - h;
						g2d.fillRect(x1, y1, pixelsPerBp, h);
					}
				}

				// draw group line
				g2d.setColor(defaultColor);
				g2d.drawLine(x0 + marginX, y0, x0 + marginX + nrPixelsX, y0);

				// print group name
				if (groupnames != null) {
					g2d.drawString(groupnames[dataset][group], x0 + marginX + nrPixelsX + 10, y0);
				}
				groupctr++;
			}
		}

		Range r = new Range(region.getStart(), 0, region.getStop(), 1);
		if (highlight != null) {
			g2d.setColor(highlightColor);
			for (int bp : highlight) {
				double perc = r.getRelativePositionX(bp);
				int y0 = this.y0 + marginY - trackheight;
				int y1 = this.y0 + marginY + (groupctr * trackheight) + (groupctr * trackmargin) - trackheight;
				int x1 = this.x0 + marginX + (int) Math.floor((perc * nrPixelsX));
				g2d.drawLine(x1, y0, x1, y1);

			}
		}

	}
}
