package umcg.genetica.graphics.panels;


import umcg.genetica.features.BedGraphFeature;
import umcg.genetica.features.Feature;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;

/**
 * Created by hwestra on 12/3/15.
 */
public class GraphAnnotationPanel extends Panel {
	Feature region;
	ArrayList<ArrayList<BedGraphFeature>> data;
	String[] filenames;
	double maxHeight = 100;
	int trackheight = 50;
	private Range range;
	private int marginBetween;

	public GraphAnnotationPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(Feature region, ArrayList<ArrayList<BedGraphFeature>> data, String[] bfiles) {
		this.region = region;
		this.data = data;
		this.filenames = bfiles;
	}

	public void setRange(Range r) {
		this.range = r;
	}

	@Override
	public void draw(DefaultGraphics g) {

		Graphics2D g2d = g.getG2d();

		int figureWidth = width;
		int regionSize = region.getStop() - region.getStart();
		int nrPixelsX = figureWidth - (2 * marginX);

		int marginBetween = 5;

		Color defaultLightGrey = new Color(175, 175, 175);
		Color defaultColor = new Color(90, 90, 90);
		Color highlight = new Color(208, 83, 77);
		g2d.setColor(defaultColor);

		if (range == null) {
			range = new Range(0, 0, maxHeight, maxHeight);
		}


		for (int i = 0; i < data.size(); i++) {

			// get subset within region
			ArrayList<BedGraphFeature> t = data.get(i);

			// check whether the data is stranded
			boolean isStranded = false;
			for (int q = 0; q < t.size(); q++) {
				BedGraphFeature f = t.get(q);
				if (f.isStranded()) {
					isStranded = true;
				}
			}

			int trackYpos = marginY + y0 + (i * trackheight) + (i * marginBetween);

			int startX = x0 + marginX;
			g2d.setColor(defaultLightGrey);

			if (!isStranded) {
				g2d.drawLine(startX, trackYpos + trackheight, startX + width, trackYpos + trackheight);
			} else {
				g2d.drawLine(startX, trackYpos + (trackheight / 2), startX + width, trackYpos + (trackheight / 2));
			}
			g2d.setColor(defaultColor);

			// plot the name of the annotation
			String name = filenames[i];
			File file = new File(name);
			String filename = file.getName();
			g2d.setFont(new Font("default", Font.BOLD, 8));

			int namePixelStart = x0 + marginX + 5 + nrPixelsX;
			if (!isStranded) {
				g2d.drawString(filename, namePixelStart, trackYpos + trackheight);
			} else {
				g2d.drawString(filename, namePixelStart, trackYpos + (trackheight / 2));
			}

			// plot the annotations
			for (int q = 0; q < t.size(); q++) {
				BedGraphFeature f = t.get(q);
				int start = f.getStart();
				int stop = f.getStop();
				if (start > region.getStart() && stop < region.getStop()) {
					int featurewidth = f.getStop() - start;
					int relativeStart = start - region.getStart();
					if (relativeStart < 0) {
						featurewidth -= Math.abs(relativeStart);
						relativeStart = 0;
					}

					int relativeStop = relativeStart + featurewidth;
					if (relativeStop > regionSize) {
						relativeStop = regionSize;
					}


					double percStart = (double) relativeStart / regionSize;
					double percStop = (double) relativeStop / regionSize;

					int pixelStart = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);
					int pixelStop = x0 + marginX + (int) Math.ceil(percStop * nrPixelsX);
					int boxwidth = pixelStop - pixelStart;
					if (boxwidth <= 0) {
						boxwidth = 1;
					}

					if (isStranded) {

						double halfTrackHeight = (double) trackheight / 2;
						double halfMaxY = range.getMaxY() / 2;
						double posVal = f.getPos();
						double negVal = f.getNeg();


						if (posVal > halfMaxY) {
							posVal = halfMaxY;
						}
						if (negVal > halfMaxY) {
							negVal = halfMaxY;
						}

						double percPos = posVal / halfMaxY;
						double percNeg = negVal / halfMaxY;

						double nrPixelsPos = percPos * halfTrackHeight;
						double nrPixelsNeg = percNeg * halfTrackHeight;

						int boxHeight = (int) Math.ceil(nrPixelsPos + nrPixelsNeg);

						int yStartPos = (int) Math.ceil(trackYpos + (trackheight / 2) - nrPixelsPos);

						g2d.fillRect(pixelStart, yStartPos, boxwidth, boxHeight);

					} else {
						int y1 = trackYpos;
						double v = f.getValue();
						int boxHeight = 0;

						if (v > range.getMaxY()) {
							v = range.getMaxY();
						}
						boxHeight = (int) ((Math.ceil((v / range.getMaxY()) * trackheight)));
						g2d.fillRect(pixelStart, y1 + trackheight - boxHeight, boxwidth, boxHeight);
					}
				}
			}

			g2d.setColor(defaultColor);


		}

	}


	public int getTrackheight() {
		return trackheight;
	}

	public int getMarginBetween() {
		return marginBetween;
	}
}
