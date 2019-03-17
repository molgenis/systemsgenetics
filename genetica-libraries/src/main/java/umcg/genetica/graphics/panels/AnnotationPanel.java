package umcg.genetica.graphics.panels;


import umcg.genetica.features.Feature;
import umcg.genetica.features.FeatureTree;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.graphics.DefaultGraphics;

import java.awt.*;
import java.util.ArrayList;
import java.util.NavigableSet;

/**
 * Created by hwestra on 12/3/15.
 */
public class AnnotationPanel extends Panel {


	private Feature region;
	private ArrayList<FeatureTree> annotations;
	private ArrayList<SNPFeature> testOverlapWith;
	int trackheight = 10;
	int marginBetween = 2;

	public AnnotationPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(Feature region, ArrayList<FeatureTree> annotations) {
		this.region = region;
		this.annotations = annotations;
	}

	public void setOverlappingFeatures(ArrayList<SNPFeature> f) {
		testOverlapWith = f;
	}

	public int getTrackheight() {
		return trackheight;
	}

	public int getMarginBetween() {
		return marginBetween;
	}

	@Override
	public void draw(DefaultGraphics g) {
		Graphics2D g2d = g.getG2d();

		int figureWidth = width;
		int regionSize = region.getStop() - region.getStart();
		int nrPixelsX = figureWidth - (2 * marginX);


		Color defaultLightGrey = new Color(175, 175, 175);
		Color defaultColor = new Color(90, 90, 90);
		Color highlight = new Color(208, 83, 77);
		g2d.setColor(defaultColor);

		for (int i = 0; i < annotations.size(); i++) {


			int trackYpos = marginY + y0 + (i * trackheight) + (i * marginBetween);

			int startX = x0 + marginX;
			g2d.setColor(defaultLightGrey);
			g2d.drawLine(startX, trackYpos + trackheight, startX + width, trackYpos + trackheight);
			g2d.setColor(defaultLightGrey);
			// get subset within region
			FeatureTree t = annotations.get(i);
			NavigableSet<Feature> set = t.getFeatureSet(region);
			ArrayList<Feature> list = new ArrayList<Feature>();
			if (set != null) {
				for (Feature f : set) {
					if (f != null) {
						list.add(f);
					}
				}

			}


			boolean[] mark = new boolean[list.size()];

			if (testOverlapWith != null) {
				mark = testOverlap(list);
			}

			for (int q = 0; q < list.size(); q++) {
				Feature f = list.get(q);
				int start = f.getStart();
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

				int y1 = trackYpos;

				int boxwidth = pixelStop - pixelStart;
				if (boxwidth <= 0) {
					boxwidth = 1;
				}

				Color col = g2d.getColor();
				if (mark[q]) {
					g2d.setColor(highlight);
				} else {
					g2d.setColor(defaultColor);

				}
				g2d.fillRect(pixelStart, y1, boxwidth, trackheight);
			}

			g2d.setColor(defaultColor);

			// plot the name of the annotation

			g2d.setFont(new Font("default", Font.PLAIN, 10));

			int pixelStart = x0 + marginX + 5 + nrPixelsX;
			g2d.drawString(t.getName(), pixelStart, trackYpos + trackheight);


		}


	}


	private boolean[] testOverlap(ArrayList<Feature> list) {
		boolean[] output = new boolean[list.size()];

		for (int i = 0; i < list.size(); i++) {
			Feature f = list.get(i);
			boolean overlap = false;
			for (int q = 0; q < testOverlapWith.size(); q++) {
				Feature f2 = testOverlapWith.get(q);
				if (f.overlaps(f2)) {
					overlap = true;
					break;
				}
			}
			output[i] = overlap;
		}
		return output;
	}


}
