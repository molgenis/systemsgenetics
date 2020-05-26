package umcg.genetica.graphics.panels;


import umcg.genetica.containers.Pair;
import umcg.genetica.features.Feature;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * Created by hwestra on 5/5/16.
 */
public class LDPanel extends Panel {

	private ArrayList<Pair<Integer, Integer>> positions;
	private ArrayList<Double> vals;
	private Feature region;
	private boolean plottopvalues = false;

	public LDPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}


	public void setData(Feature region, ArrayList<Pair<Integer, Integer>> positions, ArrayList<Double> vals) {
		this.positions = positions;
		this.vals = vals;
		this.region = region;
	}

	private boolean useQuadScale = false;

	public void scaleQuadratic(boolean b) {
		useQuadScale = true;
	}

	@Override
	public void draw(DefaultGraphics g) {


		// sort the vals on size
		ArrayList<SortObj> objs = new ArrayList<>();
		for (int i = 0; i < positions.size(); i++) {
			objs.add(new SortObj(vals.get(i), positions.get(i)));
		}

		Collections.sort(objs, new SortObjComp());

		Range r = new Range(region.getStart(), 0, region.getStop(), 1);
		double bp = r.getRangeX();
		int nrPixelsMaxX = width - (2 * marginX);
		int nrPixelsMaxY = height - (2 * marginY);
		double ppbp = nrPixelsMaxX / bp;
		if (ppbp < 1) {
			ppbp = 1;
		}

		System.out.println(nrPixelsMaxX + " max pixels");
		System.out.println(ppbp + " pixels per bp");
		int ppbpi = (int) Math.floor(ppbp);


		Graphics2D g2d = g.getG2d();

		for (int i = 0; i < objs.size(); i++) {

			double pos1 = objs.get(i).pos.getLeft();
			double pos2 = objs.get(i).pos.getRight();
			double perc1 = r.getRelativePositionX(pos1);
			double perc2 = r.getRelativePositionX(pos2);

			// flip positions when needed
			double tmpperc1 = perc1;
			if (perc2 < perc1) {
				tmpperc1 = perc1;
				perc1 = perc2;
				perc2 = tmpperc1;
			}
			if (perc1 >= 0 && perc2 >= 0 && perc1 <= 1 && perc2 <= 1) {

				double yperc = objs.get(i).p;
				if (useQuadScale) {
					yperc *= yperc;
				}
				int color = (int) Math.floor(255 - (yperc * 255));
				if (color < 0) {
					color = 0;
				}


				double startX1 = x0 + marginX + (perc1 * nrPixelsMaxX);
				double startY1 = y0 + marginY + (perc2 * nrPixelsMaxX);


				/// System.out.println(objs.get(i).p + "\t" + pos1 + "\t" + perc1 + "\t" + startX1 + "\t" + pos2 + "\t" + perc2 + "\t" + startY1);
				g2d.setColor(new Color(color, color, color));
				g2d.fillRect((int) Math.floor(startX1), (int) Math.floor(startY1), ppbpi, ppbpi);
				if (yperc == 1 && plottopvalues) {
					g2d.setColor(new Color(255, 0, 0));
					g2d.fillOval((int) Math.floor(startX1) - 5, (int) Math.floor(startY1) - 5, 10, 10);
				}

			}
		}
	}

	class SortObjComp implements Comparator<SortObj> {

		@Override
		public int compare(SortObj o1, SortObj o2) {
			if (o1.p > o2.p) {
				return 1;
			} else if (o1.p < o2.p) {
				return -1;
			} else {
				return 0;
			}
		}
	}

	class SortObj {
		public double p;
		public Pair<Integer, Integer> pos;

		public SortObj(double p, Pair<Integer, Integer> pos) {
			this.p = p;
			this.pos = pos;
		}

		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;

			SortObj sortObj = (SortObj) o;

			return Double.compare(sortObj.p, p) == 0;

		}

		@Override
		public int hashCode() {
			long temp = Double.doubleToLongBits(p);
			return (int) (temp ^ (temp >>> 32));
		}
	}

}


