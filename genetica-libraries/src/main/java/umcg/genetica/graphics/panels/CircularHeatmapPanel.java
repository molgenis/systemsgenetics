package umcg.genetica.graphics.panels;


import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.enums.SNPClass;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.math.Goniometry;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 9/12/16.
 */
public class CircularHeatmapPanel extends Panel {


	private double[][][] data; // format [dataset][group][bin]
	private String[] rownames;

	private String[] groupnames;

	private String[][] binnames;
	private SNPClass[][][] binAnnotations; // [group][bin][binAnnotations]
	private ArrayList<Triple<Integer, Integer, String>> groups;
	private Range range;
	private boolean[][] groupAnnotations; // [dataset][group]


	public CircularHeatmapPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}


	public void setData(String[] rowNames, String[] groupnames, double[][][] data) {
		this.rownames = rowNames;
		this.groupnames = groupnames;
		this.data = data;
	}

	public void setData(String[] rowNames, String[] groupnames, double[][] data) {
		this.rownames = rowNames;
		this.groupnames = groupnames;

		this.data = new double[data.length][data[0].length][1];
		for (int r = 0; r < data.length; r++) {
			for (int c = 0; c < data[r].length; c++) {
				this.data[r][c][0] = data[r][c];
			}
		}
	}

	public void setBinnames(String[][] binnames) {
		this.binnames = binnames;
	}

	public void setBinAnnotations(SNPClass[][][] binAnnotations) {
		this.binAnnotations = binAnnotations;
	}

	public void setGroupAnnotations(boolean[][] annotations) {
		this.groupAnnotations = annotations;
	}

	@Override
	public void draw(DefaultGraphics g) {

		Range rangeValues = null;
		if (range != null) {
			rangeValues = range;
		} else {
			rangeValues = new Range(data);
		}

		int nrRows = data.length + 1;
//		int nrCols = data[0].length + 1;

		int plotWidth = width - (2 * marginX);
		int plotHeight = height - (2 * marginY);

		int maxWidth = plotWidth;
		if (plotHeight < plotWidth) {
			maxWidth = plotHeight;
		}


		Graphics2D g2d = g.getG2d();


		double widthPerDataset = maxWidth / nrRows / 2;

		double startX = marginX + x0;
		double startY = marginY + y0;

		DefaultTheme theme = new DefaultTheme();

		System.out.println(marginX);
		System.out.println(marginY);
		System.out.println(x0);
		System.out.println(y0);
		System.out.println(maxWidth);

		g2d.setStroke(theme.getThickStroke());

		double degreesForRowNames = 45d;

		double angleOffSet = 90d - (degreesForRowNames / 2);
		double degreesPerSegment = (360d - degreesForRowNames) / data[0].length;

		double degreesPerGroup = 4;
		HashMap<Integer, Integer> colToGroup = null;
		if (groups != null) {
			// this is really really dumb, but what the hack.
			colToGroup = new HashMap<Integer, Integer>();
			System.out.println("found " + groups.size() + " groups");

			double groupdeg = (groups.size() - 1) * degreesPerGroup;
			degreesPerSegment = (360d - degreesForRowNames - groupdeg) / data[0].length;

			for (int group = 0; group < groups.size(); group++) {
				int s = groups.get(group).getLeft();
				int e = groups.get(group).getMiddle();
				for (int q = s; q < e; q++) {
					colToGroup.put(q, group);
				}
			}
		}

		int minalpha = 25;
		double[] alphabins = new double[5];
		for (int d = 0; d < alphabins.length; d++) {
			double perc = (1d / alphabins.length) * (d + 1);
			alphabins[d] = minalpha + ((255 - minalpha) * perc);
		}

		// draw values per dataset.
		for (int dataset = 0; dataset < data.length; dataset++) {
			// plot a white circle first
			g2d.setColor(Color.white);
			int dsWidth = (int) Math.floor(maxWidth - (widthPerDataset * dataset));
			System.out.println(dataset + "\t" + dsWidth);
			double remainder = (maxWidth - dsWidth) / 2;
			int x0 = (int) Math.floor(startX + remainder);
			int y0 = (int) Math.floor(startY + remainder);
			g2d.fillOval(x0, y0, dsWidth, dsWidth);

			Color color = theme.getColor(dataset);

			for (int groupNum = 0; groupNum < data[dataset].length; groupNum++) {

				// TODO: need to get rid of this
				int groupId = 0;
				if (colToGroup != null) {
					groupId = colToGroup.get(groupNum);
				}

				double angle0 = angleOffSet - (degreesPerSegment * groupNum) - (groupId * degreesPerGroup);
				double degreesPerSubCol = degreesPerSegment / data[dataset][groupNum].length;
				System.out.println("colvsangle 0:\t" + groupNum + "\t" + groupId + "\t" + angle0 + "\t" + degreesPerSegment + "\t" + degreesPerGroup);
				for (int bin = 0; bin < data[dataset][groupNum].length; bin++) {
					double angle1 = (angle0 - (degreesPerSubCol * bin));
					double value = data[dataset][groupNum][bin];

					if (Double.isNaN(value)) {
						// plot a grey thingy
						g2d.setColor(new Color(0, 0, 220));
						Arc2D arc = new Arc2D.Double(x0, y0, dsWidth, dsWidth, angle1, -degreesPerSubCol, Arc2D.PIE);
						g2d.fill(arc);
						if (data[dataset].length < 200) {
							g2d.setColor(new Color(220, 220, 220));
							g2d.draw(arc);
						}
					} else {
						double perc = rangeValues.getRelativePositionY(value);

						int alphabin = (int) Math.floor(perc * alphabins.length);
						if (alphabin >= alphabins.length) {
							alphabin = alphabins.length - 1;
						}

						Color color2 = new Color(color.getRed(), color.getGreen(), color.getBlue(), (int) Math.floor(alphabins[alphabin]));

						if (value >= 0) {
							g2d.setColor(color2);

							Arc2D arc = new Arc2D.Double(x0, y0, dsWidth, dsWidth, angle1, -degreesPerSubCol, Arc2D.PIE);
							g2d.fill(arc);

							if (data[dataset].length < 200) {

								g2d.setColor(new Color(220, 220, 220));
								g2d.draw(arc);
							}
						} else {
							System.out.println("Did not expect value: " + value + "\t" + groupNum + "\t" + bin);
						}
					}


				}

			}
		}

		// fill inner circle with white
		int dsWidth = (int) Math.floor(maxWidth - (widthPerDataset * data.length));
		double remainder = (maxWidth - dsWidth) / 2;
		int x0 = (int) Math.floor(startX + remainder);
		int y0 = (int) Math.floor(startY + remainder);

		// first draw the outlines for the groups, if any
		if (groups != null) {
			g2d.setColor(theme.getDarkerColor(theme.getDarkGrey(), 0.5));
			g2d.setStroke(theme.getStroke());
			for (int group = 0; group < groups.size(); group++) {
				int col0 = groups.get(group).getLeft() - 1;
				int col1 = groups.get(group).getMiddle() - 1;
				int diff = col1 - col0;

				double angle0 = angleOffSet - (degreesPerSegment * col0) - (group * degreesPerGroup);

				double deg = degreesPerSegment; // degreesPerSubCol * diff;
				double angle1 = angle0 - deg;

				int dsWidth2 = (int) Math.floor(maxWidth - (widthPerDataset * 0));
				double remainder2 = (maxWidth - dsWidth2) / 2;
				int x02 = (int) Math.floor(startX + remainder2);
				int y02 = (int) Math.floor(startY + remainder2);
				Arc2D arc2 = new Arc2D.Double(x02, y02, maxWidth, maxWidth, angle1, -deg, Arc2D.PIE);
				g2d.draw(arc2);

				Arc2D arc = new Arc2D.Double(x0, y0, dsWidth, dsWidth, angle1, -deg, Arc2D.PIE);
				g2d.draw(arc);

			}

		}

		g2d.setColor(Color.white);
		g2d.fillOval(x0 + 1, y0 + 1, dsWidth - 2, dsWidth - 2);


		// draw dataset names
		g2d.setFont(theme.getLargeFont());
		FontMetrics metrics = g2d.getFontMetrics(theme.getLargeFont());

		g2d.setColor(theme.getDarkGrey());
		for (int d = 0; d < rownames.length; d++) {
			int strlen = metrics.stringWidth(rownames[d]) / 2;
			int midpointX = (int) Math.floor(startX + (maxWidth / 2)) - strlen;
			int dsWidth1 = (int) Math.floor(maxWidth - (widthPerDataset * d));
			int dsWidth2 = (int) Math.floor(maxWidth - (widthPerDataset * (d + 1)));
			int remainder1 = (maxWidth - dsWidth1) / 2;
			int remainder2 = (maxWidth - dsWidth2) / 2;
			int midpointY1 = (int) Math.floor(startY + (remainder1));
			int midpointY2 = (int) Math.floor(startY + (remainder2));
			int actualMidPointY1 = (midpointY1 + midpointY2) / 2;
			g2d.drawString(rownames[d], midpointX, actualMidPointY1);
		}

		// draw group names
		for (int group = 0; group < data[0].length; group++) {

			double angle0 = degreesPerSegment * group;
			angle0 -= (degreesForRowNames) + (degreesForRowNames / 2) - (degreesPerSegment / 2);

			int groupId = 0;
			if (colToGroup != null) {
				groupId = colToGroup.get(group);
			}

			angle0 += (groupId * degreesPerGroup);

			// use the outer edge
			double originX = startX + (maxWidth / 2);
			double originY = startX + (maxWidth / 2);

			double radius = (double) (maxWidth + 120) / 2;
			Pair<Integer, Integer> xy0 = Goniometry.calcPosOnCircle(radius, originX, originY, angle0);
			drawRotate(g2d, xy0.getLeft(), xy0.getRight(), angle0, groupnames[group]);
			System.out.println(group + "\t" + groupnames[group] + "\t" + angle0);
		}

		// draw binAnnotations, if any....
		if (binAnnotations != null) {
			double radius = (double) (maxWidth + 30) / 2;
			// [group][bin][binAnnotations]
			for (int groupNum = 0; groupNum < binAnnotations.length; groupNum++) {

				// TODO: need to get rid of this
				int groupId = 0;
				if (colToGroup != null) {
					groupId = colToGroup.get(groupNum);
				}

				double angle0 = angleOffSet - (degreesPerSegment * groupNum) - (groupId * degreesPerGroup);
				angle0 = (degreesPerSegment * groupNum) + (groupId * degreesPerGroup) + 90 + (degreesForRowNames / 2);
				double degreesPerSubCol = degreesPerSegment / binAnnotations[groupNum].length;

				// iterate snps
				for (int bin = 0; bin < binAnnotations[groupNum].length; bin++) {
					double angle1 = (angle0 + (degreesPerSubCol * bin)) + (degreesPerSubCol / 2);
					System.out.println("colvsangle 0:\t" + groupNum + "\t" + groupId + "\t1: " + angle0 + "\t2: " + angle1 + "\t" + degreesPerSegment + "\t" + degreesPerGroup);
					for (int a = 0; a < binAnnotations[groupNum][bin].length; a++) {
						if (bin == 3 && groupNum == 5) {
							System.out.println("Found it");
						}

						SNPClass c = binAnnotations[groupNum][bin][a];
						double radius2 = radius + (a * 15);
						// use the outer edge
						double originX = startX + (maxWidth / 2);
						double originY = startX + (maxWidth / 2);
						Pair<Integer, Integer> xy0 = Goniometry.calcPosOnCircle(radius2, originX, originY, angle1 - 180);

						int circlesize = 8;
						int halfcircle = 4;
						g2d.setColor(new Color(200, 200, 200));
						g2d.drawOval(xy0.getLeft() - halfcircle, xy0.getRight() - halfcircle, circlesize, circlesize);

						if (c != null && !c.equals(SNPClass.NONCODING) && !c.equals(SNPClass.INTRONIC)) {
							Color col = theme.getColor(c.getNumber());
							g2d.setColor(col);
							g2d.fillOval(xy0.getLeft() - halfcircle, xy0.getRight() - halfcircle, circlesize, circlesize);
						}
					}
				}
			}
		}

		// draw group annotations (if any)
		if (groupAnnotations != null) {
			dsWidth = (int) Math.floor(maxWidth - 2 * (widthPerDataset * (data.length)));
			double radius = dsWidth + (widthPerDataset / 2) - 15;

			// [dataset][group]
			for (int groupNum = 0; groupNum < groupAnnotations[0].length; groupNum++) {
				// TODO: need to get rid of this
				int groupId = 0;
				if (colToGroup != null) {
					groupId = colToGroup.get(groupNum);
				}

				double angle0 = (degreesPerSegment * groupNum) + (groupId * degreesPerGroup) + 90 + (degreesForRowNames / 2) + (degreesPerSegment / 2);
				double degreesPerSubCol = degreesPerSegment / binAnnotations[groupNum].length;
				for (int ds = 0; ds < groupAnnotations.length; ds++) {

					double radius2 = radius - (ds * 15);
					// use the outer edge
					double originX = startX + (maxWidth / 2);
					double originY = startX + (maxWidth / 2);

					Pair<Integer, Integer> xy0 = Goniometry.calcPosOnCircle(radius2, originX, originY, angle0 - 180);

					int circlesize = 8;
					int halfcircle = 4;
					g2d.setColor(new Color(200, 200, 200));
					g2d.drawOval(xy0.getLeft() - halfcircle, xy0.getRight() - halfcircle, circlesize, circlesize);

					if (groupAnnotations[ds][groupNum]) {
						Color col = theme.getColor(0);
						g2d.setColor(col);
						g2d.fillOval(xy0.getLeft() - halfcircle, xy0.getRight() - halfcircle, circlesize, circlesize);
					}
				}
			}
		}

		// draw a small legend in the top left corner
		for (int d = 0; d < data.length; d++) {
			Color color = theme.getColor(d);
			int y = (int) Math.floor(startY + (d * 15));
			for (int q = 0; q < alphabins.length; q++) {
				int x = (int) Math.floor(startX + (q * 10));
				Color color2 = new Color(color.getRed(), color.getGreen(), color.getBlue(), (int) Math.floor(alphabins[q]));
				g2d.setColor(color2);
				g2d.fillRect(x, y, 10, 10);
			}
		}


	}

	public void setGroups(ArrayList<Triple<Integer, Integer, String>> groups) {
		this.groups = groups;
	}

	public void setRange(Range range) {
		this.range = range;
	}
}
