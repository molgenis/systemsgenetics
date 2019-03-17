package umcg.genetica.graphics.panels;

import umcg.genetica.containers.Pair;
import umcg.genetica.features.Feature;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.themes.DefaultTheme;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/27/15.
 */
public class AssociationPanel extends Panel {

	Feature region;
	HashSet<Feature> sequencedRegions;
	ArrayList<ArrayList<Pair<Integer, Double>>> allPValues;
	String[] datasetNames;
	double maxPval = Double.NaN;
	boolean plotGWASSignificance = true;
	private ArrayList<Pair<Integer, Double>> ld;
	private boolean[][] markDifferentShape;
	private double threshold;
	private double[] LDData;
	private boolean roundUpYAxis = true;

	public ArrayList<ArrayList<Pair<Integer, Double>>> getAllPValues() {
		return allPValues;
	}

	public AssociationPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setDataSingleDs(Feature region,
								HashSet<Feature> sequencedRegions,
								ArrayList<Pair<Integer, Double>> allPValues,
								String datasetName) {
		this.region = region;
		this.sequencedRegions = sequencedRegions;
		ArrayList<ArrayList<Pair<Integer, Double>>> tmp = new ArrayList<ArrayList<Pair<Integer, Double>>>();
		tmp.add(allPValues);
		this.allPValues = tmp;
		this.datasetNames = new String[]{datasetName};
	}

	public void setData(Feature region, HashSet<Feature> sequencedRegions, ArrayList<ArrayList<Pair<Integer, Double>>> allPValues, String[] datasetNames) {
		this.region = region;
		this.sequencedRegions = sequencedRegions;
		this.allPValues = allPValues;
		this.datasetNames = datasetNames;
	}

	public void setMaxPVal(double d) {
		maxPval = d;
	}

	public void setPlotGWASSignificance(boolean b, double threshold) {
		plotGWASSignificance = b;
		this.threshold = threshold;
	}

	public double getMaxP() {
		return maxPval;
	}

	public String getTitle() {
		return title;
	}

	@Override
	public void draw(DefaultGraphics g) {

		Graphics2D g2d = g.getG2d();

		int figureWidth = width;

		int figureHeight = height;
		int nrPixelsY = figureHeight - (2 * marginY);
		int nrDatasets = allPValues.size();

		int minDotSize = 2;
		int regionSize = region.getStop() - region.getStart();
		int nrPixelsX = figureWidth - (2 * marginX);

		int plotStarty = y0 + marginY + nrPixelsY;

		System.out.println(nrPixelsX);

		// draw sequenced regions
		Color red = new Color(231, 79, 19);
		g2d.setColor(red);
		Color highlight = red;

		// plot sequenced regions
		Font defaultfont = theme.getMediumFont();
		g2d.setFont(defaultfont);

		if (sequencedRegions != null) {
			g2d.drawString("Targeted regions in sequencing", x0 + marginX, y0 + marginY - 20);
			for (Feature f : sequencedRegions) {
				if (f.overlaps(region)) {
					int start = f.getStart();
					int featurewidth = f.getStop() - start;
					int relativeStart = start - region.getStart();
					if (relativeStart < 0) {
						featurewidth -= Math.abs(relativeStart);
						relativeStart = 0;
					}

					int relativeStop = relativeStart + featurewidth;
					if (relativeStop > region.getStop()) {
						relativeStop = regionSize;
					}

					double percStart = (double) relativeStart / regionSize;
					double percStop = (double) relativeStop / regionSize;

					int pixelStart = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);
					int pixelStop = x0 + marginX + (int) Math.ceil(percStop * nrPixelsX);

					int y1 = marginY + y0 - 20;

					int boxwidth = pixelStop - pixelStart;
					if (boxwidth <= 0) {
						boxwidth = 1;
					}

					g2d.fillRect(pixelStart, y1, boxwidth, 10);
				}
			}
		}

		// determine max Y
		if (Double.isNaN(maxPval)) {
			System.out.println("Determining max P");
			maxPval = -Double.MAX_VALUE;
			for (int d = 0; d < allPValues.size(); d++) {
				ArrayList<Pair<Integer, Double>> pvals = allPValues.get(d);
				for (Pair<Integer, Double> p : pvals) {
					if (p.getRight() > maxPval) {
						maxPval = p.getRight();
					}
				}
			}
		}


		Color[] colors = new Color[nrDatasets];
		for (int i = 0; i < colors.length; i++) {
			switch (i) {
				case 0:
					colors[i] = new Color(70, 67, 58);
					break;
				case 1:
					colors[i] = new Color(174, 164, 140);
					break;
				case 2:
					colors[i] = red;
					break;
				case 3:
					colors[i] = new Color(98, 182, 177);
					break;
				case 4:
					colors[i] = new Color(116, 156, 80);
					break;

			}
		}


		Range range = new Range(0, 0, 0, 0);
		double unit = range.determineUnit(maxPval);
		double remainder = maxPval % unit;
		System.out.println(maxPval + " maxp");
		System.out.println(unit + " unit");
		System.out.println(remainder + " remainder");

		if (roundUpYAxis) {
			maxPval += (unit - remainder); // round off using unit
			System.out.println(maxPval + " maxp");
		}

		g2d.setFont(defaultfont);
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

		if (ld != null) {
			DefaultTheme theme = new DefaultTheme();

			Color currentColor = g2d.getColor();
			g2d.setColor(theme.getLightGrey());

			for (int q = 0; q < ld.size(); q++) {
				Pair<Integer, Double> d = ld.get(q);
				Integer pos = d.getLeft();
				Double val = d.getRight();
				// x-coord
				int relativeStart = pos - region.getStart();
				double percStart = (double) relativeStart / regionSize;
				int pixelStartX = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);

				// y-coord


				int pixelY = (int) Math.ceil(val * nrPixelsY);

				int dotsize = 2;
				g2d.fillOval(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);// x-coord
			}
			g2d.setColor(currentColor);
		}


		// determine unit
		double steps = maxPval / 1;
		// draw red line near 5E-8)
		if (plotGWASSignificance) {

			double gwas = -Math.log10(threshold);
			if (maxPval >= gwas) {
				double yperc = gwas / maxPval;
				int pixelY = (int) Math.ceil(yperc * nrPixelsY);
				g2d.setColor(red);
				g2d.drawLine(x0 + marginX, plotStarty - pixelY, x0 + marginX + nrPixelsX, plotStarty - pixelY);


				g2d.setFont(theme.getSmallFont());
				int strwidth = getStringWidth(g2d, "Significance");


//				g2d.setColor(Color.white);
//				g2d.fillRect(x0 + marginX + 10 - 2, plotStarty - pixelY - 5, strwidth + 6, 10);
//				g2d.setColor(red);
//				g2d.drawString("" + threshold, x0 + marginX + nrPixelsX + 5, plotStarty - pixelY + 3);

			}
			g2d.setFont(defaultfont);
		}

		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);

		g2d.setFont(defaultfont);
		g2d.setColor(new Color(70, 67, 58));

		// y-axis
		g2d.drawLine(x0 + marginX - 5,
				plotStarty,
				x0 + marginX - 5,
				y0 + marginY);

		// tick lines
		for (double i = 0; i < maxPval + steps; i += steps) {
			if (i <= maxPval) {
				int plusY = (int) Math.ceil((i / maxPval) * nrPixelsY);
				g2d.drawLine(x0 + marginX - 10, plotStarty - plusY, x0 + marginX - 5, plotStarty - plusY);
				String formattedStr = decimalFormat.format(i);
				int adv = metrics.stringWidth(formattedStr);
				int hgt = metrics.getHeight();
				Dimension size = new Dimension(adv + 10, hgt + 10);
				g2d.drawString(formattedStr, x0 + marginX - (int) size.getWidth() - 10, plotStarty - plusY + 5);
			}
		}

		// x-axis
		g2d.drawLine(x0 + marginX, plotStarty + 5, x0 + marginX + nrPixelsX, plotStarty + 5);

		int xunit = (int) Math.ceil(range.determineUnit(regionSize));
		while (regionSize / xunit < 10 && xunit > 1) {
			xunit /= 2;
		}
		while (regionSize / xunit > 10) {
			xunit *= 2;
		}

		int[] tickpos = new int[]{
				region.getStart() + (xunit - (region.getStart() % xunit)),
				region.getStop() - ((region.getStop() % xunit))
		};

		for (int i = 0; i < tickpos.length; i++) {

			int relativeStart = tickpos[i] - region.getStart();
			double percStart = (double) relativeStart / regionSize;
			int pixelStart = (int) Math.ceil(percStart * nrPixelsX);
			g2d.drawLine(x0 + marginX + pixelStart, plotStarty + 5, x0 + marginX + pixelStart, plotStarty + 10);
			String formattedString = decimalFormat.format(tickpos[i]);
			int adv = metrics.stringWidth(formattedString);
			int hgt = metrics.getHeight();
			g2d.drawString(formattedString, x0 + marginX + pixelStart - (adv / 2), plotStarty + 25);

		}


		if (LDData != null) {
			for (int i = 0; i < 101; i++) {
				double ld = i / 100d;
				Color c2 = highlight;
				Color c1 = new Color(153, 153, 153, 125); // colors[0];
				Color interpolate = g.interpolateColor(c2, c1, ld);
				g2d.setColor(interpolate);
				g2d.fillRect(10 + i, 10, 1, 10);
			}
			g2d.setColor(colors[0]);

			g2d.drawString("0", 10, 35);
			g2d.drawString("1", 105, 35);
		}

		// plot the values
		for (int z = allPValues.size() - 1; z > -1; z--) {
			ArrayList<Pair<Integer, Double>> toPlot = allPValues.get(z);
			boolean[] mark = null;
			if (markDifferentShape != null) {
				mark = markDifferentShape[z];
			}
			if (toPlot != null) {
				g2d.setColor(colors[z]);

				for (int v = 0; v < toPlot.size(); v++) {

					if (mark == null || !mark[v]) {
						dot(g2d, toPlot, v, z, regionSize, nrPixelsX, nrPixelsY, colors, highlight, g, plotStarty, false);
					}
				}

				for (int v = 0; v < toPlot.size(); v++) {
					if (mark != null && mark[v]) {
						dot(g2d, toPlot, v, z, regionSize, nrPixelsX, nrPixelsY, colors, highlight, g, plotStarty, true);
					}
				}

				int adv = metrics.stringWidth(datasetNames[z]);
				int hgt = metrics.getHeight();
				Dimension size = new Dimension(adv + 10, hgt + 10);
				g2d.drawString(datasetNames[z], x0 + marginX, y0 + marginY - 30);
			}
		}

	}


	private void dot(Graphics2D g2d,
					 ArrayList<Pair<Integer, Double>> toPlot,
					 int v,
					 int z,
					 int regionSize,
					 int nrPixelsX,
					 int nrPixelsY,
					 Color[] colors,
					 Color highlight,
					 DefaultGraphics g,
					 int plotStarty,
					 boolean mark) {
		Pair<Integer, Double> p = toPlot.get(v);
		Integer pos = p.getLeft();
		Double pval = p.getRight();

		// x-coord
		int relativeStart = pos - region.getStart();
		double percStart = (double) relativeStart / regionSize;
		int pixelStartX = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);

		// y-coord
		double yperc = pval / maxPval;
		int pixelY = (int) Math.ceil(yperc * nrPixelsY);


		g2d.setStroke(new BasicStroke(0.5f));
		g2d.setColor(colors[z]);
		Color interpolate = null;
		if (LDData != null) {
			Color c2 = highlight;
			Color c1 = new Color(153, 153, 153, 125); // colors[0];
			double LD = LDData[v];
			if (Double.isNaN(LD) || LD < 0 && LD > 1) {
				LD = 0;
			}
			interpolate = g.interpolateColor(c2, c1, LD);
			int b = interpolate.getBlue();
			int r = interpolate.getRed();
			int gr = interpolate.getGreen();
			interpolate = new Color(r, gr, b, 190);
			g2d.setColor(interpolate);

		}
		if (mark) {
			int dotsize = 2 + (int) Math.ceil(yperc * 10);
			g2d.fillRect(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
//			g2d.setColor(Color.white);
//			g2d.drawRect(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
		} else {
			int dotsize = 2 + (int) Math.ceil(yperc * 10);
			g2d.fillOval(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
//			g2d.setColor(Color.white);
//			g2d.drawRect(pixelStartX - (dotsize / 2), plotStarty - pixelY - (dotsize / 2), dotsize, dotsize);
		}
	}

	public void setMarkDifferentShape(boolean[] markDifferentShape) {
		this.markDifferentShape = new boolean[1][0];
		this.markDifferentShape[0] = markDifferentShape;
	}

	public void setLDData(double[] LDData) {
		this.LDData = LDData;
	}

	public void setRoundUpYAxis(boolean roundUpYAxis) {
		this.roundUpYAxis = roundUpYAxis;
	}
}
