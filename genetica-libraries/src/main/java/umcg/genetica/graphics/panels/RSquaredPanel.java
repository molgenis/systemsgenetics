package umcg.genetica.graphics.panels;

import umcg.genetica.containers.Pair;
import umcg.genetica.features.Feature;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/27/15.
 */
public class RSquaredPanel extends Panel {
	public RSquaredPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	Feature region;
	HashSet<Feature> sequencedRegions;
	String[] datasetNames;
	ArrayList<ArrayList<Feature>> inputBeforeImputation;
	ArrayList<ArrayList<Pair<Feature, Double>>> allRSquaredValues;

	public void setData(Feature region, HashSet<Feature> sequencedRegions, String[] datasetNames,
						ArrayList<ArrayList<Feature>> inputBeforeImputation, ArrayList<ArrayList<Pair<Feature, Double>>> allRSquaredValues) {
		this.region = region;
		this.sequencedRegions = sequencedRegions;
		this.datasetNames = datasetNames;
		this.inputBeforeImputation = inputBeforeImputation;
		this.allRSquaredValues = allRSquaredValues;
	}

	public void setData(Feature region,
						HashSet<Feature> sequencedRegions,
						String datasetName,
						ArrayList<Feature> inputBeforeImp,
						ArrayList<Pair<Feature, Double>> rsquareds) {
		this.region = region;
		this.sequencedRegions = sequencedRegions;
		this.datasetNames = new String[]{datasetName};

		this.inputBeforeImputation = new ArrayList<ArrayList<Feature>>();
		inputBeforeImputation.add(inputBeforeImp);

		this.allRSquaredValues = new ArrayList<ArrayList<Pair<Feature, Double>>>();
		allRSquaredValues.add(rsquareds);
	}

	@Override
	public void draw(DefaultGraphics g) {
		Graphics2D g2d = g.getG2d();

		int figureWidth = width;
		int regionSize = region.getStop() - region.getStart();
		int nrPixelsX = figureWidth - (2 * marginX);
		int nrPixelsY = height - (2 * marginY);
		int nrDatasets = allRSquaredValues.size();
		int starty = y0 + marginY + nrPixelsY;

		g2d.setColor(new Color(208, 83, 77));

		// plot sequenced regions
		Font defaultfont = g2d.getFont();
		g2d.setFont(new Font("default", Font.BOLD, 16));

//		g2d.drawString("Targeted regions in sequencing", margin, margin - 20);

		g2d.setFont(defaultfont);

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

				System.out.println(f.toString());
				System.out.println(pixelStart + "\t" + pixelStop);
				int y1 = marginY + y0 - 20;
				int boxwidth = pixelStop - pixelStart;
				if (boxwidth <= 0) {
					boxwidth = 1;
				}

				g2d.fillRect(pixelStart, y1, boxwidth, 10);
			}
		}

		Color grey = new Color(0, 0, 0, 28);
		g2d.setColor(grey);


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
					colors[i] = new Color(208, 83, 77);
					break;
				case 3:
					colors[i] = new Color(98, 182, 177);
					break;
				case 4:
					colors[i] = new Color(116, 156, 80);
					break;
				case 5:
					colors[i] = new Color(156, 67, 109);
					break;
				case 6:
					colors[i] = new Color(81, 65, 156);
					break;

			}
		}

		Range range = new Range(0, 0, 0, 0);

		double unit = 0.1;
		Font originalfont = g2d.getFont();
		g2d.setFont(new Font("default", Font.BOLD, 16));
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

		for (int dataset = 0; dataset < nrDatasets; dataset++) {
			ArrayList<Pair<Feature, Double>> rsquareds = allRSquaredValues.get(dataset);
			ArrayList<Feature> beforeImputation = inputBeforeImputation.get(dataset);

			if (rsquareds != null) {

				g2d.setColor(new Color(98, 182, 177));

				// plot input before imputation (if any)
				if (beforeImputation != null) {

					for (Feature f : beforeImputation) {
						if (region.overlaps(f)) {
							int pos = f.getStart();
							int relativeStart = pos - region.getStart();
							double percStart = (double) relativeStart / regionSize;
							int pixelStartX = x0 + marginX + (int) Math.ceil(percStart * nrPixelsX);
							int pixelY = 20;
							int dotsize = 4;
							g2d.fillOval(pixelStartX - (dotsize / 2), starty - pixelY - (dotsize / 2), dotsize, dotsize);
						}
					}
				}

				g2d.setColor(colors[dataset]);

				for (Pair<Feature, Double> pair : rsquareds) {
					Feature f = pair.getLeft();
					int pos = f.getStart();
					Double ar2 = pair.getRight();
					if (region.overlaps(f)) {
						// plot
						int relativeStart = pos - region.getStart();
						double percStart = (double) relativeStart / regionSize;
						int pixelStartX = x0 + marginX + +(int) Math.ceil(percStart * nrPixelsX);

						int pixelY = (int) Math.ceil(ar2 * nrPixelsY);

						int dotsize = 2 + (int) Math.ceil(ar2 * 10);
						g2d.fillOval(pixelStartX - (dotsize / 2), starty - pixelY - (dotsize / 2), dotsize, dotsize);
					}
				}
				int adv = metrics.stringWidth(datasetNames[dataset]);
				int hgt = metrics.getHeight();
				g2d.drawString(datasetNames[dataset], x0 + marginX, y0 + marginY - 30);
			}
		}


		g2d.setFont(originalfont);


		// draw line near R-squared == 0.8
		double yperc = 0.8;
		int pixelY = (int) Math.ceil(yperc * nrPixelsY);
		g2d.setColor(new Color(208, 83, 77));
		g2d.drawLine(x0 + marginX, starty - pixelY, x0 + marginX + nrPixelsX, starty - pixelY);

		// plotVariantsUniqueIneachDataset coordinates
		g2d.setColor(new Color(70, 67, 58));

		// y-axis

		g2d.drawLine(x0 + marginX - 10, starty, x0 + marginX - 10, starty - nrPixelsY);
		double steps = 0.1;

		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);

		for (double i = 0; i < 1 + steps; i += steps) {
			if (i <= 1) {
				int plusY = (int) Math.ceil(((double) i / 1) * nrPixelsY);
				g2d.drawLine(x0 + marginX - 13, starty - plusY, x0 + marginX - 7, starty - plusY);
				String formattedStr = decimalFormat.format(i);
				int adv = metrics.stringWidth(formattedStr);
				int hgt = metrics.getHeight();
				Dimension size = new Dimension(adv + 10, hgt + 10);
				g2d.drawString(formattedStr, x0 + marginX - (int) Math.ceil(size.getWidth()) - 10, starty - plusY + 5);
			}
		}

		// x-axis

		g2d.drawLine(x0 + marginX - 5, starty + 5, x0 + marginX + nrPixelsX + 5, starty + 5);

		int xunit = (int) Math.ceil(range.determineUnit(regionSize));
		while (regionSize / xunit < 10 && xunit > 1) {
			xunit /= 2;
		}
		while (regionSize / xunit > 10) {
			xunit *= 2;
		}


		for (int i = region.getStart(); i < region.getStop(); i++) {
			if (i % xunit == 0) {
				int relativeStart = i - region.getStart();
				double percStart = (double) relativeStart / regionSize;
				int pixelStart = (int) Math.ceil(percStart * nrPixelsX);
				g2d.drawLine(x0 + marginX + pixelStart, starty, x0 + marginX + pixelStart, starty + 10);
				String formattedString = decimalFormat.format(i);
				int adv = metrics.stringWidth(formattedString);
				int hgt = metrics.getHeight();

				g2d.drawString(formattedString, x0 + marginX + pixelStart - (int) Math.ceil((double) adv / 2), starty + 5 + 20);
			}
		}

	}


}
