package umcg.genetica.graphics.panels;


import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;

import java.awt.*;
import java.text.DecimalFormat;

/**
 * Created by hwestra on 9/12/15.
 */
public class HeatmapPanel extends Panel {
	private double[][] data;
	private Range range;
	private String[][] labels;
	private String[] rowLabels;
	private String[] colLabels;
	private MODE plotMode;
	private Range rangeLower;
	private Range rangeUpper;

	int baseOpacity = 25;
	Color[] colors = new Color[]{
			new Color(80, 131, 127),
			new Color(124, 89, 148)
	};

	public void setRange(Range range) {
		this.range = range;
	}

	public void setPlotMode(MODE plotMode) {
		this.plotMode = plotMode;
	}

	public void setRangeLower(Range rangeLower) {
		this.rangeLower = rangeLower;
	}

	public void setRangeUpper(Range rangeUpper) {
		this.rangeUpper = rangeUpper;
	}

	public enum MODE {
		UPPER,
		LOWER,
		TWODS, FULL
	}

	public HeatmapPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	public void setData(double[][] data) {
		this.data = data;
	}

	public void setData(double[][] data, String[] rowLabels, String[] colLabels) {
		this.data = data;
		this.rowLabels = rowLabels;
		this.colLabels = colLabels;
	}

	private Color changeOpacity(Color c, int v) {
		try {
			Color o = new Color(c.getRed(), c.getGreen(), c.getBlue(), v);
			return o;
		} catch (IllegalArgumentException e) {
			System.out.println(c.getRed());
			System.out.println(c.getGreen());
			System.out.println(c.getBlue());
			System.out.println(v);
			System.exit(-1);

		}
		return null;
	}

	@Override
	public void draw(DefaultGraphics g) {

		if (range == null) {
			// determine min and max
			range = determineRange(data);
		}


		Graphics2D g2d = g.getG2d();

		// plot boxes
		int plotWidthX = width - (2 * marginX);
		int plotWidthY = height - (2 * marginY);
		int boxHeight = plotWidthY / data.length;
		int boxWidth = plotWidthX / data[0].length;

		int startY = marginY + y0;
		int startX = marginX + x0;

		if (plotMode.equals(MODE.FULL)) {
			for (int i = 0; i < data.length; i++) {
				for (int j = 0; j < data[i].length; j++) {
					plotBox(startX, startY, i, j, boxWidth, boxHeight, range, g2d, plotMode);
				}
			}
		} else if (plotMode.equals(MODE.UPPER)) {
			for (int i = 0; i < data.length; i++) {
				for (int j = i + 1; j < data[i].length; j++) {
					plotBox(startX, startY, i, j, boxWidth, boxHeight, rangeUpper, g2d, plotMode);
				}
			}
		} else if (plotMode.equals(MODE.LOWER)) {
			for (int i = 0; i < data.length; i++) {
				for (int j = 0; j < i + 1; j++) {
					plotBox(startX, startY, i, j, boxWidth, boxHeight, rangeLower, g2d, plotMode);
				}
			}
		} else {
			// upper
			for (int i = 0; i < data.length; i++) {
				for (int j = i + 1; j < data[i].length; j++) {
					plotBox(startX, startY, i, j, boxWidth, boxHeight, rangeUpper, g2d, MODE.UPPER);
				}
			}
			// lower
			for (int i = 0; i < data.length; i++) {
				for (int j = 0; j < i + 1; j++) {
					plotBox(startX, startY, i, j, boxWidth, boxHeight, rangeLower, g2d, MODE.LOWER);
				}
			}
		}

		// draw a small legend
		if (plotMode.equals(MODE.TWODS) && rangeLower != null && rangeUpper != null) {
			// plot two legends
			for (int r = 0; r < 2; r++) {

				Range tmprange = rangeLower;

				if (r == 1) {
					tmprange = rangeUpper;
				}

//				double rangeTicks = tmprange.getRangeY() / 5;
//				int tickNr = 0;
				String pattern = "###,###,###.##";
				DecimalFormat decimalFormat = new DecimalFormat(pattern);
				double stepsize = tmprange.getRangeY() / 5;
				int dx = x0 + plotWidthX + 10 + (r * 100);
				for (int tickNr = 0; tickNr < 5; tickNr++) {
					double v = tmprange.getMinY() + (tickNr * stepsize);
					int dy = y0 - (boxHeight * (tickNr - 1));

					double yPerc = tmprange.getRelativePositionY(v);
					if (yPerc > 1d) {
						yPerc = 1;
						System.out.println("Weird value: " + v + "\t" + tmprange.getMaxY() + "\t" + yPerc);
					} else if (yPerc < 0d) {

						System.out.println("Weird value: " + v + "\t" + tmprange.getMinY() + "\t" + yPerc);
						yPerc = 0;
					}
					int op = baseOpacity + (int) Math.ceil((255 - baseOpacity) * yPerc);

					Color c = changeOpacity(colors[r], op);
					g2d.setColor(c);


					g2d.fillRect(dx, dy, boxWidth, boxHeight);

					g2d.setColor(theme.getDarkGrey());
					g2d.drawString(decimalFormat.format(v), dx + boxWidth + 5, dy + boxHeight);
				}

				double v = tmprange.getMinY() + (5 * stepsize);
				int dy = y0 - (boxHeight * (5 - 1));
				int op = baseOpacity + (int) Math.ceil((255 - baseOpacity) * 1);
				Color c = changeOpacity(colors[r], op);
				
				g2d.setColor(c);
				g2d.fillRect(dx, dy, boxWidth, boxHeight);
				g2d.setColor(theme.getDarkGrey());

				g2d.drawString(decimalFormat.format(v), dx + boxWidth + 5, dy + boxHeight);
			}

		} else {
//			double rangeTicks = range.getRangeY() / 5;
//			int tickNr = 0;
//			String pattern = "###,###,###.##";
//			DecimalFormat decimalFormat = new DecimalFormat(pattern);
//			for (double i = range.getMinY(); i < range.getMaxY() + rangeTicks; i += rangeTicks) {
//
//				int dy = y0 - (boxHeight * tickNr);
//
//				double yPerc = range.getRelativePositionY(i);
//				if (yPerc > 1) {
//					yPerc = 1;
//					System.out.println("Weird value: " + i + "\t" + range.getMaxY());
//				} else if (yPerc < 0) {
//
//					System.out.println("Weird value: " + i + "\t" + range.getMinY());
//					yPerc = 0;
//				}
//				int op = baseOpacity + (int) Math.ceil((255 - baseOpacity) * yPerc);
//
//
//				Color c = new Color(0, 128, 255, op);
//
//				g2d.setColor(c);
//
//				int dx = x0 + plotWidthX + 10;
//				g2d.fillRect(dx, dy, boxWidth, boxHeight);
//
//				g2d.setColor(theme.getDarkGrey());
//				g2d.drawString(decimalFormat.format(i), dx + boxWidth + 5, dy + 10);
//				tickNr++;
//			}
		}


		// draw the labels, if any..
		if (colLabels != null || rowLabels != null) {
			g2d.setColor(theme.getDarkGrey());
			g2d.setFont(theme.getSmallFont());
			FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());
			int fontheight = metrics.getHeight();
			if (colLabels != null) {
				int y = startY - 5;
				for (int col = 0; col < colLabels.length; col++) {
					int x = startX + (col * boxWidth);
					String str = colLabels[col];
					int widthOfStr = metrics.stringWidth(str);
					drawRotate(g2d, x + fontheight, y, -90, str);
				}
			}

			if (rowLabels != null) {
				int x = startX - 5;
				for (int row = 0; row < colLabels.length; row++) {
					int y = startY + (row * boxHeight) + fontheight;
					String str = rowLabels[row];
					int widthOfStr = metrics.stringWidth(str);
					g2d.drawString(str, x - widthOfStr, y);
				}
			}

		}
	}


	private void plotBox(int startX, int startY, int i, int j, int boxWidth, int boxHeight, Range tmprange, Graphics2D g2d, MODE m) {
		int dy = (i * boxHeight) + startY;
		int dx = startX + (j * boxWidth);
		double v = data[i][j];

		if (!Double.isNaN(v)) {
			double yPerc = tmprange.getRelativePositionY(v);
			// generate an opacity level
			int op = baseOpacity + (int) Math.ceil((255 - baseOpacity) * yPerc);

			// range opacity between 50 and 255

			if (op > 255 || op < 0) {
				System.out.println("ERROR: value out of range: " + v + "\t" + tmprange.getMinX() + " - " + tmprange.getMaxY() + "\tcoords: " + i + " x " + j);
				System.exit(-1);
			}

			Color c = changeOpacity(colors[0], op);
			if (m.equals(MODE.UPPER)) {
				c = changeOpacity(colors[1], op);
			}

			g2d.setColor(c);
			g2d.fillRect(dx, dy, boxWidth, boxHeight);
		}

	}

	public void drawRotate(Graphics2D g2d, double x, double y, int angle, String text) {
		g2d.translate((float) x, (float) y);
		g2d.rotate(Math.toRadians(angle));
		g2d.drawString(text, 0, 0);
		g2d.rotate(-Math.toRadians(angle));
		g2d.translate(-(float) x, -(float) y);
	}

	private Range determineRange(double[][] data) {
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				if (data[i][j] > max) {
					max = data[i][j];
				}
				if (data[i][j] < min) {
					min = data[i][j];
				}
			}
		}

		return new Range(0, min, 0, max);
	}

}
