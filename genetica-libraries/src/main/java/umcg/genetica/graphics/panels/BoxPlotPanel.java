package umcg.genetica.graphics.panels;


import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.graphics.themes.Theme;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by hwestra on 9/12/15.
 */
public class BoxPlotPanel extends Panel {

	private double[][][] data;
	private Range plotRange;
	private boolean drawDataPoints;
	private boolean useMeanAndSd;
	private String outputIQRS;
	private String[] binLabels;
	private boolean useTukeysDefault;

	public BoxPlotPanel(int nrRows, int nrCols) {
		super(nrRows, nrCols);
	}

	@Override
	public void draw(DefaultGraphics g) {

		int nrDatasets = data.length;
		int nrBinClusters = data[0].length;

		double[][][] iqrs = new double[data.length][data[0].length][0];
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (int ds = 0; ds < nrDatasets; ds++) {
			for (int bin = 0; bin < nrBinClusters; bin++) {
				iqrs[ds][bin] = collectStats(data[ds][bin]);
				double localmin = iqrs[ds][bin][0];
				double localmax = iqrs[ds][bin][4];
//			System.out.println(localmin + "\t" + localmax);
				if (localmin < min) {
					min = localmin;
				}
				if (localmax > max) {
					max = localmax;
				}
			}
		}

		if (plotRange == null) {
			plotRange = new Range(0, min, 0, max);
		}

		int nrPixelsMaxX = width - (2 * marginX);
		int nrPixelsMaxY = height - (2 * marginY);

//		if (plottype == PLOTTYPE.CLUSTERED) {
		int marginBetweenBinClusters = 10;
		int marginBetweenDatasets = 2;

		// calculate bar width
		int widthPerBinCluster = ((nrPixelsMaxX) / nrBinClusters);
		int widthPerBinClusterBin = widthPerBinCluster - marginBetweenBinClusters;
		int widthPerDataset = (widthPerBinClusterBin / nrDatasets);
		int widthPerDatasetBin = widthPerDataset - marginBetweenDatasets;

		int starty = y0 + marginY;

		Theme theme = new DefaultTheme();

		Color c = theme.getDarkGrey();
		Graphics2D g2d = g.getG2d();
		int halfWidthOfDataset = widthPerDatasetBin / 2;

		try {
			TextFile iqrout = null;
			if (outputIQRS != null) {
				iqrout = new TextFile(outputIQRS, TextFile.W);
			}


			for (int binCluster = 0; binCluster < nrBinClusters; binCluster++) {
				for (int dataset = 0; dataset < nrDatasets; dataset++) {
					int startX = x0 + marginX + (binCluster * widthPerBinCluster) + (dataset * widthPerDataset);

					// retrieve IQRS
					double datamin = iqrs[dataset][binCluster][0];
					double dataq1 = iqrs[dataset][binCluster][1];
					double datamedian = iqrs[dataset][binCluster][2];
					double dataq3 = iqrs[dataset][binCluster][3];
					double datamax = iqrs[dataset][binCluster][4];
					double mean = iqrs[dataset][binCluster][5];
					double sd = iqrs[dataset][binCluster][6];

					double[] datasetDataForBin = data[dataset][binCluster];

					if (datasetDataForBin.length > 2) {
						if (drawDataPoints) {

							Pair<double[][], double[]> densityData = determineDensity(datasetDataForBin, plotRange);

							double[][] values = densityData.getLeft();
							double[] frequencies = densityData.getRight();

							for (int j = 0; j < frequencies.length; j++) {
								double[] dataInFreqBin = values[j];
								double maxPercJitter = frequencies[j];
								for (int v = 0; v < dataInFreqBin.length; v++) {
									double d = dataInFreqBin[v];

									if (d > plotRange.getMaxY()) {
										d = plotRange.getMaxY();
									}
									if (d < plotRange.getMinY()) {
										d = plotRange.getMinY();
									}

									// determine where this point falls in the plotRange
									double yperc = plotRange.getRelativePositionY(d);
									int relY = (int) Math.ceil(yperc * nrPixelsMaxY);
									int plotY = starty + nrPixelsMaxY - relY;

									// add some jitter to the x-direction
									// should be dependent on density

									int direction = 1;
									if (Math.random() > 0.5) {
										direction = -1;
									}

									int jittersize = (int) Math.ceil(Math.random() * Math.random() * (maxPercJitter * halfWidthOfDataset));
									int jitter = jittersize * direction;
									int plotX = startX + halfWidthOfDataset + jitter;
									g2d.setColor(c);
									if (d > dataq1 && d < dataq3) {
										// falls within IQR
										// set opacity to 30%;
										Color c2 = new Color(70, 67, 58, 128);
										g2d.setColor(c2);
									} else {
										g2d.setColor(c);
									}
									g2d.fillOval(plotX - 2, plotY - 2, 4, 4);
								}
							}
						} else {
							// just draw the box plot
							g2d.setColor(theme.getLightGrey());
							g2d.setStroke(theme.getStroke());
							// draw line from min to q1


							if (useMeanAndSd) {

							} else if (useTukeysDefault) {

								double iqr = (dataq3 - dataq1);
								double iqrmax = dataq3 + (1.5 * iqr);
								double iqrmin = dataq1 - (1.5 * iqr);

								// plot the outliers
								g2d.setColor(theme.getLightGrey());
								for (int i = 0; i < datasetDataForBin.length; i++) {
									double d = datasetDataForBin[i];
									if (d > iqrmax || d < iqrmin) {
										double pos = plotRange.getRelativePositionY(d);
										int outlierY = (int) Math.ceil(pos * nrPixelsMaxY);
										int outlierYPx = starty + nrPixelsMaxY - outlierY;
										g2d.fillOval(startX + halfWidthOfDataset - 2, outlierYPx - 2, 4, 4);
									}
								}

								g2d.setColor(theme.getColor(dataset));

								boolean clippingbottom = false;
								boolean clippingtop = false;

								if (iqrmin < plotRange.getMinY()) {
									iqrmin = plotRange.getMinY();
									clippingbottom = true;
								}

								if (iqrmax > plotRange.getMaxY()) {
									iqrmax = plotRange.getMaxY();
									clippingtop = true;
								}


								double iqrminy = plotRange.getRelativePositionY(iqrmin);
								double q1y = plotRange.getRelativePositionY(dataq1);
								double m2y = plotRange.getRelativePositionY(datamedian);
								double q3y = plotRange.getRelativePositionY(dataq3);
								double iqrmaxy = plotRange.getRelativePositionY(iqrmax);
								double meany = plotRange.getRelativePositionY(mean);

								int iqrMinYPx = (int) Math.ceil(iqrminy * nrPixelsMaxY);
								int q1yPx = (int) Math.ceil(q1y * nrPixelsMaxY);
								int m2yPx = (int) Math.ceil(m2y * nrPixelsMaxY);
								int q3yPx = (int) Math.ceil(q3y * nrPixelsMaxY);
								int iqrMaxYPx = (int) Math.ceil(iqrmaxy * nrPixelsMaxY);
								int meanYPx = (int) Math.ceil(meany * nrPixelsMaxY);

								int plotYiqrMinY = starty + nrPixelsMaxY - iqrMinYPx;
								int plotYq1yPx = starty + nrPixelsMaxY - q1yPx;
								int plotYm2yPx = starty + nrPixelsMaxY - m2yPx;
								int plotYq3yPx = starty + nrPixelsMaxY - q3yPx;
								int plotYiqrMaxY = starty + nrPixelsMaxY - iqrMaxYPx;
								int plotYmeanY = starty + nrPixelsMaxY - meanYPx;


								// draw the box
								g2d.fillRect(startX, plotYq3yPx, widthPerDatasetBin, plotYq1yPx - plotYq3yPx);

//								if (clippingbottom) {
//									g2d.setStroke(theme.getStrokeDashed());
//									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty + pixelsY, startX + halfWidthOfDataset + halfWidthOfDataset, starty + pixelsY);
//								}
//								g2d.setStroke(theme.getStroke());
//
//								// draw horizontal leg of the whisker
//								g2d.drawLine(startX + halfWidthOfDataset, plotYiqrMinY, startX + halfWidthOfDataset, plotYiqrMinY);
//
//								if (clippingtop) {
//									g2d.setStroke(theme.getStrokeDashed());
//									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty, startX + halfWidthOfDataset + halfWidthOfDataset, starty);
//								}

								// draw vertical parts of leg
								g2d.setStroke(theme.getStroke());
								g2d.drawLine(startX + halfWidthOfDataset, plotYiqrMinY, startX + halfWidthOfDataset, plotYq1yPx);
								g2d.drawLine(startX + halfWidthOfDataset, plotYq3yPx, startX + halfWidthOfDataset, plotYiqrMaxY);

								// draw the median
								g2d.setColor(Color.white);
								g2d.drawLine(startX, plotYm2yPx, startX + widthPerDatasetBin, plotYm2yPx);
								// draw the mean

								g2d.setColor(Color.white);
								g2d.fillOval(startX + halfWidthOfDataset - 5, plotYmeanY - 5, 10, 10);
								g2d.setColor(theme.getColor(dataset));
								g2d.fillOval(startX + halfWidthOfDataset - 3, plotYmeanY - 3, 6, 6);


								// System.out.println(startX + "\t" + plotYq1y + "\t" + widthPerDataset + "\t" + (plotYq1y - plotYq3y));

							} else {

								boolean clippingbottom = false;
								boolean clippingtop = false;

								if (datamin < plotRange.getMinY()) {
									datamin = plotRange.getMinY();
									clippingbottom = true;
								}

								if (datamax > plotRange.getMaxY()) {
									datamax = plotRange.getMaxY();
									clippingtop = true;
								}

								double m1y = plotRange.getRelativePositionY(datamin);
								double q1y = plotRange.getRelativePositionY(dataq1);
								double m2y = plotRange.getRelativePositionY(datamax);
								double q3y = plotRange.getRelativePositionY(dataq3);

								int m1yPx = (int) Math.ceil(m1y * nrPixelsMaxY);
								int q1yPx = (int) Math.ceil(q1y * nrPixelsMaxY);


								int m2yPx = (int) Math.ceil(m2y * nrPixelsMaxY);
								int q3yPx = (int) Math.ceil(q3y * nrPixelsMaxY);

								int plotYm1y = starty + nrPixelsMaxY - m1yPx;
								int plotYq1y = starty + nrPixelsMaxY - q1yPx;


								if (clippingbottom) {
									g2d.setStroke(theme.getStrokeDashed());
									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty + nrPixelsMaxY, startX + halfWidthOfDataset + halfWidthOfDataset, starty + nrPixelsMaxY);
								}
								g2d.setStroke(theme.getStroke());

								g2d.drawLine(startX + halfWidthOfDataset, plotYm1y, startX + halfWidthOfDataset, plotYq1y);


								int plotYm2y = starty + nrPixelsMaxY - m2yPx;
								int plotYq3y = starty + nrPixelsMaxY - q3yPx;

								if (clippingtop) {
									g2d.setStroke(theme.getStrokeDashed());
									g2d.drawLine(startX + halfWidthOfDataset - halfWidthOfDataset, starty, startX + halfWidthOfDataset + halfWidthOfDataset, starty);
								}
								g2d.setStroke(theme.getStroke());
								g2d.drawLine(startX + halfWidthOfDataset, plotYm2y, startX + halfWidthOfDataset, plotYq3y);


								// draw the median
								double q2y = plotRange.getRelativePositionY(datamedian);
								int q2yPx = (int) Math.ceil(q2y * nrPixelsMaxY);
								g2d.setColor(theme.getColor(1));
								int plotYq2y = starty + nrPixelsMaxY - q2yPx;
								g2d.fillOval(startX + halfWidthOfDataset - 3, plotYq2y - 3, 6, 6);
								g2d.setColor(theme.getColor(0));

								// draw a lightgrey box
								g2d.setColor(theme.getLightGrey());

								g2d.drawRect(startX, plotYq1y, widthPerDatasetBin, plotYq3y - plotYq1y);
							}


//						//g2d.drawRect(startX, plotYq3y, widthPerBox, plotYq1y - plotYq3y);
//						if (outputIQRS != null) {
//							iqrout.writeln(datamin + "\t" + dataq1 + "\t" + datamedian + "\t" + dataq3 + "\t" + datamax
//									+ "\t" + plotYm1y + "\t" + plotYq1y + "\t" + plotYq2y + "\t" + plotYq3y + "\t" + plotYm2y);
//						}
						}
					}


					// plot an axis
				}
			}
			if (outputIQRS != null) {
				iqrout.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		g2d.setColor(theme.getDarkGrey());

		// plot y-axis
		double tickUnitY = plotRange.getRangeY() / 10;
		String pattern = "###,###,###.##";
		DecimalFormat decimalFormat = new DecimalFormat(pattern);

		g2d.setFont(theme.getMediumFont());
		FontMetrics metrics = g2d.getFontMetrics(g2d.getFont());

		int xPosYAxis = x0 + marginX - 10;
		int yPosYAxis = y0 + marginY;
		int pixelsY = height - (2 * marginY);
		g2d.drawLine(xPosYAxis, yPosYAxis, xPosYAxis, yPosYAxis + pixelsY);

		int maxlen = 0;
		for (double y = plotRange.getMinY(); y < plotRange.getMaxY() + (tickUnitY / 2); y += tickUnitY) {
			double yPerc = plotRange.getRelativePositionY(y);

			int ypos = y0 + marginY + (int) Math.ceil((1 - yPerc) * pixelsY);
			int startx = xPosYAxis - 5;
			int stopx = xPosYAxis;
			g2d.drawLine(startx, ypos, stopx, ypos);
			String formattedStr = decimalFormat.format(y);
			int adv = metrics.stringWidth(formattedStr);
			if (adv > maxlen) {
				maxlen = adv;
			}
			g2d.setFont(theme.getMediumFont());
			g2d.drawString(formattedStr, startx - adv - 5, ypos);
		}


		// draw an x-axis

		// plot x-axis

		int yPosXAxis = y0 + marginY + pixelsY + 10;
		int xPosXAxis = x0 + marginX;
		int nrPixelsX = width - (2 * marginX);
		g2d.drawLine(xPosXAxis, yPosXAxis, xPosXAxis + nrPixelsX, yPosXAxis);
		for (int binCluster = 0; binCluster < nrBinClusters; binCluster++) {

			int startX = x0 + marginX + (binCluster * widthPerBinCluster);

			int halfBinClusterWidth = widthPerBinCluster / 2;
			startX += halfBinClusterWidth - marginBetweenBinClusters;
			g2d.drawLine(startX, yPosXAxis, startX, yPosXAxis + 5);
		}

		if (binLabels != null) {
			g2d.setColor(theme.getDarkGrey());
			g2d.setFont(theme.getMediumFont());
			metrics = g2d.getFontMetrics(g2d.getFont());
			int fontheight = metrics.getHeight();

			int y = yPosXAxis + 5;
			for (int binCluster = 0; binCluster < binLabels.length; binCluster++) {
				int startX = x0 + marginX + (binCluster * widthPerBinCluster);

				int halfBinClusterWidth = widthPerBinCluster / 2;
				startX += halfBinClusterWidth - marginBetweenBinClusters;

				String str = binLabels[binCluster];
				int widthOfStr = metrics.stringWidth(str);
				drawRotate(g2d, startX + (fontheight / 2), y + widthOfStr + 10, -90, str);
			}
		}


	}

	private Pair<double[][], double[]> determineDensity(double[] doubles, Range yrange) {
		ArrayList<ArrayList<Double>> binobj = new ArrayList<ArrayList<Double>>();
		int nrBins = 100;
		for (int i = 0; i < nrBins; i++) {
			binobj.add(new ArrayList<Double>());
		}
		for (int i = 0; i < doubles.length; i++) {
			double d = doubles[i];
			double perc = yrange.getRelativePositionY(d);
			double bin = perc * nrBins;
			int ibin = (int) Math.ceil(bin);
			if (ibin < 0) {
				ibin = 0;
			} else if (ibin >= nrBins) {
				ibin = nrBins - 1;
			}
			binobj.get(ibin).add(d);
		}

		// make 100 bins
		double[] freqs = new double[nrBins];
		double[][] binnedData = new double[nrBins][];

		// determine number of values in each bin
		for (int i = 0; i < nrBins; i++) {
			binnedData[i] = Primitives.toPrimitiveArr(binobj.get(i).toArray(new Double[0]));
			double percOfTotalData = (double) binobj.get(i).size() / doubles.length;
			freqs[i] = percOfTotalData;
		}

		return new Pair<double[][], double[]>(binnedData, freqs);

	}

	public void drawRotate(Graphics2D g2d, double x, double y, int angle, String text) {
		g2d.translate((float) x, (float) y);
		g2d.rotate(Math.toRadians(angle));
		g2d.drawString(text, 0, 0);
		g2d.rotate(-Math.toRadians(angle));
		g2d.translate(-(float) x, -(float) y);
	}

	public void setPlotRange(Range plotRange) {
		this.plotRange = plotRange;
	}

	public void setData(double[][] data) {
		this.data = new double[1][data.length][];
		for (int i = 0; i < data.length; i++) {
			this.data[0][i] = data[i];
		}
	}

	private double[] collectStats(double[] dataset) {

		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		double firstquartile = 0;
		double thirdquartile = 0;
		double median = 0;

		double[] datacopy = new double[dataset.length];
		System.arraycopy(dataset, 0, datacopy, 0, dataset.length);

		double mean = Descriptives.mean(datacopy);
		double sd = Math.sqrt(Descriptives.variance(datacopy));

		boolean even = false;
		if (dataset.length % 2 == 0) {
			even = true;
		}

		if (datacopy.length < 2) {
			if (datacopy.length == 1) {
				return new double[]{
						datacopy[0], datacopy[0], datacopy[0], datacopy[0], datacopy[0], datacopy[0], 0
				};
			} else {
				return new double[]{
						0, 0, 0, 0, 0, 0, 0
				};
			}

		} else {
			Arrays.sort(datacopy);
			if (even) {
				int middle = (int) Math.ceil((double) dataset.length / 2);
				double mid1 = datacopy[middle - 1];
				double mid2 = datacopy[middle];
				median = (mid1 + mid2) / 2;
				int firstqpos = (int) Math.floor((double) dataset.length / 4);
				int thirdqpos = (int) Math.floor((double) dataset.length * .75);
				firstquartile = datacopy[firstqpos];
				thirdquartile = datacopy[thirdqpos];


			} else {
				int middle = (int) Math.floor((double) dataset.length / 2);
				int firstqpos = (int) Math.floor((double) dataset.length / 4);
				int thirdqpos = middle + firstqpos;
				median = datacopy[middle];
				firstquartile = datacopy[firstqpos];
				thirdquartile = datacopy[thirdqpos];
			}


			min = datacopy[0];
			max = datacopy[datacopy.length - 1];

			return new double[]{
					min, firstquartile, median, thirdquartile, max, mean, sd
			};
		}


	}

	public void setDrawDataPoints(boolean drawDataPoints) {
		this.drawDataPoints = drawDataPoints;
	}


	public void setUseMeanAndSd(boolean useMeanAndSd) {
		this.useMeanAndSd = useMeanAndSd;
	}

	public void setOutputIQRS(String outputIQRS) {
		this.outputIQRS = outputIQRS;
	}

	public void setBinLabels(String[] binLabels) {
		this.binLabels = binLabels;
	}

	public void setData(double[][][] bins) {
		data = bins;
	}

	public void useTukeysDefault(boolean b) {
		useTukeysDefault = b;
	}
}
