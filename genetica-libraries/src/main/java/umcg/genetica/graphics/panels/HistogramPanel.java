package umcg.genetica.graphics.panels;

import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.graphics.themes.Theme;

import java.awt.*;
import java.text.DecimalFormat;

/**
 * Created by hwestra on 7/16/15.
 */
public class HistogramPanel extends Panel {

    private double[][] histogram;
    private String[] datasetLabels;
    private String xAxisLabel;
    private String yAxisLabel;
    private Range dataRange;
    private Theme theme = new DefaultTheme();
    private LOG xAxisLog = LOG.NONE;
    private LOG yAxisLog = LOG.NONE;
    private String[] binLabels;
    private Double threshold;
    private String thresholdvalue;

    public void setMarginBetweenBinClusters(int marginBetweenBinClusters) {
        this.marginBetweenBinClusters = marginBetweenBinClusters;
    }

    private int marginBetweenBinClusters;


    public HistogramPanel(int nrRows, int nrCols) {
        super(nrRows, nrCols);
        marginBetweenBinClusters = 5;
    }

    public void setData(int[] dataset) {
        double[][] tmpData = new double[1][dataset.length];

        for (int j = 0; j < dataset.length; j++) {
            tmpData[0][j] = dataset[j];
        }
        this.histogram = tmpData;
    }

    public void setAxisLog(LOG xAxis, LOG yAxis) {
        this.xAxisLog = xAxis;
        this.yAxisLog = yAxis;

    }

    public void setThreshold(double threshold, String thresholdvalue) {
        this.threshold = threshold;
        this.thresholdvalue = thresholdvalue;
    }

    public enum PLOTTYPE {
        STACKED,
        CLUSTERED
    }

    public enum DATASETPLOTTYPE {
        BAR,
        POLY
    }

    public enum LOG {
        NONE,
        TWO,
        TEN
    }

    DATASETPLOTTYPE[] datasetplottypes = null;
    PLOTTYPE plottype = PLOTTYPE.CLUSTERED;

    public void setData(int[][] histogram) {

        double[][] tmpData = new double[histogram.length][histogram[0].length];
        for (int i = 0; i < histogram.length; i++) {
            for (int j = 0; j < histogram[i].length; j++) {
                tmpData[i][j] = histogram[i][j];
            }
        }
        this.histogram = tmpData;
    }

    public void setDatasetPlotTypes(DATASETPLOTTYPE[] types) {
        this.datasetplottypes = types;
    }


    public void setData(double[] dataset) {
        double[][] tmpData = new double[1][dataset.length];

        for (int j = 0; j < dataset.length; j++) {
            tmpData[0][j] = dataset[j];
        }
        this.histogram = tmpData;
    }

    public void setData(double[][] histogram) {
        this.histogram = histogram;
    }

    public void setRange(Range range) {
        this.dataRange = range;
    }

    public void setAxisLabels(String xAxis, String yAxis) {
        this.xAxisLabel = xAxis;
        this.yAxisLabel = yAxis;
    }

    public void setBinLabels(String[] labels) {
        this.binLabels = labels;
    }

    public void setDatasetLabels(String[] datasetLabels) {
        this.datasetLabels = datasetLabels;
    }

    public void draw(DefaultGraphics g) {

        Graphics2D g2d = g.getG2d();

        if (dataRange == null) {
            System.out.println("Getting range from data");
            dataRange = new Range(histogram);
            System.out.println(dataRange);
        }

        System.out.println("Range: " + dataRange);
        Range plotRange = new Range(dataRange.getMinX(), dataRange.getMinY(), dataRange.getMaxX(), dataRange.getMaxY());
        plotRange.round();


        int nrDatasets = histogram.length;
        int nrBinClusters = histogram[0].length;

        System.out.println(nrDatasets + " datasets, " + nrBinClusters + " bin clusters");

        int nrPixelsMaxX = width - (2 * marginX);
        int nrPixelsMaxY = height - (2 * marginY);

//		if (plottype == PLOTTYPE.CLUSTERED) {

        int marginBetweenDatasets = 2;

        // calculate bar width
        int widthPerBinCluster = ((nrPixelsMaxX) / nrBinClusters);
        int widthPerBinClusterBin = widthPerBinCluster - marginBetweenBinClusters;
        int widthPerDataset = (widthPerBinClusterBin / nrDatasets);
        int widthPerDatasetBin = widthPerDataset - marginBetweenDatasets;

        System.out.println(plotRange);
        // plot bins
        for (int binCluster = 0; binCluster < nrBinClusters; binCluster++) {
            for (int dataset = 0; dataset < nrDatasets; dataset++) {
                g2d.setColor(theme.getColor(dataset));
                int startX = x0 + marginX + (binCluster * widthPerBinCluster) + (dataset * widthPerDataset);

                double percY1 = plotRange.getRelativePositionY(histogram[dataset][binCluster]);
                int yHeight1 = (int) Math.ceil(percY1 * nrPixelsMaxY);
                int yPos1 = y0 + marginY + nrPixelsMaxY - yHeight1;
                g2d.fillRect(startX, yPos1, widthPerDatasetBin, yHeight1);
                System.out.println(binCluster + "\t" + histogram[dataset][binCluster] + "\t" + percY1);
            }
        }
//		}


        // draw axes
        // find out where axes intersect
        // leave the x/y intersection for now
//			double yx0Perc = plotRange.getRelativePositionX(0d);
//			double xy0Perc = plotRange.getRelativePositionY(0d);


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

        // draw title
        if (title != null) {
            int titlePosX = x0 + marginX;
            int titlePosY = y0 + marginY - theme.getLargeFont().getSize() - 5;
            g2d.setColor(theme.getDarkGrey());
            g2d.drawString(title, titlePosX, titlePosY);

        }

        // draw legend
        if (datasetLabels != null) {

        }


        // plot threshold
        if (threshold != null) {
            double perc = dataRange.getRelativePositionX(threshold);
            int startX = x0 + marginX + (int) Math.floor(nrPixelsX * perc);

            g2d.setColor(Color.red);
            pixelsY = height - (2 * marginY);
            g2d.drawLine(startX, yPosYAxis, startX, yPosYAxis + pixelsY);
            drawString(g2d, thresholdvalue, startX + 5, yPosYAxis);
        }
    }

    private void drawString(Graphics g, String text, int x, int y) {
        int n = 0;
        for (String line : text.split("\n")) {
            if (n == 0) {
                g.drawString(line, x, y);
            } else {
                g.drawString(line, x, y += g.getFontMetrics().getHeight());
            }
            n++;
        }

    }


}
