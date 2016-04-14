package deconvolution;

import java.io.File;
import java.io.IOException;
import java.util.Random;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

public class ExampleHistogram {
	public static void main(String[] args) {
		double[] value = new double[100];
		Random generator = new Random();
		for (int i=1; i < 10; i++) {
			value[i] = generator.nextDouble();
			System.out.println(value[i]);
		}
		
		int bins = 20;
		HistogramDataset dataset = new HistogramDataset();
		dataset.setType(HistogramType.RELATIVE_FREQUENCY);
		dataset.addSeries("Histogram",value,bins);
		String plotTitle = "Histogram"; 
		String xaxis = "number";
		String yaxis = "value"; 
		PlotOrientation orientation = PlotOrientation.VERTICAL; 
		boolean show = false; 
		boolean toolTips = false;
		boolean urls = false; 
		JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis, 
				dataset, orientation, show, toolTips, urls);
		int width = 500;
		int height = 300; 
		try {
			ChartUtilities.saveChartAsPNG(new File("histogram.PNG"), chart, width, height);
		} catch (IOException e) {}
	}
}
