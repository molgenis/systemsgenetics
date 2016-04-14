package deconvolution;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.BoxAndWhiskerToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

public class Plots {
	public static void boxPlot(OLSMultipleLinearRegression regr, InteractionModel model, String filename) throws IOException, IllegalAccessException{
		/* Make a boxplot of the beta * values of linear regression model
		 * 
		 */
		DefaultBoxAndWhiskerCategoryDataset dataset = new DefaultBoxAndWhiskerCategoryDataset();
		int x = 0;
		// loop through each parameter (b0,b1,b2,b3 etc... <- b0 not there when intersect is removed)
		for (Double parameter : regr.estimateRegressionParameters()){
			List<Double> betaTimesValue = new ArrayList<Double>();
			for (int i = 0; i < model.GetObservedValues().length; i++){
				// log modulus transformation, so that negative numbers are preserved
				Boolean negative = false;
				double value = model.GetObservedValues()[i][x]*parameter;
				// check if original value is negative
				if(value < 0){
					negative = true;
				}
				double logValue = Math.log10(Math.abs(value)+1);
				if (negative){
					// if original value was negative, put sign b ack
					logValue *= -1;
				}
				betaTimesValue.add(logValue);
			}
				dataset.add(betaTimesValue,"B"+Integer.toString(x+1)+" * "+ model.GetIndependentVariables().get(x), "Parameter * explanatory variables");
				x++;
		}
		final CategoryAxis xAxis = new CategoryAxis("Explanatory variable");
        final NumberAxis yAxis = new NumberAxis("Log10 Modulus");
        yAxis.setAutoRangeIncludesZero(false);
        final BoxAndWhiskerRenderer renderer = new BoxAndWhiskerRenderer();
        renderer.setFillBox(true);
        renderer.setSeriesToolTipGenerator(1, new BoxAndWhiskerToolTipGenerator());
        final CategoryPlot plot = new CategoryPlot(dataset, xAxis, yAxis, renderer);

        final JFreeChart chart = new JFreeChart(
            "Beta * explanatory variable",
            plot
        );
        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(3000,1800));
        ChartUtilities.saveChartAsPNG(new File(filename), chart, 1000, 600);
       // setContentPane(chartPanel);
	}
	
	public static void drawHistogram(PermutationResult permutationResults, int numPermutations, String outfolder) throws IOException {
		// This assumes that all deconvolution reslts in deconResults have the same cell types
		for (String celltype : permutationResults.getCelltypes()){
				double[] pvalueVector = new double[permutationResults.getPvalues(celltype).size()];
				for (int j = 0; j < pvalueVector.length; j++) {
					pvalueVector[j] = permutationResults.getPvalues(celltype).get(j);
				}
				HistogramDataset dataset = new HistogramDataset();
				dataset.setType(HistogramType.FREQUENCY);
				// 20 = amount of bins
				dataset.addSeries("Histogram", pvalueVector, 10);
				String plotTitle = "Permuation pvalue distribution ("+numPermutations+" permuations): \n"+celltype;
				String xaxis = "pvalue";
				String yaxis = "count";
				PlotOrientation orientation = PlotOrientation.VERTICAL;
				boolean show = false;
				boolean toolTips = false;
				boolean urls = false;
				JFreeChart chart = ChartFactory.createHistogram(plotTitle, xaxis, yaxis, dataset, orientation, show,
						toolTips, urls);
				int width = 500;
				int height = 300;
	
				ChartUtilities.saveChartAsPNG(new File(outfolder+'/'+celltype+"_pvalueDistribution.PNG"), chart, width, height);	
	
				System.out.printf("Permutation distrbution written to: %s/%s_pvalueDistribution.PNG\n",outfolder, celltype);
			}
	}
}
