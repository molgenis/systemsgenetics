package deconvolution;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
	private static double betaTimesValueTempFunc (double[] estimatedRegressionParameters, int[] celltypeIndices, InteractionModel model) throws IllegalAccessException{
		List<Double> betaTimesInteraction = new ArrayList<Double>();
		double logValue = 0;
		for (int i = 0; i < model.GetObservedValues().length; i++){
			double beta = estimatedRegressionParameters[celltypeIndices[0]];
			double betaInteraction = 0;
			if(celltypeIndices.length == 2){
				//betaInteraction = estimatedRegressionParameters[celltypeIndices[1]];
			}
			beta += betaInteraction;
			betaTimesInteraction.add(beta);
			//betas.add(betaTimesValueTempFunc(estimatedRegressionParameters, celltypeIndices, model, i));

			// Calculate Beta1 * celltype
			double valueCelltype = estimatedRegressionParameters[celltypeIndices[0]] * model.GetObservedValues()[i][celltypeIndices[0]];
			// instantiate valueBetaInteraction here so that it is 0  
			double valueBetaInteraction = 0;
			/* If there is an interaciton term included for this celltype
			 * Calculate Beta2 * celltype:genotype */
			if(celltypeIndices.length == 2){
				valueBetaInteraction = estimatedRegressionParameters[celltypeIndices[1]] * model.GetObservedValues()[i][celltypeIndices[1]];
			}
			// Calculate the beta*cell% + beta*cell%*GT. If there is not cell%:GT for current celltype, valueBetaInteraction = 0
			// TODO: More informative variable name, need to find official terminology
			valueCelltype += valueBetaInteraction;
			logValue = Statistics.logmodulus(valueCelltype);
		}
		return(logValue);
	}
	public static void boxPlot(OLSMultipleLinearRegression regression, InteractionModel model, String filename) throws IOException, IllegalAccessException{
		/* Make a boxplot of the beta * values of linear regression model
		 * 
		 */
		DefaultBoxAndWhiskerCategoryDataset dataset = new DefaultBoxAndWhiskerCategoryDataset();
		// get each parameter (b0,b1,b2,b3 etc... <- b0 not there when intersect is removed)
		double[] estimatedRegressionParameters = regression.estimateRegressionParameters();
		double[] parameterErrors = regression.estimateRegressionParametersStandardErrors();
		List<Double> betasFull = new ArrayList<Double>();
		//for(double a : estimatedRegressionParameters){
		//	System.out.println(a);
		//}
		//System.exit(0);
		String xAxisLabel = "";
		// celltypeIndices = [[0, 4], [1,5], [2,6], [3,7]] where first index is cellcount% and second index is cellcount%:GT
		if(estimatedRegressionParameters.length == 8){
			for (int[] celltypeIndices : model.GetCelltypeVariablesIndex()){
				List<Double> betas = new ArrayList<Double>();
				//	betaTimesValueTempFunc
				for (int i = 0; i < model.GetObservedValues().length; i++){
					double beta = estimatedRegressionParameters[celltypeIndices[0]] + parameterErrors[celltypeIndices[0]];
					double betaInteraction = 0;
					if(celltypeIndices.length == 2){
						betaInteraction = estimatedRegressionParameters[celltypeIndices[1]] + parameterErrors[celltypeIndices[0]];
					}
					beta += 2*betaInteraction;
					beta *= model.GetObservedValues()[i][celltypeIndices[0]];
					betas.add(beta);
					betasFull.add(beta);
				}
			
				//String xAxisLabel = "B" + Integer.toString(celltypeIndices[0]+1) + " * " + model.GetIndependentVariableNames().get(celltypeIndices[0]);
				xAxisLabel = "(B" + Integer.toString(celltypeIndices[0]+1) + " + " +  "2*B" + Integer.toString(celltypeIndices[1]+1) + ") * cellcount";
				if(celltypeIndices.length == 2){ 
					//xAxisLabel += " + B" + Integer.toString(celltypeIndices[1]+1);
					//xAxisLabel += " + B" + Integer.toString(celltypeIndices[1]+1) + " * " + model.GetIndependentVariableNames().get(celltypeIndices[1]);
				}
				
				
				xAxisLabel += " of "+model.GetIndependentVariableNames().get(celltypeIndices[0]);
				dataset.add(betas, xAxisLabel, "Beta of celltype");// + beta of interactionterm");
			}

			// only plot if for at least 1 sample Beta1*celltype + beta2 * celltype:GT is negative
			if(Collections.min(betasFull) < 0){
				final CategoryAxis xAxis = new CategoryAxis("Beta");
				final NumberAxis yAxis = new NumberAxis("Value (lowest = "+Double.toString(Collections.min(betasFull))+")");
				yAxis.setAutoRangeIncludesZero(false);
				final BoxAndWhiskerRenderer renderer = new BoxAndWhiskerRenderer();
				renderer.setFillBox(true);
				renderer.setSeriesToolTipGenerator(1, new BoxAndWhiskerToolTipGenerator());
				renderer.setMeanVisible(false);

				final CategoryPlot plot = new CategoryPlot(dataset, xAxis, yAxis, renderer);

				final JFreeChart chart = new JFreeChart(
						"Beta of celltype",// + beta of interaction term",
						plot
						);
				final ChartPanel chartPanel = new ChartPanel(chart);
				chartPanel.setPreferredSize(new java.awt.Dimension(3000,1800));
				ChartUtilities.saveChartAsPNG(new File(filename), chart, 1000, 600);
				//System.out.printf("Plot written to: %s\n", filename);
			}
		}
		// setContentPane(chartPanel);
	}

	public static void drawHistogram(PermutationResult permutationResults, String outfolder) throws IOException {
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
			String plotTitle = "Permuation pvalue distribution ("+permutationResults.getNumberOfPermutations()+" permuations): \n"+celltype;
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

			ChartUtilities.saveChartAsPNG(new File(outfolder+'/'+celltype+"_" +
					Integer.toString(permutationResults.getNumberOfPermutations()) +
					"_pvalueDistribution.PNG"), chart, width, height);	

			System.out.printf("Permutation distrbution written to: %s/%s_%d_pvalueDistribution.PNG\n",outfolder, celltype,permutationResults.getNumberOfPermutations());
		}
	}
}
