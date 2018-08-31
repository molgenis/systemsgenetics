package main.java.decon_eQTL_simple;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;

public class ExpressionData {
	private ArrayList<String> sampleNames;
	private HashMap<String, double[]> geneExpression = new HashMap<String, double[]>();
	
	public ExpressionData(){}
	public ExpressionData(String expressionFile) throws IOException{
		LineIterator expressionIterator = FileUtils.lineIterator(new File(expressionFile), "UTF-8");
		this.sampleNames = new ArrayList<String>( Arrays.asList(expressionIterator.next().split("\t")) );
		this.sampleNames.removeAll(Arrays.asList("", null));
		int rowNumber = 0;
		while (expressionIterator.hasNext()) {
			rowNumber++;
			if(rowNumber % 5000 == 0){
				DeconvolutionLogger.log.info(String.format("Processed %d lines", rowNumber));
			}
			String[] expressionStringVector =  expressionIterator.next().split("\t");
			String geneName = expressionStringVector[0];
			if(this.sampleNames.size() != expressionStringVector.length-1){
				DeconvolutionLogger.log.info(String.format("Expression table %s does not have the same number of columns as there are in the header at row %d",expressionFile,rowNumber));
				DeconvolutionLogger.log.info(String.format("Number of header columns: %d",this.sampleNames.size()));
				DeconvolutionLogger.log.info(String.format("Number of columns at row %d: %d", rowNumber, expressionStringVector.length-1));
				throw new RuntimeException(String.format("Expressione table does not have the same number of columns as there are celltypes at row %d",rowNumber));
			}
			
			
			double[] expressionValues = null;
			try{
				expressionValues = Utils.StringVectorToDoubleArrayList(expressionStringVector, 1);
			}catch(NumberFormatException e){
				DeconvolutionLogger.log.warning(String.format("Gene %s contains expression values that can not be converted to Double, SKIPPING!", geneName));
			}
			
			geneExpression.put(geneName, expressionValues);
		}
	}

	public ArrayList<String> getSampleNames(){
		return this.sampleNames;
	}

	public HashMap<String, double[]>  getGeneExpression() throws IllegalAccessException{
		if(this.geneExpression == null){
			throw new IllegalAccessException("geneExpression not set ExpressionData");
		}
		return this.geneExpression;
	}

}
