package deconvolution;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import deconvolution.Utils;

public class CellCount {

	private static List<String> celltypes = new ArrayList<String> ();
	private static List<String> samplenames = new ArrayList<String> ();
	private static double[][] cellcountPercentages;
	private static List<List<String>> cellcountTable;
	private static int numberOfCelltypes;
	private static int numberOfSamples;
	public CellCount() {};
	/**
	 * Read in cellcount file 
	 * 
	 * @param cellcountFile File with cellcount percentages, 
	 * 		   with columns = celltype, rows is samples (includes column headers and row names)   
	 */
	public CellCount( String cellCountFile) throws IOException{
		// the cell type names are the first row of cellcount file, extract for
		// later printing
		// is now saved as table of strings, should be changed to table of doubles so we only have to convert them ones
		cellcountTable = Utils.readTabDelimitedColumns(cellCountFile);
		numberOfCelltypes = cellcountTable.size();
		DeconvolutionLogger.log.info(String.format("Celltypes to use:"));
		for(int i = 0; i < numberOfCelltypes; i++){
			celltypes.add(cellcountTable.get(i).get(0));
			DeconvolutionLogger.log.info(cellcountTable.get(i).get(0));
		}
		
		DeconvolutionLogger.log.info(String.format("Number of celltypes: %d", numberOfCelltypes));
		
		// minus one because the size includes the celltype header
		numberOfSamples = cellcountTable.get(0).size()-1;
		DeconvolutionLogger.log.info(String.format("Number of samples: %d", numberOfSamples));
		for(String sampleName : cellcountTable.get(0)){
			samplenames.add(sampleName);
		}
		cellcountPercentages = new double[numberOfSamples][numberOfCelltypes];
		for (int j = 0; j <= numberOfSamples-1; j++) {
			for (int i = 0; i < numberOfCelltypes; i++) {
				cellcountPercentages[j][i] =  Double.parseDouble(cellcountTable.get(i).get(j+1));
			}
		}
	}

	public List<List<String>> getCellcountTable(){
		return(cellcountTable);
	}	
	public List<String> getCelltypes(){
		return(celltypes);
	}	
	public int getNumberOfCelltypes(){
		return(numberOfCelltypes);
	}	
	public int getNumberOfSamples(){
		return(numberOfSamples);
	}
	public double[][] getCellcountPercentages(){
		return(cellcountPercentages);
	}

	public void makeCellcountsRelative(String outputFolder) throws IOException{
		for (int i = 0; i < numberOfCelltypes; i++) {
			double average = 0;
			for (int j = 0; j <= numberOfSamples-1; j++) {
				average +=  cellcountPercentages[j][i];
			}
			average = average/cellcountPercentages.length;
			for (int j = 0; j <= numberOfSamples-1; j++) {
				cellcountPercentages[j][i] =  cellcountPercentages[j][i]/average;
			}
		}
		/*
		 * Write the new cellcounts to output folder so that it can be used after running
		 */
		List<String> newCCoutput = new ArrayList<String>();
		String header = new String(); 
		for (String celltype : getCelltypes()){
			header += celltype+"\t";
		}
		header = header.substring(0, header.length() - 2);
		newCCoutput.add(header);
		for (int j = 0; j <= getNumberOfSamples()-1; j++) {
			String line = "";
			for (int i = 0; i < getNumberOfCelltypes(); i++) {
				line += getCellcountPercentages()[j][i]+"\t";
			}
			line = line.substring(0, line.length() - 2);
			newCCoutput.add(line);
		}
		Path file = Paths.get(outputFolder+"normalizedCellcounts.csv");
		Files.write(file, newCCoutput, Charset.forName("UTF-8"));
		DeconvolutionLogger.log.info(String.format("New cellcount (cause -cc option used) written to  %s", file));
	}
}
