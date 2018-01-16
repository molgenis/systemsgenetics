package deconvolution;

import java.io.IOException;
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
	@SuppressWarnings("unchecked")
	public CellCount( String cellCountFile) throws IOException{
		// the cell type names are the first row of cellcount file, extract for
		// later printing
		// is now saved as table of strings
		Object[] cellCountData = Utils.readTabDelimitedColumns(cellCountFile);
		samplenames = (ArrayList<String>) cellCountData[0];
		cellcountTable = (List<List<String>>) cellCountData[1];
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
		cellcountPercentages = new double[numberOfSamples][numberOfCelltypes];
		for (int j = 0; j <= numberOfSamples-1; j++) {
			for (int i = 0; i < numberOfCelltypes; i++) {
				cellcountPercentages[j][i] =  Double.parseDouble(cellcountTable.get(i).get(j+1));
			}
		}
	}

	public void emptyCellcountPercentages(){
		cellcountPercentages = null;
	}
	
	public List<String> getSampleNames(){
		return samplenames;
	}
	
	public List<List<String>> getCellcountTable(){
		return(cellcountTable);
	}	
	public List<String> getAllCelltypes(){
		return(celltypes);
	}	
	public String getCelltype(int index){
		return(celltypes.get(index));
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
}
