package decon_eQTL;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import decon_eQTL.Utils;

public class CellCount {

	private List<String> cellTypes = new ArrayList<String> ();
	private List<String> sampleNames = new ArrayList<String> ();
	private double[][] cellCountPercentages;
	private List<List<String>> cellCountTable;
	private int numberOfCelltypes;
	private int numberOfSamples;
	public CellCount() {};
	/**
	 * Read in cell count file 
	 * 
	 * @param cellCountFile File with cell count percentages, 
	 * 		   with columns = cell type, rows is samples (includes column headers and row names)
	 *
	 * @throws IOException	If cell count file can not be read
	 */
	@SuppressWarnings("unchecked")
	public CellCount( String cellCountFile) throws IOException{
		// the cell type names are the first row of cell count file, extract for
		// later printing is now saved as table of strings
		Object[] cellCountData = Utils.readTabDelimitedColumns(cellCountFile);
		sampleNames = (ArrayList<String>) cellCountData[0];
		cellCountTable = (List<List<String>>) cellCountData[1];
		numberOfCelltypes = cellCountTable.size();
		DeconvolutionLogger.log.info(String.format("Cell types to use:"));
		for(int i = 0; i < numberOfCelltypes; i++){
			cellTypes.add(cellCountTable.get(i).get(0));
			DeconvolutionLogger.log.info(cellCountTable.get(i).get(0));
		}
		DeconvolutionLogger.log.info(String.format("Number of cell types: %d", numberOfCelltypes));
		
		// minus one because the size includes the cell type header
		numberOfSamples = cellCountTable.get(0).size()-1;
		DeconvolutionLogger.log.info(String.format("Number of samples: %d", numberOfSamples));
		cellCountPercentages = new double[numberOfSamples][numberOfCelltypes];
		for (int j = 0; j <= numberOfSamples-1; j++) {
			for (int i = 0; i < numberOfCelltypes; i++) {
				cellCountPercentages[j][i] =  Double.parseDouble(cellCountTable.get(i).get(j+1));
			}
		}
	}

	public void emptyCellCountPercentages(){
		cellCountPercentages = null;
	}
	
	public List<String> getSampleNames(){
		return sampleNames;
	}
	
	public List<List<String>> getCellCountTable(){
		return cellCountTable;
	}	
	public List<String> getAllCelltypes(){
		return cellTypes;
	}	
	public String getCelltype(int index){
		return cellTypes.get(index);
	}	
	public int getNumberOfCelltypes(){
		return numberOfCelltypes;
	}	
	public int getNumberOfSamples(){
		return numberOfSamples;
	}
	public double[][] getCellCountPercentages(){
		return cellCountPercentages;
	}
}
