package decon_eQTL;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;

public class GenotypeData {
	private ArrayList<String> sampleNames;
	private HashMap<String, double[]> genotypes = new HashMap<String, double[]>();;
	
	public GenotypeData(){};
	public GenotypeData(String genotypeFile) throws IOException{
		LineIterator genotypeIterator = FileUtils.lineIterator(new File(genotypeFile), "UTF-8");
		this.sampleNames = new ArrayList<String>( Arrays.asList(genotypeIterator.next().split("\t")) );
		this.sampleNames.removeAll(Arrays.asList("", null));
		int rowNumber = 0;
		while (genotypeIterator.hasNext()) {
			++rowNumber;
			if(rowNumber % 5000 == 0){
				DeconvolutionLogger.log.info(String.format("Processed %d lines", rowNumber));
			}
			String[] genotypeStringVector = genotypeIterator.next().split("\t");
			String snpName = genotypeStringVector[0];
			if(this.sampleNames.size() != genotypeStringVector.length-1){
				DeconvolutionLogger.log.info(String.format("Genotype table %s does not have the same number of columns as there are in the header at row %d",genotypeFile,rowNumber));
				DeconvolutionLogger.log.info(String.format("Number of header columns: %d",this.sampleNames.size()));
				DeconvolutionLogger.log.info(String.format("Number of columns at row %d: %d", rowNumber, genotypeStringVector.length-1));
				throw new RuntimeException(String.format("Expressione table does not have the same number of columns as there are celltypes at row %d",rowNumber));
			}
			
			double[] genotypeValues = null;
			try{
				genotypeValues = Utils.StringVectorToDoubleArrayList(genotypeStringVector, 1);
			}catch(NumberFormatException e){
				DeconvolutionLogger.log.warning(String.format("SNP %s contains genotype values that can not be converted to Double, SKIPPING!", snpName));
			}
			
			genotypes.put(snpName, genotypeValues);
		}
	}
	public void setSampleNames(ArrayList<String> sampleNames){
		this.sampleNames = sampleNames;
	}
	public ArrayList<String> getSampleNames(){
		return this.sampleNames;
	}
	public void setGenotypes(HashMap<String, double[]> genotypes){
		this.genotypes = genotypes;
	}
	public HashMap<String, double[]>  getGenotypes() throws IllegalAccessException{
		if(this.genotypes == null){
			throw new IllegalAccessException("genotypes not set GenotypesData");
		}
		return this.genotypes;
	}
}

