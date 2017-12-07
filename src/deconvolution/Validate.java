package deconvolution;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;

/**
 * Validate the deconvolution results using a tab-separated text file with correlation values of cell-type specific QTLs
 * 
 */
public class Validate {
	private File validationFile;
	private List<DeconvolutionResult> deconvolutionResults = new ArrayList<DeconvolutionResult>();
	private int numberOfTopCorrelations = 100; 
	// topCorrelationColumnNames is used to retrieve the correct columns from the validation file
	private static Map<String, String> validationColumnNames;
	static
	{
		validationColumnNames = new HashMap<String, String>();
		validationColumnNames.put("CD4..T.Cells.Correlation", "CD4");
		validationColumnNames.put("CD8..T.Cells.Correlation", "CD8");
		validationColumnNames.put("Lymphoblastoid.cell.lines.Correlation", "Lymph");
		validationColumnNames.put("B.Cells.Correlation", "B");
		validationColumnNames.put("Monocytes.Correlation", "Mono");
		validationColumnNames.put("Neutrophils.Correlation", "Neut");
	}
	private static Map<String, ArrayList<String>> deconvolutionSignificant;
	static
	{
		deconvolutionSignificant = new HashMap<String, ArrayList<String>>();
		deconvolutionSignificant.put("Neut_Perc", new ArrayList<String>());
		deconvolutionSignificant.put("Lymph_Perc", new ArrayList<String>());
		deconvolutionSignificant.put("Eos_Perc", new ArrayList<String>());
		deconvolutionSignificant.put("Mono_Perc", new ArrayList<String>());
	}
	private static Map<String, ArrayList<String>> deconvolutionSignificantUnique;
	static
	{
		deconvolutionSignificantUnique = new HashMap<String, ArrayList<String>>();
		deconvolutionSignificantUnique.put("Neut_Perc_unique", new ArrayList<String>());
		deconvolutionSignificantUnique.put("Lymph_Perc_unique", new ArrayList<String>());
		deconvolutionSignificantUnique.put("Eos_Perc_unique", new ArrayList<String>());
		deconvolutionSignificantUnique.put("Mono_Perc_unique", new ArrayList<String>());
	}
	private static Map<String, ArrayList<String>> deconvolutionInsignificant;
	static
	{
		deconvolutionInsignificant = new HashMap<String, ArrayList<String>>();
		deconvolutionInsignificant.put("Neut_Perc", new ArrayList<String>());
		deconvolutionInsignificant.put("Lymph_Perc", new ArrayList<String>());
		deconvolutionInsignificant.put("Eos_Perc", new ArrayList<String>());
		deconvolutionInsignificant.put("Mono_Perc", new ArrayList<String>());
	}
	private static Map<String, ArrayList<String>> deconvolutionSignificantCorrected;
	static
	{
		deconvolutionSignificantCorrected = new HashMap<String, ArrayList<String>>();
		deconvolutionSignificantCorrected.put("Neut_Perc", new ArrayList<String>());
		deconvolutionSignificantCorrected.put("Lymph_Perc", new ArrayList<String>());
		deconvolutionSignificantCorrected.put("Eos_Perc", new ArrayList<String>());
		deconvolutionSignificantCorrected.put("Mono_Perc", new ArrayList<String>());
	}
	private static Map<String, ArrayList<String>> deconvolutionSignificantCorrectedUnique;
	static
	{
		deconvolutionSignificantCorrectedUnique = new HashMap<String, ArrayList<String>>();
		deconvolutionSignificantCorrectedUnique.put("Neut_Perc_unique", new ArrayList<String>());
		deconvolutionSignificantCorrectedUnique.put("Lymph_Perc_unique", new ArrayList<String>());
		deconvolutionSignificantCorrectedUnique.put("Eos_Perc_unique", new ArrayList<String>());
		deconvolutionSignificantCorrectedUnique.put("Mono_Perc_unique", new ArrayList<String>());
	}
	private static Map<String, ArrayList<String>> deconvolutionInsignificantCorrected;
	static
	{
		deconvolutionInsignificantCorrected = new HashMap<String, ArrayList<String>>();
		deconvolutionInsignificantCorrected.put("Neut_Perc", new ArrayList<String>());
		deconvolutionInsignificantCorrected.put("Lymph_Perc", new ArrayList<String>());
		deconvolutionInsignificantCorrected.put("Eos_Perc", new ArrayList<String>());
		deconvolutionInsignificantCorrected.put("Mono_Perc", new ArrayList<String>());
	}
	private Map<Integer, String> columnIndexes = new HashMap<Integer, String>();
	private Map<String, HashMap<String, Double>> correlationsPerQTL = new HashMap<String, HashMap<String, Double>>();
	private ArrayList<ArrayList<String>> validationData = new ArrayList<ArrayList<String>>();

	public Validate(){};
	public Validate( List<DeconvolutionResult> deconvolutionResults , String validationFile, String figuresOutfolder) throws IOException, IllegalAccessException{
		/*
		 * Do the validation
		 * 
		 * @param deconvolutionResults Results from deconvoluting QTLs
		 * @param validationFile Text file containing the correlations of cell-type specific QTLs
		 */
		this.validationFile = new File(validationFile);
		this.deconvolutionResults = deconvolutionResults;
		new File(figuresOutfolder).mkdirs();
		extractCorrelations();
		significantDeconResults();
		getEffectSizeSignificants(); 
	}

	private void extractCorrelations() throws IOException {
		/*
		 * Get the highest correlation values from the validation file
		 * 
		 * 1. Get the indexes of the columns that contain the correlation values per celltype
		 * 2. Per line in validation file, for each celltype column, extract the correlation value
		 */
		LineIterator validationFileIterator = FileUtils.lineIterator(this.validationFile, "UTF-8");
		ArrayList<String> header = new ArrayList<String>(Arrays.asList(validationFileIterator.next().split("\t")));
		// 1. Get the indexes of the columns that contain the correlation values per celltype
		for(int i = 0; i < header.size(); i++){
			String columnName = header.get(i);
			String shortCelltypeName = validationColumnNames.get(columnName);
			// check if columnName is an existing key, if so save index
			if (shortCelltypeName != null) {
				columnIndexes.put(i, shortCelltypeName);
			}
		}
		// 2. Per line in validation file, for each celltype column, extract the correlation value
		while(validationFileIterator.hasNext()){
			ArrayList<String> validateStringVector = new ArrayList<String>(Arrays.asList(validationFileIterator.next().split("\t")));
			// save the validation data per line in an array, so that when we have the indexes of the top correlations we can extract them from the file
			validationData.add(validateStringVector);
			for (int i = 0; i < validateStringVector.size(); i++){
				String shortCelltypeName = columnIndexes.get(i);
				String QTL = validateStringVector.get(0);
				if (shortCelltypeName != null) {
					String correlationString = validateStringVector.get(i); 
					if (correlationString.equals("-")){
						// have to keep the element in so that when we search for the top 100 highest values we still have the indeex
						//correlationString = "-333";
						continue;
					}
					if (correlationsPerQTL.get(QTL) == null){
						correlationsPerQTL.put(validateStringVector.get(0), new HashMap <String, Double>());
					}
					correlationsPerQTL.get(QTL).put(shortCelltypeName, Double.parseDouble(correlationString));
				}
			}
		}
	}

	private void significantDeconResults() throws IllegalAccessException{
		/*
		 * Get list of significant QTLs for exclusive and non-exclusive decon results
		 */
		for (DeconvolutionResult deconResult : deconvolutionResults){
			// significantCount, significantCelltype and significantQTL to find is it is significant for only 1 celltype
			int significantCount = 0;
			String qtl = deconResult.getQtlName();
			//if(qtl.contains("ENSG00000240720")){
				//System.out.println(qtl);
			//}
			String significantCelltype = null;
			for (String celltype : deconResult.getPvaluePerCelltype().keySet()){
				double pval = deconResult.getPvaluePerCelltype().get(celltype);
				if (pval < 0.05){
					significantCount++;
					significantCelltype = celltype;
					deconvolutionSignificant.get(celltype).add(qtl);
				}
				else{
					deconvolutionInsignificant.get(celltype).add(qtl);
				}
			}
			if(significantCount == 1){
				deconvolutionSignificantUnique.get(significantCelltype+"_unique").add(qtl);
			}
		}
	}

	/**
	 * For the decon results, calculate the effect size
	 */
	private void getEffectSizeSignificants() throws IllegalAccessException, IOException{
		for (String celltypeDecon : deconvolutionSignificant.keySet()){
			// too many nested hashmaps should be able to do this in a much simpler method
			HashMap<String, ArrayList<Double>> correlationPerSignificantDecon = new HashMap<String, ArrayList<Double>>();
			HashMap<String, ArrayList<Double>> correlationPerInsignificantDecon = new HashMap<String, ArrayList<Double>>();
			HashMap<String, ArrayList<Double>> correlationPerSignificantUniqueDecon = new HashMap<String, ArrayList<Double>>();
			HashMap<String, ArrayList<Double>> correlationPerSignificantCorrectedDecon = new HashMap<String, ArrayList<Double>>();
			HashMap<String, ArrayList<Double>> correlationPerInsignificantCorrectedDecon = new HashMap<String, ArrayList<Double>>();
			HashMap<String, ArrayList<Double>> correlationPerSignificantCorrectedUniqueDecon = new HashMap<String, ArrayList<Double>>();
			for (String significantQTL : deconvolutionSignificant.get(celltypeDecon)){
				if(!celltypeDecon.equals("Neut_Perc")){
					continue;
				}
					//System.out.printf(significantQTL+";");
				//System.out.println(celltypeDecon);
				// correlationsPerQTL contains for every QTL the correlation values per celltype
				if (correlationsPerQTL.get(significantQTL) != null) {
					System.out.println(significantQTL);
					for (String celltypeCorrelations : correlationsPerQTL.get(significantQTL).keySet()){
						if (correlationPerSignificantDecon.get(celltypeCorrelations) == null){
							correlationPerSignificantDecon.put(celltypeCorrelations, new ArrayList<Double>());
						}
						correlationPerSignificantDecon.get(celltypeCorrelations).add(correlationsPerQTL.get(significantQTL).get(celltypeCorrelations));
					}
				}
			}
			System.out.println(correlationPerSignificantDecon.get("Neut").size());
			System.exit(1);

			for (String insignificantQTL : deconvolutionInsignificant.get(celltypeDecon)){
				if (correlationsPerQTL.get(insignificantQTL) != null) {
					for (String celltypeCorrelations : correlationsPerQTL.get(insignificantQTL).keySet()){
						if (correlationPerInsignificantDecon.get(celltypeCorrelations) == null){
							correlationPerInsignificantDecon.put(celltypeCorrelations, new ArrayList<Double>());
						}
						correlationPerInsignificantDecon.get(celltypeCorrelations).add(correlationsPerQTL.get(insignificantQTL).get(celltypeCorrelations));
					}
				}
			}

			for (String significantUniqueQTL : deconvolutionSignificantUnique.get(celltypeDecon+"_unique")){
				if (correlationsPerQTL.get(significantUniqueQTL) != null) {
					for (String celltypeCorrelations : correlationsPerQTL.get(significantUniqueQTL).keySet()){
						if (correlationPerSignificantUniqueDecon.get(celltypeCorrelations) == null){
							correlationPerSignificantUniqueDecon.put(celltypeCorrelations, new ArrayList<Double>());
						}
						correlationPerSignificantUniqueDecon.get(celltypeCorrelations).add(correlationsPerQTL.get(significantUniqueQTL).get(celltypeCorrelations));
					}
				}
			}

			for (String significantCorrectedQTL : deconvolutionSignificantCorrected.get(celltypeDecon)){
				if (correlationsPerQTL.get(significantCorrectedQTL) != null) {
					for (String celltypeCorrelations : correlationsPerQTL.get(significantCorrectedQTL).keySet()){
						if (correlationPerSignificantCorrectedDecon.get(celltypeCorrelations) == null){
							correlationPerSignificantCorrectedDecon.put(celltypeCorrelations, new ArrayList<Double>());
						}
						correlationPerSignificantCorrectedDecon.get(celltypeCorrelations).add(correlationsPerQTL.get(significantCorrectedQTL).get(celltypeCorrelations));
					}
				}
			}

			for (String insignificantCorrectedQTL : deconvolutionInsignificantCorrected.get(celltypeDecon)){
				if (correlationsPerQTL.get(insignificantCorrectedQTL) != null) {
					for (String celltypeCorrelations : correlationsPerQTL.get(insignificantCorrectedQTL).keySet()){
						if (correlationPerInsignificantCorrectedDecon.get(celltypeCorrelations) == null){
							correlationPerInsignificantCorrectedDecon.put(celltypeCorrelations, new ArrayList<Double>());
						}
						correlationPerInsignificantCorrectedDecon.get(celltypeCorrelations).add(correlationsPerQTL.get(insignificantCorrectedQTL).get(celltypeCorrelations));
					}
				}
			}

			for (String significantCorrectedUniqueQTL : deconvolutionSignificantCorrectedUnique.get(celltypeDecon+"_unique")){
				if (correlationsPerQTL.get(significantCorrectedUniqueQTL) != null) {
					for (String celltypeCorrelations : correlationsPerQTL.get(significantCorrectedUniqueQTL).keySet()){
						if (correlationPerSignificantCorrectedUniqueDecon.get(celltypeCorrelations) == null){
							correlationPerSignificantCorrectedUniqueDecon.put(celltypeCorrelations, new ArrayList<Double>());
						}
						correlationPerSignificantCorrectedUniqueDecon.get(celltypeCorrelations).add(correlationsPerQTL.get(significantCorrectedUniqueQTL).get(celltypeCorrelations));
					}
				}
			}	
			//System.out.println(celltypeDecon);
			String celltype = "Neut";
			System.out.println(celltype);
			System.out.println(correlationPerSignificantDecon.get(celltype).size());
			System.out.println(correlationPerInsignificantDecon.get(celltype).size());
			System.out.println(correlationPerSignificantUniqueDecon.get(celltype).size());
			System.exit(1);
			//System.out.println("-------------");
			//System.out.println("===================");
		}
	}

	// getter and setter functions

	public int getNumberOfTopCorrelations(){
		return this.numberOfTopCorrelations;
	}

	public void setNumberOfTopCorrelations(int numberOfTopCorrelations){
		this.numberOfTopCorrelations = numberOfTopCorrelations;
	}

	public Map<String, String> getCorrelationColumnNames(){
		return Validate.validationColumnNames;
	}

	public void setNumberOfTopCorrelations(Map<String, String> topCorrelationColumnNames){
		Validate.validationColumnNames = topCorrelationColumnNames;
	}

	public Map<String, ArrayList<String>> getDeconvolutionSignificant(){
		return Validate.deconvolutionSignificant;
	}

	public void setDeconvolutionSignificant(Map<String, ArrayList<String>> deconvolutionCelltypeNames){
		Validate.deconvolutionSignificant = deconvolutionCelltypeNames;
	}

	public Map<String, ArrayList<String>> getDeconvolutionSignificantCorrected(){
		return Validate.deconvolutionSignificantCorrected;
	}

	public void setDeconvolutionSignificantCorrected(Map<String, ArrayList<String>> deconvolutionCelltypeNames){
		Validate.deconvolutionSignificantCorrected = deconvolutionCelltypeNames;
	}
}





