package decon_eQTL;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;

public class Utils {
	/**
	 * Add a value to one of the arrays in a 2D array
	 * 
	 * @param allColumns 2D array to add a value to
	 * 
	 * @param arrayIndex Index of the array to add a value to
	 * 
	 * @param value Value to add to the array
	 * 
	 * First try to add an element from row to the *i*th list, if
	 * the *i*th list does not exist yet
	 * catch the IndexOutOfBoundsException, make a new list and it
	 * at the *i*th position of the 2D array allColumns
	 */
	private static List<List<String>> addSingleValueTo2DArray(List<List<String>> allColumns, int arrayIndex, String value){

		try {
			allColumns.get(arrayIndex).add(value);
		} catch (IndexOutOfBoundsException e) {
			List<String> newColumn = new ArrayList<String>();
			newColumn.add(value);
			allColumns.add(newColumn);	
		}
		return allColumns;
	}
	
	
	/**
	 * Reads tab delimited file and returns them as list of list, with [x] =
	 * colummn and [x][y] is value in column. Needed for reading counts
	 * file, as there the rows are the samples, as opposed to expression and
	 * genotype file where the columns are the samples. Needs to be read in
	 * memory, so minimal memory requirement is larger than the size of the
	 * counts file.
	 * 
	 * 
	 * @param filepath The path to a tab delimited file to read
	 * 
	 * @return A 2D array with each array being one column from filepath except first column
	 * 		   and a 1D array with the first column (without header)
	 * 
	 * @throws IOException	If file at filepath can not be read
	 */
	public static Object[] readTabDelimitedColumns(String filepath) throws IOException {
		List<List<String>> allColumns = new ArrayList<List<String>>();
		// parses file on tabs
		CSVParser parser = new CSVParser(new FileReader(filepath), CSVFormat.newFormat('\t'));
		Boolean header = true;
		int rowNumber = 0;
		int columnIndexHeader = 0;
		List<String> firstColumn = new ArrayList<String>();
		for (CSVRecord row : parser) {
			rowNumber++;
			// starts at 1 because 1st element of column is the samplename, unless its the header row
			int columnStart = 1;
			if(header){
				columnStart = 0;
			}
			for (int columnIndex = columnStart; columnIndex < row.size(); columnIndex++) {
				// header can start from 0 if it is R styled, so check if element 0 has a value
				// R style is e.g.
				// colNameA	colNameB
				// rowNameA	AAValue	AAvalue
				// rownameB ABValue BAvalue
				// while csv style has a tab before colNameA
				if(header){
					String columnValue = row.get(columnIndex);
					if(columnValue.length() == 0){
						continue;
					}
					allColumns = addSingleValueTo2DArray(allColumns, columnIndexHeader,columnValue);
					columnIndexHeader++;
					continue;
				}
				else{
					// This changes the allColumns list of list in place, e.g. for example loop -> [[]] -> [[1]] -> [[1,2]] -> [[1,2],[3]] -> etc
					allColumns = addSingleValueTo2DArray(allColumns, columnIndex - 1, row.get(columnIndex));
				}
			}
			if(!header){
				firstColumn.add(row.get(0));
				if(row.size()-1 != columnIndexHeader){
					DeconvolutionLogger.log.info(String.format("Table %s does not have the same number of columns as there are in the header at row %d",filepath,rowNumber));
					DeconvolutionLogger.log.info(String.format("Number of header columns: %d",columnIndexHeader));
					DeconvolutionLogger.log.info(String.format("Number of columns at row %d: %d", rowNumber, row.size()-1));
					DeconvolutionLogger.log.info(row.toString());
					parser.close();
					throw new RuntimeException(String.format("Cellcount percentage table does not have the same number of columns as there are celltypes at row %d",rowNumber));
				}
			}
			if(header){
				header = false;
			}

		}
		parser.close();
		return new Object[] {firstColumn, allColumns};
	}

	/**
	 * Converting a vector of string to a vector of doubles
	 * 
	 * @param vector A vector of strings
	 * 
	 * @param start of vector from where to convert string to double
	 * 
	 * @return A vector of doubles
	 */
	public static double[] StringVectorToDoubleArrayList(String[] vector, int start) {
		int vectorLength = vector.length;
		double[] doubles = new double[vectorLength-start];
		// start at 1 because first element is sampleName
		for (int i = start; i < vectorLength; i++) {
				doubles[i-start] = Double.parseDouble(vector[i]);
		}
		return doubles;
	}

	public static <T> String listToTabSeparatedString(List<T> list)
	{
		/* Turn list into tab separated string*/
		StringBuilder builder = new StringBuilder();
		for(Object o: list)
		{
			builder.append(o+"\t");
		}
		return builder.toString().trim();
	}
	public static <T> String listToTabSeparatedString(double[] list)
	{
		/* Turn list into tab separated string*/
		StringBuilder builder = new StringBuilder();
		for(Double o: list)
		{
			builder.append(Double.toString(o)+"\t");
		}
		return builder.toString().trim();
	}

	/**
	 * Turn list into tab separated string
	 * 
	 *  @param <T> List type paramaeter
	 *  
	 *  @param list	List to separate into tabs
	 *  
	 *  @param append String to be added at end of each element
	 *  
	 *  @return Tab separated string
	 */
	public static <T> String listToTabSeparatedString(List<T> list, String append)
	{
		StringBuilder builder = new StringBuilder();
		for(Object o: list)
		{
			builder.append(o+append+"\t");
		}
		return builder.toString().trim();
	}

	/**
	 * Find the index of the top n values in array orig (from http://stackoverflow.com/a/17623521/651779)
	 * 
	 * @param orig Array of doubles to find highest values in
	 * @param n Number of highest values to report
	 */
	static ArrayList<Integer> indexesOfTopElements(double[] orig, int n) {
		double[] copy = Arrays.copyOf(orig,orig.length);
		Arrays.sort(copy);
		double[] honey = Arrays.copyOfRange(copy,copy.length - n, copy.length);
		ArrayList<Integer> result = new ArrayList<Integer>();
		for(int i = 0; i < orig.length; i++) {
			double onTrial = orig[i];
			int index = Arrays.binarySearch(honey,onTrial);
			if(index < 0) continue;
			result.add(i);
		}
		return result;
	}

	public static double[] convertDoubles(List<Double> doubles)
	{
		double[] ret = new double[doubles.size()];
		Iterator<Double> iterator = doubles.iterator();
		int i = 0;
		while(iterator.hasNext())
		{
			ret[i] = iterator.next();
			i++;
		}
		return ret;
	}

	static boolean equalLists(List<String> one, List<String> two){     
		if (one == null && two == null){
			return true;
		}

		if((one == null && two != null) 
				|| one != null && two == null
				|| one.size() != two.size()){
			return false;
		}

		//to avoid messing the order of the lists we will use a copy
		//as noted in comments by A. R. S.
		one = new ArrayList<String>(one); 
		two = new ArrayList<String>(two);   

		Collections.sort(one);
		Collections.sort(two);      
		return one.equals(two);
	}
	public static HashMap<String, ArrayList<String>> parseSnpPerGeneFile(String snpsToTestFile) throws IOException {
		LineIterator snpGenePairIterator = FileUtils.lineIterator(new File(snpsToTestFile), "UTF-8");
		HashMap<String, ArrayList<String>> geneSnpPairs = new HashMap<String, ArrayList<String>>();
		int totalSnpsToTest = 0;
		snpGenePairIterator.next();
		while (snpGenePairIterator.hasNext()) {
			ArrayList<String> snpGeneStringVector = new ArrayList<String>(Arrays.asList(snpGenePairIterator.next().split("\t")));
			String gene = snpGeneStringVector.get(0);
			String snp = snpGeneStringVector.get(1);
			ArrayList<String> snps = geneSnpPairs.get(gene);
			if (snps==null) {
				snps = new ArrayList<String>();
				geneSnpPairs.put(gene, snps);
			}
			snps.add(snp);
			totalSnpsToTest++;
		}
		DeconvolutionLogger.log.info(String.format("SNPs to deconvolute: %d", totalSnpsToTest));

		return geneSnpPairs;
	}
	
	/**
	 * Create permutations of binary numbers of length iterations. E.g. if iterations == 2, would give
	 * 00, 10, 01, 11
	 * 
	 * @param soFar The string that will contain the binary permutation. Should be empty when given to the function (used in the recursion)
	 * 
	 * @param iterations The length of the binary string
	 * 
	 * @param permutations List that contains all the permutations of the binary string. Should be given empty (used to save the recursion results)
	 * 
	 * @return List of permutations of 0's and 1's of given length
	 */
		public static ArrayList<String> binaryPermutations(String soFar, double iterations, ArrayList<String> permutations) {
		    if(iterations == 0) {
		    	permutations.add(soFar);
		    }
		    else {
		    	binaryPermutations(soFar + "0", iterations - 1, permutations);
		        binaryPermutations(soFar + "1", iterations - 1, permutations);
		    }
			return permutations;
		}
}
