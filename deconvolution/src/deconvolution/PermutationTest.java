package deconvolution;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang3.time.DurationFormatUtils;

public class PermutationTest {
	private Map<String, List<Double>> pvaluePerCelltype = new HashMap<String, List<Double>>();
	private Set<String> celltypes = new HashSet<String>();
	private int numberOfPermutations = 0;
	public PermutationTest(){};
	
	private void addPerCelltype( String celltype, Double pvalues){
		celltypes.add(celltype);
		if (pvaluePerCelltype.get(celltype) == null) { //gets the value for an id)
			pvaluePerCelltype.put(celltype, new ArrayList<Double>()); //no ArrayList assigned, create new ArrayList
		}
		pvaluePerCelltype.get(celltype).add(pvalues); //adds value to list.
   }
	
	public void add( DeconvolutionResult deconvolutionResult, int numberOfPermutations) throws IllegalAccessException{
		for(String celltype : deconvolutionResult.getCelltypes()){
			celltypes.add(celltype);
			this.numberOfPermutations = numberOfPermutations;
			addPerCelltype(celltype, deconvolutionResult.getPvaluePerCelltype().get(celltype));
		}
   }
	
	public Set<String> getCelltypes(){
		return(celltypes);
   }

	public List<Double> getPvalues(String celltype){
		return(pvaluePerCelltype.get(celltype));
   }
	
	public int getNumberOfPermutations(){
		return(numberOfPermutations);
	}
	
	public void permutationTest(String expressionFile, String genotypeFile, int numberOfPermutation, String permutationType) throws Exception {
		/*
		 * Scramble the expression values for all the samples n number of times
		 * and calculate p-value, for comparison to p-value from unscrambled
		 * p-value
		 * 
		 * @param expressionIterator Iterator containing all rows of expression file 
		 * 
		 * @param genotypeIterator Iterator containing all rows of genotype file, same order as expression iterator 
		 * 
		 * @param numberOfPermutations The number of permutations to do
		 * 
		 * @param permutationType Type of permuation to do, e.g. permute genotype or permute expression
		 * 
		 * @return PermutationResult, containing hasmap of per cell type all pvalues
		 */
		System.out.println("Starting permutation test...");
		long time = System.currentTimeMillis();
		// make iterator of files so that we can loop over 2 at the same time
		for(int i = 0; i < numberOfPermutation; i++){
			LineIterator expressionIterator = FileUtils.lineIterator(new File(expressionFile), "UTF-8");
			LineIterator genotypeIterator = FileUtils.lineIterator(new File(genotypeFile), "UTF-8");
			// Skip over the header
			expressionIterator.next();
			genotypeIterator.next();
			// set a seed so that all rows in the expression file can be shuffled identically
			long seed = System.nanoTime();
			// change number after % to get more/less prints of progress of permutations
			if (i % 1 == 0) {
				long completedIn = System.currentTimeMillis() - time;
				DeconvolutionLogger.log.info(String.format("Processed %d/%d permutations - %s\n", i, numberOfPermutation,
						DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS")));

			}
			while (expressionIterator.hasNext() && genotypeIterator.hasNext()) {

				ArrayList<String> expressionStringVector = new ArrayList<String>(Arrays.asList(expressionIterator.next().split("\t")));

				ArrayList<String> genotypeStringVector = new ArrayList<String>(Arrays.asList(genotypeIterator.next().split("\t")));
				// shuffle on expression or genotype, because at the moment not sure which one should be used for permutation
				// and using this variable is easier than commenting in/out the 3 related lines
				String qtlName = "";
				if (permutationType.equals("expression")){
					qtlName = expressionStringVector.remove(0);
					// shuffle the expression vector with the previously set seed
					Collections.shuffle(expressionStringVector, new Random(seed));
					// add back the qtl name at the front
					expressionStringVector.add(0, qtlName);
				}
				else if (permutationType.equals("genotype")){
					qtlName = genotypeStringVector.remove(0);
					// shuffle the genotype vector with the previously set seed
					Collections.shuffle(genotypeStringVector, new Random(seed));
					// add back the qtl name at the front
					genotypeStringVector.add(0, qtlName);
				}
				else{
					throw new Exception("shuffle should be either expression or genotype, not "+permutationType);
				}

				// deconvolution() gives for the current SNP-Gene QTL pair for each celltype in cellount table a pvalue. deconResultsPermutation.add() adds
				// the pvalue to a hashmap of pvalues per celltype. After permutation for loop is finished, deconResultsPermutation will contain a hashmap
				// with for each celltype all pvalues found during permutaiton testing
					throw new java.lang.UnsupportedOperationException("Not implemented yet.");
				//try{
					//deconResultsPermutation.add(deconvolution(expressionStringVector, genotypeStringVector), numberOfPermutations);
				//}
				//catch(NotEnoughGenotypesException e){
					// If there are not enough samples per genotype for this QTL, exclude it from testing
				//}
			}
		}

		throw new java.lang.UnsupportedOperationException("Not implemented yet.");
		/*commandLineOptions.setOutfolder(outputFolder+"/shuffled_" + commandLineOptions.getPermutationType());
		if(commandLineOptions.getForceNormalExpression()){
			commandLineOptions.setOutfolder(outputFolder+"_forcedNormalExpression");
			commandLineOptions.setOutfolder(outputFolder+"_"+commandLineOptions.getNormalizationType());
		}
		if(commandLineOptions.getForceNormalCellcount()){
			commandLineOptions.setOutfolder(outputFolder+"_forceNormalCellcount");
		}
		new File(outputFolder).mkdirs();
		Plots.drawHistogram(deconResultsPermutation, outputFolder);*/
	}
}
