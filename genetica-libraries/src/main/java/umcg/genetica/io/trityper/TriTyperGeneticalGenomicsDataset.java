/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;

/**
 * @author harmjan
 */
public final class TriTyperGeneticalGenomicsDataset implements Comparable<TriTyperGeneticalGenomicsDataset> {
	
	private TriTyperGenotypeData genotypeData;
	private TriTyperExpressionData expressionData;
	private THashMap<String, String> genotypeToExpressionCouplings;
	private TriTyperGeneticalGenomicsDatasetSettings settings;
	private short[] expressionToGenotypeIdArray;
	private short totalGGSamples;
	private boolean expressionDataLoadedCorrectly = true;
	private short[] genotypeToExpressionIdArray;
	DoubleMatrixDataset<String, String> covariates = null;
	
	public TriTyperGeneticalGenomicsDataset(TriTyperGeneticalGenomicsDatasetSettings settings, Pair<List<String>, List<List<String>>> pathwayDefinitions, boolean displayWarnings) throws IOException, Exception {
		
		this.settings = settings;
		
		settings.genotypeLocation = Gpio.formatAsDirectory(settings.genotypeLocation);
		
		if (settings.expressionLocation == null) {
			settings.expressionLocation = settings.genotypeLocation + "ExpressionData.txt";
		}
		
		// load the genotype metadata
		genotypeData = new TriTyperGenotypeData();
		genotypeData.displayWarnings = displayWarnings;
		genotypeData.load(settings.genotypeLocation);
		THashSet<String> includedExpressionIndividuals = new THashSet<String>();
		Boolean[] isIncluded = genotypeData.getIsIncluded();
		
		// preload the sample coupling file
		loadCouplings();
		
		// determine which expression samples to include
		Set<Entry<String, String>> entries = genotypeToExpressionCouplings.entrySet();
		for (Entry<String, String> entry : entries) {
			String genotypeIndividual = entry.getKey();
			Integer genotypeIndividualId = genotypeData.getIndividualId(genotypeIndividual);
			
			if (genotypeIndividualId != -9 && isIncluded[genotypeIndividualId] != null && isIncluded[genotypeIndividualId]) {
				includedExpressionIndividuals.add(entry.getValue());
			}
		}
		
		if (includedExpressionIndividuals.isEmpty()) {
			System.err.println("ERROR: none of the expression samples will be included with your current settings.\nPlease check the links between genotype and gene expression samples and/or your PhenotypeInformation.txt");
			System.exit(-1);
		}
		
		// load the expression data
		expressionData = new TriTyperExpressionData();
		expressionData.confineToProbes(settings.tsProbesConfine);
		expressionData.setConfineToProbesThatMapToAnyChromosome(settings.confineProbesToProbesMappingToAnyChromosome);
		expressionData.setConfineToProbesThatMapToChromosome(settings.confineProbesToProbesThatMapToChromosome);
		expressionData.setIncludeIndividuals(includedExpressionIndividuals);
		expressionData.setPathwayDefinitions(pathwayDefinitions);
		expressionDataLoadedCorrectly = expressionData.load(settings.expressionLocation, settings.probeannotation, settings.expressionplatform, (settings.cisAnalysis && settings.transAnalysis));
		
		pruneGenotypeToExpressionCouplings();
		
		if (settings.quantilenormalize) {
			QuantileNormalization.quantilenormalize(expressionData.getMatrix());
		}
		
		if (settings.logtransform) {
			Log2Transform.log2transform(expressionData.getMatrix());
		}
		
		if (settings.covariateFile != null && Gpio.exists(settings.covariateFile)) {
			// load covariates..
			System.out.println("Loading covariates: " + settings.covariateFile);
			HashSet<String> individualSet = new HashSet<String>();
			individualSet.addAll(Arrays.asList(expressionData.getIndividuals()));
			covariates = new DoubleMatrixDataset<String, String>(settings.covariateFile, null, individualSet);
			
			if (covariates.colObjects.isEmpty()) {
				// try the transpose
				System.out.println("Could not find matching sample identifiers between covariate file and expression file.\nTransposing your covariate file.");
				covariates = new DoubleMatrixDataset<String, String>(settings.covariateFile, individualSet);
				if (covariates.rowObjects.isEmpty()) {
					System.err.println("Could not find matching samples between expression data and covariate data.");
					System.exit(-1);
				} else {
					covariates.transposeDataset(); // put the covariates on the rows, samples on the columns
					covariates.recalculateHashMaps();
				}
			}
			
			covariates.removeColumnsWithNaNs();
			covariates.recalculateHashMaps();
			if (covariates.colObjects.isEmpty()) {
				System.err.println("ERROR: after removing samples with NaN values, no covariates remain");
				System.exit(-1);
			}
			
			System.out.println(covariates.rowObjects.size() + " covariates loaded for " + covariates.colObjects.size() + " samples");
			
			// remove expression samples without covariates, and reorder expression data
			expressionData.pruneAndReorderSamples(covariates.colObjects);
			
			// prune expression dataset to samples having covariates
			loadCouplings();
			pruneGenotypeToExpressionCouplings();
		}
		
	}
	
	public TriTyperGeneticalGenomicsDataset(TriTyperGeneticalGenomicsDatasetSettings settings) throws IOException, Exception {
		this(settings, null, true);
	}
	
	public TriTyperGeneticalGenomicsDataset(TriTyperGeneticalGenomicsDatasetSettings triTyperGeneticalGenomicsDatasetSettings, Pair<List<String>, List<List<String>>> pathwayDefinitions) throws Exception {
		this(triTyperGeneticalGenomicsDatasetSettings, pathwayDefinitions, true);
	}
	
	/**
	 * @return the genotypeData
	 */
	public TriTyperGenotypeData getGenotypeData() {
		return genotypeData;
	}
	
	/**
	 * @param genotypeData the genotypeData to set
	 */
	public void setGenotypeData(TriTyperGenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}
	
	/**
	 * @return the expressionData
	 */
	public TriTyperExpressionData getExpressionData() {
		return expressionData;
	}
	
	/**
	 * @param expressionData the expressionData to set
	 */
	public void setExpressionData(TriTyperExpressionData expressionData) {
		this.expressionData = expressionData;
	}
	
	/**
	 * @return the genotypeToExpressionCouplings
	 */
	public THashMap<String, String> getGenotypeToExpressionCouplings() {
		return genotypeToExpressionCouplings;
	}
	
	/**
	 * @param genotypeToExpressionCouplings the genotypeToExpressionCouplings to
	 *                                      set
	 */
	public void setGenotypeToExpressionCouplings(THashMap<String, String> genotypeToExpressionCouplings) {
		this.genotypeToExpressionCouplings = genotypeToExpressionCouplings;
	}
	
	/**
	 * @return the settings
	 */
	public TriTyperGeneticalGenomicsDatasetSettings getSettings() {
		return settings;
	}
	
	/**
	 * @param settings the settings to set
	 */
	public void setSettings(TriTyperGeneticalGenomicsDatasetSettings settings) {
		this.settings = settings;
	}
	
	public int getTotalGGSamples() {
		return totalGGSamples;
	}
	
	private void loadCouplings() throws IOException {
		genotypeToExpressionCouplings = new THashMap<String, String>();
		String genotypeToExpressionCoupling = settings.genotypeToExpressionCoupling;
		if (genotypeToExpressionCoupling != null && genotypeToExpressionCoupling.trim().length() > 0) {
			if (!Gpio.exists(genotypeToExpressionCoupling)) {
				throw new IOException("Error: genotype to expression coupling file: " + genotypeToExpressionCoupling + " does not exist.");
			}
			TextFile in = new TextFile(genotypeToExpressionCoupling, TextFile.R);
			String[] elems = in.readLineElemsReturnReference(TextFile.tab);
			
			while (elems != null) {
				if (elems.length > 1) {
					String key = new String(elems[0].getBytes("UTF-8"));
					String value = new String(elems[1].getBytes("UTF-8"));
					if (genotypeToExpressionCouplings.get(key) != null) {
						System.out.println("ERROR: your genotype to expression coupling file contains duplicate entries for individual: " + key);
						System.exit(0);
					} else {
						genotypeToExpressionCouplings.put(key, value);
					}
					
				}
				elems = in.readLineElemsReturnReference(TextFile.tab);
			}
			
			in.close();
		} else {
			Boolean[] isIncluded = genotypeData.getIsIncluded();
			int i = 0;
			String[] individuals = genotypeData.getIndividuals();
			for (String ind : individuals) {
				if (isIncluded[i] != null && isIncluded[i]) {
					if (genotypeToExpressionCouplings.get(ind) != null) {
						System.out.println("ERROR: your genotype data contains duplicate individuals: " + ind);
						System.exit(0);
					} else {
						genotypeToExpressionCouplings.put(ind, ind);
					}
				}
				i++;
			}
		}
	}
	
	public int[] getExpressionToGenotypeIdArray() {
		int[] intExpToGArr = new int[expressionToGenotypeIdArray.length];
		for (int i = 0; i < intExpToGArr.length; i++) {
			intExpToGArr[i] = expressionToGenotypeIdArray[i];
		}
		return intExpToGArr;
	}
	
	public short[] getExpressionToGenotypeIdArrayShort() {
		return expressionToGenotypeIdArray;
	}
	
	@Override
	public int compareTo(TriTyperGeneticalGenomicsDataset o) {
		int numIndsOther = o.getGenotypeData().getIndividuals().length;
		return genotypeData.getIndividuals().length - numIndsOther;
	}
	
	public boolean equals(TriTyperGeneticalGenomicsDataset o) {
		int numIndsOther = o.getGenotypeData().getIndividuals().length;
		if (genotypeData.getIndividuals().length == numIndsOther) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Permutes the mapping between genotype and gene expression samples
	 */
	public void permuteSampleLables(Random r) {
		ArrayList<Short> alIndWGA = new ArrayList<Short>();
		int numSamples = expressionToGenotypeIdArray.length;
		for (int i = 0; i < numSamples; i++) {
			if (expressionToGenotypeIdArray[i] != -1) {
				alIndWGA.add(expressionToGenotypeIdArray[i]);
			}
		}
		
		short[] indWGANew = new short[numSamples];
		
		genotypeToExpressionIdArray = new short[genotypeData.getIndividuals().length];
		for (int i = 0; i < numSamples; i++) {
			if (expressionToGenotypeIdArray[i] == -1) {
				indWGANew[i] = -1;
			} else {
				short genotypeId = alIndWGA.remove((int) (r.nextDouble() * (double) alIndWGA.size()));
				indWGANew[i] = genotypeId;
				
				genotypeToExpressionIdArray[genotypeId] = (short) i;
			}
		}
		expressionToGenotypeIdArray = indWGANew;
		
	}
	
	public void permuteCovariates(Random r) {
		// shuffle covariate, if any
		if (covariates != null) {
			System.out.println("Randomizing covariates");
			for (int covariate = 0; covariate < covariates.nrRows; covariate++) {
				ArrayList<Double> covariateData = new ArrayList<Double>();
				for (int sample = 0; sample < covariates.nrRows; sample++) {
					covariateData.add(covariates.rawData[covariate][sample]);
				}
				Collections.shuffle(covariateData, r);
				for (int sample = 0; sample < covariates.nrRows; sample++) {
					covariates.rawData[covariate][sample] = covariateData.get(sample);
				}
			}
		}
	}
	
	public void resetGenotypeToExpressionCouplings() throws IOException {
		loadCouplings();
	}
	
	public void pruneGenotypeToExpressionCouplings() {
		// now check whether each genotype is actually linked to an expression individual...
		String[] individuals = genotypeData.getIndividuals();
		
		Boolean[] isReallyIncluded = new Boolean[individuals.length];
		THashMap<String, String> realGenotypeToExpressionCouplings = new THashMap<String, String>();
		totalGGSamples = 0;
		for (int i = 0; i < isReallyIncluded.length; i++) {
			String genotypeInd = individuals[i];
			if (!genotypeToExpressionCouplings.containsKey(genotypeInd)) {
				isReallyIncluded[i] = false;
			} else {
				String coupledExpressionSample = genotypeToExpressionCouplings.get(genotypeInd);
				if (coupledExpressionSample != null) {
					Integer expressionSampleId = expressionData.getIndividualId(coupledExpressionSample);
					if (expressionSampleId == -9) {
						isReallyIncluded[i] = false;
					} else {
						isReallyIncluded[i] = true;
						realGenotypeToExpressionCouplings.put(genotypeInd, coupledExpressionSample);
						totalGGSamples++;
					}
				}
			}
		}
		
		// exclude genotypes for which no expression data is available.
		genotypeData.setIsIncluded(isReallyIncluded);
		genotypeToExpressionCouplings = realGenotypeToExpressionCouplings;
		
		// couple expression IDs to genotype IDs for quick reference
		Set<Entry<String, String>> entries = realGenotypeToExpressionCouplings.entrySet();
		expressionToGenotypeIdArray = new short[totalGGSamples];
		HashSet<Integer> visitedNumbers = new HashSet<Integer>();
		for (Entry<String, String> entry : entries) {
			Integer expressionIndId = expressionData.getIndividualId(entry.getValue());
			Integer genotypeIndId = genotypeData.getIndividualId(entry.getKey());
			if (expressionIndId != -9 && genotypeIndId != -9) {
				if (visitedNumbers.contains(expressionIndId)) {
					System.out.println("ERROR: your dataset contains duplicate samples!");
				} else {
					expressionToGenotypeIdArray[expressionIndId] = genotypeIndId.shortValue();
					visitedNumbers.add(expressionIndId);
				}
			}
		}
	}
	
	public HashMap<Integer, Integer> getGenotypeToExpressionIdHash() {
		HashMap<Integer, Integer> gte = new HashMap<Integer, Integer>();
		int expressionIndId = 0;
		for (int genotypeIndId : expressionToGenotypeIdArray) {
			gte.put(genotypeIndId, expressionIndId);
			expressionIndId++;
		}
		return gte;
	}
	
	public HashMap<Integer, Integer> getExpressionToGenotypeIdHash() {
		HashMap<Integer, Integer> etg = new HashMap<Integer, Integer>();
		int expressionIndId = 0;
		for (int genotypeIndId : expressionToGenotypeIdArray) {
			etg.put(expressionIndId, genotypeIndId);
			expressionIndId++;
		}
		return etg;
	}
	
	/**
	 * @return the expressionDataLoadedCorrectly
	 */
	public boolean isExpressionDataLoadedCorrectly() {
		return expressionDataLoadedCorrectly;
	}
	
	/**
	 * @param expressionDataLoadedCorrectly the expressionDataLoadedCorrectly to
	 *                                      set
	 */
	public void setExpressionDataLoadedCorrectly(boolean expressionDataLoadedCorrectly) {
		this.expressionDataLoadedCorrectly = expressionDataLoadedCorrectly;
	}
	
	public DoubleMatrixDataset<String, String> getCovariateData() {
		return covariates;
	}
}
