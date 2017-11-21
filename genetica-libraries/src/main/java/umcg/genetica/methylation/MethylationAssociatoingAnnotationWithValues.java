package umcg.genetica.methylation;

import JSci.maths.ArrayMath;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import org.apache.commons.collections.primitives.ArrayDoubleList;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.inference.OneWayAnova;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Heterogeneity;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.ZScores;

import static umcg.genetica.methylation.ParseTcgaMethylationFile.ENCODING;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author MarcJan
 */
public class MethylationAssociatoingAnnotationWithValues {
	
	private static Pattern SPLIT_ON_TAB = Pattern.compile("\\t");
	private static Pattern SPLIT_PARTS = Pattern.compile("-");
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		
		String fileWithAnnotation = "D:\\UMCG\\Methylation_GPL8490\\TCGA+GEO_14112012\\Annotation_AllSamples.txt";
		String dataFile = "D:\\UMCG\\Methylation_GPL8490\\TCGA+GEO_14112012\\methylation_Matrix_SexFiltered.QuantileNormalized.txt";
		//String dataFile = "D:\\UMCG\\Methylation_GPL8490\\TCGA+GEO_14112012\\NeverQN-methylation_Matrix_SexFiltered.txt";
		
		System.out.print("Read annotation file .... ");
		HashMap<String, SoftfileAnnotation> sampleAnnotation = readAnnotationFile(fileWithAnnotation);
		System.out.println("done");
		
		System.out.print("Read eigenvector file .... ");
		DoubleMatrixDataset<String, String> eigenVectors = readDoubleMatrixFile(dataFile);
		System.out.println("done");
		
		ArrayList<String> setSelection = new ArrayList<String>();
		
		//Only Blood
		//setSelection.addAll(Arrays.asList("GSE20236", "GSE23638","GSE19711","GSE20067","GSE41037"));
		//No Cancer
		setSelection.addAll(Arrays.asList("GSE20236", "GSE23638", "GSE19711", "GSE20067", "GSE15745 // GSE36194", "GSE15745", "GSE41037", "GSE32393", "GSE31979", "GSE20242",
				"GSE20080", "GSE36194", "GSE22595", "GSE29661", "GSE21232", "GSE30653 // GSE30654", "GSE27097", "GSE37988", "GSE32861 // GSE32867", "GSE17448", "GSE33422",
				"GSE25033", "GSE34035", "GSE28746", "GSE32396"));
		//No Blood
//        setSelection.addAll(Arrays.asList("GSE15745 // GSE36194","GSE15745","GSE32393","TCGA.brca","GSE31979","GSE30758 // GSE30760","GSE30759 // GSE30760",
//                "GSE20080","TCGA.coad","TCGA.gbm","GSE36194","GSE22595","GSE29661","GSE21232","GSE30653 // GSE30654","TCGA.kirc","TCGA.kirp","GSE37988",
//                "GSE32861 // GSE32867","TCGA.lusc","TCGA.luad","GSE17448","GSE33422","TCGA.ov","GSE25033","TCGA.read","GSE34035","GSE28746","TCGA.stad",
//                "GSE30844 // GSE31337","TCGA.ucec"));
//        
		//Analysis type 1) Age
		String infoKey = "Age";
		
		//Analysis type 2) Gender
		//String infoKey = "Gender";
		//ArrayList<String> entries = new ArrayList<String>();
		//entries.addAll(Arrays.asList("Male", "Female"));
		
		//Analysis type 3) GSE         //Infokey is not used
		//String infoKey = "";
		
		String nameSeriesInfoColumn = "series id";
		
		LinkedHashMap<String, HashMap<String, String>> interestSets;
		
		//Depending on the data we need to transpone the data
		eigenVectors = eigenVectors.getTransposedDataset();
		
		interestSets = selectSamplesWithInformationOfInterest(sampleAnnotation, nameSeriesInfoColumn, infoKey, eigenVectors, 25);
		
		if (setSelection.size() > 0) {
			ArrayList<String> removeEntry = new ArrayList<String>();
			for (Entry<String, HashMap<String, String>> e : interestSets.entrySet()) {
				if (!(setSelection.contains(e.getKey()))) {
					removeEntry.add(e.getKey());
				}
			}
			for (String s : removeEntry) {
				interestSets.remove(s);
			}
		}
		
		//interestSets = selectSamplesWithSeriesInformation(sampleAnnotation, eigenVectors);//This needs a refresh
		
		//Willen we hier wat mee?
		//int maxSize = 150;
		//interestSets = splitInterstingSetInPortions(interestSets, maxSize);
		//System.exit(0);
		
		for (Entry<String, HashMap<String, String>> tmp : interestSets.entrySet()) {
			System.out.println(tmp.getKey() + "\t" + tmp.getValue().size());
		}
		
		System.out.println("Number of interest sets: " + interestSets.size());
		
		//Analysis type 1) Age
		correlateScoreAndItemOfInterest(eigenVectors, interestSets);
		
		//Analysis type 2) Gender
		//associateTTestScoreAndItemOfInterest(eigenVectors, interestSets, entries);
		
		//Analysis type 3) GSE
		//associateAnovaScoreAndItemOfInterest(eigenVectors, interestSets);
		
	}
	
	/**
	 * Read annotation file Tab separated file containing sample annotation
	 *
	 * @param fileWithAnnotation
	 * @return Sample annotation
	 */
	private static HashMap<String, SoftfileAnnotation> readAnnotationFile(String fileWithAnnotation) {
		
		HashMap<String, SoftfileAnnotation> sampleInfo = new HashMap<String, SoftfileAnnotation>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileWithAnnotation)), ENCODING), 8096);
			String str = in.readLine();
			
			String[] headers = SPLIT_ON_TAB.split(str);
			
			int meshInfoIndex = -1;
			
			for (int i = 1; i < headers.length; ++i) {
				if (headers[i].toLowerCase().contains("mesh")) {
					meshInfoIndex = i;
					break;
				}
			}
			
			while ((str = in.readLine()) != null) {
				String[] entries = SPLIT_ON_TAB.split(str);
				SoftfileAnnotation tmp = new SoftfileAnnotation();
				
				tmp.setAccession(entries[0]);
				
				if (!(meshInfoIndex < 0)) {
					tmp.setMeshTerms(entries[meshInfoIndex]);
				}
				
				for (int i = 1; i < entries.length; ++i) {
					tmp.putAnnotationInformation(headers[i], entries[i]);
				}
				
				sampleInfo.put(entries[0], tmp);
			}
			in.close();
		} catch (IOException e) {
			System.out.println(e.getMessage());
			System.exit(-1);
		}
		
		return (sampleInfo);
	}
	
	/**
	 * Read double matrix file Eigenvector file / pc file / probe matrix
	 *
	 * @param eigenVectorFile
	 * @return
	 */
	private static DoubleMatrixDataset<String, String> readDoubleMatrixFile(String eigenVectorFile) {
		
		DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDataset<String, String>();
		try {
			tmp = new DoubleMatrixDataset<String, String>(eigenVectorFile, "\t");
		} catch (IOException ex) {
			Logger.getLogger(MethylationAssociatoingAnnotationWithValues.class.getName()).log(Level.SEVERE, null, ex);
		}
		
		return (tmp);
	}
	
	/**
	 * Filter out interest GSE sets. Sets need to have at least 2 different
	 * values for the infoKey of interest. Automagicaly it checks if samples are
	 * in the double matrix dataset
	 *
	 * @param sampleAnnotation Annotation information
	 * @param infoKey          Key for annotation of interest
	 * @param doubleMatrix     double matrix dataset
	 * @return
	 */
	private static LinkedHashMap<String, HashMap<String, String>> selectSamplesWithInformationOfInterest(HashMap<String, SoftfileAnnotation> sampleAnnotation, String nameSeriesInfoColumn, String infoKey, DoubleMatrixDataset<String, String> eigenVectors, int minimalNumberSamplesInSeries) {
		LinkedHashMap<String, HashMap<String, String>> gseSets = new LinkedHashMap<String, HashMap<String, String>>();
		
		ArrayList<String> removeSamples = new ArrayList<String>();
		
		for (Entry<String, SoftfileAnnotation> tmp : sampleAnnotation.entrySet()) {
			if (!(tmp.getValue().getAnnotationInformation().containsKey(infoKey))) {
				System.out.print("No " + infoKey + " information");
				System.exit(0);
			}
			break;
		}
		
		
		for (String sampleName : eigenVectors.rowObjects) {
			if (sampleAnnotation.containsKey(sampleName)) {
				SoftfileAnnotation sampleAnnot = sampleAnnotation.get(sampleName);
				if (!(sampleAnnot.getAnnotationInformation().get(infoKey).isEmpty()) || !(sampleAnnot.getAnnotationInformation().get(infoKey).equals(""))) {
					String seriesId = sampleAnnot.getAnnotationInformation().get(nameSeriesInfoColumn);
					if (gseSets.containsKey(seriesId)) {
						System.out.println(sampleName + "\t" + sampleAnnot.getAnnotationInformation().get(infoKey));
						gseSets.get(seriesId).put(sampleName, sampleAnnot.getAnnotationInformation().get(infoKey));
					} else {
						HashMap<String, String> tmp = new HashMap<String, String>();
						tmp.put(sampleName, sampleAnnot.getAnnotationInformation().get(infoKey));
						
						gseSets.put(seriesId, tmp);
					}
				} else {
					System.out.println("No age info: " + sampleName);
					removeSamples.add(sampleName);
				}
				
				
			} else if (sampleName.startsWith("TCGA-")) {
				String[] sampleIdParts = SPLIT_PARTS.split(sampleName);
				String newSampleName = sampleIdParts[0] + "-" + sampleIdParts[1] + "-" + sampleIdParts[2];
				
				if (sampleAnnotation.containsKey(newSampleName)) {
					SoftfileAnnotation sampleAnnot = sampleAnnotation.get(newSampleName);
					if (!(sampleAnnot.getAnnotationInformation().get(infoKey).isEmpty()) || !(sampleAnnot.getAnnotationInformation().get(infoKey).equals(""))) {
						String seriesId = sampleAnnot.getAnnotationInformation().get(nameSeriesInfoColumn);
						if (gseSets.containsKey(seriesId)) {
							System.out.println(sampleName + "\t" + sampleAnnot.getAnnotationInformation().get(infoKey));
							gseSets.get(seriesId).put(sampleName, sampleAnnot.getAnnotationInformation().get(infoKey));
						} else {
							HashMap<String, String> tmp = new HashMap<String, String>();
							tmp.put(sampleName, sampleAnnot.getAnnotationInformation().get(infoKey));
							
							gseSets.put(seriesId, tmp);
						}
					} else {
						System.out.println("No age info: " + sampleName);
						removeSamples.add(sampleName);
					}
					
					
				} else {
					System.out.println("Not in matrix: " + sampleName);
					removeSamples.add(sampleName);
				}
			} else {
				System.out.println("Not in matrix: " + sampleName);
				removeSamples.add(sampleName);
			}
		}
		
		ArrayList<String> removeGseSets = new ArrayList<String>();
		
		int numberOfInterestSets = 0;
		int numberOfInterestSamples = 0;
		if (gseSets.size() > 0) {
			for (Entry<String, HashMap<String, String>> gse : gseSets.entrySet()) {
				
				ArrayList<String> uniqueValues = new ArrayList<String>();
				for (Entry<String, String> sample : gse.getValue().entrySet()) {
					if (!uniqueValues.contains(sample.getValue())) {
						uniqueValues.add(sample.getValue());
					}
				}
				
				if (uniqueValues.size() > 1 && gse.getValue().size() >= minimalNumberSamplesInSeries) {
					numberOfInterestSets++;
					// System.out.println(gse.getKey());
					for (Entry<String, String> sample : gse.getValue().entrySet()) {
						numberOfInterestSamples++;
						//System.out.println("\t" + sample.getKey() + "\t" + sample.getValue());
					}
				} else {
					removeGseSets.add(gse.getKey());
					for (Entry<String, String> sample : gse.getValue().entrySet()) {
						removeSamples.add(sample.getKey());
					}
				}
			}
			
		} else {
			System.out.println("Unforeseen error check Key and code");
			System.exit(0);
		}
		
		
		System.out.println("Number of sets: " + numberOfInterestSets);
		System.out.println("Total samples of interest: " + numberOfInterestSamples);
		
		for (String removeEntry : removeGseSets) {
			gseSets.remove(removeEntry);
		}
		
		for (String removeEntry : removeSamples) {
			sampleAnnotation.remove(removeEntry);
		}
		
		return (gseSets);
		
	}
	
	/**
	 * Filter out interest GSE sets. Samples need to have a GSE id. Automagicaly
	 * it checks if samples are in the double matrix dataset
	 *
	 * @param sampleAnnotation Annotation information
	 * @param infoKey          Key for annotation of interest
	 * @param doubleMatrix     double matrix dataset
	 * @return
	 */
	private static HashMap<String, HashMap<String, String>> selectSamplesWithSeriesInformation(HashMap<String, SoftfileAnnotation> sampleAnnotation, DoubleMatrixDataset<String, String> eigenVectors) {
		HashMap<String, HashMap<String, String>> gseSets = new HashMap<String, HashMap<String, String>>();
		
		ArrayList<String> removeSamples = new ArrayList<String>();
		
		for (Entry<String, SoftfileAnnotation> sample : sampleAnnotation.entrySet()) {
			if (!(sample.getValue().getAnnotationInformation().get("series id").isEmpty()) || !(sample.getValue().getAnnotationInformation().get("series id").equals(""))) {
				if (eigenVectors.rowObjects.contains(sample.getKey())) {
					String seriesId = sample.getValue().getAnnotationInformation().get("series id");
					if (gseSets.containsKey(seriesId)) {
						gseSets.get(seriesId).put(sample.getKey(), sample.getValue().getAnnotationInformation().get("series id"));
					} else {
						HashMap<String, String> tmp = new HashMap<String, String>();
						tmp.put(sample.getKey(), sample.getValue().getAnnotationInformation().get("series id"));
						gseSets.put(seriesId, tmp);
					}
				}
			} else {
				removeSamples.add(sample.getKey());
			}
			
		}
		
		for (String removeEntry : removeSamples) {
			sampleAnnotation.remove(removeEntry);
		}
		
		return (gseSets);
		
	}
	
	/**
	 * Test for difference between 2 groups with T-test.
	 *
	 * @param doubleMatrix
	 * @param interestSets
	 * @param entries      names of the two groups
	 */
	private static void associateTTestScoreAndItemOfInterest(DoubleMatrixDataset<String, String> doubleMatrix, HashMap<String, HashMap<String, String>> interestSets, ArrayList<String> entries) {
		//System.out.println(eigenVectors.nrCols);
		//System.out.println(eigenVectors.nrRows);
		HashMap<String, Double> scorePerGse = new HashMap<String, Double>();
		HashMap<String, Integer> indeces = new HashMap<String, Integer>();
		
		for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
			for (Entry<String, String> sample : set.getValue().entrySet()) {
				if (doubleMatrix.rowObjects.contains(sample.getKey())) {
					int index = doubleMatrix.rowObjects.indexOf(sample.getKey());
					indeces.put(sample.getKey(), index);
				} else {
					System.out.println("Potential mismatch between annotation and samples");
					System.out.println(sample.getKey() + " is not in value matrix");
					System.out.println("\n However :" + indeces.size() + " are in the matrix");
					System.exit(0);
				}
			}
		}
		
		//System.out.println(indeces.size());
		
		for (int i = 0; i < doubleMatrix.nrCols; ++i) {
			//System.out.println(doubleMatrix.colObjects.get(i));
			for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
				
				ArrayDoubleList valueSet1 = new ArrayDoubleList();
				ArrayDoubleList valueSet2 = new ArrayDoubleList();
				
				//ArrayList<String> keySet1 = new ArrayList<String>();
				//ArrayList<String> keySet2 = new ArrayList<String>();
				
				for (Entry<String, String> sample : set.getValue().entrySet()) {
					if (sample.getValue().equals(entries.get(0))) {
						//keySet1.add(sample.getKey());
						valueSet1.add(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
					} else if (sample.getValue().equals(entries.get(1))) {
						//keySet2.add(sample.getKey());
						valueSet2.add(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
					}
				}
				double[] set1 = valueSet1.toArray(new double[0]);
				double[] set2 = valueSet2.toArray(new double[0]);
				
				if (set1.length > 2 && set2.length > 2) {
					double zScore = TTest.testZscore(set1, set2);
					//System.out.println(doubleMatrix.colObjects.get(i)+"_"+set.getKey()+"\t"+pValue);
					scorePerGse.put(doubleMatrix.colObjects.get(i) + "_" + set.getKey(), zScore);
				}
			}
		}
	}
	
	/**
	 * Correlate values of interest to age
	 *
	 * @param doubleMatrix
	 * @param interestSets
	 */
	private static void correlateScoreAndItemOfInterest(DoubleMatrixDataset<String, String> doubleMatrix, LinkedHashMap<String, HashMap<String, String>> interestSets) {
		//System.out.println(eigenVectors.nrCols);
		//System.out.println(eigenVectors.nrRows);
		HashMap<String, Double> scorePerGse = new HashMap<String, Double>();
		HashMap<String, Integer> indeces = new HashMap<String, Integer>();
		int largestSet = 0;
		for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
			for (Entry<String, String> sample : set.getValue().entrySet()) {
				if (doubleMatrix.rowObjects.contains(sample.getKey())) {
					int index = doubleMatrix.rowObjects.indexOf(sample.getKey());
					indeces.put(sample.getKey(), index);
				} else {
					System.out.println("Potential mismatch between annotation and samples");
					System.out.println(sample.getKey() + " is not in value matrix");
					System.out.println("\n However :" + indeces.size() + " are in the matrix");
					System.exit(0);
				}
			}
			if (largestSet < set.getValue().size()) {
				largestSet = set.getValue().size();
			}
		}
		
		Correlation.correlationToZScore(largestSet);
		//System.out.println(indeces.size());
		
		double[] metaZ = new double[doubleMatrix.nrCols];
		System.out.println("Z-scores");
		System.out.print("\tMeta Z\tpValue\tHeterogeneity\tHeterogeneity pValue");
		
		for (String t : interestSets.keySet()) {
			System.out.print("\t" + t);
		}
		System.out.print("\t");
		for (String t : interestSets.keySet()) {
			System.out.print("\t" + t);
		}
		
		System.out.println("");
		SpearmansCorrelation sc = new SpearmansCorrelation();
		for (int i = 0; i < doubleMatrix.nrCols; ++i) {
			//System.out.println(doubleMatrix.colObjects.get(i));
			double[] zScores = new double[interestSets.size()];
			double[] correlations = new double[interestSets.size()];
			int[] setSizes = new int[zScores.length];
			int index = 0;
			for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
				int sizeOfGseSet = set.getValue().size();
				setSizes[index] = sizeOfGseSet;
				ArrayDoubleList valueSet = new ArrayDoubleList();
				ArrayDoubleList ageSet = new ArrayDoubleList();
				//ArrayList<String> keySet1 = new ArrayList<String>();
				
				for (Entry<String, String> sample : set.getValue().entrySet()) {
					//keySet1.add(sample.getKey());
					//System.out.println(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
					//System.out.println(Double.parseDouble(sample.getValue()));
					
					valueSet.add(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
					try {
						ageSet.add(Double.parseDouble(sample.getValue()));
					} catch (NumberFormatException ex) { // data not numerical, assume gender here as a quick fix
						ageSet.add("male".equals(sample.getValue().toLowerCase()) ? 1 : 2);
					}
				}
				double[] setValues = valueSet.toArray(new double[0]);
				double[] setAges = ageSet.toArray(new double[0]);
				
				if (setValues.length > 2) {
					double spearman = sc.correlation(setValues, setAges);
//                    double correlation = JSci.maths.ArrayMath.correlation(setValues, setAges);
//                    double zScore = Correlation.convertCorrelationToZScore(sizeOfGseSet, correlation);
					double zScore = Correlation.convertCorrelationToZScore(sizeOfGseSet, spearman);
					zScores[index] = zScore;
					correlations[index] = spearman;

//                    System.out.println(doubleMatrix.colObjects.get(i) + "_" + set.getKey() + "\t" + correlation);
					scorePerGse.put(doubleMatrix.colObjects.get(i) + "_" + set.getKey(), zScore);
				} else {
					zScores[index] = Double.NaN;
				}
				index++;
			}
			double zSum = 0;
			double sampleSizeSum = 0;
			for (int j = 0; j < zScores.length; j++) {
				if (!Double.isNaN(zScores[j])) {
					zSum += Math.sqrt(setSizes[j]) * zScores[j];
					sampleSizeSum += setSizes[j];
				}
			}
			zSum = zSum / Math.sqrt(sampleSizeSum);
			double p = ZScores.zToP(zSum);
			
			Triple<Double, Double, Integer> hg = Heterogeneity.getISq(zScores, setSizes);
			//System.out.print("\tMeta Z\tpValue\tHeterogeneity\tHeterogeneity pValue");
			System.out.print(doubleMatrix.colObjects.get(i) + "\t" + zSum + "\t" + p + "\t" + hg.getLeft() + "\t" + hg.getMiddle());
			for (double z : zScores) {
				System.out.print("\t" + z);
			}
			System.out.print("\t");
			for (double r : correlations) {
				System.out.print("\t" + r);
			}
			System.out.println("");
			
			metaZ[i] = zSum;
		}
		System.out.println("");
		for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
			for (Entry<String, String> e : set.getValue().entrySet()) {
				Integer sampleIndex = indeces.get(e.getKey());
				double correlation = ArrayMath.correlation(doubleMatrix.rawData[sampleIndex], metaZ);
				System.out.println(e.getKey() + "\t" + e.getValue() + "\t" + correlation);
			}
		}

//        RankDoubleArray ra = new RankDoubleArray();
//        double[] rank = ra.rank(metaZ);
//        System.out.println(rank[0]);
	}
	
	/**
	 * Test for difference between all GSE groups with Anova.
	 *
	 * @param doubleMatrix
	 * @param interestSets
	 * @param entries      names of the two groups
	 */
	private static void associateAnovaScoreAndItemOfInterest(DoubleMatrixDataset<String, String> doubleMatrix, HashMap<String, HashMap<String, String>> interestSets) {
		//System.out.println(eigenVectors.nrCols);
		//System.out.println(eigenVectors.nrRows);
		//HashMap<String, Double> scorePerGse = new HashMap<String, Double>();
		HashMap<String, Integer> indeces = new HashMap<String, Integer>();
		
		for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
			for (Entry<String, String> sample : set.getValue().entrySet()) {
				if (doubleMatrix.rowObjects.contains(sample.getKey())) {
					int index = doubleMatrix.rowObjects.indexOf(sample.getKey());
					indeces.put(sample.getKey(), index);
				} else {
					System.out.println("Potential mismatch between annotation and samples");
					System.out.println(sample.getKey() + " is not in value matrix");
					System.out.println("\n However :" + indeces.size() + " are in the matrix");
					System.exit(0);
				}
			}
		}
		
		//System.out.println(indeces.size());
		
		for (int i = 0; i < doubleMatrix.nrCols; ++i) {
			//System.out.println(doubleMatrix.colObjects.get(i));
			//ArrayList<String> keyList = new ArrayList<String>();
			ArrayList<double[]> valueSets = new ArrayList<double[]>();
			
			for (Entry<String, HashMap<String, String>> set : interestSets.entrySet()) {
				ArrayDoubleList valueSet = new ArrayDoubleList();
				
				for (Entry<String, String> sample : set.getValue().entrySet()) {
					valueSet.add(doubleMatrix.rawData[indeces.get(sample.getKey())][i]);
				}
				
				if (valueSet.size() > 2) {
					//keyList.add(set.getKey());
					valueSets.add(valueSet.toArray(new double[0]));
				}
				
			}
			
			OneWayAnova anova = new OneWayAnova();
			double pValue = -1;
			
			try {
				pValue = anova.anovaPValue(valueSets);
			} catch (IllegalArgumentException ex) {
				Logger.getLogger(MethylationAssociatoingAnnotationWithValues.class.getName()).log(Level.SEVERE, null, ex);
			}
			
			System.out.println("Component: " + doubleMatrix.colObjects.get(i) + " capable of discriminating between the sets, with p-value: " + pValue);
		}
	}
	
	private static LinkedHashMap<String, HashMap<String, String>> splitInterstingSetInPortions(LinkedHashMap<String, HashMap<String, String>> interestSets, int maxSize) {
		//HashMap<String, Integer> setsToSplit = new HashMap<String, Integer>();
		
		LinkedHashMap<String, HashMap<String, String>> newInterestingSets = new LinkedHashMap<String, HashMap<String, String>>();
		
		for (Entry<String, HashMap<String, String>> series : interestSets.entrySet()) {
			if (series.getValue().size() > maxSize) {
				System.out.println(series.getValue().size() % maxSize);
			}
		}
		
		return (newInterestingSets);
	}
}
