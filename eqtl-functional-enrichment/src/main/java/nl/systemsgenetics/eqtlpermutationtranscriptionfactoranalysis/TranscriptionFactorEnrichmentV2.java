/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.io.regulomedb.RegulomeDbEntry;
import umcg.genetica.io.regulomedb.RegulomeDbFile;
import umcg.genetica.io.regulomedb.RegulomeDbFiles;
import umcg.genetica.io.regulomedb.RegulomeDbSupportingData;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author Matthieu
 */
public class TranscriptionFactorEnrichmentV2 {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public static void main(String[] args)throws IOException{
		/*
		 * args[0]: eQTL file.
		 * args[1]: permutation files locations.
		 * args[2]: genotype matrix location.
		 * args[3]: regulomedb location.
		 * args[4]: window size.
		 * args[5]: r2 cutoff.
		 */
		//EQtlPermutationTranscriptionFactorAnalysisV3 eqptfa3 = new EQtlPermutationTranscriptionFactorAnalysisV3(args[0], args[1], args[2], args[3], Integer.parseInt(args[4]), Double.valueOf(args[5]));
		EQtlPermutationTranscriptionFactorAnalysisV3 eqptfa3 =
				new EQtlPermutationTranscriptionFactorAnalysisV3("", "C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\eQTLsFDR0.05.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\PermutedEQTLsPermutationRound8.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\BloodHT12Combined\\", "C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\regulomeDb\\", 0.8, "");
	}
	
	
	public TranscriptionFactorEnrichmentV2(String eQtlProbeLevelFile, String eQtlFile, String permutationFile, String genotypeLocation, String regulomeDbLocation, double r2Cutoff, String outputFile)throws IOException{
		//STEP 1.: READ REGULOMEDB DATA.
		ArrayList<RegulomeDbFile> regulomeDbFiles = new ArrayList<RegulomeDbFile>();
		regulomeDbFiles.add(new RegulomeDbFile(new File(regulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category1.txt")));
		regulomeDbFiles.add(new RegulomeDbFile(new File(regulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category2.txt")));
		regulomeDbFiles.add(new RegulomeDbFile(new File(regulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category3.txt")));
		regulomeDbFiles.add(new RegulomeDbFile(new File(regulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category4.txt")));
		regulomeDbFiles.add(new RegulomeDbFile(new File(regulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category5.txt")));
		HashMap<String, TreeMap<Integer, String[]>> regulomeDbData = readRegulomeDbData(regulomeDbFiles);
		
		
		//STEP 2.: GET A LIST OF PROBES SHARED BY ALL DATASETS.
		HashSet<String> sharedProbes = getSharedProbes(eQtlProbeLevelFile, permutationFile);
		
		//STEP 3.: READ THE EQTL FDR 0.05 DATA AND FILTER ON TOP EFFECTS AND NON TOP EFFECTS FOR SHARED PROBES.	
		EQTL[] eqtls = readEQtlData(eQtlFile);
		Set<String> rsIdList = makeRsIdList(eqtls);
		HashMap<String, EQTL> topEqtlEffects = getTopEffects(eqtls, sharedProbes);
		HashMap<String, HashSet<Integer>> nonTopEqtlEffects = getNonTopEffects(eqtls, topEqtlEffects);
		
		
		//STEP 4.: READ THE GENOTYPE MATRIX DATA.
		RandomAccessGenotypeData genotypeMatrixData = readEQtlGenotypeDataV2(genotypeLocation, rsIdList);
		HashMap<String, GeneticVariant> genotypeVariantIdMap = genotypeMatrixData.getVariantIdMap();
		System.out.println(genotypeVariantIdMap.size());
		
		//STEP 5.: PERFORM ANALYSIS FOR EQTL DATA.
		HashMap<String, Integer> eQtlCounts = new HashMap<String, Integer>();
		getTranscriptionFactorCounts(topEqtlEffects, nonTopEqtlEffects, eQtlCounts, genotypeMatrixData, regulomeDbData, r2Cutoff);
		
		
		//STEP 6.: PERFORM ANALYSIS FOR PERMUTATION DATA.
		HashMap<String, Integer> permutationCounts = new HashMap<String, Integer>();
		for(int n=1;n<=100;n++){
			EQTL[] permutationData;
			HashMap<String, EQTL> topPermutationEffects;
			HashMap<String, HashSet<Integer>> nonTopPermutationEffects;
			
			permutationData = readPermutationData(permutationFile + "PermutedEQTLsPermutationRound" + n + ".txt.gz", genotypeVariantIdMap);
			topPermutationEffects = getTopEffects(permutationData, sharedProbes);
			nonTopPermutationEffects = getNonTopEffects(permutationData, topPermutationEffects);
			getTranscriptionFactorCounts(topPermutationEffects, nonTopPermutationEffects, permutationCounts, genotypeMatrixData, regulomeDbData, r2Cutoff);
		}
		
		
		//STEP 7.: PERFORM THE FISHER EXACT TEST.
		getFisherPvalues(eQtlCounts, permutationCounts, outputFile);
	}
	
	
	/*
	 * =========================================================================
	 * = START: REGULOMEDB READING CODE.
	 * =========================================================================
	 */
	public HashMap<String, TreeMap<Integer, String[]>> readRegulomeDbData(ArrayList<RegulomeDbFile> regulomeDbFiles){
		HashMap<String, TreeMap<Integer, String[]>> regulomeDbData = new HashMap<String, TreeMap<Integer, String[]>>();
		TreeMap<Integer, String[]> tmp;
		
		RegulomeDbFiles regulomeDbFilesData = new RegulomeDbFiles(regulomeDbFiles);
		Iterator<RegulomeDbEntry> regulomeDbDataIterator = regulomeDbFilesData.iterator();
		int n = 0;
		
		while(regulomeDbDataIterator.hasNext()){
			RegulomeDbEntry rdbe = regulomeDbDataIterator.next();
			String rdbeChr = rdbe.getChr();
			int rdbePos = rdbe.getChrPos();
			String[] transcriptionFactors = getTranscriptionFactors(rdbe);
			
			if(transcriptionFactors.length > 0){
				if(regulomeDbData.containsKey(rdbeChr)){
					tmp = regulomeDbData.get(rdbeChr);
					tmp.put(rdbePos, transcriptionFactors);
				}
				else{
					tmp = new TreeMap<Integer, String[]>();
					tmp.put(rdbePos, transcriptionFactors);
					regulomeDbData.put(rdbeChr, tmp);
				}
				n++;
			}
		}
		return regulomeDbData;
	}
	
	
	public String[] getTranscriptionFactors(RegulomeDbEntry rdbe){
		Map<String, List<RegulomeDbSupportingData>> supportData = rdbe.getSupportData();
		ArrayList<String> tfs = new ArrayList<String>();
		
		Iterator it = supportData.entrySet().iterator();
		while(it.hasNext()){
			Map.Entry pairs = (Map.Entry) it.next();

			List<RegulomeDbSupportingData> aap = (List<RegulomeDbSupportingData>) pairs.getValue();
			for(RegulomeDbSupportingData rdbsd : aap){

				//Check if the annotation is protein_binding.
				if(rdbsd.getSupportClass().equalsIgnoreCase("Protein_Binding")){
					tfs.add(rdbsd.getSupportValue());
				}
			}
		}
		return tfs.toArray( new String[tfs.size()] );
	}
	/*
	 * =========================================================================
	 * = END: REGULOMEDB READING CODE.
	 * =========================================================================
	 */
	
	
	
	
	/*
	 * =========================================================================
	 * = START: EQTL PROCESSING METHODS.
	 * =========================================================================
	 */
	public EQTL[] readEQtlData(String eqtlFileLocation) throws IOException{
		QTLTextFile eqtlData = new QTLTextFile(eqtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	public HashMap<String, EQTL> getTopEffects(EQTL[] eqtls, HashSet<String> sharedProbes) throws IOException{
		HashMap<String, EQTL> eqtlMap = new HashMap<String, EQTL>();
		for(EQTL eqtl : eqtls){
			
			String eqtlProbe = eqtl.getProbe();
			if(sharedProbes.contains(eqtlProbe)){
				
				//Check if the eQTL for the shared probe has the highest absolute Z-Score.
				if( eqtlMap.containsKey(eqtl.getProbe()) ){
					EQTL tmpEqtl = eqtlMap.get(eqtl.getProbe());
					if( Math.abs(eqtl.getZscore()) > Math.abs(tmpEqtl.getZscore())){
						eqtlMap.put(eqtl.getProbe(), eqtl);
					}
				}
				else{
					eqtlMap.put(eqtl.getProbe(), eqtl);
				}
			}
		}
		return eqtlMap;
	}
	
	
	public HashMap<String, HashSet<Integer>> getNonTopEffects(EQTL[] eqtls, HashMap<String, EQTL> topEffects){
		HashMap<String, HashSet<Integer>> otherEffects = new HashMap<String, HashSet<Integer>>();
		HashSet<Integer> tmp;
		int n = 0;
		
		for(EQTL eqtl : eqtls){
			String probe = eqtl.getProbe();
			int pos = eqtl.getRsChrPos();
			
			if(topEffects.containsKey(probe)){
				EQTL topEqtl = topEffects.get(probe);
				int topPos = topEqtl.getRsChrPos();
				
				if(pos != topPos){
					if(otherEffects.containsKey(probe)){
						tmp = otherEffects.get(probe);
						tmp.add(pos);
					}
					else{
						tmp = new HashSet<Integer>();
						tmp.add(pos);
						otherEffects.put(probe, tmp);
					}
					n++;
				}
			}
		}
		return otherEffects;
	}
	
	
	public Set<String> makeRsIdList(EQTL[] eqtls){
		Set<String> rsIdList = new HashSet<String>();
		for(EQTL eqtl : eqtls){
			rsIdList.add(eqtl.getRsName());
		}
		return rsIdList;
	}
	
	
	
	public EQTL[] readPermutationData(String fileLocation, HashMap<String, GeneticVariant> rsIdMap)throws IOException{
		ArrayList<EQTL> permutations = new ArrayList<EQTL>();
		char a = '#';
		String fileLine;
		String[] fileLineData;
		
		TextFile tf = new TextFile(fileLocation, false);
		while((fileLine=tf.readLine())!=null){
			if(a != '#'){
				fileLineData = TAB_PATTERN.split(fileLine);
				double pval = Double.valueOf(fileLineData[0]);
				
				//Search the SNP back in the rsIdMap.
				if(rsIdMap.containsKey(fileLineData[1])){
					GeneticVariant gv = rsIdMap.get(fileLineData[1]);

					//Create the EQTL object
					EQTL eqtl = new EQTL();
					eqtl.setPvalue(pval);				
					eqtl.setRsChr(Byte.valueOf(gv.getSequenceName()));
					eqtl.setRsChrPos(gv.getStartPos());
					eqtl.setRsName(gv.getPrimaryVariantId());
					eqtl.setProbe(fileLineData[2]);
					
					if(pval < 1.0){
						eqtl.setZscore(ZScores.pToZ(pval));
					}
					else{
						eqtl.setZscore(0);
					}
					

					//Add the EQTL to the ArrayList
					permutations.add(eqtl);
					System.out.println(fileLineData[1] + " found.");
				}
			}
			else{
				a = '!';
			}
		}
		tf.close();
		return permutations.toArray(new EQTL[permutations.size()]);
	}
	/*
	 * =========================================================================
	 * = END: EQTL PROCESSING METHODS.
	 * =========================================================================
	 */
	
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 630000, new VariantIdIncludeFilter(variantIdFilter), new SampleIncludedFilter());
		return gonlImputedBloodGenotypeData;
	}
	
	
	public RandomAccessGenotypeData readEQtlGenotypeDataV2(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 630000, null, null);
		return gonlImputedBloodGenotypeData;
	}
	
	
	/*
	 * =========================================================================
	 * = START: TRANSCRIPTION FACTOR COUNTS CODE.
	 * =========================================================================
	 */
	public void getTranscriptionFactorCounts(HashMap<String, EQTL> topEffectData, HashMap<String, HashSet<Integer>> nonTopEffectData, HashMap<String, Integer> countsMap,
			RandomAccessGenotypeData genotypeData, HashMap<String, TreeMap<Integer, String[]>> regulomeDbData, double r2CutOff){
		Ld ld = null;
		Iterator<Map.Entry<String, EQTL>> topEffectIterator = topEffectData.entrySet().iterator();
		while(topEffectIterator.hasNext()){
			Map.Entry<String, EQTL> topEffectDataEntry = (Map.Entry) topEffectIterator.next();
			EQTL eqtl = topEffectDataEntry.getValue();
			
			//FETCH THE EQTL FROM GENOTYPE DATA.
			String rsChr = eqtl.getRsChr().toString();
			int rsChrPos = eqtl.getRsChrPos();
			String rsProbe = eqtl.getProbe();
			GeneticVariant eQtlVariant = genotypeData.getSnpVariantByPos(rsChr.toString(), rsChrPos);
			
			
			if(eQtlVariant != null){
				HashSet<Integer> nonTopEffectsPos = nonTopEffectData.get(rsProbe);
				
				if(nonTopEffectsPos != null){
					for(int eqtlPos : nonTopEffectsPos){
						GeneticVariant snpVariantByPos = genotypeData.getSnpVariantByPos(rsChr, eqtlPos);

						if(snpVariantByPos != null){
							if(eQtlVariant.isBiallelic() && snpVariantByPos.isBiallelic()){
								try {
									ld = eQtlVariant.calculateLd(snpVariantByPos);
								} catch (LdCalculatorException ex) {
									System.out.println("Error in LD calculation: " + ex.getMessage());
									System.exit(1);
								}

								if(ld.getR2() >= r2CutOff){
									//SEARCH THROUGH REGULOMEDB
									if(regulomeDbData.containsKey(rsChr)){
										TreeMap<Integer, String[]> regulomeChromosomeEntries = regulomeDbData.get(rsChr);
										
										if(regulomeChromosomeEntries.containsKey(eqtlPos)){
											String[] transcriptionFactors = regulomeChromosomeEntries.get(eqtlPos);
											for(String tf : transcriptionFactors){
												if(countsMap.containsKey(tf)){
													countsMap.put(tf, countsMap.get(tf)+1);
												}
												else{
													countsMap.put(tf, 1);
												}
											}
										}
									}

								}
							}
						}
					}
				}
				
				//CODE TO SEARCH FOR THE TOP EFFECT.
				if(regulomeDbData.containsKey(rsChr)){
					TreeMap<Integer, String[]> regulomeChromosomeEntries = regulomeDbData.get(rsChr);

					if(regulomeChromosomeEntries.containsKey(rsChrPos)){
						String[] transcriptionFactors = regulomeChromosomeEntries.get(rsChrPos);
						for(String tf : transcriptionFactors){
							if(countsMap.containsKey(tf)){
								countsMap.put(tf, countsMap.get(tf)+1);
							}
							else{
								countsMap.put(tf, 1);
							}
						}
					}
				}
				
			}
		}
	}
	/*
	 * =========================================================================
	 * = END: TRANSCRIPTION FACTOR COUNTS CODE.
	 * =========================================================================
	 */
	
	
	
	
	
	/*
	 * =========================================================================
	 * = START: FISHER EXACT TEST CODE.
	 * =========================================================================
	 */
	public int getTotalHits(HashMap<String, Integer> countsMap){
		int total = 0;
		
		Iterator<Map.Entry<String, Integer>> hitsIterator = countsMap.entrySet().iterator();
		while(hitsIterator.hasNext()){
			Map.Entry<String, Integer> next = hitsIterator.next();
			total += next.getValue();
		}
		return total;
	}
	
	
	public String getDirection(int eqtlHits, int eqtlTotalHits, int permutationHits, int permutationTotalHits){
		double eqtlRatio = ((double)eqtlHits / (double)eqtlTotalHits) * 100;
		double permutationRatio = ((double)permutationHits / (double)permutationTotalHits) * 100;
		
		if(eqtlRatio > permutationRatio){
			return "Enrichment";
		}
		
		else if( eqtlRatio == permutationRatio){
			return "None";
		}
		
		else{
			return "Impoverishment";
		}
	}
	
	
	public void getFisherPvalues(HashMap<String, Integer> eQtlCounts, HashMap<String, Integer> permutationCounts, String outputFile)throws IOException{
		int totalEQtlCounts = getTotalHits(eQtlCounts);
		int totalPermutationCounts = getTotalHits(permutationCounts);
		double bonferroniFactor = 0.05/eQtlCounts.size();
		
		PrintWriter fisherWriter = new PrintWriter(new FileWriter(outputFile));
		fisherWriter.println("#TF=Transcription Factor; FET=Fisher Exact Test P-value; BS=Bonferroni Significant?; DIR=Direction(Enrichment or not); ERA=EQTL Ratio; PRA=Permutation Ratio");
		fisherWriter.println("TF\tFET\tBS\tDIR\tERA\tPRA");
		for(Iterator<Map.Entry<String, Integer>>iter=eQtlCounts.entrySet().iterator();iter.hasNext();){
			Map.Entry<String, Integer> eQtlCountsEntry = iter.next();
			
			if( permutationCounts.containsKey(eQtlCountsEntry.getKey()) ){
				String tf = eQtlCountsEntry.getKey();
				int eQtlCount = eQtlCountsEntry.getValue();
				int permutationCount = permutationCounts.get( eQtlCountsEntry.getKey() );
				
				//Perform Fisher Exact test.		
				FisherExactTest fet = new FisherExactTest();
				double fisherPValue = fet.getFisherPValue(eQtlCount, (totalEQtlCounts - eQtlCount), permutationCount, (totalPermutationCounts - permutationCount));
				fisherWriter.println(tf + "\t" + fisherPValue + "\t" + (fisherPValue<=bonferroniFactor) + "\t"
						+ getDirection(eQtlCount, totalEQtlCounts, permutationCount, totalPermutationCounts) + "\t" + eQtlCount + "/" + totalEQtlCounts
						+ "\t" + permutationCount + "/" + totalPermutationCounts);
			}
		}
		fisherWriter.close();
	}
	/*
	 * =========================================================================
	 * = END: FISHER EXACT TEST CODE
	 * =========================================================================
	 */
	
	
	
	
	
	/*
	 * =========================================================================
	 * = START: SHARED PROBES CODE.
	 * =========================================================================
	 */
	public HashSet<String> makeProbesList(String dataFile)throws IOException{
		String fileLine;
		String[] fileLineData;
		char a = '#';
		HashSet<String> probesList = new HashSet<String>();
		
		TextFile tf = new TextFile(dataFile, false);
		while((fileLine=tf.readLine())!=null){
			if(a != '#'){
				fileLineData = TAB_PATTERN.split(fileLine);
				probesList.add(new String(fileLineData[4]));
			}
			else{
				a = '!';
			}
		}
		return probesList;
	}
	
	
	public HashSet<String> makeProbesListV2(String dataFile)throws IOException{
		String fileLine;
		String[] fileLineData;
		char a = '#';
		HashSet<String> probesList = new HashSet<String>();
		
		TextFile tf = new TextFile(dataFile, false);
		while((fileLine=tf.readLine())!=null){
			if(a != '#'){
				fileLineData = TAB_PATTERN.split(fileLine);
				probesList.add(new String(fileLineData[2]));
			}
			else{
				a = '!';
			}
		}
		return probesList;
	}
	
	
	public HashSet<String> getSharedProbes(String eQtlProbeFile, String permutationDataLocation)throws IOException{
		HashSet<String> eqtlProbes = makeProbesList(eQtlProbeFile);
		HashSet<String> permutationProbes = new HashSet<String>();
		for(int n=1;n<=100;n++){
			HashSet<String> tmpProbes = makeProbesListV2(permutationDataLocation + "PermutedEQTLsPermutationRound" + n + ".txt.gz");
			if(n > 1){
				permutationProbes.retainAll(tmpProbes);
			}
			else{
				permutationProbes = tmpProbes;
			}
		}
		
		if(eqtlProbes.size() > permutationProbes.size()){
			eqtlProbes.retainAll(permutationProbes);
			return eqtlProbes;
		}
		else{
			permutationProbes.retainAll(eqtlProbes);
			return permutationProbes;
		}
	}
	/*
	 * =========================================================================
	 * = END: SHARED PROBES CODE.
	 * =========================================================================
	 */
}
