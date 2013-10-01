/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.io.regulomedb.RegulomeDbEntry;
import umcg.genetica.io.regulomedb.RegulomeDbFile;
import umcg.genetica.io.regulomedb.RegulomeDbFiles;
import umcg.genetica.io.regulomedb.RegulomeDbSupportingData;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author Matthieu
 */
public class EQtlPermutationTranscriptionFactorAnalysis {
	
	public static void main(String[] args)throws IOException{
		/*
		 * args[0]: eQTL file.
		 * args[1]: permutation files locations.
		 * args[2]: genotype matrix location.
		 * args[3]: regulomedb location.
		 * args[4]: window size.
		 * args[5]: r2 cutoff.
		 */
		//EQtlPermutationTranscriptionFactorAnalysisV3 eqptfa3 = new EQtlPermutationTranscriptionFactorAnalysisV3(args[1], args[2], args[3], args[4], Double.valueOf(args[5]));
//		EQtlPermutationTranscriptionFactorAnalysisV3 eqptfa3 =
//				new EQtlPermutationTranscriptionFactorAnalysisV3("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\eQTLsFDR0.05.txt",
//				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\PermutedEQTLsPermutationRound8.txt",
//				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\BloodHT12Combined\\","C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\regulomeDb\\",
//				250000, 0.8);
		if(args[0].equals("tf")){
			EQtlPermutationTranscriptionFactorAnalysisV3 eqptfa3 = new EQtlPermutationTranscriptionFactorAnalysisV3(args[1], args[2], args[3], args[4], Double.valueOf(args[5]));
		}
		else if(args[0].equals("gencode")){
			//Make call to gencode analysis.
		}
		else if(args[0].equals("dbrip")){
			eQtlsInAluRipData eqiard = new eQtlsInAluRipData(args[1], args[2], args[3], args[4], Double.valueOf(args[5]));
		}
		else if(args[0].equals("repeats")){
			//Make call to perform analysis on repeat data from dasha.
		}
		else{
			System.out.println("Use one of the following commands: ");
		}
	}
	
	
	
	public EQtlPermutationTranscriptionFactorAnalysis(String eQtlFile, String permutationFile, String genotypeLocation, String regulomeDbLocation, int windowSize, double r2Cutoff)throws IOException{
		ArrayList<RegulomeDbFile> regulomeDbFiles = new ArrayList<RegulomeDbFile>();
		regulomeDbFiles.add(new RegulomeDbFile( new File(regulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category1.txt") ));
		
		
		//STEP 1.: READ REGULOMEDB DATA.
		System.out.println("[O]:> Read RegulomeDB data.");
		//HashMap<String, TreeMap<Integer, ArrayList<RegulomeDbEntry>>> regulomeDbData = readRegulomeDbData(regulomeDbFiles);
		HashMap<String, TreeMap<Integer, String[]>> regulomeDbData = readRegulomeDbData(regulomeDbFiles);
		
		
		//STEP 2.: READ THE EQTL DATA.
		System.out.println("[E]:> Read the eQTL data.");
		EQTL[] eqtls = readEQtlData(eQtlFile);
		
		System.out.println("[E]:> Make a list of all eQTL rsIDs.");
		Set<String> rsIdList = makeRsIdList(eqtls);
		
		System.out.println("[E]:> Filter out the top effects per probe.");
		HashMap<String, EQTL> topEQtlEffects = getTopEQtlEffects(eqtls);
		
		System.out.println("[E]:> Get the non top effects.");
		HashMap<String, HashSet<Integer>> nonTopEqtlEffects = getNonTopEqtlEffects(eqtls, topEQtlEffects);
		
		
		//STEP 3.: READ THE GENOTYPE DATA.
		System.out.println("[O]:> Read Genotype matrix data.");
		RandomAccessGenotypeData genotypeData = readEQtlGenotypeData(genotypeLocation, rsIdList);
		
		
		//STEP 4.: FIND LD SNPS FOR THE TOP EQTL EFFECTS AND FIND REGULOMEDB HITS.
		HashMap<String, Integer> eqtlCounts = new HashMap<String, Integer>();
		System.out.println("[E]:> Made the counts map of size " + eqtlCounts.size() + " and will search through the data.");
		getTranscriptionFactorCounts(topEQtlEffects, nonTopEqtlEffects, eqtlCounts, genotypeData, regulomeDbData, windowSize, r2Cutoff);
		System.out.println("Performed the main analysis and the counts map is now size " + eqtlCounts.size() + ".");
		
		
		
		//STEP 5.: PERFORM THE OPERATION FOR THE PERMUTATION DATA.
		System.out.println("[P]:> Reading the permutation data.");
		EQTL[] permutationData = readEQtlData(permutationFile);
		System.out.println("[P]:> Read " + permutationData.length + " permutation data elements.");
		
		System.out.println("[P]:> Filter top permutation data.");
		HashMap<String, EQTL> topPermutationData = getTopEQtlEffects(permutationData);
		System.out.println("[P]:> Filtered data contains " + topPermutationData.size() + " elements.");
		
		System.out.println("[P]:> Filter the non top permutation data.");
		HashMap<String, HashSet<Integer>> nonTopPermutationEffects = getNonTopEqtlEffects(permutationData, topPermutationData);
		System.out.println("[P]:> Done filtering. Contains " + nonTopPermutationEffects.size() + " probes.");
		
		System.out.println("[P]:> Get the transcription factors.");
		HashMap<String, Integer> permutationCounts = new HashMap<String, Integer>();
		getTranscriptionFactorCounts(topPermutationData, nonTopPermutationEffects, permutationCounts, genotypeData, regulomeDbData, 250000, 0.8);
		
		
		//STEP 6.: GET THE FISHER EXACT VALUES.
		System.out.println("[O]:> Get the Fisher exact pvalues.");
		getFisherPvalues(eqtlCounts, permutationCounts);
		
	}
	
	
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
				//System.out.println("First TF: " + transcriptionFactors[0]);
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
		//String[] tfsArray = new String[tfs.size()];
		//tfsArray = tfs.toArray(tfsArray);
		//return tfsArray;
		return tfs.toArray( new String[tfs.size()] );
	}
	
	
	
	
	
	/*
	 * =========================================================================
	 * = METHODS DONE HERE WORK CORRECT AND ACCORDING PATRICK PLAN.
	 * =========================================================================
	 */
	
	
	/*
	 * =========================================================================
	 * = START: EQTL PROCESSING METHODS.
	 * =========================================================================
	 */
	public EQTL[] readEQtlData(String eqtlFileLocation) throws IOException{
		eQTLTextFile eqtlData = new eQTLTextFile(eqtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	public HashMap<String, EQTL> getTopEQtlEffects(EQTL[] eqtls) throws IOException{
		HashMap<String, EQTL> eqtlMap = new HashMap<String, EQTL>();
		for(EQTL eqtl : eqtls){
			
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
		return eqtlMap;
	}
	
	public Set<String> makeRsIdList(EQTL[] eqtls){
		Set<String> rsIdList = new HashSet<String>();
		for(EQTL eqtl : eqtls){
			rsIdList.add(eqtl.getRsName());
		}
		return rsIdList;
	}
	
	
	public HashMap<String, HashSet<Integer>> getNonTopEqtlEffects(EQTL[] eqtls, HashMap<String, EQTL> topEffects){
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
		System.out.println("Total hits: " + n);
		return otherEffects;
	}
	/*
	 * =========================================================================
	 * = END: EQTL PROCESSING METHODS.
	 * =========================================================================
	 */
	
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 1000, new VariantIdIncludeFilter(variantIdFilter), new SampleIncludedFilter());
		return gonlImputedBloodGenotypeData;
	}
	
	
	//Provide empty countsMap as parameter for eQTLs.
	public void getTranscriptionFactorCounts(HashMap<String, EQTL> topEffectData, HashMap<String, HashSet<Integer>> nonTopEffectData, HashMap<String, Integer> countsMap,
			RandomAccessGenotypeData genotypeData, HashMap<String, TreeMap<Integer, String[]>> regulomeDbData, int windowSize, double r2CutOff){
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
			
			
			
			//CHECK IF THE GOTTEN VARIANT IS NOT NULL.
			if(eQtlVariant != null){
				int boundaryL = eQtlVariant.getStartPos() - windowSize;
				int boundaryR = eQtlVariant.getStartPos() + windowSize;
				
				
				//GET THE VARIANTS BY RANGE AND LOOP THROUGH THEM.
				Iterable<GeneticVariant> variantsByRange = genotypeData.getVariantsByRange(rsChr, boundaryL, boundaryR);
				for(GeneticVariant gv : variantsByRange){
					
					if(nonTopEffectData.containsKey(rsProbe)){
						HashSet<Integer> get = nonTopEffectData.get(rsProbe);
						//System.out.println("LD SNP " + gv.getPrimaryVariantId() + " says: le francais can suck my dick!");
						
						if(get.contains(gv.getStartPos())){
							
							if(eQtlVariant.isBiallelic() && gv.isBiallelic()){
								try {
									ld = eQtlVariant.calculateLd(gv);
								} catch (LdCalculatorException ex) {
									System.out.println("Error in LD calculation: " + ex.getMessage());
									System.exit(1);
								}


								//CHECK IF R2 IS BIGGER THAN THE SET CUTOFF
								if(ld.getR2() >= r2CutOff){
									//SEARCH THROUGH REGULOMEDB
									if(regulomeDbData.containsKey(rsChr)){
										TreeMap<Integer, String[]> regulomeChromosomeEntries = regulomeDbData.get(rsChr);
										Iterator<Entry<Integer, String[]>> regulomeChromosomeIterator = regulomeChromosomeEntries.entrySet().iterator();
										
										//PERFORM A COUNT OPERATION FOR EACH TF.
										while(regulomeChromosomeIterator.hasNext()){
											Entry<Integer, String[]> next = regulomeChromosomeIterator.next();
											for(String aap : next.getValue()){
												if(countsMap.containsKey(aap)){
													countsMap.put(aap, countsMap.get(aap)+1);
												}
												else{
													countsMap.put(aap, 1);
												}
											}
										}
									}
									
								}
								//PERFORMED THE CHECK IF THE R2 IS BIGGER THAN THE SET CUTOFF.
							}
							
						}
					}
					
					
				}
				
			}
			
		}
	}
	
	
	public void printCountsMap(HashMap<String, Integer> counts){
		Iterator<Map.Entry<String, Integer>> countsIter = counts.entrySet().iterator();
		int n = 1;
		System.out.println("Entries in count map of size " + counts.values().size() + ".");
		while(countsIter.hasNext()){
			Map.Entry<String, Integer> next = countsIter.next();
			System.out.println(n + ". " + next.getKey() + ": " + next.getValue());
			n++;
		}
	}
	
	
	public int getTotalHits(HashMap<String, Integer> countsMap){
		int total = 0;
		
		Iterator<Map.Entry<String, Integer>> hitsIterator = countsMap.entrySet().iterator();
		while(hitsIterator.hasNext()){
			Map.Entry<String, Integer> next = hitsIterator.next();
			total += next.getValue();
		}
		return total;
	}
	
	
	public void getFisherPvalues(HashMap<String, Integer> eQtlCounts, HashMap<String, Integer> permutationCounts){
		//HashMap<String, Double> fisherPValues = new HashMap<String, Double>();
		int totalEQtlCounts = getTotalHits(eQtlCounts);
		int totalPermutationCounts = getTotalHits(permutationCounts);
		
		for(Iterator<Map.Entry<String, Integer>>iter=eQtlCounts.entrySet().iterator();iter.hasNext();){
			Map.Entry<String, Integer> eQtlCountsEntry = iter.next();
			
			if( permutationCounts.containsKey(eQtlCountsEntry.getKey()) ){
				String tf = eQtlCountsEntry.getKey();
				int eQtlCount = eQtlCountsEntry.getValue();
				int permutationCount = permutationCounts.get( eQtlCountsEntry.getKey() );
				
				//Perform Fisher Exact test.
				FisherExactTest fet = new FisherExactTest();
				double fisherPValue = fet.getFisherPValue(eQtlCount, totalEQtlCounts, permutationCount, totalPermutationCounts);
				//fisherPValues.put(tf, fisherPValue);
				System.out.println(tf + "\t" + fisherPValue);
			}
		}
		
	}
	
}
