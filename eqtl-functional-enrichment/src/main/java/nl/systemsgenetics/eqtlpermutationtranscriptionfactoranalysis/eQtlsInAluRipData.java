/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.genomicboundaries.GenomicBoundary;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author Matthieu
 */
public class eQtlsInAluRipData {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public static void main(String[] args)throws IOException{
		eQtlsInAluRipData eqiard = new eQtlsInAluRipData("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\eQTLsFDR0.05.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\PermutedEQTLsPermutationRound8.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\BloodHT12Combined\\", "C:\\Users\\Matthieu\\Downloads\\dbrip\\Alu_hg18_v2h.txt", 0.8);
	}
	
	public eQtlsInAluRipData(String eQtlFile, String permutationLocation, String genotypeLocation, String aluRipData, double r2cutoff) throws IOException{
		//STEP 1.: READ THE ALU RIP DATA.
		System.out.println("[O]:> Reading the Alu DBRIP data.");
		GenomicBoundaries<Object> aluRipRegions = readAluRipData(aluRipData);
		
		
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
		
		
		//STEP 4.: FIND SNPS IN LD AND COUNT THE HITS.
		System.out.println("[E]:> Find LD SNPs and collect counts.");
		HashMap<String, Integer> aluCounts = new HashMap<String, Integer>();
		getAluCounts(topEQtlEffects, nonTopEqtlEffects, aluCounts, genotypeData, aluRipRegions, 0.8);
		System.out.println("[E]:> Found " + aluCounts.size() + " hits.");
		
		
		
		//STEP 5.: PERFORM THE OPERATION FOR THE PERMUTATION DATA.
		System.out.println("[O]:> Start the analysis for the permutation data.");
		HashMap<String, Integer> permutationCounts = new HashMap<String, Integer>();
		for(int n=1;n<=100;n++){
			EQTL[] permutationData;
			HashMap<String, EQTL> topPermutationData;
			HashMap<String, HashSet<Integer>> nonTopPermutationEffects;
			
			permutationData = readEQtlData(permutationLocation + "PermutedEQTLsPermutationRound" + n + ".txt.gz");
			topPermutationData = getTopEQtlEffects(permutationData);
			nonTopPermutationEffects = getNonTopEqtlEffects(permutationData, topPermutationData);
			getAluCounts(topPermutationData, nonTopPermutationEffects, permutationCounts, genotypeData, aluRipRegions, r2cutoff);
		}
		
		
		//STEP 6.: GET THE FISHER EXACT VALUES.
		System.out.println("[O]:> Get the Fisher exact pvalues.");
		getFisherPvalues(aluCounts, permutationCounts);
	}
	
	
	/*
	 * =========================================================================
	 * = START: ALU DBRIP DATA PROCESSING METHODS.
	 * =========================================================================
	 */
	public GenomicBoundaries<Object> readAluRipData(String aluFile)throws IOException{
		GenomicBoundaries<Object> ripRegions = new GenomicBoundaries();
		
		String fileLine;
		String[] fileLineData;
		
		BufferedReader aluRipReader = new BufferedReader(new FileReader(aluFile));
		while((fileLine=aluRipReader.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			String chr = new String(fileLineData[1]);
			int startPos = Integer.parseInt(fileLineData[2]);
			int stopPos = Integer.parseInt(fileLineData[3]);
			String id = new String(fileLineData[7]);
			String elementClass = new String(fileLineData[10]);
			String elementFamily = new String(fileLineData[11]);
			String elementSubFamily = new String(fileLineData[12]);
			String disease = new String(fileLineData[21]);
			
			//String supportData = id+"|"+elementClass+"|"+elementFamily+"|"+elementSubFamily+"|"+disease;
			ripRegions.addBoundary(chr, startPos, stopPos, elementClass);
		}
		aluRipReader.close();
		return ripRegions;
	}
	/*
	 * =========================================================================
	 * = END: ALU DBRIP DATA PROCESSING METHODS.
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
		return otherEffects;
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
	
	
	/*
	 * =========================================================================
	 * = START: ANALYSIS CODE.
	 * =========================================================================
	 */
	public void getAluCounts(HashMap<String, EQTL> topEffectData, HashMap<String, HashSet<Integer>> nonTopEffectData, HashMap<String, Integer> countsMap,
			RandomAccessGenotypeData genotypeData, GenomicBoundaries<Object> aluRipRegions, double r2CutOff){
		
		Ld ld = null;
		Iterator<Map.Entry<String, EQTL>> topEffectIterator = topEffectData.entrySet().iterator();
		while(topEffectIterator.hasNext()){
			Map.Entry<String, EQTL> topEffectDataEntry = (Map.Entry) topEffectIterator.next();
			EQTL eqtl = topEffectDataEntry.getValue();
			
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
									if(aluRipRegions.isInBoundary(rsChr, snpVariantByPos.getStartPos(), 0)){
										GenomicBoundary boundary = aluRipRegions.getBoundary(rsChr, snpVariantByPos.getStartPos(), 0);
										String annotation = (String) boundary.getAnnotation();
										if(countsMap.containsKey(annotation)){
											countsMap.put(annotation, countsMap.get(annotation)+1);
										}
										else{
											countsMap.put(annotation, 1);
										}
									}

								}
							}
						}
					}
				}
			}
		}
	}
	
	
	
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
				double fisherPValue = fet.getFisherPValue(eQtlCount, (totalEQtlCounts - eQtlCount), permutationCount, (totalPermutationCounts - permutationCount));
				//fisherPValues.put(tf, fisherPValue);
				System.out.println(tf + "\t" + fisherPValue);
			}
		}
	}
	/*
	 * =========================================================================
	 * = END: FISHER EXACT TEST CODE.
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
	
	
	public HashSet<String> getSharedProbes(String eQtlProbeFile, String permutationDataLocation)throws IOException{
		HashSet<String> eqtlProbes = makeProbesList(eQtlProbeFile);
		HashSet<String> permutationProbes = new HashSet<String>();
		for(int n=1;n<=100;n++){
			HashSet<String> tmpProbes = makeProbesList(permutationDataLocation + "PermutedEQTLsPermutationRound" + n + ".txt.gz");
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
