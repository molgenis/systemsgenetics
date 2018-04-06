/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
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
import umcg.genetica.io.gtf.GffElement;
import umcg.genetica.io.gtf.GtfReader;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author Matthieu
 */
public class eQtlsInEncodeAnnotationData {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public eQtlsInEncodeAnnotationData(String eQtlProbeLevelFile, String eQtlFile, String permutationLocation, String genotypeLocation, String gencodeFile, double r2cutoff, String outputFile) throws IOException{
		//STEP 1.: READ THE ENCODE DATA.
		GenomicBoundaries<Object> encodeData = readEncodeAnnotationData(gencodeFile);
		
		//STEP 2.: GET THE SHARED PROBES.
		HashSet<String> sharedProbes = getSharedProbes(eQtlProbeLevelFile, permutationLocation);
		
		//STEP 3.: READING AND FILTERING EQTL DATA.
		EQTL[] eqtls = readEQtlData(eQtlFile);
		Set<String> rsIdList = makeRsIdList(eqtls);
		HashMap<String, EQTL> topEQtlEffects = getTopEffects(eqtls, sharedProbes);
		HashMap<String, HashSet<Integer>> nonTopEqtlEffects = getNonTopEffects(eqtls, topEQtlEffects);
		
		//STEP 4.: READING GENOTYPE MATRIX DATA AND ENCODE DATA SOURCES.
		RandomAccessGenotypeData genotypeData = readEQtlGenotypeDataV2(genotypeLocation, rsIdList);
		
		//STEP 5.: GET THE ENCODE REGION COUNTS FOR THE REAL EQTL DATA.
		HashMap<String, Integer> eqtlCounts = new HashMap<String, Integer>();
		getEncodeRegionCounts(topEQtlEffects, nonTopEqtlEffects, eqtlCounts, genotypeData, encodeData, r2cutoff);
		
		
		//STEP 6.: PERFORM READING, FILTERING AND ANALYSIS FOR PERMUTATION DATA.
		HashMap<String, Integer> permutationCounts = new HashMap<String, Integer>();
		for(int n=1;n<=100;n++){
			EQTL[] permutationData;
			HashMap<String, EQTL> topPermutationData;
			HashMap<String, HashSet<Integer>> nonTopPermutationEffects;
			
			permutationData = readEQtlData(permutationLocation + "PermutedEQTLsPermutationRound" + n + ".txt.gz");
			topPermutationData = getTopEffects(permutationData, sharedProbes);
			nonTopPermutationEffects = getNonTopEffects(permutationData, topPermutationData);
			
			getEncodeRegionCounts(topPermutationData, nonTopPermutationEffects, permutationCounts, genotypeData, encodeData, r2cutoff);
		}
		
		//STEP 7.: PERFORM FISHER EXACT TEST.
		getFisherPvalues(eqtlCounts, permutationCounts, outputFile);
	}
	
	
	
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
	/*
	 * =========================================================================
	 * = END: EQTL PROCESSING METHODS.
	 * =========================================================================
	 */
	
	
	
	
	
	/*
	 * =========================================================================
	 * = START: GENOTYPE DATA PROCESSING METHODS.
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
	 * = END: GENOTYPE DATA PROCESSING METHODS.
	 * =========================================================================
	 */
	
	
	
	
	/*
	 * =========================================================================
	 * = START: ENCODE ANNOTATION DATA PROCESSING METHODS.
	 * =========================================================================
	 */
	public GenomicBoundaries<Object> readEncodeAnnotationData(String encodeFile)throws IOException{
		GenomicBoundaries<Object> encodeBoundaries = new GenomicBoundaries();
		
		GtfReader gtfr = new GtfReader(new File(encodeFile));
		Iterator<GffElement> encodeIterator = gtfr.iterator();
		while(encodeIterator.hasNext()){
			GffElement encodeEntry = encodeIterator.next();
			
			encodeBoundaries.addBoundary(encodeEntry.getSeqname(), encodeEntry.getStart(), encodeEntry.getEnd(),
					encodeEntry.getAttributeValue("transcript_type"));
		}
		gtfr.close();
		return encodeBoundaries;
	}
	/*
	 * =========================================================================
	 * = END: ENCODE ANNOTATION DATA PROCESSING METHODS.
	 * =========================================================================
	 */
	
	
	
	
	
	/*
	 * =========================================================================
	 * = START: ENCODE REGION COUNTS CODE.
	 * =========================================================================
	 */
	public void getEncodeRegionCounts(HashMap<String, EQTL> topEffectData, HashMap<String, HashSet<Integer>> nonTopEffectData, HashMap<String, Integer> countsMap,
			RandomAccessGenotypeData genotypeData, GenomicBoundaries<Object> boundaries, double r2CutOff){
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
									if(boundaries.isInBoundary(rsChr, eqtlPos, 0)){
										GenomicBoundary boundary = boundaries.getBoundary(rsChr, eqtlPos, 0);
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
				
				
				//CODE TO SEARCH FOR THE TOP EFFECT.
				if(boundaries.isInBoundary(rsChr, rsChrPos, 0)){
					GenomicBoundary boundary = boundaries.getBoundary(rsChr, rsChrPos, 0);
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
	/*
	 * =========================================================================
	 * = END: ENCODE REGIONS COUNTS CODE.
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
		double eqtlRatio = ((double) eqtlHits / (double) eqtlTotalHits) * 100;
		double permutationRatio = ((double) permutationHits / (double) permutationTotalHits) * 100;
		
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
		fisherWriter.println("#TF=Transcription Factor; FET=Fisher Exact Test P-value; BS=Bonferroni Significant?; DIR=Direction(Enrichment or not); ERA=EQTL Ratio; PRA+Permutation Ratio");
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
						+ getDirection(eQtlCount, totalEQtlCounts, permutationCount, totalPermutationCounts) + "\t" + eQtlCount + "/" + totalEQtlCounts +
						"\t" + permutationCount + "/" + totalPermutationCounts);
			}
		}
		fisherWriter.close();
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
