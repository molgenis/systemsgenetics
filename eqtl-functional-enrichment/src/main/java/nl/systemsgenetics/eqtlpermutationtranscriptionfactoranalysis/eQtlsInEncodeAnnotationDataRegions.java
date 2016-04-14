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
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.genomicboundaries.GenomicBoundary;
import umcg.genetica.io.gtf.GffElement;
import umcg.genetica.io.gtf.GtfReader;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author Matthieu
 */
public class eQtlsInEncodeAnnotationDataRegions {
	/*
	 * This class is ment to check which eQTLs are located in which regions.
	 * For this analysis LD SNPs with eQTL effects will be searched for each
	 *	top eQTL effect.
	 * Each LD SNP will be checked in which encode region it's located. The hits
	 * will be reported. Also counts for all the different classes will be kept
	 * as well as the amount and percentage of total LD SNPs that mapped in any
	 * encode annotation region.
	 */
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public static void main(String[] args)throws IOException {
		eQtlsInEncodeAnnotationDataRegions aap = new eQtlsInEncodeAnnotationDataRegions("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\eQTLsFDR0.05.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\gencode\\gencode.v18.annotation.b36.gtf");
	}
	
	public eQtlsInEncodeAnnotationDataRegions(String eQtlFile, String gencodeFile) throws IOException{
		EQTL[] eqtls = readEQtlData(eQtlFile);
		Set<String> rsIdList = makeRsIdList(eqtls);
		HashMap<String, EQTL> topEQtlEffects = getTopEQtlEffects(eqtls);
		HashMap<String, HashSet<Integer>> nonTopEqtlEffects = getNonTopEqtlEffects(eqtls, topEQtlEffects);
		
		//RandomAccessGenotypeData genotypeData = readEQtlGenotypeData(genotypeLocation, rsIdList);
		
		GenomicBoundaries<Object> encodeData = readEncodeAnnotationData(gencodeFile);
		HashMap<String, Integer> eQtlInEncodeRegions = eQtlInEncodeRegions(eqtls, encodeData);
		printCountsMap(eQtlInEncodeRegions);
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
	
	
	public HashMap<String, Integer> eQtlInEncodeRegions(EQTL[] eqtls, GenomicBoundaries<Object> boundaries){
		HashMap<String, Integer> encodeCounts = new HashMap<String, Integer>();
		for(EQTL eqtl : eqtls){
			if(boundaries.isInBoundary(eqtl.getRsChr().toString(), eqtl.getRsChrPos(), 0)){
				GenomicBoundary boundary = boundaries.getBoundary(eqtl.getRsChr().toString(), eqtl.getRsChrPos(), 0);
				
				String transcriptType = (String) boundary.getAnnotation();
				if(encodeCounts.containsKey(transcriptType)){
					encodeCounts.put(transcriptType, encodeCounts.get(transcriptType)+1);
				}
				else{
					encodeCounts.put(transcriptType, 1);
				}
			}
		}
		return encodeCounts;
	}
	
	
	
	
	public void printCountsMap(HashMap<String, Integer> counts){
		Iterator<Map.Entry<String, Integer>> countsIter = counts.entrySet().iterator();
		int n = 1;
		int totalHits = 0;
		System.out.println("Entries in count map of size " + counts.values().size() + ".");
		while(countsIter.hasNext()){
			Map.Entry<String, Integer> next = countsIter.next();
			System.out.println(n + ". " + next.getKey() + ": " + next.getValue());
			totalHits += next.getValue();
			n++;
		}
		System.out.println("Total hits: " + totalHits);
	}
}
