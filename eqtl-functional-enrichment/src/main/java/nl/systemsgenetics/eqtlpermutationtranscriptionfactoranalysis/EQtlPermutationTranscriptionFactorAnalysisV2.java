/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.File;
import java.io.FileNotFoundException;
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
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author Matthieu
 */
public class EQtlPermutationTranscriptionFactorAnalysisV2 {
	
	
	public static void main(String[] args)throws IOException{
		EQtlPermutationTranscriptionFactorAnalysisV2 eqptfa2 = new EQtlPermutationTranscriptionFactorAnalysisV2(args[0], args[1], args[2], args[3], Integer.parseInt(args[4]), Double.valueOf(args[5]));
	}
	
	
	public EQtlPermutationTranscriptionFactorAnalysisV2(String eQtlFile, String permutationFile, String genotypeData, String regulomeDbLocation, int windowSize, double r2CutOff)throws IOException{
		ArrayList<RegulomeDbFile> regulomeDbFiles = new ArrayList<RegulomeDbFile>();
		regulomeDbFiles.add(new RegulomeDbFile( new File(regulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category1.txt") ));
		
		/*
		 * args[0]: regulomeDbFile
		 * args[1]: eQtlFile
		 * args[2]: location of genotype data
		 * args[3]: location of permutation files
		 * 
		 * args[0]: eQTL file
		 * args[1]: permutation file
		 * args[2]: GenotypeMatrix location
		 * args[3]: RegulomeDB location
		 * args[4]: window size
		 * args[5]: r2 cutoff
		 */
		
		//Step 1.: Read RegulomeDB data.
		System.out.println("[O]:> Read RegulomeDB data.");
		HashMap<String, TreeMap<Integer, RegulomeDbEntry>> regulomeDbData = readRegulomeDbData(regulomeDbFiles);
		
		/*
		 * =====================================================================
		 * = READ AND DO THINGS FOR THE REAL EQTL DATA.
		 * =====================================================================
		 */
		//Step 2.: Read and filter eQTL data.
		System.out.println("[E]:> Reading and filtering the eQTL data.");
		EQTL[] eqtls = readEQtlResultData(eQtlFile);
		HashMap<String, EQTL> topEqtls = filterEQtlSet(eqtls);
		System.out.println("[E]:> Reading and filtering of the eQTL data succesfull.");
		
		
		//Step 3.: Read and filter permutation data.
		System.out.println("[P]:> Reading and filtering the Permutation data.");
		EQTL[] permutationData = readEQtlResultData(permutationFile);
		HashMap<String, EQTL> topPermutations = filterEQtlSet(permutationData);
		System.out.println("[P]:> Reading and filtering of the permutation data succesfull.");
		
		
		//Step 4.: Make rsID list of the eQTLs and read the genotype data.
		System.out.println("[O]:> Reading GenotypeMatrix data.");
		Set<String> rsIdList = makeRsIdList(eqtls);
		RandomAccessGenotypeData readEQtlGenotypeData = readEQtlGenotypeData(genotypeData, rsIdList);
		System.out.println("[O]:> Reading of Genotype data succesfull.");
		
		
		//Step 5.: Filter both sets.
		System.out.println("[O]:> Filter both sets on shared probes.");
		filterOnSharedProbes(topEqtls, topPermutations);
		filterOnSharedProbes(topPermutations, topEqtls);
		System.out.println("[O]:> Filtering of both sets on shared probes complete.");
		
		
		//Step 6.: Search for LD SNPs and get Counts for each eQTL.
		System.out.println("[O]:> Collect the counts for both sets.");
		HashMap<String, Integer> eQtlCounts = getEQtlCounts(topEqtls, readEQtlGenotypeData, regulomeDbData, windowSize, r2CutOff);
		HashMap<String, Integer> permutationCounts = getEQtlCounts(topPermutations, readEQtlGenotypeData, regulomeDbData, 250000, 0.8);
		System.out.println("[O]:> Collect the counts complete..");
		
		
		//Step 7.: Get Fisher Exact Counts.
		System.out.println("[O]:> Fisher Exact P-values:");
		getFisherPvalues(eQtlCounts, permutationCounts);
		
		//printCountsMap(eQtlCounts);
		/*
		 * =====================================================================
		 * = READ AND DO THINGS FOR THE REAL EQTL DATA.
		 * =====================================================================
		 */
		
		
		//Step 5.: Perform analysis for permutation data.
//		HashMap<String, Integer> permutationCounts = new HashMap<String, Integer>();
//		for(int n=1;n<=101;n++){
//			//EQTL[] permutationEQtls = readEQtlResultData("");
//			//HashMap<String, EQTL> topPermutations = filterEQtlSet(permutationEQtls);
//			//filterOnSharedProbes(topEqtls, topPermutations);
//		}
	}
	
	
	public HashMap<String, TreeMap<Integer, RegulomeDbEntry>> readRegulomeDbData(ArrayList<RegulomeDbFile> regulomeDbFiles)throws IOException{
		HashMap<String, TreeMap<Integer, RegulomeDbEntry>> regulomeDbData = new HashMap<String, TreeMap<Integer, RegulomeDbEntry>>();
		TreeMap<Integer, RegulomeDbEntry> tmp; 
		
		RegulomeDbFiles regulomeDbFilesData = new RegulomeDbFiles(regulomeDbFiles);
		Iterator<RegulomeDbEntry> regulomeDbDataIterator = regulomeDbFilesData.iterator();
		
		while(regulomeDbDataIterator.hasNext()){
			RegulomeDbEntry rdbe = regulomeDbDataIterator.next();
			String rdbeChr = rdbe.getChr();
			int rdbePos = rdbe.getChrPos();
			
			if(regulomeDbData.containsKey(rdbeChr)){
				tmp = regulomeDbData.get(rdbeChr);
				tmp.put(rdbePos, rdbe);
			}
			
			else{
				tmp = new TreeMap<Integer, RegulomeDbEntry>();
				tmp.put(rdbePos, rdbe);
				regulomeDbData.put(rdbeChr, tmp);
			}
		}
		
		return regulomeDbData;
	}
	
	public EQTL[] readEQtlResultData(String eqtlFileLocation) throws IOException{
		QTLTextFile eqtlData = new QTLTextFile(eqtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 1000, new VariantIdIncludeFilter(variantIdFilter), new SampleIncludedFilter());
		return gonlImputedBloodGenotypeData;
	}
	
	
	public HashMap<String, Integer> getEQtlCounts(HashMap<String, EQTL> topEqtls, RandomAccessGenotypeData genotypeData, HashMap<String, TreeMap<Integer, RegulomeDbEntry>> regulomeDbData, int windowSize, double r2Cutoff){
		HashMap<String, Integer> eQtlCounts = new HashMap<String, Integer>();
		Ld ld = null;
		
		//System.out.println("EQTL\tLD_SNP\tTF");
		Iterator<Entry<String, EQTL>> eQtlIterator = topEqtls.entrySet().iterator();
		while(eQtlIterator.hasNext()){
			Map.Entry<String, EQTL> pairs = (Map.Entry) eQtlIterator.next();
			EQTL eqtl = pairs.getValue();
			
			//FETCH THE EQTL FROM THE GENOTYPE DATA.
			Byte rsChr = eqtl.getRsChr();
			int rsChrPos = eqtl.getRsChrPos();
			GeneticVariant eQtlVariant = genotypeData.getSnpVariantByPos(rsChr.toString(), rsChrPos);
			
			
			if(eQtlVariant != null){
				int boundaryL = eQtlVariant.getStartPos() - windowSize;
				int boundaryR = eQtlVariant.getStartPos() + windowSize;

				//GET THE VARIANTS BY RANGE
				Iterable<GeneticVariant> variantsByRange = genotypeData.getVariantsByRange(rsChr.toString(), boundaryL, boundaryR);
				for(GeneticVariant gv : variantsByRange){
					try {
						ld = eQtlVariant.calculateLd(gv);
					} catch (LdCalculatorException ex) {
						System.out.println("Error in LD calculation: " + ex.getMessage());
						System.exit(1);
					}


					//CHECK IF R2 IS BIGGER THAN THE SET CUTOFF
					if(ld.getR2() >= r2Cutoff){
						/*
						 * TODO: Check if the LD variant is in the GENCODE
						 * 1. Check if the chr is there.
						 * 2. Get the treemap.
						 * 3. Look if the position of the LD snp is there as a key in the regulomedb data.
						 * 4. if so get the element, get the Protein_Binding and add it to the counts map.
						 */

						if(regulomeDbData.containsKey(gv.getSequenceName())){
							TreeMap<Integer, RegulomeDbEntry> regulomeDbChrEntries = regulomeDbData.get(gv.getSequenceName());

							if(regulomeDbChrEntries.containsKey(gv.getStartPos())){
								//GET REGULOMEDB ELEMENT AND ADD TO COUNTS MAP.
								RegulomeDbEntry rdbe = regulomeDbChrEntries.get(gv.getStartPos());
								
								String supportData = filterSupportData(rdbe);
								if(!supportData.equals("")){
									if(eQtlCounts.containsKey(supportData)){
										eQtlCounts.put(supportData, eQtlCounts.get(supportData)+1);
									}
									
									else{
										eQtlCounts.put(supportData, 1);
									}
								}
							}
						}
					}
				}
			}
		}
		return eQtlCounts;
	}
	
	
	public String filterSupportData(RegulomeDbEntry rdbe){
		Map<String, List<RegulomeDbSupportingData>> supportData = rdbe.getSupportData();
		
		Iterator it = supportData.entrySet().iterator();
		while(it.hasNext()){
			Map.Entry pairs = (Map.Entry) it.next();

			List<RegulomeDbSupportingData> aap = (List<RegulomeDbSupportingData>) pairs.getValue();
			for(RegulomeDbSupportingData rdbsd : aap){

				//Check if the annotation is protein_binding.
				if(rdbsd.getSupportClass().equalsIgnoreCase("Protein_Binding")){
					return rdbsd.getSupportValue();
				}
			}
		}
		
		return "";
	}
	
	
	
	
	/*
	 * =========================================================================
	 * = OLD BUT WORKING METHODS.
	 * =========================================================================
	 */
	public HashMap<String, EQTL> filterEQtlSet(EQTL[] eqtls) throws IOException{
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
	
	
	
	public void filterOnSharedProbes(HashMap<String, EQTL> eQtlData, HashMap<String, EQTL> permutationData){
		
		for(Iterator<Map.Entry<String, EQTL>>iter=permutationData.entrySet().iterator();iter.hasNext();){
			Map.Entry<String, EQTL> permutationEntry = iter.next();
			
			if( !(eQtlData.containsKey(permutationEntry.getKey())) ){
				iter.remove();
			}
		}
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
	
	
	
	
	
	public void printCountsMap(HashMap<String, Integer> counts){
		Iterator<Entry<String, Integer>> countsIter = counts.entrySet().iterator();
		int n = 1;
		System.out.println("Entries in count map of size " + counts.values().size() + ".");
		while(countsIter.hasNext()){
			Entry<String, Integer> next = countsIter.next();
			System.out.println(n + ". " + next.getKey() + ": " + next.getValue());
			n++;
		}
	}
	
	
	public int getTotalHits(HashMap<String, Integer> countsMap){
		int total = 0;
		
		Iterator<Entry<String, Integer>> hitsIterator = countsMap.entrySet().iterator();
		while(hitsIterator.hasNext()){
			Entry<String, Integer> next = hitsIterator.next();
			total += next.getValue();
		}
		return total;
	}
}
