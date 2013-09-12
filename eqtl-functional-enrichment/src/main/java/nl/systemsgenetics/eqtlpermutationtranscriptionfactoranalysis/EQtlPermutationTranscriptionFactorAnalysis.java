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
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;
import umcg.genetica.io.regulomedb.RegulomeDbEntry;
import umcg.genetica.io.regulomedb.RegulomeDbFile;
import umcg.genetica.io.regulomedb.RegulomeDbFiles;
import umcg.genetica.io.regulomedb.RegulomeDbSupportingData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author Matthieu
 */
public class EQtlPermutationTranscriptionFactorAnalysis {
	
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private RandomAccessGenotypeData eQtlGenotypeData;
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, LdCalculatorException {
		// TODO code application logic here
		
		/*
		 * Order of commandline arguments:
		 * 1. String: Path to eQTL or permutation data file 
		 * 2. String: Path to folder containing genotypeMatrix data (GenotypeMatrix.dat)
		 * 3. int: LD Search window size
		 * 4. double: R2 cutoff value (cutoff will be used as >= cutoff value)
		 * 5. String: Path to folder containing regulomeDB (make sure that the files are named RegulomeDB.dbSNP132.Category#.txt)
		 * 6. String: Path for output file location.
		 */
		//EQtlPermutationTranscriptionFactorAnalysis eqptfa = new EQtlPermutationTranscriptionFactorAnalysis();
		EQtlPermutationTranscriptionFactorAnalysis eqptfa =
				new EQtlPermutationTranscriptionFactorAnalysis(args[0], args[1], Integer.parseInt(args[2]), Double.valueOf(args[3]), args[4]);
		
		//EQtlPermutationTranscriptionFactorAnalysis eqptfa =
		//		new EQtlPermutationTranscriptionFactorAnalysis("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\2.eQtlFunctionalEnrichment\\analysis\\data\\eQTLsFDR0.05.txt",
		//		"E:\\GroningenBloodData\\BloodHT12Eur1000GImputed\\", 250000, 0.8,
		//		"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\regulomeDb\\");
		//C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\BloodHT12Combined\\
	}
	
	
	/*
	 * Constructor that calls the other methods. The current code might change later.
	 */
	public EQtlPermutationTranscriptionFactorAnalysis(String eQtlFile, String genotypeData, int ldWindow, double r2cutoff, String regulomeDbFilesFolder)throws IOException, LdCalculatorException{
		
		//Step 1.: Read the eQTL data.
		System.out.println("Read the eQTL data and apply a filtering step.");
		EQTL[] eQtlResultData = this.readEQtlResultData(eQtlFile);
		HashMap<String, EQTL> topEQtlData = filterEQtlSet(eQtlResultData);
		Set<String> eqtlRsIdList = makeRsIdList(eQtlResultData);
		
		//Step 2.: Read the RandonAccess GenotypeData.
		System.out.println("Starting to read the GenotypeMatrix data.");
		RandomAccessGenotypeData eQtlGenotypeData = readEQtlGenotypeData(genotypeData, eqtlRsIdList);
		
		//Step 3.: Calculate LD
		System.out.println("Start identifying LD Snps with eQTL effects and calculate the R2.");
		HashMap<String, TreeMap<Integer, ArrayList<Ld>>> eQtlLdData = this.takeALookAtFindLdSnps(topEQtlData, ldWindow, r2cutoff, eQtlGenotypeData);
		
		//Step 4.: Find SNPs In regulomeDB
		System.out.println("Start searching through RegulomeDB with the LD SNPs.");
		ArrayList<RegulomeDbFile> regulomeDbFiles = new ArrayList<RegulomeDbFile>();
		regulomeDbFiles.add(new RegulomeDbFile( new File(regulomeDbFilesFolder + "RegulomeDB.dbSNP132.b36.txt") ));
		findSnpsInRegulomeDb(regulomeDbFiles, eQtlLdData);
	}
	
	
	
	/*
	 * ===========================================================================
	 * = START OF DATA READING CODE.
	 * ===========================================================================
	 */
	public EQTL[] readEQtlResultData(String eqtlFileLocation) throws IOException{
		eQTLTextFile eqtlData = new eQTLTextFile(eqtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 1000, new VariantIdIncludeFilter(variantIdFilter), new SampleIncludedFilter());
		return gonlImputedBloodGenotypeData;
	}
	/*
	 * ===========================================================================
	 * = END OF DATA READING CODE.
	 * ===========================================================================
	 */
	
	
	
	/*
	 * ===========================================================================
	 * = START OF LD CALCULATION CODE.
	 * ===========================================================================
	 */
	
	/*
	 * IMPORTANT NOTE: I've chosen for the current return structure because I wanted to easily search through it when working
	 * with the regulomeDB data.
	 */
	public HashMap<String, TreeMap<Integer, ArrayList<Ld>>> calculateLd(HashMap<String, EQTL> eqtlData, int windowSize, double r2CutOff, RandomAccessGenotypeData genotypeData)throws IOException{
		//Use a window size of 250k: eQTL pos - 250k and eQTL pos + 250k
		Ld ld = null;
		HashMap<String, TreeMap<Integer, ArrayList<Ld>>> ldResults = new HashMap<String, TreeMap<Integer, ArrayList<Ld>>>();
		
		Iterator<Entry<String, EQTL>> eqtlIterator = eqtlData.entrySet().iterator();
		while(eqtlIterator.hasNext()){
			Map.Entry pairs = (Map.Entry) eqtlIterator.next();
			EQTL eqtl = (EQTL) pairs.getValue();
			
			GeneticVariant eQtlSnp = genotypeData.getSnpVariantByPos(eqtl.getRsChr().toString(), eqtl.getRsChrPos());
			
			if(eQtlSnp != null){
				
				for(GeneticVariant gv : genotypeData.getVariantsByRange(eqtl.getRsChr().toString(), 0, windowSize)){
					
					if(eQtlSnp.isBiallelic() && gv.isBiallelic()){
						try {
							ld = eQtlSnp.calculateLd(gv);
						} catch (LdCalculatorException ex) {
							System.out.println("Error in LD calculation: " + ex.getMessage());
							System.exit(1);
						}
						GeneticVariant variant1 = ld.getVariant1();
						GeneticVariant variant2 = ld.getVariant2();


						if(ld.getR2() >= r2CutOff){

							//Place results in a convenient structure for later.
							TreeMap<Integer, ArrayList<Ld>> tmp;
							ArrayList<Ld> ldList;
							if(ldResults.containsKey(variant2.getSequenceName())){
								
								tmp = ldResults.get(variant2.getSequenceName());

								if(tmp.containsKey(variant2.getStartPos())){
									ldList = tmp.get(variant2.getStartPos());
									ldList.add(ld);
								}

								else{
									ldList = new ArrayList<Ld>();
									ldList.add(ld);
									tmp.put(variant2.getStartPos(), ldList);
								}
							}

							else{
								tmp = new TreeMap<Integer, ArrayList<Ld>>();
								ldList = new ArrayList<Ld>();
								ldList.add(ld);
								tmp.put(variant2.getStartPos(), ldList);
								ldResults.put(variant2.getSequenceName(), tmp);
							}
						}
					}
					
				}
				
			}
		}
		System.out.println("The map size is: " + ldResults.values().size());
		return ldResults;
	}
	/*
	 * ===========================================================================
	 * = END OF LD CALCULATION CODE.
	 * ===========================================================================
	 */
	
	
	public void findSnpsInRegulomeDb(ArrayList<RegulomeDbFile> regulomeDbFileLocations, HashMap<String, TreeMap<Integer, ArrayList<Ld>>> ldData) throws IOException{
		HashMap<String, Integer> countsMap = new HashMap<String, Integer>();
		HashSet<String> eQtls = new HashSet<String>();
		int n = 1;
		
		RegulomeDbFiles regulomeDbData = new RegulomeDbFiles(regulomeDbFileLocations);
		Iterator<RegulomeDbEntry> regulomeDbDataIterator = regulomeDbData.iterator();
		
		while(regulomeDbDataIterator.hasNext()){
			RegulomeDbEntry rdbe = regulomeDbDataIterator.next();
			String chr = rdbe.getChr();
			int pos = rdbe.getChrPos();
			
			if(ldData.containsKey(chr)){
				TreeMap<Integer, ArrayList<Ld>> ldChrData = ldData.get(chr);
				
				if(ldChrData.containsKey(pos)){
					ArrayList<Ld> ldList = ldChrData.get(pos);
					
					for(Ld ld : ldList){
						//Add eQTL id to list :)
						String eQtlId = ld.getVariant1().getPrimaryVariantId();
						Map<String, List<RegulomeDbSupportingData>> supportData = rdbe.getSupportData();
						
						Iterator it = supportData.entrySet().iterator();
						while(it.hasNext()){
							Map.Entry pairs = (Map.Entry) it.next();
							
							List<RegulomeDbSupportingData> aap = (List<RegulomeDbSupportingData>) pairs.getValue();
							for(RegulomeDbSupportingData rdbsd : aap){
								
								//Check if the annotation is protein_binding.
								if(rdbsd.getSupportClass().equalsIgnoreCase("Protein_Binding")){
									
									String tfName = rdbsd.getSupportValue();
									if(eQtls.contains(eQtlId)){
										if( !(countsMap.containsKey(tfName)) ){
											countsMap.put(tfName, 1);
										}
									}
									
									else{
										if(countsMap.containsKey(tfName)){
											int counter = countsMap.get(tfName);
											counter++;
											countsMap.put(tfName, counter);
										}

										else{
											countsMap.put(tfName, 1);
										}
									}
								}
								eQtls.add(eQtlId);
							}
						}
					}
					n++;
				}
			}
			
		}
		printCountsMap(countsMap);
		//System.out.println("eQTLs: " + eQtls.size());
		//System.out.println("counts: " + countsMap.size());
	}
	
	
	
	public HashMap<String, TreeMap<Integer, ArrayList<Ld>>> takeALookAtFindLdSnps(HashMap<String, EQTL> eqtlData, int windowSize, double r2CutOff, RandomAccessGenotypeData genotypeData) throws IOException{
		Ld ld = null;
		HashMap<String, TreeMap<Integer, ArrayList<Ld>>> ldResults = new HashMap<String, TreeMap<Integer, ArrayList<Ld>>>();
		int n = 1;
		
		Iterator<Entry<String, EQTL>> eQtlIterator = eqtlData.entrySet().iterator();
		
		while(eQtlIterator.hasNext()){
			Map.Entry pairs = (Map.Entry) eQtlIterator.next();
			EQTL eqtl = (EQTL) pairs.getValue();
			
			int rsPos = eqtl.getRsChrPos().intValue();
			String rsName = eqtl.getRsName();
			Iterable<GeneticVariant> variantsByPos = genotypeData.getVariantsByPos(eqtl.getRsChr().toString(), rsPos);
			
			for(GeneticVariant eqtlSnp : variantsByPos){
				int boundaryL = eqtlSnp.getStartPos() - 250000;
				int boundaryR = eqtlSnp.getStartPos() + 250000;
				
				if(eqtlSnp.isBiallelic()){
					Iterable<GeneticVariant> variantsByRangeForward = genotypeData.getVariantsByRange(eqtl.getRsChr().toString(), rsPos-windowSize, rsPos+windowSize);
					
					for(GeneticVariant gv : variantsByRangeForward){
						
						if(gv.isBiallelic()){
							try {
							ld = eqtlSnp.calculateLd(gv);
							} catch (LdCalculatorException ex) {
								System.out.println("Error in LD calculation: " + ex.getMessage());
								System.exit(1);
							}
							
							
							//Check if the R2 equals or is greater than the set cutoff.
							if(ld.getR2() >= r2CutOff){
								GeneticVariant variant1 = ld.getVariant1();
								GeneticVariant variant2 = ld.getVariant2();
								
								TreeMap<Integer, ArrayList<Ld>> tmp;
								ArrayList<Ld> ldList;
								
								if(ldResults.containsKey(variant2.getSequenceName())){
									tmp = ldResults.get(variant2.getSequenceName());
									if(tmp.containsKey(variant2.getStartPos())){
										ldList = tmp.get(variant2.getStartPos());
										ldList.add(ld);
									}
									
									else{
										ldList = new ArrayList<Ld>();
										ldList.add(ld);
										tmp.put(variant2.getStartPos(), ldList);
									}
								}
								
								else{
									tmp = new TreeMap<Integer, ArrayList<Ld>>();
									ldList = new ArrayList<Ld>();
									ldList.add(ld);
									tmp.put(variant2.getStartPos(), ldList);
									ldResults.put(variant2.getSequenceName(), tmp);
								}
								n++;
							}
							
							
						}
					}
				}
			}
		}
		//System.out.println("Amount of valid things saved: " + n);
		//System.out.println("Size of LD map: " + ldResults.size());
		
		return ldResults;
	}
	
	
	
	
	
	
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
	
	
	public void printCountsMap(HashMap<String, Integer> counts){
		Iterator<Entry<String, Integer>> countsIter = counts.entrySet().iterator();
		int n = 1;
		System.out.println("Entries in count map.");
		while(countsIter.hasNext()){
			Entry<String, Integer> next = countsIter.next();
			System.out.println(n + ". " + next.getKey() + ": " + next.getValue());
			n++;
		}
	}
	
}
