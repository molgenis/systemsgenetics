/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
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

/**
 *
 * @author Matthieu
 */
public class eQtlAndLdInRepeatRegions {
	private static final Pattern CHR_PATTERN = Pattern.compile("\\d{1,2}");
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	private static final Pattern SUPPORTDATA_PATTERN = Pattern.compile("|");
	
	public static void main(String[] args) throws IOException, Exception{
		eQtlAndLdInRepeatRegions eqalirr = new eQtlAndLdInRepeatRegions();
	}
	
	public eQtlAndLdInRepeatRegions()throws IOException, Exception{
		System.out.println("[1]: Read the repeat data.");
		GenomicBoundaries<Object> readRepeatData = readRepeatData("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\4.eQtlsInRepeats\\testFile.txt");
		System.out.println("[1]: Read " + readRepeatData.getBoudaryCount());
		
		System.out.println("\n[2]: Read the eQTL data.");
		EQTL[] readEQtlData = readEQtlData("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\4.eQtlsInRepeats\\testEqtlFile.txt");
		System.out.println("[2]: Amount of eQTLs read: " + readEQtlData.length);
		
		System.out.println("\n[3]: Get a list of rsID's.");
		Set<String> listOfRsIds = makeRsIdList(readEQtlData);
		System.out.println("[3]: Got " + listOfRsIds.size() + " unique IDs.");
		
		System.out.println("\n[4]: Read some genotype matrix data with a filter.");
		RandomAccessGenotypeData genotypeMatrixData = readEQtlGenotypeData("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\BloodHT12Combined\\", listOfRsIds);
		System.out.println("[4]: Done reading the genotype matrix data.");
		
		System.out.println("\n[5]: Start searching LD SNPs for each eQTL.");
		HashMap<String, ArrayList<Ld>> ldSnps = findLdSnps(readEQtlData, genotypeMatrixData, 250000, 0.8);
		System.out.println("[5]: Done searching LD SNPs for each eQTL.");
		
		//printLdData(ldSnps);
		
		System.out.println("\n[6]: Start searching with the LD SNPs through the repeat sequences.");
		analysis(ldSnps, readRepeatData);
		System.out.println("[6]: Done searching through the repeat sequences.");
		
	}
	
	
	/*
	 * Method to read the repeat data into Genomic Boundaries. Hurray :)
	 */
	public GenomicBoundaries<Object> readRepeatData(String fileLocation) throws IOException{
		GenomicBoundaries<Object> repeatBoundaries = new GenomicBoundaries();
		String fileLine;
		String[] fileLineData;
		
		TextFile repeatTextFile = new TextFile(fileLocation, false);
		while( (fileLine = repeatTextFile.readLine())!=null ){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			//Add a filter to exclude chromosome X and other non chromosome specific entries.
			Matcher matcher = CHR_PATTERN.matcher(fileLineData[0]);
			if(matcher.matches()){
				String chr = fileLineData[0];
				int startPos = Integer.parseInt(fileLineData[1]);
				int stopPos = Integer.parseInt(fileLineData[2]);
				String supportData = fileLineData[3]+"|"+fileLineData[4]+"|"+fileLineData[5];
				
				repeatBoundaries.addBoundary(fileLineData[0], startPos, stopPos, supportData);
			}
			
		}
		repeatTextFile.close();
		return repeatBoundaries;
	}
	
	
	/*
	 * Just read the eQTL data.
	 */
	public EQTL[] readEQtlData(String eQtlFileLocation) throws IOException{
		QTLTextFile eqtlData = new QTLTextFile(eQtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData, Set<String> variantIdFilter) throws IOException{
		//Provide a Set<String> containing rsID of all significant eQTLs.
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData), 1000, new VariantIdIncludeFilter(variantIdFilter), new SampleIncludedFilter());
		return gonlImputedBloodGenotypeData;
	}
	
	
	
	/*
	 * Find LD SNPs
	 */
	public HashMap<String, ArrayList<Ld>> findLdSnps(EQTL[] eqtls, RandomAccessGenotypeData genotypeData, int windowSize, double r2Cutoff){
		//TODO: Code to find ld snps.
		HashMap<String, ArrayList<Ld>> ldSnpsPerEQtl = new HashMap<String, ArrayList<Ld>>();
		Ld ld = null;
		ArrayList<Ld> tmpLdList;
		
		//
		for(EQTL eqtl : eqtls){
			String rsName = eqtl.getRsName();
			String rsChr = eqtl.getRsChr().toString();
			int rsPos = eqtl.getRsChrPos();
			
			//Search eQtl back in GenotypeMatrix data.
			Iterable<GeneticVariant> variantsByPos = genotypeData.getVariantsByPos(rsChr, rsPos);
			for(GeneticVariant eqtlSnp : variantsByPos){
				int boundaryL = eqtlSnp.getStartPos() - windowSize;
				int boundaryR = eqtlSnp.getStartPos() + windowSize;
				String chr = eqtlSnp.getSequenceName();
				
				
				if(eqtlSnp.isBiallelic()){
					Iterable<GeneticVariant> ldVariantsByRange = genotypeData.getVariantsByRange(chr, boundaryL, boundaryR);
					
					for(GeneticVariant rangeVariant : ldVariantsByRange){
						if(rangeVariant.isBiallelic()){
							try{
								ld = eqtlSnp.calculateLd(rangeVariant);
							} catch(LdCalculatorException ex){
								System.out.println("Error in LD calculation: " + ex.getMessage());
								System.exit(1);
							}
							
							
							//BLOCK START: CHECK THE R2 VALUE.
							if(ld.getR2() >= r2Cutoff){
								GeneticVariant eQtlVariant = ld.getVariant1();
								GeneticVariant ldVariant = ld.getVariant2();
								String eQtlRsId = eQtlVariant.getPrimaryVariantId();
								
								if(ldSnpsPerEQtl.containsKey(eQtlRsId)){
									tmpLdList = ldSnpsPerEQtl.get(eQtlRsId);
									tmpLdList.add(ld);
								}
								
								else{
									tmpLdList = new ArrayList<Ld>();
									tmpLdList.add(ld);
									ldSnpsPerEQtl.put(eQtlRsId, tmpLdList);
								}
							}
							//BLOCK END: CHECK THE R2 VALUE.
						}
					}
				}
			}
			//Processed one EQTL element from the EQTL[] thingie.
		}
		return ldSnpsPerEQtl;
	}
	
	
	public void analysis(HashMap<String, ArrayList<Ld>> ldData, GenomicBoundaries<Object> repeatData) throws Exception{
		Iterator<Entry<String, ArrayList<Ld>>> ldIterator = ldData.entrySet().iterator();
		while(ldIterator.hasNext()){
			Map.Entry<String, ArrayList<Ld>> next = ldIterator.next();
			ArrayList<Ld> value = next.getValue();
			
			for(Ld ld : value){
				GeneticVariant eqtlVariant = ld.getVariant1();
				GeneticVariant ldVariant = ld.getVariant2();
				System.out.println("eqtl variant seq name: " + eqtlVariant.getSequenceName() );
				
				//Some testing:
				if(repeatData.isChromosomeInBoundary(eqtlVariant.getSequenceName())){
					TreeMap<Integer, ArrayList<GenomicBoundary<Object>>> genomicBoundariesMap = repeatData.getGenomicBoundariesMap(eqtlVariant.getSequenceName());
					NavigableMap<Integer, ArrayList<GenomicBoundary<Object>>> headMap = genomicBoundariesMap.headMap(eqtlVariant.getStartPos(), true);
					System.out.println("Start Searching");
					searchThroughRepeats(ldVariant, headMap);
				}
				else{
					System.out.println("aap!");
				}
			}
		}
	}
	
	
	public void searchThroughRepeats(GeneticVariant gv, NavigableMap<Integer, ArrayList<GenomicBoundary<Object>>> selectRepeatData){
		Iterator<Entry<Integer, ArrayList<GenomicBoundary<Object>>>> iterator = selectRepeatData.entrySet().iterator();
		System.out.println("Zoek");
		while(iterator.hasNext()){
			Entry<Integer, ArrayList<GenomicBoundary<Object>>> next = iterator.next();
			ArrayList<GenomicBoundary<Object>> value = next.getValue();
			
			for(GenomicBoundary boundary : value){
				System.out.println( boundary.isInBoundarie(gv.getStartPos()) );
			}
		}
	}
	
	
	
	/*
	 * Make a list of rsIDs for the 
	 */
	public Set<String> makeRsIdList(EQTL[] eqtls){
		Set<String> rsIdList = new HashSet<String>();
		for(EQTL eqtl : eqtls){
			rsIdList.add(eqtl.getRsName());
		}
		return rsIdList;
	}
	
	
	public void printLdData(HashMap<String, ArrayList<Ld>> ldData){
		System.out.println("Print LD data.");
		
		int n = 1;
		Iterator<Entry<String, ArrayList<Ld>>> iterator = ldData.entrySet().iterator();
		while( iterator.hasNext() ){
			Entry<String, ArrayList<Ld>> next = iterator.next();
			ArrayList<Ld> value = next.getValue();
			
			for(Ld ld : value){
				System.out.println(ld.getVariant2().getSequenceName() + "\t"
						+ ld.getVariant1().getPrimaryVariantId() + "\t"
						+ ld.getVariant2().getPrimaryVariantId() + "\t" + n);
				n++;
			}
		}
	}
}
