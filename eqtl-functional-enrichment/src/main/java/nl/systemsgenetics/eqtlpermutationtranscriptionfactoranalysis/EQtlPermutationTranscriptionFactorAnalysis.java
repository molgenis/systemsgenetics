/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.regex.Pattern;
import nl.umcg.deelenp.regulomedb.RegulomeDbEntry;
import nl.umcg.deelenp.regulomedb.RegulomeDbFile;
import nl.umcg.deelenp.regulomedb.RegulomeDbFiles;
import nl.umcg.deelenp.regulomedb.RegulomeDbSupportingData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
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
		EQtlPermutationTranscriptionFactorAnalysis eqptfa = new EQtlPermutationTranscriptionFactorAnalysis();
	}
	
	
	/*
	 * Constructor that calls the other methods. The current code might change later.
	 */
	public EQtlPermutationTranscriptionFactorAnalysis()throws IOException, LdCalculatorException{
		//Read the eQTL data (preferably from eQTLProbesFDR0.05-ProbeLevel.txt file.
		EQTL[] eQtlResultData = this.readEQtlResultData("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\results\\Imputed\\eQTLProbesFDR0.05-ProbeLevel.txt");
		
		//Read the GenotypeMatrix.dat data using a RandomAccessGenotypeData object.
		RandomAccessGenotypeData eQtlGenotypeData = readEQtlGenotypeData("E:\\GroningenBloodData\\BloodHT12Combined\\");
		
		//Calculate the LD.
		HashMap<String, TreeMap<Integer, ArrayList<Ld>>> eQtlLdData = this.calculateLd(eQtlResultData, 250000, 0.8, eQtlGenotypeData);
		
		//Find SNPs In regulomeDB
		ArrayList<RegulomeDbFile> regulomeDbFiles = new ArrayList<RegulomeDbFile>();
		regulomeDbFiles.add(new RegulomeDbFile( new File("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Data\\regulomeDb\\RegulomeDB.dbSNP132.Category7.txt") ));
		
	}
	
	
	
	/*
	 * ===========================================================================
	 * = START OF DATA READING CODE.
	 * ===========================================================================
	 */
	public EQTL[] readEQtlResultData(String eqtlFileLocation) throws IOException{
		eQTLTextFile eqtlData = new eQTLTextFile(eqtlFileLocation, false);
		System.out.print("Read the eQTL data.\n\n");
		return eqtlData.read();
	}
	
	
	public RandomAccessGenotypeData readEQtlGenotypeData(String genotypeData) throws IOException{
		RandomAccessGenotypeData gonlImputedBloodGenotypeData = new TriTyperGenotypeData( new File(genotypeData) );
		System.out.println("I guess I've read some genotype data matrix things.");
		System.out.print("I contain " + gonlImputedBloodGenotypeData.getSeqNames().size() + " sequence names.\n\n");
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
	public HashMap<String, TreeMap<Integer, ArrayList<Ld>>> calculateLd(EQTL[] eqtlData, int windowSize, double r2CutOff, RandomAccessGenotypeData genotypeData)throws IOException{
		//Use a window size of 250k: eQTL pos - 250k and eQTL pos + 250k
		Ld ld = null;
		HashMap<String, ArrayList<Ld>> ldResultsMap = new HashMap<String, ArrayList<Ld>>();
		HashMap<String, TreeMap<Integer, ArrayList<Ld>>> ldResults = new HashMap<String, TreeMap<Integer, ArrayList<Ld>>>();
		int n = 1;
		
		System.out.println("I will now start with the LD calculation for each eQTL.\n");
		for(EQTL eqtl : eqtlData){
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
							System.out.println("Found LD SNP " + gv.getPrimaryVariantId() + " for eQTL " + eQtlSnp.getPrimaryVariantId());

							//Place results in a convenient structure for later.
							TreeMap<Integer, ArrayList<Ld>> tmp;
							ArrayList<Ld> ldList;
							if(ldResultsMap.containsKey(variant2.getSequenceName())){
								
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

								/*
								tmp = ldResults.get(variant2.getSequenceName());
								//tmp.put(variant2.getStartPos(), ld);
								ldList = tmp.get(variant2.getStartPos());
								ldList.add(ld);
								*/
							}

							else{
								tmp = new TreeMap<Integer, ArrayList<Ld>>();
								ldList = new ArrayList<Ld>();
								ldList.add(ld);
								tmp.put(variant2.getStartPos(), ldList);
								ldResults.put(variant2.getSequenceName(), tmp);
							}
							System.out.println("Data: " + variant2.getPrimaryVariantId() + " , " + variant2.getStartPos());
						}
					}
					else{
						System.out.println("eQTL " + eQtlSnp.getPrimaryVariantId() + " with " + gv.getPrimaryVariantId() + " are not biallelic.");
					}
					
				}
				
			}
		}
		
		return ldResults;
	}
	
	/*
	 * ===========================================================================
	 * = END OF LD CALCULATION CODE.
	 * ===========================================================================
	 */
	
	
	
	public void findSnpsInRegulomeDb(ArrayList<RegulomeDbFile> regulomeDbFileLocations, HashMap<String, TreeMap<Integer, ArrayList<Ld>>> ldData) throws IOException{
		//TextFile resultsOutputFile = new TextFile("", true);
		//resultsOutputFile.write("eQTL\tLD_SNP\tRegulomeDbScore\tRegulomeData");
		
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
						//Write results to a file. :)
						System.out.println( ld.getVariant1().getPrimaryVariantId() + "\t" + ld.getVariant2().getPrimaryVariantId() + "\t" + rdbe.getRegulomeDbScore() + "\t"
								+ rdbe.getSupportData().values().toArray().toString());
					}
				}
			}
			
		}
		//resultsOutputFile.close();
	}
}
