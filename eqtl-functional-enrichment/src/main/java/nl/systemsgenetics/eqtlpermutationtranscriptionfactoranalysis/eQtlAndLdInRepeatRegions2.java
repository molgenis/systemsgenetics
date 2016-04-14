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
import java.util.Map;
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
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author Matthieu
 */
public class eQtlAndLdInRepeatRegions2 {
	private static final Pattern CHR_PATTERN = Pattern.compile("\\d{1,2}");
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	
	public static void main(String[] args){
		
	}
	
	public GenomicBoundaries<Object> readRepeatData(String fileLocation) throws IOException{
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
			}
			
		}
		repeatTextFile.close();
		//return repeatBoundaries;
		return null;
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
	public HashMap<String, TreeMap<Integer, ArrayList<Ld>>> calculateLd(HashMap<String, EQTL> eqtlData, RandomAccessGenotypeData genotypeData, int windowSize, double r2CutOff)throws IOException{
		//Use a window size of 250k: eQTL pos - 250k and eQTL pos + 250k
		Ld ld = null;
		HashMap<String, TreeMap<Integer, ArrayList<Ld>>> ldResults = new HashMap<String, TreeMap<Integer, ArrayList<Ld>>>();
		
		Iterator<Map.Entry<String, EQTL>> eqtlIterator = eqtlData.entrySet().iterator();
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
		return ldResults;
	}
	
}
