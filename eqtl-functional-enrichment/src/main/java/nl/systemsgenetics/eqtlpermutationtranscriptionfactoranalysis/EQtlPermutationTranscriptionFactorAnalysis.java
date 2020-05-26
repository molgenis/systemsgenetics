/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import com.google.common.primitives.Doubles;
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
import umcg.genetica.io.regulomedb.RegulomeDbFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author Matthieu
 */
public class EQtlPermutationTranscriptionFactorAnalysis {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public static void main(String[] args)throws IOException{
		/*
		 * args[0]: Type of analysis
		 *			- tf: Transcription Factor Enrichment
		 *			- gencode: ENCODE Region Enrichment
		 *			- dbrip: Retrotransposon Insertion Polymorphisms Enrichment
		 *			- repeats: Repeat Region Enrichment
		 * 
		 * args[1]: eQTL ProbeLevel file.
		 * args[2]: eQTL file.
		 * args[3]: permutation files locations.
		 * args[4]: genotype matrix location.
		 * args[5]: regulomedb location
		 *			gencode annotation file
		 *			repeat data file
		 *			dbrip data file
		 * args[6]: r2 cutoff.
		 * args[7]: output file location.
		 */
		
		if(args.length == 3 || args.length == 8){
			
			//START CHECKING WHICH TYPE OF ANALYSIS TO PERFORM.
			if(args[0].equalsIgnoreCase("tf")){
				EQtlPermutationTranscriptionFactorAnalysisV3 eqptfa3 = new EQtlPermutationTranscriptionFactorAnalysisV3(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("tf2")){
				TranscriptionFactorEnrichmentV2 tfe2 = new TranscriptionFactorEnrichmentV2(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("gencode")){
				eQtlsInEncodeAnnotationData eqiead = new eQtlsInEncodeAnnotationData(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("dbrip")){
				//eQtlsInAluRipData eqiard = new eQtlsInAluRipData(args[1], args[2], args[3], args[4], Double.valueOf(args[5]));
				System.out.println("This analysis type is disabled for now. No interesting results were produced earlier.");
			}
			
			else if(args[0].equalsIgnoreCase("repeats")){
				//Make call to perform analysis on repeat data from dasha.
				eQtlsInRepeatData eqird = new eQtlsInRepeatData(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("pseudogenes")){
				eQtlsInPseudogeneOrgData eqipod = new eQtlsInPseudogeneOrgData(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("histones") || args[0].equalsIgnoreCase("dnase") || args[0].equalsIgnoreCase("methylation")){
				//GET THE PROBES SHARED BY ALL USED DATASETS.
				EQtlPermutationTranscriptionFactorAnalysis eqptfa = new EQtlPermutationTranscriptionFactorAnalysis();
				HashSet<String> sharedProbes = eqptfa.getSharedProbes(args[1], args[3]);

				//READ THE EQTL DATA AND FILTER THE SET ON TOP AND NON TOP EQTL EFFECTS.
				eQtlDataParser eqdp = new eQtlDataParser();
				EQTL[] eqtls = eqdp.readEQtlData(args[2]);
				Set<String> rsIdList = eqdp.makeRsIdList(eqtls);
				HashMap<String, EQTL> topEqtlEffects = eqdp.getTopEffects(eqtls, sharedProbes);
				HashMap<String, HashSet<Integer>> nonTopEqtlEffects = eqdp.getNonTopEffects(eqtls, topEqtlEffects);

				//READ THE GENOTYPE DATA.
				EnrichmentDataReader edr = new EnrichmentDataReader();
				RandomAccessGenotypeData genotypeMatrixData = edr.readEQtlGenotypeDataV2(args[4], rsIdList);
				
				
				ArrayList<Double> eQtlHitScores = new ArrayList<Double>();
				ArrayList<Double> permutationHitScores = new ArrayList<Double>();
				
				//READ THE NEEDED DATA.
				//GenomicBoundaries boundaries = edr.readHistoneDataFromText(args[5]);
				GenomicBoundaries boundaries = edr.readHistoneNarrowPeakFileData(args[5]);
				
				//PERFORM THE ENRICHMENT FOR THE EQTLS.
				HistoneSiteEnrichment hse = new HistoneSiteEnrichment();
				hse.performAnalysis(topEqtlEffects, nonTopEqtlEffects, eQtlHitScores, genotypeMatrixData, boundaries, Double.valueOf(args[6]));
				double[] eQtlHits = Doubles.toArray(eQtlHitScores);
				
				//PERFORM THE ENRICHMENT FOR THE PERMUTATION DATA.
				for(int n=1;n<=100;n++){
					EQTL[] permutationData;
					HashMap<String, EQTL> topPermutationEffects;
					HashMap<String, HashSet<Integer>> nonTopPermutationEffects;

					permutationData = eqdp.readEQtlData(args[3] + "PermutedEQTLsPermutationRound" + n + ".txt.gz");
					topPermutationEffects = eqdp.getTopEffects(permutationData, sharedProbes);
					nonTopPermutationEffects = eqdp.getNonTopEffects(permutationData, topPermutationEffects);
					hse.performAnalysis(topPermutationEffects, nonTopPermutationEffects, permutationHitScores, genotypeMatrixData, boundaries, Double.valueOf(args[6]));
				}
				double[] permutationHits = Doubles.toArray(permutationHitScores);
				
				//First print the values.
				System.out.println(hse.getWilcoxonPValue(eQtlHits, permutationHits));
				System.out.println(hse.getWilcoxonUac(eQtlHits, permutationHits));
				
				
				//Write result scores to a file.
				eqptfa.writeThingsToTextFile("/target/gpfs2/gcc/groups/gcc/projects/eQtlFunctionalEnrichment/results/histoneEnrichment/resultsReal.txt", eQtlHits);
				eqptfa.writeThingsToTextFile("/target/gpfs2/gcc/groups/gcc/projects/eQtlFunctionalEnrichment/results/histoneEnrichment/resultsChance.txt", permutationHits);
				
				//GET THE RESULTS.
				hse.drawBoxPlot(eQtlHits, permutationHits, "pdf", "Histone Enrichment", args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("convertWigToText")){
				WriteWigToText wwtt = new WriteWigToText();
				wwtt.writeWigToText(args[1], args[2]);
			}
			else if(args[0].equalsIgnoreCase("copyWigToText")){
				WriteWigToText wwtt = new WriteWigToText();
				wwtt.copyWigToText(args[1], args[2]);
			}
			
			
			else if(args[0].equalsIgnoreCase("tfTotal")){
				String sharedProbesProbeLevelFile = args[1];
				String eQtlFileLocation = args[2];
				String permutationLocation = args[3];
				String GenotypeDataLocation = args[4];
				String RegulomeDbLocation = args[5];
				double r2Cutoff = Double.valueOf(args[6]);
				String outputFileLocation = args[7];
				
				//Read the Genotype data.
				EnrichmentDataReader edr = new EnrichmentDataReader();
				RandomAccessGenotypeData genotypeMatrixData = edr.readEQtlGenotypeDataV2(GenotypeDataLocation, null);
				HashMap<String, GeneticVariant> genotypeVariantIdMap = genotypeMatrixData.getVariantIdMap();
				
				//Read the RegulomeDB data.
				ArrayList<RegulomeDbFile> regulomeDbFiles = new ArrayList<RegulomeDbFile>();
				regulomeDbFiles.add(new RegulomeDbFile(new File(RegulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category1.txt")));
				regulomeDbFiles.add(new RegulomeDbFile(new File(RegulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category2.txt")));
				regulomeDbFiles.add(new RegulomeDbFile(new File(RegulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category3.txt")));
				regulomeDbFiles.add(new RegulomeDbFile(new File(RegulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category4.txt")));
				regulomeDbFiles.add(new RegulomeDbFile(new File(RegulomeDbLocation + "RegulomeDB.dbSNP132.b36.Category5.txt")));
				HashMap<String, TreeMap<Integer, String[]>> regulomeDbData = edr.readRegulomeDbData(regulomeDbFiles);
				
				//Filter on Shared Probes
				EQtlPermutationTranscriptionFactorAnalysis eqptfa = new EQtlPermutationTranscriptionFactorAnalysis();
				HashSet<String> sharedProbes = eqptfa.getSharedProbes(sharedProbesProbeLevelFile, permutationLocation);
				
				//READ THE EQTL DATA AND FILTER THE SET ON TOP AND NON TOP EQTL EFFECTS.
				eQtlDataParser eqdp = new eQtlDataParser();
				EQTL[] eqtls = eqdp.readEQtlData(eQtlFileLocation);
				HashMap<String, EQTL> topEqtlEffects = eqdp.getTopEffects(eqtls, sharedProbes);
				HashMap<String, HashSet<Integer>> nonTopEqtlEffects = eqdp.getNonTopEffects(eqtls, topEqtlEffects);
				
				
				HashMap<String, Integer> eQtlCounts = new HashMap<String, Integer>();
				eqptfa.getTranscriptionFactorCounts(topEqtlEffects, nonTopEqtlEffects, eQtlCounts, genotypeMatrixData, regulomeDbData, r2Cutoff);


				//STEP 6.: PERFORM ANALYSIS FOR PERMUTATION DATA.
				HashMap<String, Integer> permutationCounts = new HashMap<String, Integer>();
				for(int n=1;n<=100;n++){
					EQTL[] permutationData;
					HashMap<String, EQTL> topPermutationEffects;
					HashMap<String, HashSet<Integer>> nonTopPermutationEffects;

					permutationData = eqdp.readEQtlData(permutationLocation + "PermutedEQTLsPermutationRound" + n + ".txt.gz");
					topPermutationEffects = eqdp.getTopEffects(permutationData, sharedProbes);
					nonTopPermutationEffects = eqdp.getNonTopEffects(permutationData, topPermutationEffects);
					eqptfa.getTranscriptionFactorCounts(topPermutationEffects, nonTopPermutationEffects, permutationCounts, genotypeMatrixData, regulomeDbData, r2Cutoff);
				}
				
				
				
				
				/*
				 * NOW PERFORM THE ANALYSIS FOR ITERATIONS 1 THROUGH 9.
				 */
				for(int n=1;n<10;n++){
					EQTL[] iterationEQTLs;
					HashMap<String, EQTL> topIterationEffects;
					HashMap<String, HashSet<Integer>> nonTopIterationEffects;
					HashSet<String> sharedIterationProbes;
					
					//1. Get shared probes.
					sharedIterationProbes = eqptfa.getSharedIterationProbes("/target/gpfs2/gcc/groups/gcc/projects/eQTLMapping_Matthieu/RegressOut_Old/Iteration" + n + "/eQTLProbesFDR0.05-ProbeLevel.txt", "/target/gpfs2/gcc/groups/gcc/projects/eQTLMapping_Matthieu/RegressOut_Old/Iteration" + n + "/");
					
					//2. Filter data for real iteration eQTLs.
					iterationEQTLs = eqdp.readEQtlData("/target/gpfs2/gcc/groups/gcc/projects/eQTLMapping_Matthieu/RegressOut_Old/Iteration" + n + "/eQTLsFDR0.05.txt.gz");
					topIterationEffects = eqdp.getTopEffects(iterationEQTLs, sharedIterationProbes);
					nonTopIterationEffects = eqdp.getNonTopEffects(iterationEQTLs, topIterationEffects);
					
					//3. Perform analysis for real iteration eQTLs.
					eqptfa.getTranscriptionFactorCounts(topIterationEffects, nonTopIterationEffects, eQtlCounts, genotypeMatrixData, regulomeDbData, r2Cutoff);
					
					//4.Perform the analysis for the permutation data.
					for(int o=1;o<=100;o++){
						EQTL[] permutationData;
						HashMap<String, EQTL> topIterationPermutationEffects;
						HashMap<String, HashSet<Integer>> nonTopIterationPermutationEffects;
						
						permutationData = eqptfa.readPermutationData("/target/gpfs2/gcc/groups/gcc/projects/eQTLMapping_Matthieu/RegressOut_Old/Iteration" + n + "/PermutedEQTLsPermutationRound" + o + ".txt.gz", genotypeVariantIdMap);
						topIterationPermutationEffects = eqdp.getTopEffects(permutationData, sharedProbes);
						nonTopIterationPermutationEffects = eqdp.getNonTopEffects(permutationData, topIterationPermutationEffects);
						eqptfa.getTranscriptionFactorCounts(topIterationPermutationEffects, nonTopIterationPermutationEffects, permutationCounts, genotypeMatrixData, regulomeDbData, r2Cutoff);
					}	
				}
				
				eqptfa.getFisherPvalues(eQtlCounts, permutationCounts, outputFileLocation);
			}
			
			else{
				EQtlPermutationTranscriptionFactorAnalysis aap = new EQtlPermutationTranscriptionFactorAnalysis();
				aap.printHelp();
			}
		}
		else{
			EQtlPermutationTranscriptionFactorAnalysis aap = new EQtlPermutationTranscriptionFactorAnalysis();
			aap.printHelp();
		}
	}
	
	public EQtlPermutationTranscriptionFactorAnalysis(){
		
	}
	
	
	public void printHelp(){
		System.out.println("Use one of the following values as the first parameter:");
		System.out.println("\t\"tf\": Enrichment analysis for RegulomeDB Transcription Factors.");
		System.out.println("\t\"gencode\": Enrichment analysis for ENCODE annotation data.");
		System.out.println("\t\"dbrip\": Enrichment analysis for DBRIP data.");
		System.out.println("\t\"repeats\": Enrichment analysis for repeat data.");
		System.out.println("\t\"pseudogenes\": Enrichment for pseudogene data from pseudogenes.org");
		System.out.println("\t\"histones\": Enrichment for ENCODE histone sites. (Can use Build36 wig or Build36 txt)");

		System.out.println("\nOther things you can do:");
		System.out.println("\t\"convertWigToText\": Data from wig file will be placed in a txt as: chr\tstart\tstop\tscore");
		System.out.println("\t\"copyWigToText\": Data from wig file will be copied to a txt file.");
	}
	
	/*
	 * =========================================================================
	 * = CODE NEEDED BY ALL CLASSES.
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
	
	
	public void writeThingsToTextFile(String outFile, double[] vals)throws IOException{
		
		PrintWriter pw = new PrintWriter(new FileWriter(outFile));
		for(double val : vals){
			pw.println(val);
		}
		pw.close();
	}
	
	
	
	
	
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
					
					if(!(gv.getSequenceName().equalsIgnoreCase("X") || gv.getSequenceName().equalsIgnoreCase("Y") || gv.getSequenceName().equalsIgnoreCase("Mt"))){
						
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
					}
				}
			}
			else{
				a = '!';
			}
		}
		tf.close();
		return permutations.toArray(new EQTL[permutations.size()]);
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
	
	
	public HashSet<String> getSharedIterationProbes(String eQtlProbeFile, String permutationDataLocation)throws IOException{
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
}
