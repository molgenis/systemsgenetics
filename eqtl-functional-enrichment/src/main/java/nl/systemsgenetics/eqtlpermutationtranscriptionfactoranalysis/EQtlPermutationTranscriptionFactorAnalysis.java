/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import com.google.common.primitives.Doubles;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;

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
				RandomAccessGenotypeData genotypeMatrixData = edr.readEQtlGenotypeData(args[4], rsIdList);
				
				
				ArrayList<Double> eQtlHitScores = new ArrayList<Double>();
				ArrayList<Double> permutationHitScores = new ArrayList<Double>();
				
				//READ THE NEEDED DATA.
				GenomicBoundaries boundaries = edr.readHistoneDataFromText(args[5]);
				
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
				
				//GET THE RESULTS.
				hse.drawBoxPlot(eQtlHits, permutationHits, "pdf", args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("histonetest")){
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
				RandomAccessGenotypeData genotypeMatrixData = edr.readEQtlGenotypeData(args[4], rsIdList);
				
				
				ArrayList<Double> eQtlHitScores = new ArrayList<Double>();
				ArrayList<Double> permutationHitScores = new ArrayList<Double>();
				
				//READ THE NEEDED DATA.
				GenomicBoundaries boundaries = edr.readHistoneDataFromText(args[5]);
				
				//PERFORM THE ENRICHMENT FOR THE EQTLS.
				HistoneSiteEnrichment hse = new HistoneSiteEnrichment();
				hse.performAnalysis(topEqtlEffects, nonTopEqtlEffects, eQtlHitScores, genotypeMatrixData, boundaries, Double.valueOf(args[6]));
				double[] eQtlHits = Doubles.toArray(eQtlHitScores);
				System.out.println("Amount of eQTL hits: " + eQtlHits.length);
				
				
				//PERFORM THE ENRICHMENT FOR THE PERMUTATION DATA.
				EQTL[] permutationData = eqdp.readEQtlData(args[3]);
				HashMap<String, EQTL> topPermutationEffects = eqdp.getTopEffects(permutationData, sharedProbes);
				HashMap<String, HashSet<Integer>> nonTopPermutationEffects = eqdp.getNonTopEffects(permutationData, topPermutationEffects);
				hse.performAnalysis(topPermutationEffects, nonTopPermutationEffects, permutationHitScores, genotypeMatrixData, boundaries, Double.valueOf(args[6]));
				
				double[] permutationHits = Doubles.toArray(permutationHitScores);
				System.out.println("Amount of permutation hits: " + permutationHits.length);
				
				//GET THE RESULTS.
				hse.drawBoxPlot(eQtlHits, permutationHits, "pdf", args[7]);
			}
			
			else if(args[0].equalsIgnoreCase("convertWigToText")){
				WriteWigToText wwtt = new WriteWigToText();
				wwtt.writeWigToText(args[1], args[2]);
			}
			else if(args[0].equalsIgnoreCase("copyWigToText")){
				WriteWigToText wwtt = new WriteWigToText();
				wwtt.copyWigToText(args[1], args[2]);
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
}
