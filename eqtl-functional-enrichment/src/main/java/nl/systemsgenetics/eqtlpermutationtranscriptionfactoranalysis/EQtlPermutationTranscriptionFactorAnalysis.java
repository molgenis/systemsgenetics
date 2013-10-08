/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
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
import umcg.genetica.io.regulomedb.RegulomeDbEntry;
import umcg.genetica.io.regulomedb.RegulomeDbFile;
import umcg.genetica.io.regulomedb.RegulomeDbFiles;
import umcg.genetica.io.regulomedb.RegulomeDbSupportingData;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.math.stats.FisherExactTest;

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
		
		if(args[0].equals("tf")){
			EQtlPermutationTranscriptionFactorAnalysisV3 eqptfa3 = new EQtlPermutationTranscriptionFactorAnalysisV3(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
		}
		else if(args[0].equals("gencode")){
			eQtlsInEncodeAnnotationData eqiead = new eQtlsInEncodeAnnotationData(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
		}
		else if(args[0].equals("dbrip")){
			//eQtlsInAluRipData eqiard = new eQtlsInAluRipData(args[1], args[2], args[3], args[4], Double.valueOf(args[5]));
			System.out.println("This analysis type is disabled for now. No interesting results were produced earlier.");
		}
		else if(args[0].equals("repeats")){
			//Make call to perform analysis on repeat data from dasha.
			eQtlsInRepeatData eqird = new eQtlsInRepeatData(args[1], args[2], args[3], args[4], args[5], Double.valueOf(args[6]), args[7]);
		}
		else{
			System.out.println("Use one of the following values as the first parameter:");
			System.out.println("\t\"tf\": Enrichment analysis for RegulomeDB Transcription Factors.");
			System.out.println("\t\"gencode\": Enrichment analysis for ENCODE annotation data.");
			System.out.println("\t\"dbrip\": Enrichment analysis for DBRIP data.");
			System.out.println("\t\"repeats\": Enrichment analysis for repeat data.");
		}
		
	}
	
	
	
	public EQtlPermutationTranscriptionFactorAnalysis(String eQtlProbeLevelFile, String eQtlFile, String permutationFile, String genotypeLocation, String regulomeDbLocation, int windowSize, double r2Cutoff)throws IOException{
		
	}
}
