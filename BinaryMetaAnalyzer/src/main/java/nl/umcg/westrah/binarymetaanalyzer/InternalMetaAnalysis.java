/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

import javax.xml.bind.annotation.adapters.HexBinaryAdapter;

import umcg.genetica.io.bin.BinaryFile;

import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import org.apache.commons.collections.primitives.ArrayIntList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.math.stats.ZScores;

/**
 * @author Marc-Jan
 */
public class InternalMetaAnalysis {
	
	//Check why it ran over the number of observed permutations!
	//What is the SNP index?
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private int[] snpIndex;
	private String[] snpList;
	private InternalMetaAnalysisDataset dataset;
	private final InternalMetaAnalysisSettings settings;
	
	THashMap<String, ArrayList<String>> traitMap;
	ArrayList<String> outputEntries = null;
	
	public static void main(String[] args) {
		
		String settingsFile = null;
		String textToReplace = null;
		String replaceTextWith = null;
		
		if (args.length == 3) {
			settingsFile = args[0];
			textToReplace = args[1];
			replaceTextWith = args[2];
			
		} else if (args.length == 1) {
			settingsFile = args[0];
			
		} else {
			System.out.println("Usage of the internal meta-analysis: settings.xml");
			
			System.exit(-1);
		}
		settingsFile = "D:\\Data\\eqtlgen\\InternalMetaFix\\settings_Illumina_PC_corrected_DILGOM_14082017v1012.xml";
		new InternalMetaAnalysis(settingsFile, textToReplace, replaceTextWith);
	}
	
	public InternalMetaAnalysis(String settingsFile, String textToReplace, String replaceTextWith) {
		// initialize settings
		settings = new InternalMetaAnalysisSettings();
		settings.parse(settingsFile, textToReplace, replaceTextWith);
		
		try {
			run();
		} catch (IOException ex) {
			Logger.getLogger(InternalMetaAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
		//ToDo, check for non overlapping sites.
		//Check allele assesed?
		
	}
	
	public void run() throws IOException {
		
		String outdir = settings.getOutput();
		
		traitMap = new THashMap<String, ArrayList<String>>();
		LinkedHashMap<String, Integer> traitLocationMap = new LinkedHashMap<String, Integer>();
		TextFile entryReader = new TextFile(settings.getProbeToFeature(), TextFile.R);
		String row;
		int index = 0;
		while ((row = entryReader.readLine()) != null) {
//            System.out.println(row);
			String[] parts = TAB_PATTERN.split(row);
			
			if (!traitLocationMap.containsKey(parts[1])) {
				traitLocationMap.put(parts[1], index);
				index++;
			}
			
			if (!traitMap.containsKey(parts[0])) {
				traitMap.put(parts[0], new ArrayList<String>());
			}
			traitMap.get(parts[0]).add(parts[1]);
		}
		outputEntries = new ArrayList<String>(traitLocationMap.keySet());
		
		if (!Gpio.exists(settings.getOutput())) {
			Gpio.createOuputDir(new File(settings.getOutput()));
		}
		
		
		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Will try to make use of " + cores + " CPU cores");
		System.out.println();
		ExecutorService threadPool = Executors.newFixedThreadPool(cores);
		int nrtasks = 0;
		for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
			nrtasks++;
		}
		MultiThreadProgressBar mtpb = new MultiThreadProgressBar(nrtasks);
		int taskid = 0;
		for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
			// run task
			InternalMetaAnalysisTask t = new InternalMetaAnalysisTask(settings,
					traitMap,
					outputEntries,
					traitLocationMap,
					permutation,
					outdir,
					mtpb,
					taskid);
			threadPool.submit(t);
			taskid++;
		}
		
		while (!mtpb.allCompleted()) {
			mtpb.display();
			try {
				Thread.sleep(2000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		
		threadPool.shutdown();
		
		
	}
	
	
	// index the probes
//    private void createProbeIndex() throws IOException {
//
//        HashSet<String> confineToTheseProbes = new HashSet<String>(traitMap.keySet());
//
//        probeIndex = new Integer[dataset.getProbeList().length];
//
//        String[] probes = dataset.getProbeList();
//
//        for (int p = 0; p < probes.length; p++) {
//
//            String t = traitMap.get(probes[p]);
////            System.out.println(probes[p]);
////            System.out.println(t);
//            Integer index = outputEntries.get(t);
//
//            if (index != null && confineToTheseProbes.contains(probes[p])) {
//                probeIndex[p] = index;
//            } else {
//                probeIndex[p] = null;
//            }
//        }
//    }
//
	private void writeBinaryResult(String snpname,
								   double hwe,
								   double cr,
								   double maf,
								   int numberCalled,
								   String alleles,
								   String minorAllele,
								   String alleleassessed,
								   double[] datasetZScores,
								   String[] probeNames,
								   BinaryFile outfile,
								   TextFile snpfile) throws IOException {
		StringBuilder sb = null;
		for (int p = 0; p < datasetZScores.length; p++) {
			float z = (float) datasetZScores[p];
			// System.out.println(p + "\t" + alleleassessed + "\t" + m_probeList[p] + "\t" + z + "\t" + currentWP.getFlipSNPAlleles()[d]);
			if (probeNames != null) {
				// take into account that not all probes have been tested..
				String probeName = probeNames[p];
				outfile.writeFloat(z);
				if (sb == null) {
					sb = new StringBuilder();
				} else {
					sb.append("\t");
				}
				sb.append(probeName);
			} else {
				outfile.writeFloat(z);
			}
		}
		
		StringBuilder buffer = new StringBuilder();
		buffer.append(snpname)
				.append("\t")
				.append(alleles)
				.append("\t")
				.append(minorAllele)
				.append("\t")
				.append(alleleassessed)
				.append("\t")
				.append(numberCalled)
				.append("\t")
				.append(maf)
				.append("\t")
				.append(hwe)
				.append("\t")
				.append(cr)
				.append("\t")
				.append(datasetZScores.length)
				.append("\t");
		if (sb != null) {
			buffer.append(sb.toString());
		} else {
			buffer.append("-");
		}
		
		snpfile.writeln(buffer.toString());
	}
	
}
