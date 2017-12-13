/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.conditionalanalysis;

import eqtlmappingpipeline.metaqtl3.EQTLRegression;
import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import gnu.trove.set.hash.THashSet;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;

/**
 * @author harm-jan
 */
public class IterativeConditionalAnalysis extends MetaQTL3 {
	
	public void run(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
					String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
					String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads) throws IOException, Exception {
		
		
		initialize(xmlSettingsFile, texttoreplace, texttoreplacewith, null, null, ingt, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, null, null, null, true, true, null, null, null);
		
		double fdrthreshold = m_settings.fdrCutOff;
		m_settings.provideBetasAndStandardErrors = true;
		m_settings.provideFoldChangeData = true;
		String origOutputDir = m_settings.outputReportsDir;
		boolean prevIterHasSignResults = true;
		int iteration = 1;
		
		EQTLRegression eqr = new EQTLRegression();
		
		while (prevIterHasSignResults) {
			m_settings.outputReportsDir = origOutputDir + "/Iteration" + iteration + "/";
			m_settings.plotOutputDirectory = origOutputDir + "/Iteration" + iteration + "/";
			Gpio.createDir(m_settings.plotOutputDirectory);
			Gpio.createDir(m_settings.outputReportsDir);
			
			System.out.println("Iteration: " + iteration);
			
			if (iteration == 1) {
				mapEQTLs();
			} else {
				// check whether there were significant results in the previous iteration
				String efilename = origOutputDir + "/Iteration" + (iteration - 1) + "/eQTLProbesFDR" + fdrthreshold + "-ProbeLevel.txt";
				if (!Gpio.exists(efilename)) {
					System.err.println("Previous iteration (" + (iteration - 1) + ") did not have any significant results.");
					System.err.println("File: " + efilename + " does not exist.");
					prevIterHasSignResults = false;
				} else {
					TextFile tf = new TextFile(efilename, TextFile.R);
					int nrlns = tf.countLines();
					tf.close();
					
					if (nrlns == 1) {
						System.err.println("Previous iteration (" + (iteration - 1) + ") did not have any significant results.");
						System.err.println("File: " + efilename + " has no entries.");
						prevIterHasSignResults = false;
					} else {
						System.err.println("Previous iteration (" + (iteration - 1) + ") yielded " + (nrlns - 1) + " significant results.");
					}
				}
				
				if (prevIterHasSignResults) {
					// get the list of eQTLs to regress out...
					ArrayList<Pair<String, String>> toRegress = collectEQTLs(origOutputDir, iteration, fdrthreshold);
					
					// get the significant probes from the previous run
					m_settings.tsProbesConfine = collectEQTLProbes(origOutputDir, iteration, fdrthreshold);
					
					// reset the datasets
					reinit();
					
					// regress significant eQTLs
					eqr.regressOutEQTLEffects(toRegress, m_gg);
					
					// recalculate mean and SD
					for (int i = 0; i < m_gg.length; i++) {
						if (!m_settings.performParametricAnalysis) {
							m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
						}
						m_gg[i].getExpressionData().calcAndSubtractMean();
						m_gg[i].getExpressionData().calcMeanAndVariance();
						numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
					}
					
					// then map eQTLs
					mapEQTLs();
				}
				
				
			}
			
			iteration++;
			
			
		}
	}
	
	private void reinit() throws IOException, Exception {
		m_gg = null;
		
		int numDatasets = m_settings.datasetSettings.size();
		m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
		numAvailableInds = 0;
		for (int i = 0; i < numDatasets; i++) {
			
			System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
			System.out.println(ConsoleGUIElems.LINE);
			m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i));
		}
		
		
		for (int i = 0; i < numDatasets; i++) {
			if (!m_settings.performParametricAnalysis) {
				m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
			}
			m_gg[i].getExpressionData().calcAndSubtractMean();
			m_gg[i].getExpressionData().calcMeanAndVariance();
			numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
		}
		
		if (m_settings.regressOutEQTLEffectFileName != null && m_settings.regressOutEQTLEffectFileName.trim().length() > 0) {
			EQTLRegression eqr = new EQTLRegression();
			eqr.regressOutEQTLEffects(m_settings.regressOutEQTLEffectFileName, false, m_gg);
			for (int i = 0; i < numDatasets; i++) {
				if (!m_settings.performParametricAnalysis) {
					m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
				}
				m_gg[i].getExpressionData().calcAndSubtractMean();
				m_gg[i].getExpressionData().calcMeanAndVariance();
				numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
			}
		}
		
		System.out.println(ConsoleGUIElems.LINE);
		System.out.println("");
		
		
		
		System.out.println("Accumulating available data...");
		System.out.print(ConsoleGUIElems.LINE);
		
		createSNPList();
		createProbeList();
		
		// create WorkPackage objects
		determineSNPProbeCombinations();
		
		if (m_workPackages == null || m_workPackages.length == 0) {
			System.err.println("Error: No work detected");
			System.exit(0);
		}
		
		// determine number of threadss
		if (m_settings.nrThreads == null) {
			m_settings.nrThreads = Runtime.getRuntime().availableProcessors();
		} else {
			int numProcs = Runtime.getRuntime().availableProcessors();
			if (m_settings.nrThreads > numProcs || m_settings.nrThreads < 1) {
				m_settings.nrThreads = numProcs;
			}
		}
		
		if (m_workPackages.length < m_settings.nrThreads) {
			m_settings.nrThreads = m_workPackages.length;
		}
		printSummary();
	}
	
	private ArrayList<Pair<String, String>> collectEQTLs(String origOutputDir, int currentIteration, double fdr) throws IOException {
		
		HashSet<Pair<String, String>> eqtls = new HashSet<Pair<String, String>>();
		for (int iteration = 1; iteration < currentIteration; iteration++) {
			String iterationFile = origOutputDir + "/Iteration" + iteration + "/eQTLProbesFDR" + fdr + "-ProbeLevel.txt";
			TextFile tf = new TextFile(iterationFile, TextFile.R);
			tf.readLineElems(TextFile.tab);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				eqtls.add(new Pair<String, String>(elems[1], elems[4]));
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}
		ArrayList<Pair<String, String>> pairs = new ArrayList<Pair<String, String>>();
		pairs.addAll(eqtls);
		return pairs;
	}
	
	private THashSet<String> collectEQTLProbes(String origOutputDir, int currentIteration, double fdr) throws IOException {
		
		THashSet<String> output = new THashSet<String>();
		String iterationFile = origOutputDir + "/Iteration" + (currentIteration - 1) + "/eQTLProbesFDR" + fdr + "-ProbeLevel.txt";
		TextFile tf = new TextFile(iterationFile, TextFile.R);
		tf.readLineElems(TextFile.tab);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			output.add(elems[4]);
			elems = tf.readLineElems(TextFile.tab);
		}
		return output;
	}
}
