/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.conditionalanalysis;

import eqtlmappingpipeline.metaqtl3.EQTLRegression;
import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import eqtlmappingpipeline.normalization.Normalizer;
import eqtlmappingpipeline.util.QTLFileMerger;
import gnu.trove.set.hash.THashSet;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * @author harm-jan
 */
public class IterativeConditionalAnalysis extends MetaQTL3 {


	public static void main(String[] args) {

//		Normalizer z = new Normalizer();
//		try {
//			z.rank("D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM.txt.gz",
//					"D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM-ranked.txt.gz");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		System.exit(0);

		IterativeConditionalAnalysis s = new IterativeConditionalAnalysis();
		try {
			s.run("D:\\TMP\\geuvadistest\\metaqtlsettings.xml", null,
					null, null, null, null, null, null, null,
					true, false, 10, true, false, null, 4);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private Integer startIter = 2;
	boolean useOLS = true;

	public void run(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
					String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
					String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads) throws IOException, Exception {


		initialize(xmlSettingsFile, texttoreplace, texttoreplacewith, ingt, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, null, null, null, true, true, null, null, null);

		double fdrthreshold = m_settings.fdrCutOff;
		m_settings.provideBetasAndStandardErrors = true;
		m_settings.provideFoldChangeData = true;
		m_settings.displayWarnings = false;
		String origOutputDir = m_settings.outputReportsDir;
		boolean prevIterHasSignResults = true;
		int iteration = startIter;

		boolean saveIntermediateResiduals = m_settings.regressOutEQTLEffectsSaveOutput;

		EQTLRegression eqr = new EQTLRegression();

		while (prevIterHasSignResults) {
			m_settings.outputReportsDir = origOutputDir + "/Iteration" + iteration + "/";
			m_settings.plotOutputDirectory = origOutputDir + "/Iteration" + iteration + "/";
			Gpio.createDir(m_settings.plotOutputDirectory);
			Gpio.createDir(m_settings.outputReportsDir);

			System.out.println("Iteration: " + iteration);

			if (iteration == 1) {
//				if (saveIntermediateResiduals) {
//					exportResidualsToDisk(origOutputDir, iteration);
//				}
				mapEQTLs();
			} else {
				// check whether there were significant results in the previous iteration
				String efilename = origOutputDir + "/Iteration" + (iteration - 1) + "/eQTLProbesFDR" + fdrthreshold + "-ProbeLevel.txt.gz";
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
					try {
						eqr.regressOutEQTLEffects(toRegress, m_gg, useOLS);
					} catch (Exception e) {
						e.printStackTrace();
						System.exit(-1);
					}

					if (saveIntermediateResiduals) {
						exportResidualsToDisk(origOutputDir, iteration);
					}

					numAvailableInds = 0;

//					 recalculate mean and SD
					for (int i = 0; i < m_gg.length; i++) {
//						if (!m_settings.performParametricAnalysis) {
//						m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
//						}
//						m_gg[i].getExpressionData().calcAndSubtractMean();
//						m_gg[i].getExpressionData().calcMeanAndVariance();
						numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
					}

					// then map eQTLs
					mapEQTLs();

				}


			}

			iteration++;
		}


		System.out.println("Done with iterations. Will now save residual expression matrix.");


		// get the list of eQTLs to regress out...
		ArrayList<Pair<String, String>> toRegress = collectEQTLs(origOutputDir, iteration - 1, fdrthreshold);

		if (toRegress.isEmpty()) {
			System.out.println("No significant eQTLs found, and thus no need to save residual gene expression matrix.");
		} else {
			// get the significant probes from the previous run
			m_settings.tsProbesConfine = null;

			// reset the datasets
			reinit();

			numAvailableInds = 0;
			// recalculate mean and SD
			for (int i = 0; i < m_gg.length; i++) {
				if (!m_settings.performParametricAnalysis) {
					m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
				}
				m_gg[i].getExpressionData().calcAndSubtractMean();
				m_gg[i].getExpressionData().calcMeanAndVariance();
				numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
			}

			// regress significant eQTLs
			eqr.regressOutEQTLEffects(toRegress, m_gg, useOLS);

			// save the output
			exportResidualsToDisk(origOutputDir, 0);
		}
	}

	private void exportResidualsToDisk(String origOutputDir, int iter) throws Exception {
		for (int d = 0; d < m_gg.length; d++) {
			TriTyperGeneticalGenomicsDataset ds = m_gg[d];
			TriTyperExpressionData dsexp = ds.getExpressionData();
			double[][] matrix = dsexp.getMatrix();
			String[] probes = dsexp.getProbes();
			String[] individuals = dsexp.getIndividuals();
			String filename = ds.getSettings().expressionLocation;
			File f = new File(filename);
			String fname = f.getName();
			DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<>();
			dsout.setRowObjects(Arrays.asList(probes));
			dsout.setColObjects(Arrays.asList(individuals));
			dsout.setMatrix(matrix);


			String foutname = origOutputDir + fname + "-EQTLEffectsRemoved-Iteration-" + iter + ".txt.gz";
			if (iter == 0) {
				foutname = origOutputDir + fname + "-EQTLEffectsRemoved-Iteration-Last.txt.gz";
			}
			System.out.println("Saving expression file after removal of eQTL effects: " + foutname);
			dsout.save(foutname);
		}

	}

	private void reinit() throws IOException, Exception {
		m_gg = null;

		int numDatasets = m_settings.datasetSettings.size();
		m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
		numAvailableInds = 0;
		IntStream.range(0, numDatasets).parallel().forEach(i -> {
			System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
			System.out.println(ConsoleGUIElems.LINE);
			try {
				m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i), null, false);
			} catch (Exception e) {
				e.printStackTrace();
			}
		});

		AtomicInteger avinds = new AtomicInteger();

		IntStream.range(0, numDatasets).parallel().forEach(i -> {
			if (!m_settings.performParametricAnalysis) {
				m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
			}
			m_gg[i].getExpressionData().calcAndSubtractMean();
			m_gg[i].getExpressionData().calcMeanAndVariance();
			avinds.getAndAdd(m_gg[i].getExpressionToGenotypeIdArray().length);
//			numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
		});
		numAvailableInds = avinds.get();

		if (m_settings.regressOutEQTLEffectFileName != null && m_settings.regressOutEQTLEffectFileName.trim().length() > 0) {
			EQTLRegression eqr = new EQTLRegression();
			eqr.regressOutEQTLEffects(m_settings.regressOutEQTLEffectFileName, false, m_gg, useOLS);
			numAvailableInds = 0;
			AtomicInteger avinds2 = new AtomicInteger();

			IntStream.range(0, numDatasets).parallel().forEach(i -> {
				if (!m_settings.performParametricAnalysis) {
					m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
				}
				m_gg[i].getExpressionData().calcAndSubtractMean();
				m_gg[i].getExpressionData().calcMeanAndVariance();
				avinds.getAndAdd(m_gg[i].getExpressionToGenotypeIdArray().length);

			});
			numAvailableInds = avinds.get();
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
			String iterationFile = origOutputDir + "/Iteration" + iteration + "/eQTLProbesFDR" + fdr + "-ProbeLevel.txt.gz";
			TextFile tf = new TextFile(iterationFile, TextFile.R);
			tf.readLineElems(TextFile.tab);
			String[] elems = tf.readLineElems(TextFile.tab);
			int ctr = 0;
			while (elems != null) {
				eqtls.add(new Pair<String, String>(elems[1], elems[4]));
				ctr++;
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			System.out.println("Iteration " + iteration + " has " + ctr + " effects. Total sofar: " + eqtls.size());
		}

		ArrayList<Pair<String, String>> pairs = new ArrayList<Pair<String, String>>();
		pairs.addAll(eqtls);
		return pairs;
	}

	private THashSet<String> collectEQTLProbes(String origOutputDir, int currentIteration, double fdr) throws IOException {

		THashSet<String> output = new THashSet<String>();
		String iterationFile = origOutputDir + "/Iteration" + (currentIteration - 1) + "/eQTLProbesFDR" + fdr + "-ProbeLevel.txt.gz";
		TextFile tf = new TextFile(iterationFile, TextFile.R);
		tf.readLineElems(TextFile.tab);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			output.add(elems[4]);
			elems = tf.readLineElems(TextFile.tab);
		}
		System.out.println("Iteration " + (currentIteration - 1) + " has " + output.size() + " significant probes.");
		return output;
	}

	public void setStartIter(Integer startiter) {

		this.startIter = startiter;
	}
}
