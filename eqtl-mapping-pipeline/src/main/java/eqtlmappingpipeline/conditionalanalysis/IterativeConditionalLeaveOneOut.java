/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.conditionalanalysis;

import eqtlmappingpipeline.metaqtl3.EQTLRegression;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import eqtlmappingpipeline.metaqtl3.containers.Settings;
import gnu.trove.set.hash.THashSet;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.Console;
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
public class IterativeConditionalLeaveOneOut extends MetaQTL3 {

    private boolean limitConseqcutiveIterationsOnSignificantGenes = true;

    private int startIter = 1;
    private int stopIter = 2;
    private int takeEQTLsUpToIter = 2;
    boolean useOLS = true;

    public static void waitForEnter(String message, Object... args) {
        Console c = System.console();
        if (c != null) {
            // printf-like arguments
            if (message != null)
                c.format(message, args);
            c.format("\nPress ENTER to proceed.\n");
            c.readLine();
        }
    }

    public void run(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
                    String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
                    String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads) throws IOException, Exception {

        System.out.println("Iterative leave one out analysis");
        if (xmlSettingsFile == null) {
            System.err.println("Sorry you have to supply a settings file.");
            System.exit(-1);
        }

        m_settings = new Settings();
        m_settings.load(xmlSettingsFile);
        String origOutputDir = m_settings.outputReportsDir;
        double fdrthreshold = m_settings.fdrCutOff;

        // check whether the needed files are there
        boolean filesOK = true;
        for (int i = startIter; i < (stopIter + 1); i++) {
            String eqtlsToTestFile = origOutputDir + "/Iteration" + i + "/eQTLProbesFDR" + fdrthreshold + "-ProbeLevel.txt.gz";
            if (m_settings.fdrType.equals(FDR.FDRMethod.FULL)) {
                eqtlsToTestFile = origOutputDir + "/Iteration" + i + "/eQTLProbesFDR" + fdrthreshold + ".txt.gz";
            } else if (m_settings.fdrType.equals(FDR.FDRMethod.SNPLEVEL)) {
                eqtlsToTestFile = origOutputDir + "/Iteration" + i + "/eQTLProbesFDR" + fdrthreshold + "-SNPLevel.txt.gz";
            } else if (m_settings.fdrType.equals(FDR.FDRMethod.GENELEVEL)) {
                eqtlsToTestFile = origOutputDir + "/Iteration" + i + "/eQTLProbesFDR" + fdrthreshold + "-GeneLevel.txt.gz";
            }
            if (!Gpio.exists(eqtlsToTestFile)) {
                System.out.println("Missing file: " + eqtlsToTestFile);
                filesOK = false;
            }
        }

        if (!filesOK) {
            System.out.println("Stopping because not all files are present for requested iterations.");
            System.exit(-1);
        }

//        initialize(xmlSettingsFile, texttoreplace, texttoreplacewith, ingt, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, null, null, null, true, true, null, null, null);

        m_settings.provideBetasAndStandardErrors = true;
        m_settings.provideFoldChangeData = true;
        m_settings.displayWarnings = false;
        boolean prevIterHasSignResults = true;
        boolean saveIntermediateResiduals = m_settings.regressOutEQTLEffectsSaveOutput;
        THashSet<String> originalProbeConfine = m_settings.tsProbesConfine;
        EQTLRegression eqr = new EQTLRegression();

        // this analysis dumps all relevant eQTLs to disk.
        m_settings.dumpeverythingtodisk = true;
        m_settings.nrPermutationsFDR = 0;
        m_settings.skipFDRCalculation = true;
        m_settings.maxNrMostSignificantEQTLs = 10;

        for (int i = startIter; i < (stopIter + 1); i++) {
            System.out.println("Starting iteration " + i);
            // determine probes to confine on; retest those genes that were significant in this iteration
            ArrayList<Integer> regressedIters = new ArrayList<>();
            m_settings.tsProbesConfine = collectEQTLProbes(origOutputDir, i, fdrthreshold);
            ArrayList<Pair<String, String>> toRegress = new ArrayList<>();
            for (int j = startIter; j < (takeEQTLsUpToIter + 1); j++) {
                if (i != j) {
                    toRegress.addAll(collectEQTLs(origOutputDir, j, fdrthreshold));
                    regressedIters.add(j);
                }
            }

            System.out.println("Iteration: " + i + "\t" + m_settings.tsProbesConfine.size() + " genes to test, will regress out " + toRegress.size() + " QTLs");
            if (m_settings.tsProbesConfine.isEmpty() || toRegress.isEmpty()) {
                System.err.println("Error: no work remaining.");
                System.exit(-1);
            }

            String regressedItersStr = Strings.concat(Primitives.toPrimitiveArr(regressedIters.toArray(new Integer[0])), Strings.dash);

            m_settings.outputReportsDir = origOutputDir + "/Iteration" + i + "-OtherIterationsRemoved/";
            m_settings.plotOutputDirectory = origOutputDir + "/Iteration" + i + "-OtherIterationsRemoved/";
            Gpio.createDir(m_settings.plotOutputDirectory);
            Gpio.createDir(m_settings.outputReportsDir);

            // reset the datasets
//            waitForEnter("Press enter to initialize");
            reinit();

            // regress significant eQTLs
            try {
                eqr.setLog(m_settings.outputReportsDir, i);
                eqr.regressOutEQTLEffects(toRegress, m_gg, useOLS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(-1);
            }

            if (saveIntermediateResiduals) {
                exportResidualsToDisk(origOutputDir, i);
            }

            numAvailableInds = 0;

//					 recalculate mean and SD
            for (int d = 0; d < m_gg.length; d++) {
//						if (!m_settings.performParametricAnalysis) {
//						m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
//						}
//						m_gg[i].getExpressionData().calcAndSubtractMean();
//						m_gg[i].getExpressionData().calcMeanAndVariance();
                numAvailableInds += m_gg[d].getExpressionToGenotypeIdArray().length;
            }

            // then map eQTLs
            mapEQTLs();

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


            String foutname = origOutputDir + ds.getSettings().name + "-EQTLEffectsRemoved-Iteration-" + iter + ".txt.gz";
            if (iter == 0) {
                foutname = origOutputDir + ds.getSettings().name + "-EQTLEffectsRemoved-Iteration-Last.txt.gz";
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

    private ArrayList<Pair<String, String>> collectEQTLs(String origOutputDir, int iteration, double fdr) throws IOException {

        HashSet<Pair<String, String>> eqtls = new HashSet<Pair<String, String>>();

        String iterationFile = origOutputDir + "/Iteration" + iteration + "/eQTLProbesFDR" + fdr + "-ProbeLevel.txt.gz";

        if (m_settings.fdrType.equals(FDR.FDRMethod.FULL)) {
            iterationFile = origOutputDir + "/Iteration" + (iteration) + "/eQTLProbesFDR" + fdr + ".txt.gz";
        } else if (m_settings.fdrType.equals(FDR.FDRMethod.SNPLEVEL)) {
            iterationFile = origOutputDir + "/Iteration" + (iteration) + "/eQTLProbesFDR" + fdr + "-SNPLevel.txt.gz";
        } else if (m_settings.fdrType.equals(FDR.FDRMethod.GENELEVEL)) {
            iterationFile = origOutputDir + "/Iteration" + (iteration) + "/eQTLProbesFDR" + fdr + "-GeneLevel.txt.gz";
        }

        int ctr = 0;
        System.out.println("Trying to collect QTLs from " + iterationFile);
        if (Gpio.exists(iterationFile)) {
            TextFile tf = new TextFile(iterationFile, TextFile.R);
            tf.readLineElems(TextFile.tab);
            String[] elems = tf.readLineElems(TextFile.tab);

            while (elems != null) {
                eqtls.add(new Pair<String, String>(elems[1], elems[4]));
                ctr++;
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
        }
        System.out.println("Iteration " + iteration + " has " + ctr + " effects. Total sofar: " + eqtls.size());


        ArrayList<Pair<String, String>> pairs = new ArrayList<Pair<String, String>>();
        pairs.addAll(eqtls);
        return pairs;
    }

    private THashSet<String> collectEQTLProbes(String origOutputDir, int currentIteration, double fdr) throws IOException {

        THashSet<String> output = new THashSet<String>();
        String iterationFile = origOutputDir + "/Iteration" + (currentIteration) + "/eQTLProbesFDR" + fdr + "-ProbeLevel.txt.gz";
        if (m_settings.fdrType.equals(FDR.FDRMethod.FULL)) {
            iterationFile = origOutputDir + "/Iteration" + (currentIteration) + "/eQTLProbesFDR" + fdr + ".txt.gz";
        } else if (m_settings.fdrType.equals(FDR.FDRMethod.SNPLEVEL)) {
            iterationFile = origOutputDir + "/Iteration" + (currentIteration) + "/eQTLProbesFDR" + fdr + "-SNPLevel.txt.gz";
        } else if (m_settings.fdrType.equals(FDR.FDRMethod.GENELEVEL)) {
            iterationFile = origOutputDir + "/Iteration" + (currentIteration) + "/eQTLProbesFDR" + fdr + "-GeneLevel.txt.gz";
        }

        System.out.println("Trying to collect genes/probes from: " + iterationFile);
        if (Gpio.exists(iterationFile)) {
            TextFile tf = new TextFile(iterationFile, TextFile.R);
            tf.readLineElems(TextFile.tab);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                output.add(elems[4]);
                elems = tf.readLineElems(TextFile.tab);
            }
            System.out.println("Iteration " + (currentIteration - 1) + " has " + output.size() + " significant probes.");
        }


        return output;
    }

    public void setStartIter(Integer startiter) {

        this.startIter = startiter;
    }

    public void setStopIter(Integer stopIter) {
        this.stopIter = stopIter;
        if (this.takeEQTLsUpToIter < this.stopIter) {
            this.takeEQTLsUpToIter = this.stopIter;
        }
    }

    public void setTakeEQTLsUpToIter(Integer takeEQTLsUpToIter) {
        this.takeEQTLsUpToIter = takeEQTLsUpToIter;
    }

    public void setLimitConsecutiveIterationsToSignificantGenes(boolean limitConseqcutiveIterationsOnSignificantGenes) {
        this.limitConseqcutiveIterationsOnSignificantGenes = limitConseqcutiveIterationsOnSignificantGenes;

    }
}
