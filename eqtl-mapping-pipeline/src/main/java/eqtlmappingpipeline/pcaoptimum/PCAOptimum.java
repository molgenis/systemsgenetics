/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.pcaoptimum;

import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import eqtlmappingpipeline.metaqtl3.MetaQTL3Settings;
import eqtlmappingpipeline.normalization.Normalizer;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class PCAOptimum extends MetaQTL3 {

    protected String inexpplatform;
    protected String inexpannot;
    protected String ingt;
    protected String gte;
    protected Integer m_threads;
    protected int permutations = 10;
    protected boolean covariatesremoved = false;
    protected String cissnps;
    protected String transsnps;
    private boolean performEigenVectorQTLMapping;

    public void setCovariatesRemoved(boolean b) {
        covariatesremoved = b;
    }

    public void setSNPSets(String cissnps, String transsnps) {
        this.cissnps = cissnps;
        this.transsnps = transsnps;
    }

    public void setPerformpcqtlNormalization(boolean performEigenvectorQTLMapping) {
        this.performEigenVectorQTLMapping = performEigenvectorQTLMapping;
    }

    @Override
    public void initialize(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
            String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
            String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads, Integer maxNrResults, String regressouteqtls, String snpprobecombofile) throws IOException, Exception {
        if (!out.endsWith("/")) {
            out += "/";
        }
        if (!Gpio.exists(out)) {
            Gpio.createDir(out);
        }

        permutations = perm;

        String origInExp = inexp;

        m_settings = new MetaQTL3Settings();
        int nrProcs = Runtime.getRuntime().availableProcessors();
        if (threads != null && threads > 0 && threads <= nrProcs) {
            //
            m_threads = threads;
        } else {
            if (threads != null && threads > nrProcs) {
                System.out.println("The number of threads you set using the command line is not correct for your system. You set " + threads + " threads, while your machine has " + nrProcs + " processors");
            }
            threads = nrProcs;
        }

        m_threads = threads;
        int round = 0;


        HashSet<String> cisSnpsToTest = new HashSet<String>();
        HashSet<String> transSnpsToTest = new HashSet<String>();



        if (cissnps != null) {
            cis = true;
        } else {
            cis = false;
        }

        if (transsnps != null) {
            trans = true;
        } else {
            trans = false;
        }

        if (cis) {
            System.out.println("Loading cis SNP set from: " + cissnps);
            TextFile tf = new TextFile(cissnps, TextFile.R);
            cisSnpsToTest.addAll(tf.readAsArrayList());
            tf.close();
        }

        if (trans) {
            System.out.println("Loading trans SNP set from: " + transsnps);
            TextFile tf = new TextFile(transsnps, TextFile.R);
            transSnpsToTest.addAll(tf.readAsArrayList());
            tf.close();
        }

        //String nextInExp = origInExp + ".QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz";
        String outputdir = out + "Cis-0PCAsRemoved";
        if (!outputdir.endsWith("/")) {
            outputdir += "/";
        }
        if (!Gpio.exists(outputdir)) {
            Gpio.createDir(outputdir);
        }

        // set standard cis-settings
        this.inexpannot = inexpannot;
        this.inexpplatform = inexpplatform;
        this.ingt = ingt;
        this.gte = gte;

        String parentDir = Gpio.getParentDir(origInExp);
        System.out.println("Looking for PCA corrected files in folder: " + parentDir);
        String[] fileList = Gpio.getListOfFiles(parentDir);
        ArrayList<Integer> pcs = new ArrayList<Integer>();
        HashMap<Integer, String> stepToFile = new HashMap<Integer, String>();
        for (String f : fileList) {
            if (f.toLowerCase().contains("pcasoversamplesremoved") && !f.toLowerCase().contains("geneticvectorsnotremoved")) {
                String[] fileParts = f.split("\\.");
                for (String p : fileParts) {
                    if (p.toLowerCase().contains("pcasoversamplesremoved")) {
                        Integer pc = Integer.parseInt(p.toLowerCase().replace("pcasoversamplesremoved", ""));
                        pcs.add(pc);
                        stepToFile.put(pc, f);
                        System.out.println("Found file for PC: " + pc + "\t" + f);
                        break;
                    }
                }
            }
        }

        if (pcs.isEmpty()) {
            System.out.println("No PCA corrected files."
                    + "\n Please first run the normalization procedure.");
            System.exit(0);
        }

        Collections.sort(pcs);

        int max = 0;
        int stepSize = 0;
        for (int i = 0; i < (pcs.size() - 1); i++) {
            if (i == 0) {
                if (pcs.get(pcs.size() - 1) > max) {
                    max = pcs.get(pcs.size() - 1);
                }
            }
            if (pcs.get(i) > max) {
                max = pcs.get(i);
            }
            stepSize += pcs.get(i + 1) - pcs.get(i);
        }

        if (pcs.size() > 1) {
            if ((((double) stepSize / (pcs.size() - 1)) % 1) != 0) {
                System.out.println("Step size is invalid."
                        + "\n Please look in to the input directory for missing files");
                System.out.println((((double) stepSize / (pcs.size() - 1)) % 1));
                System.exit(0);
            }
            stepSize = (int) ((double) stepSize / (pcs.size() - 1));
        } else {
            stepSize = pcs.get(0);
            max = pcs.get(0);
        }



        System.out.println("Determined max: " + max);
        System.out.println("Determined stepsize: " + stepSize);

        if (performEigenVectorQTLMapping) {
            performeQTLMappingOverEigenvectorMatrixAndReNormalize(origInExp, out, stepSize, max);

            fileList = Gpio.getListOfFiles(parentDir);
            pcs = new ArrayList<Integer>();
            stepToFile = new HashMap<Integer, String>();
            for (String f : fileList) {

                if (f.toLowerCase().contains("pcasoversamplesremoved-geneticvectorsnotremoved")) {
                    String[] fileParts = f.split("\\.");
                    for (String p : fileParts) {
                        if (p.toLowerCase().contains("pcasoversamplesremoved-geneticvectorsnotremoved")) {
                            Integer pc = Integer.parseInt(p.toLowerCase().replace("pcasoversamplesremoved-geneticvectorsnotremoved", ""));
                            pcs.add(pc);
                            stepToFile.put(pc, f);
                            System.out.println("Found file for PC: " + pc + "\t" + f);
                            break;
                        }
                    }
                }
            }
        }

        String nextInExp = "";
        for (int pca = 0; pca <= max; pca += stepSize) {
            if (pca == 0) {
                //nextInExp = origInExp + ".QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz";
                nextInExp = origInExp;
            } else {
//                String tmpName = origInExp.replaceAll(".txt", "");
//                tmpName = tmpName.replaceAll(".gz", "");
//                nextInExp = tmpName + "." + pca + "PCAsOverSamplesRemoved.txt.gz";
                nextInExp = parentDir + "/" + stepToFile.get(pca);
            }

            if (!Gpio.exists(nextInExp)) {
                System.err.println("Could not find file for pca: " + pca + "\t" + nextInExp);
            } else {
                if (cis) {
                    String outputDir = out + "Cis-" + pca + "PCAsRemoved/";
                    if (performEigenVectorQTLMapping && pca > 0) {
                        outputDir = out + "Cis-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/";
                    }
                    if ((pca == 0 && !Gpio.exists(outputDir + "eQTLProbesFDR0.05.txt")) || pca > 0) {
                        performeQTLMapping(true, false, nextInExp, outputDir, cisSnpsToTest, null, threads);
                        cleanup();
                    }
                }
                if (trans) {
                    String outputDir = out + "Trans-" + pca + "PCAsRemoved/";
                    if (performEigenVectorQTLMapping && pca > 0) {
                        outputDir = out + "Trans-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/";
                    }
                    if ((pca == 0 && !Gpio.exists(outputDir + "eQTLProbesFDR0.05.txt")) || pca > 0) {
                        performeQTLMapping(false, true, nextInExp, outputDir, transSnpsToTest, null, threads);
                        cleanup();
                    }
                }
            }

        }

        PCAOptimumInventorize pi = new PCAOptimumInventorize();
        if (performEigenVectorQTLMapping) {
            pi.inventorypcqtl(out, cis, trans, max, stepSize);
        } else {
            pi.inventory(out, cis, trans, max, stepSize);
        }
    }

    protected void init() throws IOException, Exception {
        int numDatasets = 1;
        m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
        numAvailableInds = 0;
        for (int i = 0; i < numDatasets; i++) {

            System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
            System.out.println(ConsoleGUIElems.LINE);
            m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i));
            if (!m_settings.performParametricAnalysis) {

                m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
            }
            numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
            System.out.println(ConsoleGUIElems.LINE);
            System.out.println("");
        }


        // sort datasets on size to increase efficiency of random reads..
        // for some reason, it is faster to load the largest dataset first.
        Arrays.sort(m_gg, Collections.reverseOrder());

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


        // create threadss
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

    protected void cleanup() {
        for (int i = 0; i < m_gg.length; i++) {
            this.m_gg[i] = null;
        }

        this.m_probeList = null;
        this.m_probeTranslationTable = null;
        this.m_snpList = null;
        this.m_snpTranslationTable = null;
        this.m_workPackages = null;
    }

    protected void performeQTLMapping(boolean cis, boolean trans, String inFile, String out, HashSet<String> snpsToTest, HashSet<String> probesToTest, int threads) throws IOException, Exception {
//
        String nextInExp = inFile;

        String outputdir = out;
        if (!Gpio.exists(outputdir)) {
            Gpio.createDir(outputdir);
        }

        // set output dir
        // set standard cis-settings
        m_settings = new MetaQTL3Settings();
        TriTyperGeneticalGenomicsDatasetSettings s = new TriTyperGeneticalGenomicsDatasetSettings();

        s.name = "Dataset";
        s.expressionLocation = nextInExp;
        s.expressionplatform = inexpplatform;
        s.probeannotation = inexpannot;
        s.genotypeLocation = ingt;
        s.genotypeToExpressionCoupling = gte;
        s.cisAnalysis = cis;
        s.transAnalysis = trans;

        if (probesToTest != null) {
            s.tsProbesConfine = probesToTest;
        }
        m_settings.datasetSettings = new ArrayList<TriTyperGeneticalGenomicsDatasetSettings>();
        m_settings.datasetSettings.add(s);

        if (cis) {
            m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance = 250000;
        } else {
            m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance = 5000000;
        }

        m_settings.nrThreads = threads;
        m_settings.cisAnalysis = cis;
        m_settings.transAnalysis = trans;

        m_settings.nrPermutationsFDR = permutations;
        m_settings.tsSNPsConfine = snpsToTest;
        if (probesToTest != null) {
            m_settings.tsProbesConfine = probesToTest;
        }
        m_settings.outputReportsDir = outputdir;
        m_settings.createTEXTOutputFiles = true;
        m_settings.createBinaryOutputFiles = false;

        init();

        // set standard trans settings
        super.mapEQTLs();
        cleanup();
    }

    protected void compareZScores(EQTL[] ciseqtls, EQTL[] transeqtls, EQTL[] originalCisEQTLs, EQTL[] originalTransEQTLs, String out, int pca) {
        HashMap<String, EQTL> origCis = new HashMap<String, EQTL>();
        HashMap<String, EQTL> origTrans = new HashMap<String, EQTL>();

        for (EQTL e : originalCisEQTLs) {
            origCis.put(e.getRsName() + "///" + e.getProbe(), e);
        }

        for (EQTL e : originalTransEQTLs) {
            origTrans.put(e.getRsName() + "///" + e.getProbe(), e);
        }

        int cisShared = 0;
        int transShared = 0;

        ArrayList<Double> cisXAL = new ArrayList<Double>();
        ArrayList<Double> cisYAL = new ArrayList<Double>();
        ArrayList<Double> transXAL = new ArrayList<Double>();
        ArrayList<Double> transYAL = new ArrayList<Double>();

        for (EQTL e : ciseqtls) {
            EQTL origE = origTrans.get(e.getRsName() + "///" + e.getProbe());
            if (origE != null) {
                cisShared++;
                cisXAL.add(e.getZscore());
                cisYAL.add(origE.getZscore());
            }
        }

        for (EQTL e : transeqtls) {
            EQTL origE = origTrans.get(e.getRsName() + "///" + e.getProbe());
            if (origE != null) {
                transShared++;
                transXAL.add(e.getZscore());
                transYAL.add(origE.getZscore());
            }
        }


        double[] cisx = toArray(cisXAL);
        double[] cisy = toArray(cisYAL);
        double[] transx = toArray(transXAL);
        double[] transy = toArray(transYAL);

        PCAOptimumPlot p = new PCAOptimumPlot(500, 1000, true, out);
        p.plot(cisx, cisy, transx, transy, originalCisEQTLs.length, ciseqtls.length, cisShared, originalTransEQTLs.length, transeqtls.length, transShared);
        p.draw(out);

    }

    protected double[] toArray(ArrayList<Double> x) {
        double[] y = new double[x.size()];
        int i = 0;
        for (Double z : x) {
            y[i] = z;
            i++;
        }
        return y;
    }

    private void performeQTLMappingOverEigenvectorMatrixAndReNormalize(String origInExp, String out, int stepSize, int max) throws IOException, Exception {
        // Eigenvector mapping
        TextFile tf = new TextFile(origInExp, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        int nrCols = header.length;
        tf.close();

        int nrToRemove = nrCols;

        HashSet<String> probesToTest = new HashSet<String>();
        for (int i = 1; i < nrToRemove; i++) {
            probesToTest.add("Comp" + i);
        }

        String parentDir = Gpio.getParentDir(origInExp);
        parentDir += "/";
        String[] files = Gpio.getListOfFiles(parentDir);
        String pcaOverSampleEigenVectorsTransposedFile = null;
        for (String file : files) {
            if (file.toLowerCase().contains("pcaoversampleseigenvectorstransposed")) {
                pcaOverSampleEigenVectorsTransposedFile = parentDir + "" + file;
                break;
            }
        }

        if (pcaOverSampleEigenVectorsTransposedFile == null) {
            System.out.println("Error! Could not find " + pcaOverSampleEigenVectorsTransposedFile);
            System.exit(0);
        }

        // ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.PCAOverSamplesEigenvectorsTransposed


        String nextInExp = pcaOverSampleEigenVectorsTransposedFile;
        performeQTLMapping(true, true, nextInExp, out + "CisTrans-PCAEigenVectors/", null, probesToTest, m_threads);
        cleanup();

        eQTLTextFile etf = new eQTLTextFile(out + "CisTrans-PCAEigenVectors/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
        EQTL[] eigenvectorEQTLs = etf.read();
        etf.close();

        HashSet<Integer> geneticEigenVectors = new HashSet<Integer>();
        for (EQTL e : eigenvectorEQTLs) {
            Double fdr = e.getFDR();
            if (fdr == null) {
                System.out.println("Error with eQTL file!: FDR == null: " + out + "CisTrans-EigenVectors/eQTLProbesFDR0.05.txt");
                System.exit(0);
            }
            if (fdr > 0) {
                break;
            } else {
                String probe = e.getProbe();
                Integer compId = Integer.parseInt(probe.replace("Comp", ""));
                // quick hack: component 1 captures population stratification information...
                if (compId > 1) {
                    geneticEigenVectors.add(compId);
                }
            }
        }

        System.out.println("Repeating PCA analysis, without removal of genetically controlled PCs");
        System.out.println("These PCs are under genetic control: " + Strings.concat(geneticEigenVectors.toArray(new Integer[0]), Strings.comma));
        System.out.println();
        if (geneticEigenVectors.size() > 0) {
            Normalizer n = new Normalizer();
            n.repeatPCAOmitCertainPCAs(geneticEigenVectors, origInExp, max, stepSize);
        } else {
            System.out.println("No PCA vectors seem to be genetically associated.\n"
                    + "To find the optimum number of PCs to use, rerun this command without --pcqtl, if you have not done so already.");
            System.exit(0);
        }

    }
}
