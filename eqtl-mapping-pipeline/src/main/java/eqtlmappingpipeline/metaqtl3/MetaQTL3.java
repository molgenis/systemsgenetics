/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import com.lowagie.text.DocumentException;
import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import eqtlmappingpipeline.metaqtl3.containers.Result;
import umcg.genetica.math.stats.Descriptives;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLDotPlot;
import eqtlmappingpipeline.metaqtl3.graphics.EQTLPlotter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.util.RunTimer;

/**
 *
 * @author harmjan
 */
public class MetaQTL3 {

    protected MetaQTL3Settings m_settings;
    protected TriTyperGeneticalGenomicsDataset[] m_gg = null;
    protected String[] m_snpList;
    protected String[] m_probeList;
    protected Integer[][] m_probeTranslationTable;
    protected Integer[][] m_snpTranslationTable;
    protected int numAvailableInds;
    protected WorkPackage[] m_workPackages;

    public MetaQTL3() {
    }

    public MetaQTL3(MetaQTL3Settings settings) throws IOException, Exception {
        m_settings = settings;
        initialize(null, null, null, null, null, null, null, null, null, true, true, 0, true, false, null, null, null, null, null);
    }

    public void setOutputPlotThreshold(double d) {
        m_settings.plotOutputPValueCutOff = d;
        m_settings.plotOutputDirectory = m_settings.outputReportsDir;

    }

    public void initialize(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
            String ingt, String inexp, String inexpplatform, String inexpannot,
            String gte, String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads, Integer maxNrResults,
            String regressouteqtls, String snpprobecombofile) throws IOException, Exception {


        if (m_settings == null && xmlSettingsFile == null && ingt != null) {

            // check the settings
            boolean settingsOk = true;
            if (inexp == null || inexp.trim().length() == 0) {
                System.err.println("ERROR: you did not specify a gene expression file.");
                settingsOk = false;
            }

            if (inexpannot != null && inexpannot.trim().length() != 0) {
                if (inexpplatform == null || inexpplatform.trim().length() == 0) {
                    System.err.println("ERROR: you specified " + inexpannot + " but you did not specify the platform (using --inexpplatform)!");
                    settingsOk = false;
                }
            }

            if (out == null || out.trim().length() == 0) {
                System.err.println("ERROR: you did not specify an output directory.");
                settingsOk = false;
            }

            if (!settingsOk) {
                System.out.println();
                System.exit(0);
            }

            m_settings = new MetaQTL3Settings();
            TriTyperGeneticalGenomicsDatasetSettings s = new TriTyperGeneticalGenomicsDatasetSettings();

            s.name = "Dataset";
            s.expressionLocation = inexp;
            s.expressionplatform = inexpplatform;
            s.probeannotation = inexpannot;
            s.genotypeLocation = ingt;
            s.genotypeToExpressionCoupling = gte;
            s.cisAnalysis = cis;
            s.transAnalysis = trans;

            m_settings.cisAnalysis = cis;
            m_settings.transAnalysis = trans;

            boolean cistrans = false;
            if (m_settings.cisAnalysis && m_settings.transAnalysis) {
                m_settings.confineProbesToProbesMappingToAnyChromosome = true;
            }

            m_settings.datasetSettings = new ArrayList<TriTyperGeneticalGenomicsDatasetSettings>();



            m_settings.regressOutEQTLEffectFileName = regressouteqtls;
            m_settings.datasetSettings.add(s);
            m_settings.nrThreads = threads;
            m_settings.cisAnalysis = cis;
            m_settings.transAnalysis = trans;
            m_settings.nrPermutationsFDR = perm;
            if (!out.endsWith("/")) {
                out += "/";
            }
            if (!Gpio.exists(out)) {
                Gpio.createDir(out);
            }
            if (snpfile != null) {
                m_settings.tsSNPsConfine = new HashSet<String>();
                TextFile ts = new TextFile(snpfile, TextFile.R);
                m_settings.strConfineSNP = snpfile;
                m_settings.tsSNPsConfine.addAll(ts.readAsArrayList());
                ts.close();
            }

            if (snpprobecombofile != null) {
                m_settings.loadSNPProbeConfinement(snpprobecombofile);
            }

            m_settings.outputReportsDir = out;
            m_settings.createTEXTOutputFiles = textout;
            m_settings.createBinaryOutputFiles = binout;
            if (maxNrResults != null && maxNrResults > 0) {
                m_settings.maxNrMostSignificantEQTLs = maxNrResults;
            }

        } else if (m_settings == null && xmlSettingsFile != null) {
            // parse settings
            m_settings = new MetaQTL3Settings();
            m_settings.settingsTextReplaceWith = texttoreplacewith;
            m_settings.settingsTextToReplace = texttoreplace;
            m_settings.load(xmlSettingsFile);
        } else if (m_settings == null) {
            System.out.println("ERROR: No input specified");
            System.exit(0);
        }

        // initialize dataset
        if (!m_settings.cisAnalysis && !m_settings.transAnalysis) {
            System.err.println("! WARNING: defaulting to CIS analysis (override with --trans or --trans and --cis))");
            m_settings.cisAnalysis = true;
        }

        m_settings.writeSettingsToDisk();

        int numDatasets = m_settings.datasetSettings.size();
        m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
        numAvailableInds = 0;
        int nrOfDatasetsWithGeneExpressionData = 0;
        for (int i = 0; i < numDatasets; i++) {

            System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
            m_settings.datasetSettings.get(i).confineProbesToProbesMappingToAnyChromosome = m_settings.confineProbesToProbesMappingToAnyChromosome;
            System.out.println(ConsoleGUIElems.LINE);
            m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i));

            if (m_gg[i].isExpressionDataLoadedCorrectly()) {
                nrOfDatasetsWithGeneExpressionData++;
            }

        }

        if (nrOfDatasetsWithGeneExpressionData == 0) {
            System.out.println("Error: none of your datasets contain any gene expression data for the settings you have specified");
            System.exit(0);
        }

        if (nrOfDatasetsWithGeneExpressionData != m_gg.length) {
            System.out.println("WARNING: was able to load gene expression data for " + nrOfDatasetsWithGeneExpressionData + " while you specified " + m_gg.length + " datasets in the settings.");

            // remove the datasets without expression data.

            TriTyperGeneticalGenomicsDataset[] tmp_gg = new TriTyperGeneticalGenomicsDataset[nrOfDatasetsWithGeneExpressionData];
            int ctr = 0;
            for (TriTyperGeneticalGenomicsDataset d : m_gg) {
                if (d.isExpressionDataLoadedCorrectly()) {
                    tmp_gg[ctr] = d;
                }
            }
            m_gg = tmp_gg;
        }



        for (int i = 0; i < numDatasets; i++) {
            if (!m_settings.performParametricAnalysis) {
                m_gg[i].getExpressionData().rankAllExpressionData(m_settings.equalRankForTies);
            }
            numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;

        }

        if (m_settings.regressOutEQTLEffectFileName != null && m_settings.regressOutEQTLEffectFileName.trim().length() > 0) {
            if (!Gpio.exists(m_settings.regressOutEQTLEffectFileName)) {
                System.err.println("ERROR: you have specified an eQTL file to regress out, but the file was not found " + m_settings.regressOutEQTLEffectFileName);
                System.exit(0);
            }
            EQTLRegression eqr = new EQTLRegression();
            eqr.regressOutEQTLEffects(m_settings.regressOutEQTLEffectFileName, m_gg);
        }

        System.out.println(ConsoleGUIElems.LINE);
        System.out.println("");

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

    protected void createSNPList() throws IOException {
        ArrayList<String> availableSNPs = new ArrayList<String>();
        HashSet<String> snptree = new HashSet<String>();
        int chrYSNPs = 0;
        int unknownPos = 0;
        int unknownchr = 0;

        TextFile excludedSNPs = new TextFile(m_settings.outputReportsDir + "excludedSNPsBySNPFilter.txt.gz", TextFile.W);

        if ((m_settings.confineSNPsToSNPsPresentInAllDatasets && m_gg.length > 1)) {
            System.out.println("- Confining to SNPs shared by all datasets.");
            String[] snps = m_gg[0].getGenotypeData().getSNPs();

            int absentInSomeDatasets = 0;
            int nonIdenticalMapping = 0;
            for (int s = 0; s < snps.length; s++) {
                String reason = snps[s];
                boolean presentInAllDatasets = true;
                boolean identicalMapping = true;
                boolean snpOk = false;
                String d1SNP = snps[s];
                if ((m_settings.tsSNPsConfine == null || m_settings.tsSNPsConfine.contains(d1SNP)) && (m_settings.confineToSNPsThatMapToChromosome == null || m_gg[0].getGenotypeData().getChr(s) == m_settings.confineToSNPsThatMapToChromosome)) {
                    for (int d = 1; d < m_gg.length; d++) {
                        Integer d2SNPId = m_gg[d].getGenotypeData().getSnpToSNPId().get(d1SNP);
                        if (d2SNPId == null) {
                            presentInAllDatasets = false;
                        } else if (m_gg[d].getGenotypeData().getChr(d2SNPId) != m_gg[0].getGenotypeData().getChr(s) || m_gg[d].getGenotypeData().getChrPos(d2SNPId) != m_gg[0].getGenotypeData().getChrPos(s)) {
                            identicalMapping = false;
                        }
                    }
                    if (!identicalMapping) {
                        nonIdenticalMapping++;
                        reason += "\tSNPs map to different datasets";
                    }

                    if (!presentInAllDatasets) {
                        absentInSomeDatasets++;
                        reason += "\tSNPs are not present in all datasets";
                    }

                    if (presentInAllDatasets && identicalMapping) {
                        if (!snptree.contains(d1SNP)) {
                            availableSNPs.add(d1SNP);
                            snptree.add(d1SNP);
                            snpOk = true;
                        }
                    }
                } else {
                    if (!m_settings.tsSNPsConfine.contains(d1SNP)) {
                        reason += "\tSNP not confined to";
                    } else {
                        reason += "\tSNP chromosome does not match requested chromosome.";
                    }

                }

                if (!snpOk) {
                    excludedSNPs.write(reason + "\n");
                }
            }
            System.out.println("- " + absentInSomeDatasets + " snps not present in all datasets, " + nonIdenticalMapping + " without identical mapping. " + availableSNPs.size() + " included.");
        } else {
            System.out.println("- Not confining to SNPs shared by all datasets.");
            if (m_gg.length == 1) {
                int d = 0;
                String[] snps = m_gg[d].getGenotypeData().getSNPs();
                int snpNum = 0;
                for (String d1SNP : snps) {
                    String reason = d1SNP;
                    Byte chr = m_gg[d].getGenotypeData().getChr(snpNum);
                    Integer chrPos = m_gg[d].getGenotypeData().getChrPos(snpNum);
                    boolean snpOk = false;
                    if (chr == null) {
                        unknownchr++;
                        reason += "\tSNP has unknown chromosome: " + unknownchr;
                    } else if (chr >= 24) {
                        chrYSNPs++;
                        reason += "\tSNP is located on Y, MT, XY chromosome";
                    } else if (chrPos == null || chrPos <= 0) {
                        unknownPos++;
                        reason += "\tSNP has unknown mapping";
                    } else if ((m_settings.tsSNPsConfine == null || m_settings.tsSNPsConfine.contains(d1SNP))
                            && (m_settings.confineToSNPsThatMapToChromosome == null || m_gg[d].getGenotypeData().getChr(snpNum) == m_settings.confineToSNPsThatMapToChromosome)) {
                        if (!snptree.contains(d1SNP)) {
                            availableSNPs.add(d1SNP);
                            snpOk = true;
                        } else {
                            System.err.println("\t!- WARNING! Duplicate SNP detected:\t" + d1SNP);
                            reason += "\t - SNP is duplicate!";
                        }
                    } else {
                        if (!m_settings.tsSNPsConfine.contains(d1SNP)) {
                            reason += "\tSNP not confined to";
                        } else {
                            reason += "\tSNP chromosome does not match requested chromosome: " + m_gg[d].getGenotypeData().getChr(snpNum);
                        }

                    }

                    if (!snpOk) {
                        excludedSNPs.write(reason + "\n");
                    }

                    snptree.add(d1SNP);
                    snpNum++;
                }

            } else {
                for (int d = 0; d < m_gg.length; d++) {
                    String[] snps = m_gg[d].getGenotypeData().getSNPs();
                    int snpNum = 0;
                    for (String d1SNP : snps) {
                        boolean snpOk = false;
                        String reason = d1SNP;
                        if (!snptree.contains(d1SNP)) {
                            int d2 = d + 1;
                            if ((m_settings.tsSNPsConfine == null || m_settings.tsSNPsConfine.contains(d1SNP))
                                    && (m_settings.confineToSNPsThatMapToChromosome == null || m_gg[d].getGenotypeData().getChr(snpNum) == m_settings.confineToSNPsThatMapToChromosome)) {
                                // if this is the last set, and this SNP is unique to this dataset,

                                Byte chr = m_gg[d].getGenotypeData().getChr(snpNum);
                                Integer chrPos = m_gg[d].getGenotypeData().getChrPos(snpNum);

                                if (chr == null) {
                                    unknownchr++;
                                    reason += "\tSNP has unknown chromosome: " + unknownchr;
                                } else if (chr >= 24) {
                                    chrYSNPs++;
                                    reason += "\tSNP is located on Y, MT, XY chromosome";
                                } else if (chrPos == null || chrPos <= 0) {
                                    unknownPos++;
                                    reason += "\tSNP has unknown mapping";
                                } else if (d2 > m_gg.length - 1) { // if this is the last dataset, and the SNP has not been visited, it must be unique to this dataset.
                                    if (!snptree.contains(d1SNP)) {
                                        availableSNPs.add(d1SNP);
                                    }
                                } else { // else, check whether the mappings are identical accross the other datasets.
                                    boolean identicalMapping = true;
                                    for (d2 = d + 1; d2 < m_gg.length; d2++) {
                                        Integer d2SNPId = m_gg[d2].getGenotypeData().getSnpToSNPId().get(d1SNP);
                                        if (d2SNPId != null) {
                                            if (m_gg[d2].getGenotypeData().getChr(d2SNPId) != m_gg[d].getGenotypeData().getChr(snpNum) || m_gg[d2].getGenotypeData().getChrPos(d2SNPId) != m_gg[d].getGenotypeData().getChrPos(snpNum)) {
                                                identicalMapping = false;
                                                reason += "\tSNP has different mappings between datasets";
                                            }
                                        }
                                    }

                                    if (identicalMapping) {
                                        if (!snptree.contains(d1SNP)) {
                                            availableSNPs.add(d1SNP);
                                            snpOk = true;
                                        }
                                    }

                                }
                            } else {
                                if (!m_settings.tsSNPsConfine.contains(d1SNP)) {
                                    reason += "\tSNP not confined to";
                                } else {
                                    reason += "\tSNP chromosome does not match requested chromosome: " + m_gg[d].getGenotypeData().getChr(snpNum);
                                }
                            }

                            snptree.add(d1SNP); // add it to the search list anyway, so we don't query this SNP twice.
                        } else {
                            reason += "\t SNP is duplicate!";

                        }
                        snpNum++;

                        if (!snpOk) {
                            excludedSNPs.write(reason + "\n");
                        }
                    }
                }

            }
            System.out.println("- " + chrYSNPs + " chromosome Y SNPs, " + unknownPos + " SNPS with unknown position, " + unknownchr + " with unknown chromosome.");
            System.out.println("- Remaining SNPs: " + availableSNPs.size());
        }

        int i = 0;
        m_snpList = new String[availableSNPs.size()];
        for (String s : availableSNPs) {
            m_snpList[i] = s;
            i++;
        }



        // create snp translation table..
        m_snpTranslationTable = new Integer[m_gg.length][m_snpList.length];

        for (int p = 0; p < m_snpList.length; p++) {
            String snp = m_snpList[p];
            for (int d = 0; d < m_gg.length; d++) {
                m_snpTranslationTable[d][p] = m_gg[d].getGenotypeData().getSnpToSNPId().get(snp);
            }
        }
        excludedSNPs.close();

        // now determine which of the SNPs that was queried for does not exist in any of the datasets.
        if (m_settings.tsSNPsConfine != null) {
            Iterator<String> it = m_settings.tsSNPsConfine.iterator();
            if (m_settings.tsSNPsConfine.isEmpty()) {
                System.err.println("ERROR: a SNP confinement file is specified in the settings, but it is apparently empty? " + m_settings.strConfineSNP);
            } else {
                String next = it.next();

                TextFile querySNPNotPresent = new TextFile(m_settings.outputReportsDir + "querySNPsNotPresentInDataset.txt.gz", TextFile.W);
                while (it.hasNext()) {
                    boolean isPresentInAnyDataset = false;
                    for (int d = 0; d < m_gg.length; d++) {
                        Integer id = m_gg[d].getGenotypeData().getSnpToSNPId().get(next);
                        if (id != null) {
                            isPresentInAnyDataset = true;
                        }
                    }
                    if (!isPresentInAnyDataset) {
                        querySNPNotPresent.writeln(next);
                    }
                    next = it.next();
                }
                querySNPNotPresent.close();
            }
        }

    }

    protected void createProbeList() throws IOException {
        TextFile probeLog = new TextFile(m_settings.outputReportsDir + "ProbeQCLog.txt.gz", TextFile.W);

        System.out.println("- Determining available probes.");

        System.out.println("\t- Saving logfile to: " + m_settings.outputReportsDir + "ProbeQCLog.txt.gz");
        if (m_settings.confineProbesToProbesPresentInAllDatasets) {
            System.out.println("\t- Confining to probes present in all datasets");
        } else {
            System.out.println("\t- Not confining to probes present in all datasets");
        }

        if (m_settings.confineProbesToProbesMappingToAnyChromosome) {
            System.out.println("\t- Confining to probes that map to any chromosome (including probes without a valid position)");
        } else {
            System.out.println("\t- Confining to probes that map to autosomes, X, Y and MT chromosomes");
        }

        if (m_settings.confineToProbesThatMapToChromosome != null) {
            System.out.println("\t- Confining to probes that map to chromosome: " + m_settings.confineToProbesThatMapToChromosome);
        }


        int numExcluded = 0;
        // first check all datasets for duplicate probes...

        HashSet<String> duplicateProbes = new HashSet<String>();
        HashSet<String> allAvailableProbes = new HashSet<String>();
        for (int d = 0; d < m_gg.length; d++) {
            HashSet<String> availableProbes = new HashSet<String>();
            String[] probes = m_gg[d].getExpressionData().getProbes();
            for (String probe : probes) {
                if (availableProbes.contains(probe)) {
                    duplicateProbes.add(probe);
                    probeLog.writeln("Removing probe:\t" + probe + "\tis a duplicate");
                }
                availableProbes.add(probe);
                allAvailableProbes.add(probe);
            }
        }

        System.out.println("\t- " + allAvailableProbes.size() + " available probes for all datasets. Will now look for duplicate probes.");
        // remove duplicates from the analysis...

        String[] duplicates = duplicateProbes.toArray(new String[duplicateProbes.size()]);
        for (String duplicate : duplicates) {
            allAvailableProbes.remove(duplicate);
        }
        System.out.println("\t- " + allAvailableProbes.size() + " available probes for all datasets after removing " + duplicateProbes.size() + " duplicates");

        // determine in how many datasets the probes are present
        String[] availableProbeArray = allAvailableProbes.toArray(new String[0]);



        if (m_settings.tsProbesConfine == null) {
        } else {
            System.out.println("Probe confinement list has " + m_settings.tsProbesConfine.size() + " probes");
            availableProbeArray = m_settings.tsProbesConfine.toArray(new String[0]);
        }

        int mappingToDifferentPositionsAcrossDatasets = 0;
        int mapToWrongChromosome = 0;
        int invalidMappingPosition = 0;
        int nrProbesNotInAllDatasets = 0;
        HashSet<String> finalProbeList = new HashSet<String>();
        for (String probe : availableProbeArray) {
            int presence = 0;

            Byte chr = null;
            Integer chrPosStart = null;
            Integer chrPosEnd = null;

            boolean hasIdenticalMappingAcrossDatasets = true;
            String mappingOutput = "";
            for (int d = 0; d < m_gg.length; d++) {
                Integer probeId = m_gg[d].getExpressionData().getProbeToId().get(probe);
                if (probeId != null) {
                    presence++;
                    if (chr == null) {
                        chr = m_gg[d].getExpressionData().getChr()[probeId];
                        if (chr == -1) {
                            chr = null;
                        }
                        chrPosStart = m_gg[d].getExpressionData().getChrStart()[probeId];
                        chrPosEnd = m_gg[d].getExpressionData().getChrStop()[probeId];

                        mappingOutput += "\t" + m_gg[d].getSettings().name + ": Chr " + chr + "; Pos " + chrPosStart + "-" + chrPosEnd;

                    } else {

                        Byte chr2 = m_gg[d].getExpressionData().getChr()[probeId];
                        if (chr2 == -1) {
                            chr2 = null;
                        }
                        Integer chrPosStart2 = m_gg[d].getExpressionData().getChrStart()[probeId];
                        Integer chrPosEnd2 = m_gg[d].getExpressionData().getChrStop()[probeId];

                        if (chr2 == null && chr != null) {
                            hasIdenticalMappingAcrossDatasets = false;
                        } else if (chrPosStart2 == null && chrPosStart != null) {
                            hasIdenticalMappingAcrossDatasets = false;
                        } else if (chrPosEnd2 == null && chrPosEnd != null) {
                            hasIdenticalMappingAcrossDatasets = false;
                        } else if (!chr.equals(chr2)) {
                            hasIdenticalMappingAcrossDatasets = false;
                        } else if (!chrPosStart.equals(chrPosStart2)) {
                            hasIdenticalMappingAcrossDatasets = false;
                        } else if (!chrPosEnd.equals(chrPosEnd2)) {
                            hasIdenticalMappingAcrossDatasets = false;
                        }
                        mappingOutput += "\t" + m_gg[d].getSettings().name + ": Chr " + chr2 + "; Pos " + chrPosStart2 + "-" + chrPosEnd2;
                    }
                }
            }

            // check whether the probe has a proper mapping
            boolean includeProbe = true;
            // check whether we would want to include this probe on the basis of its presence in other datasets
//	    if (m_settings.cisAnalysis && m_settings.transAnalysis && m_settings.confineProbesToProbesPresentInAllDatasets) { // ?

            if (m_settings.tsProbesConfine != null && m_settings.tsProbesConfine.contains(probe)) {
                includeProbe = true;
            } else if (m_settings.cisAnalysis && m_settings.transAnalysis) { // ?
                includeProbe = true;
            } else if (m_settings.confineProbesToProbesPresentInAllDatasets && presence < m_gg.length) {
                includeProbe = false;
                nrProbesNotInAllDatasets++;
                probeLog.writeln("Removing probe:\t" + probe + "\tpresent in " + presence + " / " + m_gg.length + "\tdatasets");
            } else if (chr == null || chr >= 25 || chr <= 0) {
                if (!m_settings.confineProbesToProbesMappingToAnyChromosome) {
                    includeProbe = false;
                    invalidMappingPosition++;
                    probeLog.writeln("Removing probe:\t" + probe + "\t has no valid mapping position in any dataset: " + mappingOutput);
                }
            } else if (m_settings.confineToProbesThatMapToChromosome != null && chr != m_settings.confineToProbesThatMapToChromosome) {
                // check whether this chromosome was requested to be analysed
                includeProbe = false;
                mapToWrongChromosome++;
                probeLog.writeln("Removing probe:\t" + probe + "\tmaps to wrong chromosome: " + mappingOutput);
            } else if (!hasIdenticalMappingAcrossDatasets) {
                // exclude the probe if it does not have identical mappings across datasets
                includeProbe = false;
                mappingToDifferentPositionsAcrossDatasets++;
                probeLog.writeln("Removing probe:\t" + probe + "\tmaps to different positions in datasets: " + mappingOutput);
            }

            if (includeProbe) {
                finalProbeList.add(probe);
            }
        }

        System.out.println("\t- " + finalProbeList.size() + "\tprobes finally included: ");
        if (m_settings.confineProbesToProbesPresentInAllDatasets) {
            System.out.println("\t\t- " + nrProbesNotInAllDatasets + " are not present in all datasets");
        }
        if (!m_settings.confineProbesToProbesMappingToAnyChromosome) {
            System.out.println("\t\t- " + invalidMappingPosition + " have an invalid mapping position");
        }
        System.out.println("\t\t- " + mappingToDifferentPositionsAcrossDatasets + " probes map to different positions in one or more datasets");
        if (m_settings.confineToProbesThatMapToChromosome != null) {
            System.out.println("\t\t- " + mapToWrongChromosome + " probes map to a different chromosome than the one selected (Chr: " + m_settings.confineToProbesThatMapToChromosome + ")");
        }
        if (m_settings.confineProbesToProbesPresentInAllDatasets) {
        }

        if (finalProbeList.isEmpty()) {
            System.err.println("Error: no probes remaining after filter. Are your settings correct?");
            probeLog.close();
            System.exit(0);
        }
        m_probeList = finalProbeList.toArray(new String[finalProbeList.size()]);

        // create probe translation table..
        m_probeTranslationTable = new Integer[m_gg.length][m_probeList.length];

        for (int p = 0; p < m_probeList.length; p++) {
            String probe = m_probeList[p];
            for (int d = 0; d < m_gg.length; d++) {
                m_probeTranslationTable[d][p] = m_gg[d].getExpressionData().getProbeToId().get(probe);
            }
        }

        probeLog.close();
    }

    public void mapEQTLs() throws IOException {

        // create work packages
        RunTimer t = new RunTimer();

        CalculationThread[] pool = new CalculationThread[m_settings.nrThreads];

        SNPLoader[] snploaders = new SNPLoader[m_gg.length];
        TriTyperExpressionData[] expressiondata = new TriTyperExpressionData[m_gg.length];
        for (int d = 0; d < snploaders.length; d++) {
            snploaders[d] = m_gg[d].getGenotypeData().createSNPLoader();
            expressiondata[d] = m_gg[d].getExpressionData();
        }

        // initialize lookup tables
        int maxNrSamples = 0;
        for (int d = 0; d < m_gg.length; d++) {
            if (m_gg[d].getExpressionToGenotypeIdArray().length > maxNrSamples) {
                maxNrSamples = m_gg[d].getExpressionToGenotypeIdArray().length;
            }
        }
        Correlation.correlationToZScore(maxNrSamples);
        Descriptives.lookupSqrt(numAvailableInds);            // pre-calculate a square root lookup table
        Descriptives.zScoreToPValue();

        boolean permuting = false;

        System.gc();
        System.gc();
        System.gc();

        System.out.println("Will write output to dir: " + m_settings.outputReportsDir);

        int permStart = 0;
        int permEnd = m_settings.nrPermutationsFDR + 1;
        if (m_settings.startWithPermutation != null) {
            permStart = m_settings.startWithPermutation;
            permEnd = permStart + m_settings.nrPermutationsFDR + 1;
        }

        boolean hasResults = true;

        for (int permutationRound = permStart; permutationRound < permEnd; permutationRound++) {
            RunTimer permtime = new RunTimer();

            if (permutationRound > 0) {
                System.out.print("Permuting data, round: " + permutationRound + " of " + m_settings.nrPermutationsFDR + "\n" + ConsoleGUIElems.LINE);

                for (int d = 0; d < m_gg.length; d++) {
                    int[] indWGAOriginal = m_gg[d].getExpressionToGenotypeIdArray();
                    m_gg[d].permuteSampleLables();

                    int[] indWGAPerm = m_gg[d].getExpressionToGenotypeIdArray();
                    int identical = 0;
                    for (int i = 0; i < indWGAPerm.length; i++) {
                        if (indWGAOriginal[i] == indWGAPerm[i]) {
                            identical++;
                        }
                    }

//                    System.out.println("After permuting: " + identical + " unchanged, " + indWGAOriginal.length + " total");
                }
                permuting = true;
            } else {
                System.out.print("Running real eQTL analysis\n" + ConsoleGUIElems.LINE);
            }

            int[][] expressionToGenotypeIds = new int[m_gg.length][0];
            for (int d = 0; d < m_gg.length; d++) {
                expressionToGenotypeIds[d] = m_gg[d].getExpressionToGenotypeIdArray();
            }

            LinkedBlockingQueue<WorkPackage> resultQueue = new LinkedBlockingQueue<WorkPackage>(100);
            ResultProcessorThread resultthread = new ResultProcessorThread(m_settings.nrThreads, resultQueue, m_settings.createBinaryOutputFiles,
                    m_gg, m_settings, m_probeTranslationTable, permuting, permutationRound, m_snpList, m_probeList, m_workPackages);
            resultthread.setName("ResultProcessorThread");
            resultthread.start();

            // start production in advance
            LinkedBlockingQueue<WorkPackage> packageQueue = new LinkedBlockingQueue<WorkPackage>(100000);
            WorkPackageProducer producer = new WorkPackageProducer(packageQueue, m_workPackages, m_snpList, m_probeList, m_probeTranslationTable, m_snpTranslationTable, m_gg, snploaders, m_settings, permuting);
            producer.setName("WorkPackageProducerThread");
            producer.start();

            // run calculations
            for (int tnum = 0; tnum < pool.length; tnum++) {
                EQTLPlotter plotter = null;
                if (!permuting) {
                    plotter = new EQTLPlotter(m_gg, m_settings, m_probeList, m_probeTranslationTable);
                }

                pool[tnum] = new CalculationThread(permutationRound, packageQueue, resultQueue, expressiondata, m_probeTranslationTable,
                        expressionToGenotypeIds, m_settings, plotter, m_settings.createBinaryOutputFiles);
                pool[tnum].setName("CalcThread-" + tnum);
                pool[tnum].start();

            }

            // kill the threads
            try {

                producer.join();



//                System.out.println("Joining calculation threads");
                for (int threadNum = 0; threadNum < m_settings.nrThreads; threadNum++) {
                    pool[threadNum].join();
                }

                WorkPackage poison = new WorkPackage();
                poison.results = new Result(true);
                resultQueue.put(poison);
                resultthread.join();

            } catch (InterruptedException e) {
                System.err.println("Exception: Main Thread interrupted.");
            }
            System.out.print(ConsoleGUIElems.LINE);
            System.out.println("Round done. Elapsed time:\t" + permtime.getTimeDesc());
            System.out.println("");
            resultQueue.clear();
            packageQueue.clear();

            resultQueue = null;
            packageQueue = null;
            permtime = null;
            producer = null;
            resultthread = null;
            expressionToGenotypeIds = null;
            for (int i = 0; i < pool.length; i++) {
                pool[i] = null;
            }

            // check whether there were results..
            String fileName;
            if (permutationRound > 0) {
                fileName = m_settings.outputReportsDir + "PermutedEQTLsPermutationRound" + permutationRound + ".txt.gz";
            } else {
                fileName = m_settings.outputReportsDir + "eQTLs.txt.gz";
            }
            TextFile tf = new TextFile(fileName, TextFile.R);
            tf.readLine(); // skip header
            int lnCounter = 0;
            String line = tf.readLine();
            while (line != null) {
                lnCounter++;
                if (lnCounter > 1) {
                    break;
                }
                line = tf.readLine();
            }
            tf.close();
            if (lnCounter == 0) {
                System.err.println("WARNING: QTL Mapping did not yield any results.");
                hasResults = false;
            }



        }

        for (int d = 0; d < snploaders.length; d++) {
            snploaders[d].close();
        }

        if (!m_settings.runOnlyPermutations && hasResults) {
            if (m_settings.createTEXTOutputFiles && m_settings.nrPermutationsFDR > 0) {
                System.out.println("Calculating FDR:\n" + ConsoleGUIElems.LINE);
                FDR.calculateFDR(m_settings.outputReportsDir, m_settings.nrPermutationsFDR, m_settings.maxNrMostSignificantEQTLs, m_settings.fdrCutOff, m_settings.createQQPlot, null, null);

                if (m_settings.createDotPlot) {
                    EQTLDotPlot edp = new EQTLDotPlot();
                    try {
                        edp.draw(m_settings.outputReportsDir + "/eQTLsFDR" + m_settings.fdrCutOff + ".txt", m_settings.outputReportsDir + "/DotPlot-FDR" + m_settings.fdrCutOff + ".pdf", EQTLDotPlot.Output.PDF); // "/eQTLsFDR" + fdrCutOff + ".txt", outputReportsDir + "/eQTLsFDR" + fdrCutOff + "DotPlot.png"
                    } catch (DocumentException ex) {
                        Logger.getLogger(MetaQTL3.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    edp = null;
                }

            }
        }







        System.out.print(ConsoleGUIElems.DOUBLELINE);

        System.out.println("eQTL mapping elapsed:\t" + t.getTimeDesc() + "\n");
    }

    protected void determineSNPProbeCombinations() throws IOException {

        String loc = m_settings.outputReportsDir + "excludedSNPsBySNPProbeCombinationFilter.txt.gz";
        TextFile excludedSNPs = new TextFile(loc, TextFile.W);
        long maxNrTestsToPerform = 0;
        int[] midpoint = new int[m_probeList.length];
        byte[] chr = new byte[m_probeList.length];
        HashMap<Byte, ArrayList<Integer>> chrToProbe = new HashMap<Byte, ArrayList<Integer>>();
        System.out.println("- Calculating probe midpoint positions");
        HashSet<String> visitedProbes = new HashSet<String>();

        for (int p = 0; p < m_probeList.length; p++) {
            for (int d = 0; d < m_gg.length; d++) {
                if (m_probeTranslationTable[d][p] != null && !visitedProbes.contains(m_probeList[p])) {
                    int pid = m_probeTranslationTable[d][p];
                    int start = m_gg[d].getExpressionData().getChrStart()[pid];
                    int stop = m_gg[d].getExpressionData().getChrStop()[pid];
                    midpoint[p] = (int) Math.floor((double) (stop + start) / 2);
                    chr[p] = m_gg[d].getExpressionData().getChr()[pid];

                    ArrayList<Integer> probes = chrToProbe.get(chr[p]);
                    if (probes == null) {
                        probes = new ArrayList<Integer>();
                    }
                    probes.add(p);

                    chrToProbe.put(chr[p], probes);
                    visitedProbes.add(m_probeList[p]);
                }
            }
        }


        WorkPackage[] workPackages = new WorkPackage[m_snpList.length];
        // improve performance by sorting here, and then breaking later..
        int numWorkPackages = 0;

        System.out.println("- Determining SNP-Probe combinations to test");

        boolean cisOnly = false;
        boolean cisTrans = false;
        boolean transOnly = false;

        if (m_settings.cisAnalysis && !m_settings.transAnalysis) {
            cisOnly = true;
        } else if (!m_settings.cisAnalysis && m_settings.transAnalysis) {
            transOnly = true;
        } else if (m_settings.cisAnalysis && m_settings.transAnalysis) {
            cisTrans = true;
        }

        HashMap<String, Integer> probeNameToId = null;
        if (m_settings.tsSNPProbeCombinationsConfine != null) {
            m_settings.transAnalysis = false;
            m_settings.cisAnalysis = true;
            for (int i = 0; i < m_gg.length; i++) {
                m_gg[i].getSettings().cisAnalysis = true;
                m_gg[i].getSettings().transAnalysis = false;
            }
            cisTrans = false;
            cisOnly = true;
            transOnly = false;

            probeNameToId = new HashMap<String, Integer>();
            for (int i = 0; i < m_probeList.length; i++) {
                probeNameToId.put(m_probeList[i], i);
            }
        }



        int prevProc = 0;
        ProgressBar pb = new ProgressBar(m_snpList.length);
        for (int s = 0; s < m_snpList.length; s++) {
            WorkPackage output = null;

            SNP[] snps = new SNP[m_gg.length];

            byte snpchr = -1;
            int snppos = -1;
            int dreq = 0;

            String snpname = "";

            // load the genomic positions for this snp
            for (int d = 0; d < m_gg.length; d++) {
                Integer snpId = m_snpTranslationTable[d][s];
                if (snpId != null) {
                    snps[d] = m_gg[d].getGenotypeData().getSNPObject(snpId);
                    snpchr = snps[d].getChr();
                    snppos = snps[d].getChrPos();
                    snpname = snps[d].getName();
                    dreq = d;
                }
            }


            ArrayList<Integer> probeOnChr = chrToProbe.get(snpchr);

            // cis trans
            if (cisTrans) {
                output = new WorkPackage();
                output.setMetaSNPId(s);
                output.setSnps(snps);
                output.setProbes(null);
                numWorkPackages++;
                maxNrTestsToPerform += m_probeList.length;
                workPackages[s] = output;


                // cis or trans
            } else {
                ArrayList<Integer> probeToTest = null;
                if (m_settings.tsSNPProbeCombinationsConfine != null) {
                    HashSet<String> probesSelected = m_settings.tsSNPProbeCombinationsConfine.get(snpname);
                    if (probesSelected != null) {
                        probeToTest = new ArrayList<Integer>();
                        for (String probe : probesSelected) {
                            Integer probeId = probeNameToId.get(probe);
                            if (probeId == null) {
                                System.err.println("You selected the following SNP-Probe combination, but probe not present in dataset!?\t" + snpname + "\t-\t" + probe);
                            } else {
                                probeToTest.add(probeId);
                            }
                        }
                    }
                } else {
                    if (probeOnChr != null && !probeOnChr.isEmpty()) {
                        probeToTest = new ArrayList<Integer>();
                        for (int e = 0; e < probeOnChr.size(); e++) {
                            int p = probeOnChr.get(e);

                            if (Math.abs(midpoint[p] - snppos) < m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance) {
                                probeToTest.add(p); // depending if the test is cis or trans, we test these probes, or not.
                                // System.out.println(snps[dreq].getName()+"\t"+p+"\t"+probeList[p]+"\t"+Math.abs(midpoint[p] - snppos));
                            }
                        }
                    }
                }

                // don't add the cis-workpackage when there are no probes to test.
                if (cisOnly && (probeToTest == null || probeToTest.isEmpty())) {
                    workPackages[s] = null;
                    excludedSNPs.write(snpname + "\tNo probes within " + m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance + "bp\n");
                } else {
                    int[] testprobes = null;
                    if (probeToTest != null) {
                        testprobes = new int[probeToTest.size()];
                        for (int p = 0; p < testprobes.length; p++) {
                            testprobes[p] = probeToTest.get(p);
                        }
                    }

                    output = new WorkPackage();
                    output.setMetaSNPId(s);
                    output.setSnps(snps);
                    output.setProbes(testprobes);

                    workPackages[s] = output;

                    numWorkPackages++;
                    if (cisOnly) {
                        maxNrTestsToPerform += testprobes.length;
                    } else {
                        if (testprobes != null) {
                            maxNrTestsToPerform += (m_probeList.length - testprobes.length);
                        } else {
                            maxNrTestsToPerform += (m_probeList.length);
                        }

                    }
                }
            }

            pb.iterate();
//		if (numWorkPackages % 100000 == 0 && numWorkPackages > 0 && numWorkPackages > prevProc) {
//		    System.out.println("\t" + numWorkPackages + " SNPs Processed.");
//		    prevProc = numWorkPackages;
//		}
        }

        pb.close();
        System.out.println("");
        if (numWorkPackages != m_snpList.length) {
            m_workPackages = new WorkPackage[numWorkPackages];
            int q = 0;
            for (int i = 0; i < workPackages.length; i++) {
                if (workPackages[i] != null) {
                    m_workPackages[q] = workPackages[i];
                    m_workPackages[q].setId(q);
                    q++;
                }
            }
        } else {
            for (int i = 0; i < workPackages.length; i++) {
                workPackages[i].setId(i);
            }
            this.m_workPackages = workPackages;
        }

        excludedSNPs.close();


        System.out.println("The maximum number of SNPs to test: " + m_workPackages.length);
        System.out.println("The maximum number of SNP-Probe combinations: " + maxNrTestsToPerform);
    }

    protected void printSummary() {

        System.out.print("\nSummary\n" + ConsoleGUIElems.DOUBLELINE);
        int totalSamples = 0;
        for (int d = 0; d < m_gg.length; d++) {
            System.out.print("Dataset:\t" + m_gg[d].getSettings().name);
            System.out.print("\tprobes:\t" + m_gg[d].getExpressionData().getProbes().length);
            System.out.print("\tSNPs:\t" + m_gg[d].getGenotypeData().getSNPs().length);
            totalSamples += m_gg[d].getTotalGGSamples();
            System.out.println("\tsamples:\t" + m_gg[d].getTotalGGSamples());
        }
        System.out.println("");
        System.out.print("\nTotals\n" + ConsoleGUIElems.DOUBLELINE);
        System.out.println("Total number of datasets:\t" + m_gg.length);
        System.out.println("Total number of samples:\t" + totalSamples);
        System.out.println("Maximum number of SNPs to test:\t" + m_workPackages.length);
        System.out.println("Maximum number of Probes to test:\t" + m_probeList.length);

        if (totalSamples == 0) {
            System.err.println("ERROR!: No samples detected");
            System.exit(-1);
        }
        System.out.println("");
        System.out.print("\nAnalysis\n" + ConsoleGUIElems.DOUBLELINE);

        if (m_settings.cisAnalysis && m_settings.transAnalysis) {
            System.out.println("- cis/trans analysis");
        } else if (m_settings.cisAnalysis && !m_settings.transAnalysis) {
            System.out.println("- cis analysis");
        } else if (!m_settings.cisAnalysis && m_settings.transAnalysis) {
            System.out.println("- trans analysis");
        }
        if (!m_settings.performParametricAnalysis) {
            System.out.println("- non-parametric (Spearman ranked) correlation");
        } else {
            System.out.println("- parametric (Pearson) correlation");
        }
        System.out.println("- Mid-point distance:\t" + m_settings.ciseQTLAnalysMaxSNPProbeMidPointDistance);
        System.out.println("- FDR cutoff:\t" + m_settings.fdrCutOff);
        System.out.println("- Nr. permutations:\t" + m_settings.nrPermutationsFDR);
        System.out.println("- Nr. Threads:\t" + m_settings.nrThreads);
        System.out.println("- Max nr results:\t" + m_settings.maxNrMostSignificantEQTLs);

        if (m_settings.createBinaryOutputFiles) {
            System.out.println("- creating BINARY output");
        }

        if (m_settings.createTEXTOutputFiles) {
            System.out.println("- creating TEXT output");
        }

        System.out.println("");
    }
}
