 /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class MetaQTL4Settings {

    // Output
    private final boolean createSNPPValueSummaryStatisticsFile;               // Output SNP P-Value summary statistics
    private final boolean createSNPSummaryStatisticsFile;                     // Output SNP P-Value summary statistics
    private final boolean createEQTLPValueTable;                              // Output an eQTL p-value table (only applies for a limited number (500) of SNP)
    private final File outputReportsDir;                                            // Output directory for reports
    // SNP QC
    // Analysis settings
    private final boolean performParametricAnalysis;                          // Perform parametric analysis
    private final boolean useAbsoluteZScorePValue;                            // Use absolute Z-score? (required for finding opposite allelic effects)
    private final int ciseQTLAnalysMaxSNPProbeMidPointDistance;                       // Midpoint distance for declaring an eQTL effect CIS
    private final int maxNrMostSignificantEQTLs;                             // Max number of results stored in memory
    private final boolean performParametricAnalysisGetAccuratePValueEstimates;         // Use an accurate estimation of the P-values
    private final Integer nrThreads;                                                      // Use this number of threads
    // Multiple testing correction
    private final double fdrCutOff;                                                   // Cutoff for FDR procedure
    private final int nrPermutationsFDR;                                              // Number of permutations to determine FDR
    // confinements
    private final boolean performEQTLAnalysisOnSNPProbeCombinationSubset;             // Confine to a certain set of probe/snp combinations?
    private final Byte confineToSNPsThatMapToChromosome;                               // Confine SNP to be assessed to SNPs mapped on this chromosome
    private final boolean expressionDataLoadOnlyProbesThatMapToChromosome;    // Only load expression data for probes with a known chromosome mapping
    // plots
    private final double plotOutputPValueCutOff;                                      // Use this p-value as a cutoff for drawing plots
    private final String plotOutputDirectory;                                         // Print the plots in this directory
    private final boolean runOnlyPermutations;
    private final Integer startWithPermutation;
    private final Boolean confineSNPsToSNPsPresentInAllDatasets;
    private final boolean confineProbesToProbesPresentInAllDatasets;
    private final ArrayList<MetaQTL4DatasetSettings> datasetSettings;
    private final String regressOutEQTLEffectFileName;
    private final Double snpQCCallRateThreshold;
    private final Double snpQCHWEThreshold;
    private final Double snpQCMAFThreshold;
    private final Byte confineToProbesThatMapToChromosome;
    private final boolean createBinaryOutputFiles;
    private final boolean createTEXTOutputFiles;
    private final String strConfineSNP;
    private final String strConfineProbe;
    private final boolean provideFoldChangeData;
    private final boolean provideBetasAndStandardErrors;
    private final boolean equalRankForTies;
    private final boolean createQQPlot;
    private final boolean createDotPlot;
    private final String strSNPProbeConfine;
    private final String settingsTextToReplace;
    private final String settingsTextReplaceWith;
    private final boolean cisAnalysis;
    private final boolean transAnalysis;
    private final String probeAnnotationFile;

    public MetaQTL4Settings(String settings, String settingsTextToReplace, String settingsTextReplaceWith) throws IOException, Exception {

        this.settingsTextToReplace = settingsTextToReplace;
        this.settingsTextReplaceWith = settingsTextReplaceWith;
        XMLConfiguration config = null;
        try {
            config = new XMLConfiguration(settings);           // Use the apache XML configuration parser
        } catch (ConfigurationException e) {
            throw new Exception("Configuration not properly formatted: " + e.getMessage());
        }

        String analysisType = null;
        String correlationType = null;
        Boolean useAbsPVal = null;
        Integer nrthread = null;
        Integer cisDist = null;
        Double MAF = null;
        Double HWE = null;
        Double callrate = null;
        String correctiontype = null;
        Double mtThreshold = null;
        Integer numPermutations = null;
        String outdir = null;
        Double outputplotthreshold = null;
        String outputplotdirectory = null;

        // load the defaults:

        performEQTLAnalysisOnSNPProbeCombinationSubset = false;

        // QC defaults
        try {
            callrate = config.getDouble("defaults.qc.snpqccallratethreshold");
        } catch (Exception e) {
        }
        if (callrate != null) {
            snpQCCallRateThreshold = callrate;
        } else {
            snpQCCallRateThreshold = 0.95;
        }

        try {
            HWE = config.getDouble("defaults.qc.snpqchwethreshold");
        } catch (Exception e) {
        }

        if (HWE != null) {
            snpQCHWEThreshold = HWE;
        } else {
            snpQCHWEThreshold = 0.0001;
        }

        String probeAnnot = null;
        try{
            probeAnnot = config.getString("defaults.analysis.probeAnnotationFile");
        } catch (Exception e){
            
        }
        
        probeAnnotationFile = probeAnnot;
        
        try {
            MAF = config.getDouble("defaults.qc.snpqcmafthreshold");
        } catch (Exception e) {
        }
        if (MAF != null) {
            snpQCMAFThreshold = MAF;
        } else {
            snpQCMAFThreshold = 0.05;
        }

        // analysis settings
        try {
            analysisType = config.getString("defaults.analysis.analysistype");
        } catch (Exception e) {
        }


        createQQPlot = config.getBoolean("defaults.analysis.createqqplot", false);
        createDotPlot = config.getBoolean("defaults.analysis.createdotplot", false);
        runOnlyPermutations = config.getBoolean("defaults.analysis.onlypermutations", false);

        String strStartWithPermutation = null;
        try {
            strStartWithPermutation = config.getString("defaults.analysis.startwithpermutation", null);

            if (settingsTextToReplace != null && strStartWithPermutation.contains(settingsTextToReplace)) {
                strStartWithPermutation = strStartWithPermutation.replace(settingsTextToReplace, settingsTextReplaceWith);
            }


        } catch (Exception e) {
        }

        if (strStartWithPermutation != null) {
            startWithPermutation = Integer.parseInt(strStartWithPermutation);
        } else {
            startWithPermutation = null;
        }


        if (analysisType != null) {
            if (analysisType.toLowerCase().equals("cis")) {
                cisAnalysis = true;
                transAnalysis = false;
            } else if (analysisType.toLowerCase().equals("trans")) {
                cisAnalysis = false;
                transAnalysis = true;
            } else if (analysisType.toLowerCase().equals("cistrans")) {
                cisAnalysis = true;
                transAnalysis = true;
            } else {
                cisAnalysis = true;
                transAnalysis = false;
            }
        } else {
            cisAnalysis = true;
            transAnalysis = false;
        }

        try {
            cisDist = config.getInteger("defaults.analysis.cisanalysisprobedistance", null);
        } catch (Exception e) {
        }
        if (cisDist != null) {
            ciseQTLAnalysMaxSNPProbeMidPointDistance = cisDist;
        } else {
            ciseQTLAnalysMaxSNPProbeMidPointDistance = 250000;
        }

        try {
            correlationType = config.getString("defaults.analysis.correlationtype");
        } catch (Exception e) {
        }
        if (correlationType != null) {
            if (correlationType.toLowerCase().equals("parametric")) {
                performParametricAnalysis = true;
            } else {
                performParametricAnalysis = false;
            }
        } else {
            performParametricAnalysis = true;
        }


        Boolean useIdenticalRanksForTies = null;
        try {
            useIdenticalRanksForTies = config.getBoolean("defaults.analysis.equalrankforties");
        } catch (Exception e) {
        }
        if (useIdenticalRanksForTies != null) {
            equalRankForTies = useIdenticalRanksForTies;
        } else {
            equalRankForTies = false;
        }

        try {
            useAbsPVal = config.getBoolean("defaults.analysis.useabsolutepvalue");
        } catch (Exception e) {
        }

        if (useAbsPVal == null) {
            useAbsoluteZScorePValue = false;
        } else {
            useAbsoluteZScorePValue = useAbsPVal;
        }



        try {
            nrthread = config.getInteger("defaults.analysis.threads", null);
        } catch (Exception e) {
        }

        if (nrthread == null) {
            nrThreads = (Runtime.getRuntime().availableProcessors() - 1);
        } else {
            int numProcs = Runtime.getRuntime().availableProcessors();
            if (nrthread > numProcs || nrthread < 1) {
                nrThreads = numProcs;
            } else {
                nrThreads = nrthread;
            }
        }

        // multiple testing
        try {
            correctiontype = config.getString("defaults.multipletesting.type");
        } catch (Exception e) {
        }
        if (correctiontype != null) {
        } else {
        }

        try {
            mtThreshold = config.getDouble("defaults.multipletesting.threshold", null);
        } catch (Exception e) {
        }

        if (mtThreshold != null) {
            fdrCutOff = mtThreshold;
        } else {
            fdrCutOff = 0.05;
        }

        try {
            numPermutations = config.getInteger("defaults.multipletesting.permutations", null);
        } catch (Exception e) {
        }
        if (numPermutations != null) {
            nrPermutationsFDR = numPermutations;
        } else {
            nrPermutationsFDR = 0;
        }

        // output settings
        try {
            outdir = config.getString("defaults.output.outputdirectory");
            if (settingsTextToReplace != null && outdir.contains(settingsTextToReplace)) {
                outdir = outdir.replace(settingsTextToReplace, settingsTextReplaceWith);
            }
        } catch (Exception e) {
        }

        if (outdir != null) {
            outputReportsDir = new File(outdir);
            // check if dir exists. if it does not, create it:
            if (!Gpio.exists(outdir)) {
                Gpio.createOuputDir(outputReportsDir);
            }
            try {
				config.save(new File(outputReportsDir, "metaqtlsettings.xml"));           
            } catch (ConfigurationException e) {
                e.printStackTrace();
            }
        } else {
            outputReportsDir = null;
            System.out.println("Error: please supply an output directory.");
            System.exit(-1);
        }

//        try {
//            String probeConversionFileLoc = config.getString("defaults.analysis.probeconversion", null);
//            if (probeConversionFileLoc != null && probeConversionFileLoc.length() > 0) {
//                System.out.println(probeConversionFileLoc);
//            }
//
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(-1);
//        }


        try {
            outputplotthreshold = config.getDouble("defaults.output.outputplotthreshold");
        } catch (Exception e) {
        }


        createBinaryOutputFiles = config.getBoolean("defaults.output.binaryoutput", false);
        createTEXTOutputFiles = config.getBoolean("defaults.output.textoutput", true);

        if (outputplotthreshold != null) {
            plotOutputPValueCutOff = outputplotthreshold;
        } else {
            plotOutputPValueCutOff = Double.MAX_VALUE;
        }

        try {
            outputplotdirectory = config.getString("defaults.output.outputplotdirectory");
            if (settingsTextToReplace != null && outputplotdirectory.contains(settingsTextToReplace)) {
                outputplotdirectory = outputplotdirectory.replace(settingsTextToReplace, settingsTextReplaceWith);
            }

        } catch (Exception e) {
        }

        if (outputplotdirectory != null) {
            plotOutputDirectory = Gpio.formatAsDirectory(outputplotdirectory);
        } else {
            plotOutputDirectory = outdir + "/plots/";
        }

        // check if dir exists. if it does not, create it if plots are requested
        if (plotOutputPValueCutOff != Double.MAX_VALUE) {
            if (!Gpio.exists(plotOutputDirectory)) {
                Gpio.createDir(plotOutputDirectory);
            }
        }

        createSNPPValueSummaryStatisticsFile = config.getBoolean("defaults.output.generatesnppvaluesummarystatistics", false);
        provideFoldChangeData = config.getBoolean("defaults.output.generatefoldchangevalues", false);
        provideBetasAndStandardErrors = config.getBoolean("defaults.output.generatebetaandfoldchanges", false);
        createSNPSummaryStatisticsFile = config.getBoolean("defaults.output.generatesnpsummarystatistics", false);
        createEQTLPValueTable = config.getBoolean("defaults.output.generateeqtlpvaluetable", false);
        maxNrMostSignificantEQTLs = config.getInt("defaults.output.maxnreqtlresults", 150000);
        expressionDataLoadOnlyProbesThatMapToChromosome = config.getBoolean("defaults.confine.confineProbesThatMapToKnownChromosome", true);
        Boolean confineTpSNPsPresentInAllDs = config.getBoolean("defaults.confine.confineSNPsToSNPsPresentInAllDatasets", false);
        confineSNPsToSNPsPresentInAllDatasets = confineTpSNPsPresentInAllDs;
        confineToProbesThatMapToChromosome = null;

        // confinements on snp, probe, or snp-probe
        String confineSNP = null;
        String confineProbe = null;
        String snpProbeConfine = null;

        // confine to this list of snp
        try {
            confineSNP = config.getString("defaults.confine.snp");
            if (settingsTextToReplace != null && confineSNP.contains(settingsTextToReplace)) {
                confineSNP = confineSNP.replace(settingsTextToReplace, settingsTextReplaceWith);
            }
        } catch (Exception e) {
        }

        // confine to this list of probes
        try {
            confineProbe = config.getString("defaults.confine.probe");
            if (settingsTextToReplace != null && confineProbe.contains(settingsTextToReplace)) {
                confineProbe = confineProbe.replace(settingsTextToReplace, settingsTextReplaceWith);
            }
        } catch (Exception e) {
        }

        // confine to this list of snp-probe combinations
        try {
            snpProbeConfine = config.getString("defaults.confine.snpProbe");
            if (settingsTextToReplace != null && snpProbeConfine.contains(settingsTextToReplace)) {
                snpProbeConfine = snpProbeConfine.replace(settingsTextToReplace, settingsTextReplaceWith);
            }
        } catch (Exception e) {
        }

        strConfineSNP = confineSNP;
        strConfineProbe = confineProbe;
        strSNPProbeConfine = snpProbeConfine;

        // confine to snp present in all datasets
        performParametricAnalysisGetAccuratePValueEstimates = false;


        // confine to SNP that map to this chromosome
        String confineStr = null;
        try {
            confineStr = config.getString("defaults.confine.confineToSNPsThatMapToChromosome", null);
        } catch (Exception e) {
        }

        if (confineStr == null || confineStr.trim().length() == 0) {
            confineToSNPsThatMapToChromosome = null;
        } else {

            byte chrConfine = ChrAnnotation.parseChr(confineStr);


            if (chrConfine < 1) {
                confineToSNPsThatMapToChromosome = null;
            } else {
                confineToSNPsThatMapToChromosome = chrConfine;
            }
        }

        // confine to probes present in all datasets
        confineProbesToProbesPresentInAllDatasets = config.getBoolean("defaults.confine.confineToProbesPresentInAllDatasets", false);
        regressOutEQTLEffectFileName = config.getString("defaults.analysis.regressOutEQTLEffects");

        // dataset parameters
        int i = 0;
        String dataset = config.getString("datasets.dataset(" + i + ").name");  // see if a dataset is defined
        if (settingsTextToReplace != null && dataset.contains(settingsTextToReplace)) {
            dataset = dataset.replace(settingsTextToReplace, settingsTextReplaceWith);
        }

        datasetSettings = new ArrayList<MetaQTL4DatasetSettings>();



//            ArrayList<GeneticalGenomicsDataset> vGG = new ArrayList<GeneticalGenomicsDataset>();                    // create a dataset vector
//            GeneticalGenomicsDataset tmpDataset;                                               // create a temporary dataset object
        while (dataset != null) {

            String expressionData = null;
            String genotypeLocation = null;
            String genToExpCoupling = null;
            String probeannotation = null;
            String expressionplatform = null;
            String gte = null;



            // get the location of the expression data
            try {
                expressionData = config.getString("datasets.dataset(" + i + ").expressiondata");
                if (settingsTextToReplace != null && expressionData.contains(settingsTextToReplace)) {
                    expressionData = expressionData.replace(settingsTextToReplace, settingsTextReplaceWith);
                }


            } catch (Exception e) {
            }

            if (expressionData != null) {
                String expressionPlatform = null;
                try {
                    expressionPlatform = config.getString("datasets.dataset(" + i + ").expressionplatform");
                    if (settingsTextToReplace != null && expressionData.contains(settingsTextToReplace)) {
                        expressionPlatform = expressionPlatform.replace(settingsTextToReplace, settingsTextReplaceWith);
                    }
                } catch (Exception e) {
                }

                probeannotation = null;
                try {
                    probeannotation = config.getString("datasets.dataset(" + i + ").probeannotation");
                    if (settingsTextToReplace != null && expressionData.contains(settingsTextToReplace)) {
                        probeannotation = probeannotation.replace(settingsTextToReplace, settingsTextReplaceWith);
                    }
                    if (probeannotation.length() == 0) {
                        probeannotation = null;
                    }
                } catch (Exception e) {
                }
                expressionplatform = expressionPlatform;
            }

            // get the location of the dataset
            try {
                genotypeLocation = config.getString("datasets.dataset(" + i + ").location");
                if (settingsTextToReplace != null && genotypeLocation.contains(settingsTextToReplace)) {
                    genotypeLocation = genotypeLocation.replace(settingsTextToReplace, settingsTextReplaceWith);
                }
            } catch (Exception e) {
                System.out.println("Please provide a location on your disk where " + dataset + " is located");
                System.exit(-1);
            }

            Gpio.formatAsDirectory(genotypeLocation);

            // see if there is a genotype to expression couplings file
            try {
                genToExpCoupling = config.getString("datasets.dataset(" + i + ").genometoexpressioncoupling");
                if (settingsTextToReplace != null && genToExpCoupling.contains(settingsTextToReplace)) {
                    genToExpCoupling = genToExpCoupling.replace(settingsTextToReplace, settingsTextReplaceWith);
                }
            } catch (Exception e) {
            }

            dataset = null;
            i++;
            try {
                dataset = config.getString("datasets.dataset(" + i + ").name");
                if (settingsTextToReplace != null && dataset.contains(settingsTextToReplace)) {
                    dataset = dataset.replace(settingsTextToReplace, settingsTextReplaceWith);
                }

            } catch (Exception e) {
            }

            datasetSettings.add(new MetaQTL4DatasetSettings(dataset, genotypeLocation, expressionData, genToExpCoupling, expressionplatform));
        }
    }

    public String summarize() {
        String summary = "Following settings will be applied:\n"
                + "Settings\n----\n"
                + "settingsTextToReplace\t" + settingsTextToReplace + "\n"
                + "settingsTextReplaceWith\t" + settingsTextReplaceWith + "\n"
                + "\nOutput\n----\n"
                + "createSNPPValueSummaryStatisticsFile\t" + createSNPPValueSummaryStatisticsFile + "\n"
                + "createEQTLPValueTable\t" + createEQTLPValueTable + "\n"
                + "outputReportsDir\t" + outputReportsDir + "\n"
                + "\nAnalysis\n----\n"
                + "performCiseQTLAnalysis\t" + cisAnalysis + "\n"
                + "performTranseQTLAnalysis\t" + transAnalysis + "\n"
                + "performParametricAnalysis\t" + performParametricAnalysis + "\n"
                + "useAbsoluteZScorePValue\t" + useAbsoluteZScorePValue + "\n"
                + "ciseQTLAnalysMaxSNPProbeMidPointDistance\t" + ciseQTLAnalysMaxSNPProbeMidPointDistance + "\n"
                + "maxNrMostSignificantEQTLs\t" + maxNrMostSignificantEQTLs + "\n"
                + "performParametricAnalysisGetAccuratePValueEstimates\t" + performParametricAnalysisGetAccuratePValueEstimates + "\n"
                + "nrThreads\t" + nrThreads + "\n"
                + "fdrCutOff\t" + fdrCutOff + "\n"
                + "nrPermutationsFDR\t" + nrPermutationsFDR + "\n"
                + "regressOutEQTLEffectFileName\t" + regressOutEQTLEffectFileName + "\n"
                + "snpQCCallRateThreshold\t" + snpQCCallRateThreshold + "\n"
                + "snpQCHWEThreshold\t" + snpQCHWEThreshold + "\n"
                + "snpQCMAFThreshold\t" + snpQCMAFThreshold + "\n"
                + "\nConfinements\n----\n"
                + "performEQTLAnalysisOnSNPProbeCombinationSubset\t" + performEQTLAnalysisOnSNPProbeCombinationSubset + "\n"
                + "confineToSNPsThatMapToChromosome\t" + confineToSNPsThatMapToChromosome + "\n"
                + "expressionDataLoadOnlyProbesThatMapToChromosome\t" + expressionDataLoadOnlyProbesThatMapToChromosome + "\n"
                //                + "tsSNPsConfine\t" +tsSNPsConfine+ "\n"
                //                + "tsProbesConfine\t" +tsProbesConfine+ "\n"
                //                + "tsSNPProbeCombinationsConfine\t" +tsSNPProbeCombinationsConfine+ "\n"
                + "confineSNPsToSNPsPresentInAllDatasets\t" + confineSNPsToSNPsPresentInAllDatasets + "\n"
                + "confineProbesToProbesPresentInAllDatasets\t" + confineProbesToProbesPresentInAllDatasets + "\n"
                + "confineToProbesThatMapToChromosome\t" + confineToProbesThatMapToChromosome + "\n"
                + "expressionDataLoadOnlyProbesThatMapToChromosome\t" + expressionDataLoadOnlyProbesThatMapToChromosome + "\n"
                + "\n";

        summary += "\nDatasets\n----\n";

        for (MetaQTL4DatasetSettings settings : this.datasetSettings) {
            summary += "DatasetName\t" + settings.getName() + "\n";
            summary += "ExpressionData\t" + settings.getTraitDataLocation() + "\n";
            summary += "ExpressionPlatform\t" + settings.getAnnotationFile() + "\n";
            summary += "GenotypeData\t" + settings.getGenotypeLocation() + "\n";
            summary += "GTE\t" + settings.getGenotypeToTraitCoupling() + "\n";
        }


        return summary;

    }

    public void writeSettingsToDisk() throws Exception {
        TextFile tf = new TextFile(this.outputReportsDir + "/UsedSettings.txt", TextFile.W);
        tf.write(summarize());
        tf.close();
    }

    public String getSettingsTextToReplace() {
        return settingsTextToReplace;
    }

    public String getSettingsTextReplaceWith() {
        return settingsTextReplaceWith;
    }

    public boolean isCreateSNPPValueSummaryStatisticsFile() {
        return createSNPPValueSummaryStatisticsFile;
    }

    public boolean isCreateSNPSummaryStatisticsFile() {
        return createSNPSummaryStatisticsFile;
    }

    public boolean isCreateEQTLPValueTable() {
        return createEQTLPValueTable;
    }

    public File getOutputReportsDir() {
        return outputReportsDir;
    }

    public boolean isPerformParametricAnalysis() {
        return performParametricAnalysis;
    }

    public boolean isUseAbsoluteZScorePValue() {
        return useAbsoluteZScorePValue;
    }

    public int getCiseQTLAnalysMaxSNPProbeMidPointDistance() {
        return ciseQTLAnalysMaxSNPProbeMidPointDistance;
    }

    public int getMaxNrMostSignificantEQTLs() {
        return maxNrMostSignificantEQTLs;
    }

    public boolean isPerformParametricAnalysisGetAccuratePValueEstimates() {
        return performParametricAnalysisGetAccuratePValueEstimates;
    }

    public Integer getNrThreads() {
        return nrThreads;
    }

    public double getFdrCutOff() {
        return fdrCutOff;
    }

    public int getNrPermutationsFDR() {
        return nrPermutationsFDR;
    }

    public boolean isPerformEQTLAnalysisOnSNPProbeCombinationSubset() {
        return performEQTLAnalysisOnSNPProbeCombinationSubset;
    }

    public Byte getConfineToSNPsThatMapToChromosome() {
        return confineToSNPsThatMapToChromosome;
    }

    public boolean isExpressionDataLoadOnlyProbesThatMapToChromosome() {
        return expressionDataLoadOnlyProbesThatMapToChromosome;
    }

    public double getPlotOutputPValueCutOff() {
        return plotOutputPValueCutOff;
    }

    public String getPlotOutputDirectory() {
        return plotOutputDirectory;
    }

    public boolean isRunOnlyPermutations() {
        return runOnlyPermutations;
    }

    public Integer getStartWithPermutation() {
        return startWithPermutation;
    }

    public Boolean getConfineSNPsToSNPsPresentInAllDatasets() {
        return confineSNPsToSNPsPresentInAllDatasets;
    }

    public boolean isConfineProbesToProbesPresentInAllDatasets() {
        return confineProbesToProbesPresentInAllDatasets;
    }

    public ArrayList<MetaQTL4DatasetSettings> getDatasetSettings() {
        return datasetSettings;
    }

    public String getRegressOutEQTLEffectFileName() {
        return regressOutEQTLEffectFileName;
    }

    public Double getSnpQCCallRateThreshold() {
        return snpQCCallRateThreshold;
    }

    public Double getSnpQCHWEThreshold() {
        return snpQCHWEThreshold;
    }

    public Double getSnpQCMAFThreshold() {
        return snpQCMAFThreshold;
    }

    public Byte getConfineToProbesThatMapToChromosome() {
        return confineToProbesThatMapToChromosome;
    }

    public boolean isCreateBinaryOutputFiles() {
        return createBinaryOutputFiles;
    }

    public boolean isCreateTEXTOutputFiles() {
        return createTEXTOutputFiles;
    }

    public String getStrConfineSNP() {
        return strConfineSNP;
    }

    public String getStrConfineProbe() {
        return strConfineProbe;
    }

    public boolean isProvideFoldChangeData() {
        return provideFoldChangeData;
    }

    public boolean isProvideBetasAndStandardErrors() {
        return provideBetasAndStandardErrors;
    }

    public boolean isEqualRankForTies() {
        return equalRankForTies;
    }

    public boolean isCreateQQPlot() {
        return createQQPlot;
    }

    public boolean isCreateDotPlot() {
        return createDotPlot;
    }

    public String getStrSNPProbeConfine() {
        return strSNPProbeConfine;
    }

    public boolean isCisAnalysis() {
        return cisAnalysis;
    }

    public boolean isTransAnalysis() {
        return transAnalysis;
    }

    public int getNrDatasets() {
        return datasetSettings.size();
    }
    
    public String getProbeAnnotationFile(){
        return probeAnnotationFile;
    }
}
