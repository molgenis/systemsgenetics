/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3.containers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class Settings extends TriTyperGeneticalGenomicsDatasetSettings {

    // Output
    public String settingsTextToReplace = null;                                // Replace this text in the configuration XML file
    public String settingsTextReplaceWith = null;                              // Replace the text in settingsTextToReplace with this text
    public boolean createSNPPValueSummaryStatisticsFile = false;               // Output SNP P-Value summary statistics
    public boolean createSNPSummaryStatisticsFile = false;                     // Output SNP P-Value summary statistics
    public boolean createEQTLPValueTable = false;                              // Output an eQTL p-value table (only applies for a limited number (500) of SNP)
    public String outputReportsDir;                                            // Output directory for reports
    // SNP QC
    // Analysis settings
    public boolean performCiseQTLAnalysis;                                     // Perform Cis analysis
    public boolean performTranseQTLAnalysis;                                   // Perform Trans analysis
    public boolean performParametricAnalysis = false;                          // Perform parametric analysis
    public boolean useAbsoluteZScorePValue = false;                            // Use absolute Z-score? (required for finding opposite allelic effects)
    public int ciseQTLAnalysMaxSNPProbeMidPointDistance;                       // Midpoint distance for declaring an eQTL effect CIS
    public int maxNrMostSignificantEQTLs = 150000;                             // Max number of results stored in memory
    public boolean performParametricAnalysisGetAccuratePValueEstimates;         // Use an accurate estimation of the P-values
    public Integer nrThreads;                                                      // Use this number of threads
    // Multiple testing correction
    public double fdrCutOff;                                                   // Cutoff for FDR procedure
    public int nrPermutationsFDR;                                              // Number of permutations to determine FDR
    // confinements
    public boolean performEQTLAnalysisOnSNPProbeCombinationSubset;             // Confine to a certain set of probe/snp combinations?
    public Byte confineToSNPsThatMapToChromosome;                               // Confine SNP to be assessed to SNPs mapped on this chromosome
    public boolean expressionDataLoadOnlyProbesThatMapToChromosome = false;    // Only load expression data for probes with a known chromosome mapping
    public HashSet<String> tsSNPsConfine;                           // Confine analysis to the SNPs in this hash
    public HashMap<String, HashSet<String>> tsSNPProbeCombinationsConfine;           // Confine analysis to the combinations of SNP and Probes in this hash
    // plots
    public double plotOutputPValueCutOff;                                      // Use this p-value as a cutoff for drawing plots
    public String plotOutputDirectory;                                         // Print the plots in this directory
    public boolean runOnlyPermutations = false;
    public Integer startWithPermutation;

    
    public Boolean confineSNPsToSNPsPresentInAllDatasets;

    public boolean confineProbesToProbesPresentInAllDatasets;
    public ArrayList<TriTyperGeneticalGenomicsDatasetSettings> datasetSettings;
    public String regressOutEQTLEffectFileName;
    public Double snpQCCallRateThreshold;
    public Double snpQCHWEThreshold;
    public Double snpQCMAFThreshold;
    public Byte confineToProbesThatMapToChromosome;
    public boolean createCHARGEOutputFiles;
    public boolean createTEXTOutputFiles;

    
    public Settings() {
        
    }
    
    public void load(String settings) throws IOException, ConfigurationException {

        XMLConfiguration config = new XMLConfiguration(settings);           // Use the apache XML configuration parser

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
            snpQCHWEThreshold = 0.0000001;
        }

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

        try {
            runOnlyPermutations = config.getBoolean("defaults.analysis.onlypermutations", false);
        } catch (Exception e) {
        }

        try {
            String strStartWithPermutation = config.getString("defaults.analysis.startwithpermutation", null);

            if (settingsTextToReplace != null && strStartWithPermutation.contains(settingsTextToReplace)) {
                strStartWithPermutation = strStartWithPermutation.replace(settingsTextToReplace, settingsTextReplaceWith);
            }
            startWithPermutation = Integer.parseInt(strStartWithPermutation);

        } catch (Exception e) {
        }

        if (analysisType != null) {
            if (analysisType.toLowerCase().equals("cis")) {
                performCiseQTLAnalysis = true;
                performTranseQTLAnalysis = false;
            } else if (analysisType.toLowerCase().equals("trans")) {
                performCiseQTLAnalysis = false;
                performTranseQTLAnalysis = true;
            } else if (analysisType.toLowerCase().equals("cistrans")) {
                performCiseQTLAnalysis = true;
                performTranseQTLAnalysis = true;
            }
        } else {
            performCiseQTLAnalysis = true;
            performTranseQTLAnalysis = false;
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
            nrThreads = Runtime.getRuntime().availableProcessors();
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
            outputReportsDir = outdir;
            if (!outputReportsDir.endsWith("/")) {
                outputReportsDir += "/";
            }

            // check if dir exists. if it does not, create it:

            if (!Gpio.exists(outdir)) {
                Gpio.createDir(outdir);
            }

            config.save(outputReportsDir + "metaqtlsettings.xml");


        } else {
            System.out.println("Error: please supply an output directory.");
            System.exit(-1);
        }

        try {
            String probeConversionFileLoc = config.getString("defaults.analysis.probeconversion", null);
            if (probeConversionFileLoc != null && probeConversionFileLoc.length() > 0) {
                System.out.println(probeConversionFileLoc);


            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }


        try {
            outputplotthreshold = config.getDouble("defaults.output.outputplotthreshold");
        } catch (Exception e) {
        }
        
        try {
            createCHARGEOutputFiles = config.getBoolean("defaults.output.chargeoutput",false);
            createTEXTOutputFiles = config.getBoolean("defaults.output.textoutput",true);
        } catch (Exception e) {
        }

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

            plotOutputDirectory = outputplotdirectory;
            if (!plotOutputDirectory.endsWith("/")) {
                plotOutputDirectory += "/";
            }

            // check if dir exists. if it does not, create it:
            if (!Gpio.exists(plotOutputDirectory)) {
                Gpio.createDir(plotOutputDirectory);
            }

        } else {
            plotOutputDirectory = outdir;
        }

        try {
            createSNPPValueSummaryStatisticsFile = config.getBoolean("defaults.output.generatesnppvaluesummarystatistics", false);
        } catch (Exception e) {
        }

        try {
            createSNPSummaryStatisticsFile = config.getBoolean("defaults.output.generatesnpsummarystatistics", false);
        } catch (Exception e) {
        }

        try {
            createEQTLPValueTable = config.getBoolean("defaults.output.generateeqtlpvaluetable", false);
        } catch (Exception e) {
        }

        try {
            maxNrMostSignificantEQTLs = config.getInt("defaults.output.maxnreqtlresults", 150000);
        } catch (Exception e) {
            e.printStackTrace();
        }

        //Load only expression probes that map to a known chromosome:
        try {
            expressionDataLoadOnlyProbesThatMapToChromosome = config.getBoolean("defaults.confine.confineProbesThatMapToKnownChromosome");
        } catch (Exception e) {
        }

        
        // confinements on snp, probe, or snp-probe
        String confineSNP       = null;
        String confineProbe     = null;
        String snpProbeConfine  = null;

        confineSNPsToSNPsPresentInAllDatasets = null;

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

        if(confineSNP != null && confineSNP.trim().length() > 0 && Gpio.exists(confineSNP)){
            TextFile in = new TextFile(confineSNP, TextFile.R);
            tsSNPsConfine = new HashSet<String>();
            String[] data = in.readAsArray();
            for(String d: data){
                tsSNPsConfine.add(new String(d.getBytes("UTF-8")));
            }
           
        } 
        
        if(confineProbe != null && confineProbe.trim().length() > 0 && Gpio.exists(confineProbe)){
            TextFile in = new TextFile(confineProbe, TextFile.R);
            tsProbesConfine = new HashSet<String>();
            String[] data = in.readAsArray();
            for(String d: data){
                tsProbesConfine.add(new String(d.getBytes("UTF-8")));
            }
            
        }
        
        if(snpProbeConfine != null && snpProbeConfine.trim().length() > 0 && Gpio.exists(snpProbeConfine)){
            TextFile in = new TextFile(snpProbeConfine, TextFile.R);
            tsSNPProbeCombinationsConfine = new HashMap<String, HashSet<String>>();
            String[] elems = in.readLineElemsReturnReference(TextFile.tab);
            while(elems!=null){
                if(elems.length > 1){
                    String snp = new String(elems[0].getBytes("UTF-8"));
                    if( (tsSNPsConfine != null && tsSNPsConfine.contains(snp)) || tsSNPsConfine == null){
                        HashSet<String> probes = tsSNPProbeCombinationsConfine.get(snp);
                        if(probes == null){
                            probes = new HashSet<String>();
                        }
                        String probe = new String(elems[1].getBytes("UTF-8"));
                        if( (tsProbesConfine != null && tsProbesConfine.contains(probe)) || tsProbesConfine == null){
                            probes.add(probe);
                        }
                    } 
                }
                elems = in.readLineElemsReturnReference(TextFile.tab);
            }
            in.close();
        }

        
        // confine to snp present in all datasets
        try {
            confineSNPsToSNPsPresentInAllDatasets = config.getBoolean("defaults.confine.confineSNPsToSNPsPresentInAllDatasets");
        } catch (Exception e) {
        }

        // confine to SNP that map to this chromosome
        try {
            String confineStr = config.getString("defaults.confine.confineToSNPsThatMapToChromosome", null);
            if(confineStr == null || confineStr.trim().length() == 0){
                confineToSNPsThatMapToChromosome = null; 
            } else {
                
                confineToSNPsThatMapToChromosome = ChrAnnotation.parseChr(confineStr);
                if(confineToSNPsThatMapToChromosome < 1){
                    confineToSNPsThatMapToChromosome = null;
                }
            }
            
        } catch (Exception e) {
        }

        // confine to probes present in all datasets
        confineProbesToProbesPresentInAllDatasets = false;
        try {
            confineProbesToProbesPresentInAllDatasets = config.getBoolean("defaults.confine.confineToProbesPresentInAllDatasets", false);
        } catch (Exception e) {
        }

        

        regressOutEQTLEffectFileName = null;
        try {
            regressOutEQTLEffectFileName = config.getString("defaults.analysis.regressOutEQTLEffects");
        } catch (Exception e) {
        }

        
        // dataset parameters
        int i = 0;

        String dataset = config.getString("datasets.dataset(" + i + ").name");  // see if a dataset is defined
        if (settingsTextToReplace != null && dataset.contains(settingsTextToReplace)) {
            dataset = dataset.replace(settingsTextToReplace, settingsTextReplaceWith);
        }
      
        datasetSettings = new ArrayList<TriTyperGeneticalGenomicsDatasetSettings>();

        

//            ArrayList<GeneticalGenomicsDataset> vGG = new ArrayList<GeneticalGenomicsDataset>();                    // create a dataset vector
//            GeneticalGenomicsDataset tmpDataset;                                               // create a temporary dataset object
        while (dataset != null) {
            
            String expressionData = null;
            String dataloc = null;
            String genToExpCoupling = null;
            Boolean qnorm = false;
            Boolean logtr = false;

            TriTyperGeneticalGenomicsDatasetSettings s = new TriTyperGeneticalGenomicsDatasetSettings();
            s.name = dataset;
            datasetSettings.add(s);
            System.out.println("Loading dataset specific settings for: "+dataset);
            
            // get the location of the expression data
            try {
                expressionData = config.getString("datasets.dataset(" + i + ").expressiondata");
                if (settingsTextToReplace != null && expressionData.contains(settingsTextToReplace)) {
                    expressionData = expressionData.replace(settingsTextToReplace, settingsTextReplaceWith);
                }
            } catch (Exception e) {
            }

            s.expressionLocation = expressionData;

            // get the location of the dataset
            try {
                dataloc = config.getString("datasets.dataset(" + i + ").location");
                if (settingsTextToReplace != null && dataloc.contains(settingsTextToReplace)) {
                    dataloc = dataloc.replace(settingsTextToReplace, settingsTextReplaceWith);
                }
            } catch (Exception e) {
                System.out.println("Please provide a location on your disk where " + dataset + " is located");
                System.exit(-1);
            }

            if (!dataloc.endsWith("/")) {
                dataloc += "/";
            }

            s.genotypeLocation = dataloc;



            // see if there is a genotype to expression couplings file
            try {
                genToExpCoupling = config.getString("datasets.dataset(" + i + ").genometoexpressioncoupling");
                if (settingsTextToReplace != null && genToExpCoupling.contains(settingsTextToReplace)) {
                    genToExpCoupling = genToExpCoupling.replace(settingsTextToReplace, settingsTextReplaceWith);
                }
            } catch (Exception e) {
            }

            s.genotypeToExpressionCoupling = genToExpCoupling;

            // quantile normalize the expression data?
            try {
                qnorm = config.getBoolean("datasets.dataset(" + i + ").quantilenormalize", false);
            } catch (Exception e) {
            }

            s.quantilenormalize = qnorm;

//                if (qnorm) {
//                    tmpDataset.quantileNormalize();
//                }
//
            // log2 transform the expression data?
            try {
                logtr = config.getBoolean("datasets.dataset(" + i + ").logtranform", false);
            } catch (Exception e) {
            }
            s.logtransform = logtr;

            dataset = null;
            i++;
            try {
                dataset = config.getString("datasets.dataset(" + i + ").name");
                if (settingsTextToReplace != null && dataset.contains(settingsTextToReplace)) {
                    dataset = dataset.replace(settingsTextToReplace, settingsTextReplaceWith);
                }

            } catch (Exception e) {
            }

            s.confineProbesToProbesMappingToAnyChromosome = confineProbesToProbesMappingToAnyChromosome;
            s.confineProbesToProbesThatMapToChromosome = confineProbesToProbesThatMapToChromosome;
            s.tsProbesConfine = tsProbesConfine;
            

        }
        
        System.out.println("Done.");
        
        summarize();

        

        
    }

    private void summarize() {
        System.out.println("Following settings will be applied:\n"
                + "Settings\n----\n"
                + "settingsTextToReplace\t" +settingsTextToReplace                                                 + "\n"
                + "settingsTextReplaceWith\t" +settingsTextReplaceWith+ "\n"
                
                
                + "\nOutput\n----\n"
                + "createSNPPValueSummaryStatisticsFile\t" +createSNPPValueSummaryStatisticsFile+ "\n"
                + "createEQTLPValueTable\t" +createEQTLPValueTable+ "\n"
                + "outputReportsDir\t" +outputReportsDir+ "\n"
                
                
                + "\nAnalysis\n----\n"
                + "performCiseQTLAnalysis\t" +performCiseQTLAnalysis+ "\n"
                + "performTranseQTLAnalysis\t" +performTranseQTLAnalysis+ "\n"
                + "performParametricAnalysis\t" +performParametricAnalysis+ "\n"
                    
                + "useAbsoluteZScorePValue\t" +useAbsoluteZScorePValue+ "\n"
                + "ciseQTLAnalysMaxSNPProbeMidPointDistance\t" +ciseQTLAnalysMaxSNPProbeMidPointDistance+ "\n"
                + "maxNrMostSignificantEQTLs\t" +maxNrMostSignificantEQTLs+ "\n"
                + "performParametricAnalysisGetAccuratePValueEstimates\t" +performParametricAnalysisGetAccuratePValueEstimates+ "\n"
                + "nrThreads\t" +nrThreads+ "\n"
                + "fdrCutOff\t" +fdrCutOff+ "\n"
                + "nrPermutationsFDR\t" +nrPermutationsFDR+ "\n"
                
                + "regressOutEQTLEffectFileName\t" +regressOutEQTLEffectFileName+ "\n"
                + "snpQCCallRateThreshold\t" +snpQCCallRateThreshold+ "\n"
                + "snpQCHWEThreshold\t" +snpQCHWEThreshold+ "\n"
                + "snpQCMAFThreshold\t" +snpQCMAFThreshold+ "\n"
                
                + "\nConfinemens\n----\n"
                + "performEQTLAnalysisOnSNPProbeCombinationSubset\t" +performEQTLAnalysisOnSNPProbeCombinationSubset+ "\n"
                + "confineToSNPsThatMapToChromosome\t" +confineToSNPsThatMapToChromosome+ "\n"
                + "expressionDataLoadOnlyProbesThatMapToChromosome\t" +expressionDataLoadOnlyProbesThatMapToChromosome+ "\n"
                + "tsSNPsConfine\t" +tsSNPsConfine+ "\n"
                + "tsProbesConfine\t" +tsProbesConfine+ "\n"
                + "tsSNPProbeCombinationsConfine\t" +tsSNPProbeCombinationsConfine+ "\n"
                + "confineSNPsToSNPsPresentInAllDatasets\t" +confineSNPsToSNPsPresentInAllDatasets+ "\n"
                + "confineProbesToProbesPresentInAllDatasets\t" +confineProbesToProbesPresentInAllDatasets+ "\n"
                + "confineToProbesThatMapToChromosome\t" +confineToProbesThatMapToChromosome+ "\n"
                + "expressionDataLoadOnlyProbesThatMapToChromosome\t" +expressionDataLoadOnlyProbesThatMapToChromosome+ "\n"
                + "\n");
        
        
     

    
   
    }
}
