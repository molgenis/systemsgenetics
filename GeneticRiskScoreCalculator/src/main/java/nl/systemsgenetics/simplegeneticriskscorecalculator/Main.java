/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.simplegeneticriskscorecalculator;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.ConvertDoubleMatrixDataToTriTyper;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author MarcJan
 */
public class Main {

    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");

    public static void main(String[] args) {
        CommandLineParser parser = new GnuParser();
        Options options = new Options();

        Option FileOut = OptionBuilder.withArgName("path").hasArg().withDescription("Location (folder) for the output.").withLongOpt("OutputFolder").create("o");
        Option GenotypeTypeIn = OptionBuilder.withArgName("type").hasArg().withDescription("Type of reference data.").withLongOpt("GenotypeType").create("gt");
        Option GenotypeIn = OptionBuilder.withArgName("path").hasArg().withDescription("Location for the reference data.").withLongOpt("GenotypeLocation").create("gi");
        Option InFolder = OptionBuilder.withArgName("path").hasArg().withDescription("Location of the folder with genetic risk score information.").withLongOpt("input").create("i");
        Option rSquared = OptionBuilder.withArgName("double").hasArg().withDescription("R2 for pruning.").withLongOpt("rSquared").create("r");
        Option pValueThreshold = OptionBuilder.withArgName("double").hasArg().withDescription("P-value thresholds for genetic risk score inclusion, colon separated should be ordered from most stringent to least stringent.").withLongOpt("pValue").create("p");
        Option WindowSize = OptionBuilder.withArgName("double").hasArg().withDescription("Window size for pruning, if given two window-sizes (colon separated), a two step window approach is used.").withLongOpt("wSize").create("w");
        Option debugOpt = OptionBuilder.withArgName("boolean").withDescription("Switch on debugging.").withLongOpt("debug").create('d');
        Option sumRisksOpt = OptionBuilder.withArgName("boolean").withDescription("Sum risks.").withLongOpt("riskSum").create('s');
        Option excludeGenomicRange = OptionBuilder.withArgName("String").hasArg().withDescription("Exclude genomic range(s) from the risk score calculation. Range needs to be specified as: \"6:101-110;6:250000-350000. Warning: Chr name must be specified as expected in the genotype dataset.").withLongOpt("excludeRange").create("er");
        Option unWeightedScore = OptionBuilder.withArgName("boolean").withDescription("Use unweighted combination of risk factors.").withLongOpt("unWeighted").create('u');
        options.addOption(FileOut).addOption(GenotypeTypeIn).addOption(GenotypeIn).addOption(InFolder).addOption(rSquared).addOption(pValueThreshold).addOption(WindowSize).addOption(debugOpt).addOption(excludeGenomicRange).addOption(unWeightedScore).addOption(sumRisksOpt);

        String genotypePath = null;
        String genotypeType = null;
        String riskFolder = null;
        File outputFolder = null;
        double rSquare = 1.0d;
        double[] windowSize = null;
        double[] pValThres = null;
        String[] genomicRangesToExclude = null;
        boolean debugMode = false;
        boolean unweighted = false;
        boolean sumRisks = false;

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
            HelpFormatter formatter = new HelpFormatter();

            if (cmd.hasOption("OutputFolder") || cmd.hasOption("o")) {
                // initialise the member variable
                outputFolder = new File(cmd.getOptionValue("OutputFolder"));
            } else {
                System.out.println("Missing necessary information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("GenotypeLocation") || cmd.hasOption("gi")) {
                // initialise the member variable
                genotypePath = cmd.getOptionValue("GenotypeLocation");
            } else {
                System.out.println("Missing necessary information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("GenotypeType") || cmd.hasOption("gt")) {
                // initialise the member variable
                genotypeType = cmd.getOptionValue("GenotypeType");
            } else {
                System.out.println("Missing necessary information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("input") || cmd.hasOption("i")) {
                // initialise the member variable
                riskFolder = cmd.getOptionValue("input");
            } else {
                System.out.println("Missing necessary information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("rSquared") || cmd.hasOption("r")) {
                // initialise the member variable
                rSquare = Double.parseDouble(cmd.getOptionValue("rSquared"));
            } else {
                System.out.println("Missing necessary information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("wSize") || cmd.hasOption("w")) {
                // initialise the member variable
                String[] tempThres = cmd.getOptionValue("wSize").split(":");
                windowSize = new double[tempThres.length];
                for (int i = 0; i < tempThres.length; i++) {
                    windowSize[i] = Double.parseDouble(tempThres[i]);
                }

            } else {
                System.out.println("Missing necessary information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("pValue") || cmd.hasOption("p")) {
                // initialise the member variable
                String[] tempThres = cmd.getOptionValue("pValue").split(":");
                pValThres = new double[tempThres.length];
                for (int i = 0; i < tempThres.length; i++) {
                    pValThres[i] = Double.parseDouble(tempThres[i]);
                }
            } else {
                System.out.println("Missing necessary information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("excludeRange") || cmd.hasOption("er")) {
                // initialise the member variable
                genomicRangesToExclude = cmd.getOptionValue("excludeRange").split(";");
            }
            debugMode = cmd.hasOption('d');
            unweighted = cmd.hasOption('u');
            sumRisks = cmd.hasOption('s');

        } catch (ParseException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        try {
            if (outputFolder!=null && !(outputFolder.exists())) {
                Gpio.createDir(outputFolder.getAbsolutePath());
            }
            RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOf(genotypeType).createFilteredGenotypeData(genotypePath, 750000, null, null);
            THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks = readRiskFiles(genotypeData, riskFolder, pValThres, genomicRangesToExclude, unweighted, debugMode);
            if (windowSize.length == 1) {
                DoubleMatrixDataset<String, String> geneticRiskScoreMatrix = CalculateSimpleGeneticRiskScore.calculate(genotypeData, risks, outputFolder, rSquare, windowSize[0], debugMode, pValThres, sumRisks);
                writeMatrixToFile(geneticRiskScoreMatrix, outputFolder);
            } else if (windowSize.length == 2) {
                DoubleMatrixDataset<String, String> geneticRiskScoreMatrix = CalculateSimpleGeneticRiskScore.calculateTwoStages(genotypeData, risks, outputFolder, rSquare, windowSize, debugMode, pValThres, sumRisks);
                writeMatrixToFile(geneticRiskScoreMatrix, outputFolder);
            } else {
                System.out.println("More than two window-sizes is not supported.");
                System.exit(0);
            }

        } catch (IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> readRiskFiles(RandomAccessGenotypeData genotypeData, String riskFolder, double[] pValueThreshold, String[] genomicRangesToExclude, boolean unweighted, boolean debugMode) {
        THashMap<String, ArrayList<Pair<Integer, Integer>>> exclussionRanges = new THashMap<>();
        
        if(genomicRangesToExclude!=null){
            System.out.println("Trying to exclude genomic ranges.");
            int ranges = 0;
            for (String s : genomicRangesToExclude) {
                String[] parts = s.split(":");
                String key = parts[0];
                parts = parts[1].split("-");

                if (!exclussionRanges.contains(key)) {
                    exclussionRanges.put(key, new ArrayList<Pair<Integer, Integer>>());
                }
                exclussionRanges.get(key).add(new Pair(Integer.parseInt(parts[0]), Integer.parseInt(parts[1])));
                ranges++;
            }
            if(debugMode){
                System.out.println("Number of ranges excluded: "+ranges+" on: "+exclussionRanges.size()+" chromosomes");
            }
        }
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks = new THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>>();

        File riskFileFolder = new File(riskFolder);
        if (!riskFileFolder.exists()) {
            System.out.println("Warning: input risk folder does not exists:\n");
            System.out.println(riskFileFolder);
            System.exit(-1);
        }
        File[] riskFiles = riskFileFolder.listFiles();
        
        
        for (File f : riskFiles) {
            THashSet<String> chromosomesExcluded = new THashSet<>();
            int snpsExcluded = 0;
            try {
                TextFile readFiles = new TextFile(f.getAbsolutePath(), TextFile.R);

                System.out.println(f.getName());
                for (double p : pValueThreshold) {
                    String name = f.getName();
                    String name2 = "_P" + p;
                    if (!risks.containsKey(name)) {
                        risks.put(name, new THashMap<String, THashMap<String, ArrayList<RiskEntry>>>());
                    }
                    if (!risks.get(name).containsKey(name2)) {
                        risks.get(name).put(name2, new THashMap<String, ArrayList<RiskEntry>>());
                    }
                }

                String s = readFiles.readLine();
                while ((s = readFiles.readLine()) != null) {
                    String[] parts = TAB_PATTERN.split(s);
//                    System.out.println(s);
                    if (genotypeData.getVariantIdMap().containsKey(parts[0])) {
                        GeneticVariant snpObject = genotypeData.getVariantIdMap().get(parts[0]);
//                        System.out.print(snpObject.getSequenceName() + "\t" + snpObject.getStartPos() + "\n");
                        double currentP = Double.parseDouble(parts[3]);
                        boolean addEntry = true;
                        
                        if(unweighted){
                            if(parts[2].startsWith("-")){
                                parts[2] = "-1";
                            } else {
                                parts[2] = "1";
                            }
                        }
                        
                        if(exclussionRanges.contains(snpObject.getSequenceName())){
                            chromosomesExcluded.add(snpObject.getSequenceName());
                            for(Pair<Integer, Integer> p : exclussionRanges.get(snpObject.getSequenceName())){
                                if(p.getLeft()<=snpObject.getStartPos() && p.getRight()>=snpObject.getStartPos()){
                                    addEntry = false;
                                    snpsExcluded++;
                                }
                            }
                        }

                        if (addEntry) {
                            for (double p : pValueThreshold) {
                                if (currentP < p) {
                                    String name = f.getName();
                                    String name2 = "_P" + p;

                                    if (!risks.get(name).get(name2).containsKey(snpObject.getSequenceName())) {
                                        risks.get(name).get(name2).put(snpObject.getSequenceName(), new ArrayList<RiskEntry>());
                                    }
                                    risks.get(name).get(name2).get(snpObject.getSequenceName()).add(new RiskEntry(parts[0], snpObject.getSequenceName(), snpObject.getStartPos(), parts[1].charAt(0), parts[2], currentP));
                                }
                            }
                        }
                    }
                }
                if(debugMode){
                    System.out.println("Chromosomes where regios are excluded: " + chromosomesExcluded);
                    System.out.println("Number of SNPs excluded: " + snpsExcluded);
                }
                readFiles.close();
            } catch (IOException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            }

        }

        for (Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> e : risks.entrySet()) {

            for (Entry<String, THashMap<String, ArrayList<RiskEntry>>> e2 : e.getValue().entrySet()) {
                int entries = 0;
                for (Entry<String, ArrayList<RiskEntry>> e3 : e2.getValue().entrySet()) {
                    Collections.sort(e3.getValue());
                    entries += e3.getValue().size();
                }
//                System.out.println(e.getKey()+e2.getKey()+" has: "+entries+" entries");
            }

        }

        return risks;
    }

    private static void writeMatrixToFile(DoubleMatrixDataset<String, String> geneticRiskScoreMatrix, File outputFolder) throws IOException {

        String outputF = outputFolder + File.separator + "TT";
        if (!(new File(outputF).exists())) {
            Gpio.createDir(outputF);
        } else if (!(new File(outputF).isDirectory())) {
            System.out.println("Error file is already there but not a directory, writing to TT1");
            outputF = outputFolder + File.separator + "TT1";
            Gpio.createDir(outputF);

            System.exit(0);
        }
        geneticRiskScoreMatrix.save(outputF + File.separator + "rawScoreMatrix.txt");
        try {
            System.out.println("Writing SNPMappings.txt & SNPs.txt to file:");

            java.io.BufferedWriter outSNPMappings = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputF + File.separator + "SNPMappings.txt")));
            java.io.BufferedWriter outSNPs = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputF + File.separator + "SNPs.txt")));

            for (String var : geneticRiskScoreMatrix.getRowObjects()) {
                outSNPMappings.write("1\t1000000\t" + var + '\n');
                outSNPs.write(var + '\n');
            }

            outSNPMappings.close();
            outSNPs.close();
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
            System.exit(0);
        }

        geneticRiskScoreMatrix.setMatrix(ConvertDoubleMatrixDataToTriTyper.rankRows(geneticRiskScoreMatrix.getMatrix()));

        geneticRiskScoreMatrix.setMatrix(ConvertDoubleMatrixDataToTriTyper.rescaleValue(geneticRiskScoreMatrix.getMatrix(), 200.0d));

        try {
            System.out.println("\nWriting Individuals.txt and Phenotype.txt to file:");
            BufferedWriter outIndNew = new BufferedWriter(new FileWriter(outputF + File.separator + "Individuals.txt"));
            BufferedWriter outPhenoNew = new BufferedWriter(new FileWriter(outputF + File.separator + "PhenotypeInformation.txt"));
            for (String ind : geneticRiskScoreMatrix.getColObjects()) {
                outIndNew.write(ind + '\n');
                outPhenoNew.write(ind + "\tcontrol\tinclude\tmale\n");
            }
            outIndNew.close();
            outPhenoNew.close();
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
            System.exit(0);
        }

        int nrSNPs = geneticRiskScoreMatrix.rows();
        int nrSamples = geneticRiskScoreMatrix.columns();

        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, new File(outputF + File.separator + "GenotypeMatrix.dat"), false);
        WGAFileMatrixImputedDosage fileMatrixDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamples, new File(outputF + File.separator + "ImputedDosageMatrix.dat"), false);
        byte[] alleles = new byte[2];
        alleles[0] = 84;
        alleles[1] = 67;

        for (int snp = 0; snp < nrSNPs; snp++) {
            DoubleMatrix1D snpRow = geneticRiskScoreMatrix.getMatrix().viewRow(snp);

            byte[] allele1 = new byte[nrSamples];
            byte[] allele2 = new byte[nrSamples];
            byte[] dosageValues = new byte[nrSamples];
            for (int ind = 0; ind < nrSamples; ind++) {
                if (snpRow.get(ind) > 100) {
                    allele1[ind] = alleles[1];
                    allele2[ind] = alleles[1];
                } else {
                    allele1[ind] = alleles[0];
                    allele2[ind] = alleles[0];
                }

                int dosageInt = (int) Math.round(snpRow.get(ind));
                byte value = (byte) (Byte.MIN_VALUE + dosageInt);
                dosageValues[ind] = value;
            }
            fileMatrixGenotype.setAllele1(snp, 0, allele1);
            fileMatrixGenotype.setAllele2(snp, 0, allele2);
            fileMatrixDosage.setDosage(snp, 0, dosageValues);
        }
        fileMatrixGenotype.close();
        fileMatrixDosage.close();
    }

}
