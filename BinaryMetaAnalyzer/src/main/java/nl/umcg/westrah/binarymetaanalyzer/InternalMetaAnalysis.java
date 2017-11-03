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
import umcg.genetica.io.text.TextFile;
import javax.xml.bind.annotation.adapters.HexBinaryAdapter;
import umcg.genetica.io.bin.BinaryFile;

import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;
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
 * @author Harm-Jan
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

        InternalMetaAnalysis meta = new InternalMetaAnalysis(settingsFile, textToReplace, replaceTextWith);
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
                traitLocationMap.put(parts[1],index);
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

        for (int permutation = settings.getStartPermutations(); permutation <= settings.getNrPermutations(); permutation++) {
            boolean runningPermutation = true;
            if (permutation == 0) {
                runningPermutation = false;
            }

            //Initialize the new binaryOutput
            //Original initialize binary matrix
            String fileName = settings.getOutput() + settings.getDatasetname();

            if (runningPermutation) {
                fileName += "-PermutationRound-" + permutation;
            }
            BinaryFile zScoreBinaryFile = new BinaryFile(fileName + ".dat", BinaryFile.W);

            System.out.println("Loading dataset");
            dataset = new InternalMetaAnalysisDataset(settings.getDatasetlocation(),
                    settings.getDatasetname(),
                    settings.getDatasetPrefix(),
                    permutation);

            System.out.println("Loaded");

            // create meta-analysis SNP index. have to recreate this every permutation,
            // since the order of SNPs is generated at random.
            System.out.println("Creating SNP index");
            createSNPIndex(outdir);
            System.out.println("Total of " + snpIndex.length + " SNPs");

            // write magic number
            if (dataset.getIsCisDataset()) {
                zScoreBinaryFile.writeInt(1);
            } else {
                zScoreBinaryFile.writeInt(0);
            }

            TextFile zScoreRowNamesFile = new TextFile(fileName + "-RowNames.txt.gz", TextFile.W);
            zScoreRowNamesFile.writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate");

            TextFile tf = new TextFile(fileName + "-ColNames.txt.gz", TextFile.W);
            tf.writeList(outputEntries);
            tf.close();

            // create dataset objects
            // run analysis
            System.out.println("Starting meta-analysis");
            ProgressBar pb = new ProgressBar(snpList.length);
            for (int snp = 0; snp < snpList.length; snp++) {
                //Here now. Need to check what probes to meta-analyze, do the meta-analysis and write.
                // do we need to check if alleles are different in the same dataset?
                // get ZScores for this SNP, no matter what.
                // get list of probes to test
                int datasetSNPId = snpIndex[snp];

                if (datasetSNPId != -9) { // -9 means: snp not available
                    float[] datasetZScores = dataset.getZScores(datasetSNPId);
                    THashMap<String, DescriptiveStatistics> remappedEntries = new THashMap<String, DescriptiveStatistics>();
                    THashMap<String, SummaryStatistics> remappedEntriesAbs = new THashMap<String, SummaryStatistics>();
//                    System.out.println(dataset.getSNPs()[datasetSNPId]);
//                    System.out.println(datasetZScores.length);
                    if (dataset.getIsCisDataset()) {
                        // this requires us to retrieve the z-scores differently
                        // we need to figure out which probes match up, but their orders might be different
                        // and the number of probes tested in each dataset might differ as well

                        // get the probes tested against the SNP
                        String[] datasetCisProbes = dataset.getCisProbes(datasetSNPId);
//                        System.out.println(datasetCisProbes.length);
                        for (int i = 0; i < datasetCisProbes.length; i++) {
                            String p = datasetCisProbes[i];
                            if (traitMap.containsKey(p)) {
                                for (String feature : traitMap.get(p)) {
                                    if (!remappedEntries.containsKey(feature)) {
                                        remappedEntries.put(feature, new DescriptiveStatistics());
                                        remappedEntriesAbs.put(feature, new SummaryStatistics());
                                    }
                                    remappedEntries.get(feature).addValue(datasetZScores[i]);
                                    remappedEntriesAbs.get(feature).addValue(Math.abs(datasetZScores[i]));
//                                System.out.print(p+"\t"+datasetZScores[i]);
                                }
                            }
                        }
                        if (!remappedEntries.isEmpty()) {
                            double[] zScoresOut = new double[remappedEntries.size()];
                            String[] keysOut = new String[remappedEntries.size()];
                            int counter = 0;
                            for (Entry<String, DescriptiveStatistics> e : remappedEntries.entrySet()) {
                                keysOut[counter] = e.getKey();
                                if (e.getValue().getN() == 1) {
                                    zScoresOut[counter] = e.getValue().getValues()[0];

                                } else {
                                    //Now we need to meta-analyze it.
                                    if (settings.getzScoreMergeOption().equals("weightedzscore")) {
                                        //Need to get sample size
                                        ArrayIntList sampleSizes = new ArrayIntList();
                                        for (int sC = 0; sC < e.getValue().getN(); sC++) {
                                            sampleSizes.add(dataset.getSampleSize(datasetSNPId));
                                        }
                                        zScoresOut[counter] = ZScores.getWeightedZ(e.getValue().getValues(), sampleSizes.toArray());
                                    } else if (settings.getzScoreMergeOption().equals("mean")) {
                                        zScoresOut[counter] = e.getValue().getMean();
                                    } else if (settings.getzScoreMergeOption().equals("median")) {
                                        zScoresOut[counter] = e.getValue().getPercentile(50);
                                    }  else if (settings.getzScoreMergeOption().equals("min")) {
                                        zScoresOut[counter] = remappedEntriesAbs.get(e.getKey()).getMin();
                                    } else if (settings.getzScoreMergeOption().equals("max")) {
                                        zScoresOut[counter] = remappedEntriesAbs.get(e.getKey()).getMax();
                                    }  else {
                                        System.out.println("Not supported merging.");
                                        System.exit(0);
                                    }

                                }
                                counter++;
                            }
                            //                        writeBinaryResult(String snpname, double hwe, double cr, double maf, int numberCalled, String alleles, String minorAllele, String alleleassessed, double[] datasetZScores, String[] probeNames, BinaryFile outfile, TextFile snpfile)
                            writeBinaryResult(dataset.getSNPs()[datasetSNPId], dataset.getHwes()[datasetSNPId], dataset.getCallrates()[datasetSNPId], dataset.getMafs()[datasetSNPId], dataset.getN()[datasetSNPId], dataset.getAlleles()[datasetSNPId], dataset.getMinorAlleles()[datasetSNPId], dataset.getAlleleAssessed(datasetSNPId), zScoresOut, keysOut, zScoreBinaryFile, zScoreRowNamesFile);
                        }
                    } else { // this is not a cis dataset
                        // use the full probe index

                        // probeIndex[t.getMetaTraitId()][d] = p;
                        for (int i = 0; i < dataset.getProbeList().length; i++) {
                            String p = dataset.getProbeList()[i];
                            if (traitMap.containsKey(p)) {
                                for (String feature : traitMap.get(p)) {
                                    if (!remappedEntries.containsKey(feature)) {
                                        remappedEntries.put(feature, new DescriptiveStatistics());
                                        remappedEntriesAbs.put(feature, new SummaryStatistics());
                                    }
                                    remappedEntries.get(feature).addValue(datasetZScores[i]);
                                    remappedEntriesAbs.get(feature).addValue(datasetZScores[i]);
                                }
                            }
//                            if(dataset.getSNPs()[datasetSNPId].equals("rs2546890")){
//                                System.out.println(p+"\t"+datasetZScores[i]);
//                            }
                        }
                        if (!remappedEntries.isEmpty()) {
                            //This is almost idenitcal to the cis-data. 
                            //But here we need to make sure the remapped data is stored in the correct relative to the file description. i.e. outputEntries
                            
                            //Not observed zScores are kept at 0.
                            double[] zScoresOut = new double[outputEntries.size()];

                            for (Entry<String, DescriptiveStatistics> e : remappedEntries.entrySet()) {

                                int arrayLoc = traitLocationMap.get(e.getKey());
                               
                                if (e.getValue().getN() == 1) {
                                    zScoresOut[arrayLoc] = e.getValue().getValues()[0];

                                } else {
                                    //Now we need to meta-analyze it.
                                    if (settings.getzScoreMergeOption().equals("weightedzscore")) {
                                        //Need to get sample size
                                        ArrayIntList sampleSizes = new ArrayIntList();
                                        for (int sC = 0; sC < e.getValue().getN(); sC++) {
                                            sampleSizes.add(dataset.getSampleSize(datasetSNPId));
                                        }
                                        zScoresOut[arrayLoc] = ZScores.getWeightedZ(e.getValue().getValues(), sampleSizes.toArray());
                                    } else if (settings.getzScoreMergeOption().equals("mean")) {
                                        zScoresOut[arrayLoc] = e.getValue().getMean();
                                    } else if (settings.getzScoreMergeOption().equals("median")) {
                                        zScoresOut[arrayLoc] = e.getValue().getPercentile(50);
                                    }  else if (settings.getzScoreMergeOption().equals("min")) {
                                        zScoresOut[arrayLoc] = remappedEntriesAbs.get(e.getKey()).getMin();
                                    } else if (settings.getzScoreMergeOption().equals("max")) {
                                        zScoresOut[arrayLoc] = remappedEntriesAbs.get(e.getKey()).getMax();
                                    } else {
                                        System.out.println("Not supported merging.");
                                        System.exit(0);
                                    }

                                }
                            }
                            //                        writeBinaryResult(String snpname, double hwe, double cr, double maf, int numberCalled, String alleles, String minorAllele, String alleleassessed, double[] datasetZScores, String[] probeNames, BinaryFile outfile, TextFile snpfile)
                            writeBinaryResult(dataset.getSNPs()[datasetSNPId], dataset.getHwes()[datasetSNPId], dataset.getCallrates()[datasetSNPId], dataset.getMafs()[datasetSNPId], dataset.getN()[datasetSNPId], dataset.getAlleles()[datasetSNPId], dataset.getMinorAlleles()[datasetSNPId], dataset.getAlleleAssessed(datasetSNPId), zScoresOut, null, zScoreBinaryFile, zScoreRowNamesFile);
                        }
                    }
                }
                pb.iterate();
            }
            pb.close();

            //Close binaryOutput & Create binary check file.
            zScoreBinaryFile.close();
            fileName = "check";
            if (runningPermutation) {
                fileName += "-PermutationRound-" + permutation;
            }
            fileName += ".md5";

            HexBinaryAdapter md5Parser = new HexBinaryAdapter();
            BufferedWriter md5writer = new BufferedWriter(new FileWriter(settings.getOutput() + fileName));

            fileName = settings.getOutput() + settings.getDatasetname();
            if (runningPermutation) {
                fileName += "-PermutationRound-" + permutation;
            }
            fileName += ".dat";
            md5writer.write(md5Parser.marshal(zScoreBinaryFile.getWrittenHash()) + "  " + fileName + '\n');
            zScoreRowNamesFile.close();
            md5writer.close();

        }
    }

    private void createSNPIndex(String outdir) throws IOException {

        // create a list of all available SNPs
        HashSet<String> allSNPs = new HashSet<String>();
        String[] snps = dataset.getSNPs();
        for (String snp : snps) {
            allSNPs.add(snp);
        }

        // create a temporary map that maps each SNP to a meta-analysis position
        int ctr = 0;
        TObjectIntHashMap<String> snpMap = new TObjectIntHashMap<String>(allSNPs.size(), 0.85f, -9);
        snpList = new String[allSNPs.size()];
        for (String s : allSNPs) {
            snpMap.put(s, ctr);
            snpList[ctr] = s;
            ctr++;
        }

        // fill index
        snpIndex = new int[allSNPs.size()];
        for (int s = 0; s < allSNPs.size(); s++) {
            snpIndex[s] = -9;
        }

        for (int s = 0; s < snps.length; s++) {
            String snp = snps[s];
            int id = snpMap.get(snp);
            if (id != -9) {
                snpIndex[id] = s;
            }
        }

//        TextFile tf = new TextFile(outdir + "snpindex.txt", TextFile.W);
//        String header = "metaID";
//        header += "\t" + dataset.getName() + "-sid";
//        tf.writeln(header);
//
//        for (int s = 0; s < snpList.length; s++) {
//            String ln = snpList[s];
//            ln += "\t" + snpIndex[s];
//            tf.writeln(ln);
//        }
//        tf.close();
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
    private void writeBinaryResult(String snpname, double hwe, double cr, double maf, int numberCalled, String alleles, String minorAllele, String alleleassessed, double[] datasetZScores, String[] probeNames, BinaryFile outfile, TextFile snpfile) throws IOException {
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
