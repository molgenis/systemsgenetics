/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedWriter;
import java.io.FileWriter;
import umcg.genetica.io.text.TextFile;
import javax.xml.bind.annotation.adapters.HexBinaryAdapter;
import umcg.genetica.io.bin.BinaryFile;

import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.console.ProgressBar;

/**
 * @author Harm-Jan
 */
public class InternalMetaAnalysis {

    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
    private int[] snpIndex;
    private Integer[] probeIndex;
    private String[] snpList;
    private InternalMetaAnalysisDataset dataset;
    private final InternalMetaAnalysisSettings settings;

    HashMap<String, String> traitMap = new HashMap<String, String>();
    LinkedHashMap<String, Integer> outputEntries = null;

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

    }

    public void run() throws IOException {

        String outdir = settings.getOutput();

        outputEntries = new LinkedHashMap<String, Integer>();
        traitMap = new HashMap<String, String>();

        TextFile entryReader = new TextFile(settings.getProbeToFeature(), TextFile.R);
        String row;
        int index = 0;
        while ((row = entryReader.readLine()) != null) {
//            System.out.println(row);
            String[] parts = TAB_PATTERN.split(row);

            if (!outputEntries.containsKey(parts[1])) {
                outputEntries.put(parts[1], index);
                index++;
            }

            traitMap.put(parts[0], parts[1]);
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
            System.out.println("Creating probe index");
            createProbeIndex(outdir);
            System.out.println("Total of " + probeIndex.length + " probes");

            // write magic number
            if (dataset.getIsCisDataset()) {
                zScoreBinaryFile.writeInt(1);
            } else {
                zScoreBinaryFile.writeInt(0);
            }

            TextFile zScoreRowNamesFile = new TextFile(fileName + "-RowNames.txt.gz", TextFile.W);
            zScoreRowNamesFile.writeln("SNP\tAlleles\tMinorAllele\tAlleleAssessed\tNrCalled\tMaf\tHWE\tCallRate");

            TextFile tf = new TextFile(fileName + "-ColNames.txt.gz", TextFile.W);
            tf.writeList(new ArrayList<String>(outputEntries.keySet()));
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
                // do cis stuff
                // get all the possible traits for SNP
                HashMap<String, Integer> probeMap = new HashMap<String, Integer>();
                float[] finalZScores = new float[probeMap.size()]; // size: [possible cis-probes][nr of datasets]

                // get list of probes to test for each dataset
                //initialize z-score
                for (int p = 0; p < probeMap.size(); p++) {
                    finalZScores[p] = Float.NaN; // this is not very nice, but does prevent the metaZ method from going nuts
                }
                // load the z-scores for the dataset
                int datasetSNPId = snpIndex[snp];

                if (datasetSNPId != -9) { // -9 means: snp not available
                    float[] datasetZScores = dataset.getZScores(datasetSNPId);
                    
                    System.out.println(datasetZScores.length);
                    if (dataset.getIsCisDataset()) {
                        // this requires us to retrieve the z-scores differently
                        // we need to figure out which probes match up, but their orders might be different
                        // and the number of probes tested in each dataset might differ as well

                        // get the probes tested against the SNP
                        MetaQTL4MetaTrait[] datasetCisProbes = dataset.getCisProbes(datasetSNPId);
                        System.out.println(datasetCisProbes.length);
                        
//                        for (int i = 0; i < datasetCisProbes.length; i++) {
//                            MetaQTL4MetaTrait p = datasetCisProbes[i];
//                            if (p != null) {
//                                Integer index = cisProbeMap.get(p);
//                                if (index != null) {
//                                    float datasetZ = datasetZScores[i];
//                                    finalZScores[index] = datasetZ;
//                                }
//                            }
//                        }

                    } else { // this is not a cis dataset
                        // use the full probe index
//                        for (int probe = 0; probe < cisProbeArray.length; probe++) {
//                            MetaQTL4MetaTrait cisProbe = cisProbeArray[probe];
//                            Integer metaProbeIndex = traitMap.get(cisProbe);
//                            Integer datasetProbeId = probeIndex[metaProbeIndex][d];
//                            if (datasetProbeId != null) {
//                                finalZScores[probe][d] = datasetZScores[datasetProbeId];
//                                if (flipZScores[d]) {
//                                    finalZScores[probe][d] *= -1;
//                                }
//                            }
//                        }
                    }
                }

                //Writing from original result processing
                //Get data from a thread and write R.
                //writeBinaryResult(r);
                //meta-analyze and write output
//                for (int probe = 0; probe < finalZScores.length; probe++) {
//                    MetaQTL4MetaTrait t = cisProbeArray[probe];
//
//                    double metaZ = ZScores.getWeightedZ(finalZScores[probe], sampleSizes);
//                    double p = Descriptives.convertZscoreToPvalue(metaZ);
//
//                    if (!Double.isNaN(p) && !Double.isNaN(metaZ)) {
//                        // create output object
//                        QTL q = new QTL(p, t, snp, BaseAnnot.toByte(alleleAssessed), metaZ, BaseAnnot.toByteArray(alleles), finalZScores[probe], sampleSizes); // sort buffer if needed.
////                            System.out.println(q.getSNPId()+"\t"+q.getMetaTrait().getMetaTraitName()+"\t"+q.toString());
//                        addEQTL(q);
//                    }
//                }
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

        TextFile tf = new TextFile(outdir + "snpindex.txt", TextFile.W);
        String header = "metaID";
        header += "\t" + dataset.getName() + "-sid";
        tf.writeln(header);

        for (int s = 0; s < snpList.length; s++) {
            String ln = snpList[s];
            ln += "\t" + snpIndex[s];
            tf.writeln(ln);
        }
        tf.close();
    }

    // index the probes
    private void createProbeIndex(String outdir) throws IOException {

        HashSet<String> confineToTheseProbes = new HashSet<String>(traitMap.keySet());

        probeIndex = new Integer[dataset.getProbeList().length];

        String[] probes = dataset.getProbeList();

        for (int p = 0; p < probes.length; p++) {

            String t = traitMap.get(probes[p]);
//            System.out.println(probes[p]);
//            System.out.println(t);
            Integer index = outputEntries.get(t);

            if (index != null && confineToTheseProbes.contains(probes[p])) {
                probeIndex[p] = index;
            } else {
                probeIndex[p] = null;
            }
        }

    }

//
//    private void writeBinaryResult(Result r) throws IOException {
//
//        if (r != null) {
//            int[] numSamples = null;
//            try {
//                numSamples = r.numSamples;
//            } catch (NullPointerException e) {
//                System.out.println("ERROR: null result?");
//            }
//
//            int wpId = r.wpid;
//            WorkPackage currentWP = m_availableWorkPackages[wpId];
//            double[][] zscores = r.zscores;
//
//            if (zscores != null) {
//                SNP[] snps = currentWP.getSnps();
//                int numDatasets = zscores.length;
//                double[] finalZscores = r.finalZScore;
//                StringBuilder snpoutput = null;
//
//                // if we're doing a meta-analysis, write the meta-analysis Z to a separate binaryFile
//                if (m_gg.length > 1) {
//                    int totalSampleNr = 0;
//                    String snpname = null;
//                    for (int d = 0; d < numDatasets; d++) {
//                        if (snps[d] != null) {
//                            snpname = snps[d].getName();
//
//                            byte[] alleles = snps[d].getAlleles();
//                            byte minorAllele = snps[d].getMinorAllele();
//                            byte alleleassessed = alleles[1];
//
//                            if (currentWP.getFlipSNPAlleles()[d]) {
//                                alleleassessed = alleles[0];
//                            }
//                            if (snpoutput == null) {
//                                snpoutput = new StringBuilder();
//                                snpoutput.append(snpname);
//                                snpoutput.append("\t");
//                                snpoutput.append(BaseAnnot.getAllelesDescription(alleles));
//                                snpoutput.append("\t");
//                                snpoutput.append(BaseAnnot.toString(minorAllele));
//                                snpoutput.append("\t");
//                                snpoutput.append(BaseAnnot.toString(alleleassessed));
//                            }
//                            totalSampleNr += r.numSamples[d];
//                        }
//                    }
//
//                    StringBuilder sb = null;
//                    for (int p = 0; p < finalZscores.length; p++) {
//                        float z = (float) finalZscores[p];
//                        if (m_cisOnly) {
//                            int[] probes = currentWP.getProbes();
//                            int probeId = probes[p];
//                            String probeName = m_probeList[probeId];
//                            if (sb == null) {
//                                sb = new StringBuilder();
//                            } else {
//                                sb.append("\t");
//                            }
//                            sb.append(probeName);
//
//                            zScoreMetaAnalysisFile.writeFloat(z);
//                        } else {
//                            zScoreMetaAnalysisFile.writeFloat(z);
//                        }
//                    }
//
//                    if (snpoutput != null) {
//                        snpoutput.append("\t");
//                        snpoutput.append(totalSampleNr);
//                        snpoutput.append("\t-\t-\t-\t");
//                        snpoutput.append(finalZscores.length);
//                        snpoutput.append("\t");
//                        if (sb != null) {
//                            snpoutput.append(sb.toString());
//                        } else {
//                            snpoutput.append("-");
//                        }
//                        zScoreMetaAnalysisRowNamesFile.writeln(snpoutput.toString());
//                    }
//                }
//
//                for (int d = 0; d < numDatasets; d++) {
//                    double[] datasetZScores = zscores[d];
//                    SNP datasetSNP = snps[d];
//                    if (datasetSNP != null) {
//                        BinaryFile outfile = zScoreBinaryFile[d];
//
//                        String snpname = datasetSNP.getName();
//
//                        byte[] alleles = datasetSNP.getAlleles();
//                        byte minorAllele = datasetSNP.getMinorAllele();
//                        byte alleleassessed = alleles[1];
//                        double hwe = datasetSNP.getHWEP();
//                        double cr = datasetSNP.getCR();
//                        double maf = datasetSNP.getMAF();
//
//                        if (currentWP.getFlipSNPAlleles()[d]) {
//                            alleleassessed = alleles[0];
//                        }
//                        TextFile snpfile = zScoreRowNamesFile[d];
//                        StringBuilder sb = null;
//                        for (int p = 0; p < datasetZScores.length; p++) {
//                            float z = (float) datasetZScores[p];
//                            if (currentWP.getFlipSNPAlleles()[d]) {
//                                z *= -1;
//                            }
//                            // System.out.println(p + "\t" + alleleassessed + "\t" + m_probeList[p] + "\t" + z + "\t" + currentWP.getFlipSNPAlleles()[d]);
//                            if (m_cisOnly) {
//                                // take into account that not all probes have been tested..
//                                int[] probes = currentWP.getProbes();
//                                int probeId = probes[p];
//                                String probeName = m_probeList[probeId];
//                                outfile.writeFloat(z);
//                                if (sb == null) {
//                                    sb = new StringBuilder();
//                                } else {
//                                    sb.append("\t");
//                                }
//                                sb.append(probeName);
//                            } else {
//                                outfile.writeFloat(z);
//                            }
//                        }
//
//                        StringBuilder buffer = new StringBuilder();
//                        buffer.append(snpname)
//                                .append("\t")
//                                .append(BaseAnnot.getAllelesDescription(alleles))
//                                .append("\t")
//                                .append(BaseAnnot.toString(minorAllele))
//                                .append("\t")
//                                .append(BaseAnnot.toString(alleleassessed))
//                                .append("\t")
//                                .append(datasetSNP.getNrCalled())
//                                .append("\t")
//                                .append(maf)
//                                .append("\t")
//                                .append(hwe)
//                                .append("\t")
//                                .append(cr)
//                                .append("\t")
//                                .append(datasetZScores.length)
//                                .append("\t");
//                        if (sb != null) {
//                            buffer.append(sb.toString());
//                        } else {
//                            buffer.append("-");
//                        }
//
//                        snpfile.writeln(buffer.toString());
//
//                    }
//                }
//            }
//        }
//    }
}
