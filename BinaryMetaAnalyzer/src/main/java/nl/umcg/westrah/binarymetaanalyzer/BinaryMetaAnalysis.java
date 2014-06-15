/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author Harm-Jan
 */
public class BinaryMetaAnalysis {

    public static void main(String[] args) {

        String settingsFile = args[0];
        String textToReplace = args[1];
        String replaceTextWith = args[2];

        BinaryMetaAnalysis meta = new BinaryMetaAnalysis(settingsFile, textToReplace, replaceTextWith);

    }
    private MetaQTL4TraitAnnotation probeAnnotation;

    private BinaryMetaAnalysisDataset[] datasets = new BinaryMetaAnalysisDataset[0];
    private int[][] snpIndex;
    private String[] snpList;
    private final BinaryMetaAnalysisSettings settings;
    private String[] snpChr;
    private int[] snpPositions;
    private Integer[][] probeIndex;

    private QTL[] finalEQTLs;
    private boolean bufferHasOverFlown;
    private double maxSavedPvalue = -Double.MAX_VALUE;
    private boolean sorted;
    private int locationToStoreResult;

    public BinaryMetaAnalysis(String settingsFile, String textToReplace, String replaceTextWith) {
        // initialize settings
        settings = new BinaryMetaAnalysisSettings();
        settings.parse(settingsFile, textToReplace, replaceTextWith);
        int maxResults = settings.getFinalEQTLBufferMaxLength();
        int tmpbuffersize = (maxResults / 10);

        if (tmpbuffersize == 0) {
            tmpbuffersize = 10;
        } else if (tmpbuffersize > 250000) {
            tmpbuffersize = 250000;
        }

        finalEQTLs = new QTL[(maxResults + tmpbuffersize)];
        try {
            run();
        } catch (IOException ex) {
            Logger.getLogger(BinaryMetaAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void run() throws IOException {

        // load probe annotation and index
        // this particular probe annotation can take multiple probes for a single location into account.
        loadProbeAnnotation();

        for (int permutation = 0; permutation < settings.getNrPermutations(); permutation++) {
            // create dataset objects
            datasets = new BinaryMetaAnalysisDataset[settings.getDatasetlocations().size()];

            for (int d = 0; d < datasets.length; d++) {
                datasets[d] = new BinaryMetaAnalysisDataset(settings.getDatasetlocations().get(d), permutation, settings.getDatasetannotations().get(d), probeAnnotation);
            }

            // create meta-analysis SNP index. have to recreate this every permutation, 
            // since the order of SNPs is generated at random.
            createSNPIndex();
            createProbeIndex();

            if (snpChr == null) {
                loadSNPAnnotation();
            }

            // run analysis
            for (int snp = 0; snp < snpList.length; snp++) {

                int[] sampleSizes = new int[datasets.length];
                Boolean[] flipZScores = new Boolean[datasets.length];
                String alleles = null;
                String alleleAssessed = null;
                // get list of probes to test for each dataset
                for (int d = 0; d < datasets.length; d++) {

                    int datasetSNPId = snpIndex[snp][d];
                    if (datasetSNPId != -9) {
                        sampleSizes[d] = datasets[d].getSampleSize(datasetSNPId);

                        if (alleles == null) {
                            flipZScores[d] = false;
                            alleles = datasets[d].getAlleles(datasetSNPId);
                            alleleAssessed = datasets[d].getAlleleAssessed(datasetSNPId);
                        } else {
                            String alleles2 = datasets[d].getAlleles(datasetSNPId);
                            String alleleAssessed2 = datasets[d].getAlleleAssessed(datasetSNPId);
                            flipZScores[d] = BaseAnnot.flipalleles(alleles, alleleAssessed, alleles2, alleleAssessed2);
                        }
                    }
                }

                // get ZScores for this SNP, no matter what.
                // get list of probes to test
                if (settings.isCis() && !settings.isTrans()) {
                    // do cis stuff

                    // get all the traits near the SNP
                    Set<MetaQTL4MetaTrait> cisProbes = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
                    MetaQTL4MetaTrait[] cisProbeArray = cisProbes.toArray(new MetaQTL4MetaTrait[0]);
                    HashMap<MetaQTL4MetaTrait, Integer> cisProbeMap = new HashMap<MetaQTL4MetaTrait, Integer>();
                    int ctr = 0;
                    for (MetaQTL4MetaTrait cisProbe : cisProbes) {
                        cisProbeMap.put(cisProbe, ctr);
                        ctr++;
                    }

                    double[][] zScores = new double[cisProbeMap.size()][datasets.length];

                    // get list of probes to test for each dataset
                    for (int d = 0; d < datasets.length; d++) {
                        for (int p = 0; p < cisProbeMap.size(); p++) {
                            zScores[p][d] = Double.NaN;
                        }
                        int datasetSNPId = snpIndex[snp][d];
                        float[] z = datasets[d].getZScores(datasetSNPId);
                        if (datasetSNPId != -9) {
                            if (datasets[d].getIsCisDataset()) {
// this requires us to retrieve the z-scores differently
                                // we need to figure out which probes match up, but their orders might be different
                                // and the number of probes tested in each dataset might differ as well
                                // get the probes tested against the SNP
                                MetaQTL4MetaTrait[] probes = datasets[d].getCisProbes(datasetSNPId);

                                for (int i = 0; i < probes.length; i++) {
                                    MetaQTL4MetaTrait p = probes[i];
                                    if (p != null) {
                                        Integer index = cisProbeMap.get(probes[i]);
                                        if (index != null) {
                                            zScores[index][d] = z[i];
                                            if (flipZScores[d]) {
                                                zScores[index][d] *= -1;
                                            }
                                        }
                                    }
                                }

                            } else {
                                // use the full probe index    
                                int index = 0;

                                for (MetaQTL4MetaTrait t : cisProbes) {
                                    Integer metaId = t.getMetaTraitId();
                                    Integer probeId = probeIndex[metaId][d];
                                    if (probeId != null) {
                                        zScores[index][d] = z[probeId];
                                        if (flipZScores[d]) {
                                            zScores[index][d] *= -1;
                                        }
                                    }
                                    index++;
                                }
                            }
                        }

                    }

                    // meta-analyze!
                    for (int probe = 0; probe < zScores.length; probe++) {
                        MetaQTL4MetaTrait t = cisProbeArray[probe];

                        double metaZ = ZScores.getWeightedZ(zScores[probe], sampleSizes);
                        double p = ZScores.zToP(metaZ);

                        // create output object
                        QTL q = new QTL(p, t.getMetaTraitId(), snp, BaseAnnot.toByte(alleleAssessed), metaZ, BaseAnnot.toByteArray(alleles), zScores[probe], sampleSizes); // sort buffer if needed.
                        addEQTL(q);
                    }
                } else {
                    Set<MetaQTL4MetaTrait> cisProbes = null;
                    if (!settings.isCis()) {
                        // do not test the cis probes
                        cisProbes = probeAnnotation.getMetatraits().getTraitInWindow(snpChr[snp], snpPositions[snp], settings.getCisdistance());
                    }
                    // iterate over the probe index

                    double[][] zScores = new double[probeIndex.length][datasets.length];
                    Set<MetaQTL4MetaTrait> traits = probeAnnotation.getMetatraits();
                    for (int d = 0; d < datasets.length; d++) {
                        if(datasets[d].getIsCisDataset()){
                            System.err.println("ERROR: cannot run trans analysis on a cis dataset: "+settings.getDatasetlocations().get(d));
                            System.exit(-1);
                        }
                        int datasetSNPId = snpIndex[snp][d];
                        float[] z = datasets[d].getZScores(datasetSNPId);
                        int p = 0;

                        for (MetaQTL4MetaTrait t : traits) {
                            if (cisProbes != null && cisProbes.contains(t)) {
                                zScores[p][d] = Double.NaN;
                            } else {
                                Integer probeId = probeIndex[p][d];
                                if (probeId != null) {
                                    zScores[p][d] = z[probeId];
                                    if (flipZScores[d]) {
                                        zScores[p][d] *= -1;
                                    }
                                } else {
                                    zScores[p][d] = Double.NaN;
                                }
                                p++;
                            }
                        }
                    }

                    // meta-analyze!
                    int probe = 0;
                    for (MetaQTL4MetaTrait t : traits) {

                        double metaZ = ZScores.getWeightedZ(zScores[probe], sampleSizes);
                        double pval = ZScores.zToP(metaZ);

                        // create output object
                        QTL q = new QTL(pval, t.getMetaTraitId(), snp, BaseAnnot.toByte(alleleAssessed), metaZ, BaseAnnot.toByteArray(alleles), zScores[probe], sampleSizes); // sort buffer if needed.
                        addEQTL(q);
                        probe++;
                    }
                }
            }
            
            // write the results to disk    
        }
        
        /*
        TODO:
        - Plotting of z-scores
        - writing of output file
        - validation
        - multithreadalize
        */

        
    }

    private void createSNPIndex() throws IOException {

        HashSet<String> confineToTheseSNPs = null;
        if (settings.getSNPSelection() != null) {
            confineToTheseSNPs = new HashSet<String>();
            TextFile tf = new TextFile(settings.getSNPSelection(), TextFile.R);
            confineToTheseSNPs.addAll(tf.readAsArrayList());
            tf.close();
        }

        // create a list of all available SNPs
        HashSet<String> allSNPs = new HashSet<String>();
        for (BinaryMetaAnalysisDataset dataset : datasets) {
            String[] snps = dataset.getSNPs();
            for (String snp : snps) {
                if (confineToTheseSNPs == null || confineToTheseSNPs.contains(snp)) {
                    allSNPs.add(snp);
                }
            }

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
        snpIndex = new int[allSNPs.size()][datasets.length];
        for (int d = 0; d < datasets.length; d++) {
            for (int s = 0; s < allSNPs.size(); s++) {
                snpIndex[s][d] = -9;
            }
        }
        for (int d = 0; d < datasets.length; d++) {
            String[] snps = datasets[d].getSNPs();
            for (int s = 0; s < snps.length; s++) {
                String snp = snps[s];
                int id = snpMap.get(snp);
                if (id != -9) {
                    snpIndex[id][d] = s;
                }
            }
        }
    }

    private void loadProbeAnnotation() throws IOException {

        HashSet<String> platforms = new HashSet<String>();
        platforms.addAll(settings.getDatasetannotations());
        probeAnnotation = new MetaQTL4TraitAnnotation(new File(settings.getProbetranslationfile()), platforms);

    }

    private void loadSNPAnnotation() throws IOException {

        snpChr = new String[snpList.length];
        snpPositions = new int[snpList.length];
        TObjectIntHashMap<String> snpMap = new TObjectIntHashMap<String>(snpList.length);
        for (int s = 0; s < snpList.length; s++) {
            snpMap.put(snpList[s], s);
        }
        TextFile tf = new TextFile(settings.getSNPAnnotationFile(), TextFile.R);

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[2];
            if (snpMap.contains(snp)) {
                int id = snpMap.get(snp);
                snpChr[snpMap.get(snp)] = new String(elems[0].getBytes("UTF-8")).intern();
                snpPositions[id] = Integer.parseInt(elems[1]);
            }
        }
        tf.close();

    }

    // index the probes 
    private void createProbeIndex() {
        probeIndex = new Integer[probeAnnotation.getMetatraits().size()][datasets.length];
        for (int d = 0; d < datasets.length; d++) {
            String[] probes = datasets[d].getProbeList();
            int platformId = probeAnnotation.getPlatformId(settings.getDatasetannotations().get(d));
            for (int p = 0; p < probes.length; p++) {
                MetaQTL4MetaTrait t = probeAnnotation.getTraitForPlatformId(platformId, probes[p]);
                probeIndex[t.getMetaTraitId()][d] = p;
            }
        }
    }

    private void addEQTL(QTL q) {

        double pval = q.getPvalue();
        if (bufferHasOverFlown) {
            if (pval <= maxSavedPvalue) {

                sorted = false;

                finalEQTLs[locationToStoreResult] = q;
                locationToStoreResult++;

                if (locationToStoreResult == finalEQTLs.length) {

                    Arrays.sort(finalEQTLs);
//                    SmoothSort.sort(finalEQTLs);
//                    inplaceArrayQuickSort.sort(finalEQTLs);
                    sorted = true;
                    locationToStoreResult = settings.getFinalEQTLBufferMaxLength();
                    maxSavedPvalue = finalEQTLs[(settings.getFinalEQTLBufferMaxLength() - 1)].getPvalue();
                }
            }

        } else {
            if (pval > maxSavedPvalue) {
                maxSavedPvalue = pval;
            }

            finalEQTLs[locationToStoreResult] = q;
            locationToStoreResult++;

            if (locationToStoreResult == settings.getFinalEQTLBufferMaxLength()) {
                bufferHasOverFlown = true;
            }
        }
    }
}
