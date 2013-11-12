/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 * MetaQTL4 - Experimental Skunkworx edition
 *
 * @author harmjan
 */
public class MetaQTL4 {

    private Set<MetaQTL4MetaTrait> traitsToInclude = null;
    private Set<String> variantsToInclude = null;
    private Set<Pair<String, MetaQTL4MetaTrait>> genotypeTraitCombinations = null;
    private final MetaQTL4Settings m_settings;
    private MetaQTL4TraitAnnotation traitAnnotation;
    private MetaQTL4Dataset[] datasets;
    private Integer[][] traitIndex;
    private ArrayList<MetaQTL4MetaTrait> availableTraits;
    private HashMap<MetaQTL4MetaTrait, Integer> availableTraitsHash;
    private GeneticVariant[][] geneticVariantIndex;

    public MetaQTL4(MetaQTL4Settings settings) throws IOException, Exception {
        m_settings = settings;
        initialize();
        run();
    }

    public MetaQTL4(String settings, String replaceText, String testToReplaceWith) throws IOException, Exception {
        m_settings = new MetaQTL4Settings(settings, replaceText, testToReplaceWith);
        initialize();
        run();
    }

    public final void initialize() throws IOException {

        // load the probe annotation (only for the platforms that are loaded)..
        HashSet<String> platforms = new HashSet<String>();
        for (int d = 0; d < m_settings.getNrDatasets(); d++) {
            platforms.add(m_settings.getDatasetSettings().get(d).getTraitPlatform());
        }
        traitAnnotation = new MetaQTL4TraitAnnotation(m_settings.getProbeAnnotationFile(), platforms);

        if (!Gpio.exists(m_settings.getOutputReportsDir())) {
            Gpio.createDir(m_settings.getOutputReportsDir());
        }

        if (m_settings.getStrConfineProbe() != null && Gpio.exists(m_settings.getStrConfineProbe())) {
            traitsToInclude = loadMetaTraitList(m_settings.getStrConfineProbe());
        }

        if (m_settings.getStrConfineSNP() != null && Gpio.exists(m_settings.getStrConfineSNP())) {
            TextFile tf = new TextFile(m_settings.getStrConfineSNP(), TextFile.R);
            variantsToInclude = tf.readAsSet(0, TextFile.tab);
            tf.close();
        }

        // load list of specific SNP-trait combinations
        if (m_settings.getStrSNPProbeConfine() != null && Gpio.exists(m_settings.getStrSNPProbeConfine())) {
            Triple<HashSet<String>, HashSet<MetaQTL4MetaTrait>, HashSet<Pair<String, MetaQTL4MetaTrait>>> combinations = loadgenotypeMarkerCombinations(m_settings.getStrSNPProbeConfine());
            variantsToInclude = combinations.getLeft();
            traitsToInclude = combinations.getMiddle();
            genotypeTraitCombinations = combinations.getRight();
        }

        // PROBE FILTER...
        // TODO: there is one exception for the probe filter, 
        // which is when all the datasets are on the same platform. 
        // will implement this at a later stage, or not at all.
        System.out.println("Running probe filter");
        runProbeFilter();

        // load the datasets..
        datasets = new MetaQTL4Dataset[m_settings.getNrDatasets()];
        if (datasets.length == 0) {
            throw new IllegalStateException("Error: no dataset information provided.");
        }
        for (int d = 0; d < datasets.length; d++) {
            HashSet<String> platformSpecificProbesToInclude = null;
            if (traitsToInclude != null) {
                String platform = m_settings.getDatasetSettings().get(d).getTraitPlatform();
                Integer platformId = traitAnnotation.getPlatformId(platform);
                platformSpecificProbesToInclude = new HashSet<String>();
                for (MetaQTL4MetaTrait metaprobe : traitsToInclude) {
                    String platformProbe = metaprobe.getPlatformIds()[platformId];
                    if (platformProbe != null) {
                        platformSpecificProbesToInclude.add(platformProbe);
                    }
                }
            }

            datasets[d] = new MetaQTL4Dataset(m_settings.getDatasetSettings().get(d), platformSpecificProbesToInclude);
            if (m_settings.isPerformParametricAnalysis()) {
                datasets[d].rankTraitData(m_settings.isEqualRankForTies());
            }
        }

        // create probe index
        System.out.println("Creating trait index");
        createTraitIndex();

        // create SNP index
        System.out.println("Create SNP index");
        createSNPIndex();
    }

    public final void run() {
        // create threadpool
        System.out.println("Running software!");
        int nrPermutations = m_settings.getNrPermutationsFDR();

        // initialize random seed array
        long[] randomizationSeeds = new long[nrPermutations];
        Random rand = new Random();
        for (int permutation = 0; permutation < nrPermutations; permutation++) {
            randomizationSeeds[permutation] = rand.nextLong();
        }

        // create frequency distributions
        int distributionSize = 1000000;
        int[] realFrequencyDistribution = new int[distributionSize];
        int[] realFrequencyDistributionProbeLevel = new int[distributionSize];
        int[] permutedFrequencyDistribution = new int[distributionSize];
        int[] permutedFrequencyDistributionProbeLevel = new int[distributionSize];

        // run threads
        int nrThreads = m_settings.getNrThreads();
        ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
        CompletionService<Pair<int[], int[]>> pool = new ExecutorCompletionService<Pair<int[], int[]>>(threadPool);
        MetaQTL4CorrelationTask task = new MetaQTL4CorrelationTask(nrThreads, distributionSize, randomizationSeeds, availableTraits, availableTraitsHash, datasets, geneticVariantIndex, m_settings, traitAnnotation, traitIndex, traitsToInclude, variantsToInclude, 1);
        pool.submit(task);

        int returned = 0;
        while (returned < nrThreads) {
            try {
                Pair<int[], int[]> result = pool.take().get();
                if (result != null) {
                    int[] realDist = result.getLeft();
                    int[] permDist = result.getRight();
                    for (int i = 0; i < distributionSize; i++) {
                        realFrequencyDistribution[i] += realDist[i];
                        permutedFrequencyDistribution[i] += permDist[i];
                    }
                    returned++;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        // close threadpool
        threadPool.shutdown();
        
        // perform FDR estimation
        

        // shutdown
    }

    private Set<MetaQTL4MetaTrait> loadMetaTraitList(String location) throws IOException {
        TextFile tf = new TextFile(location, TextFile.R);
        Set<String> list = tf.readAsSet(0, TextFile.tab);

        Set<MetaQTL4MetaTrait> traits = new HashSet<MetaQTL4MetaTrait>();
        for (String s : list) {
            MetaQTL4MetaTrait t = traitAnnotation.getMetaTraitNameToObj().get(s);
            if (t != null) {
                traits.add(t);
            }
        }

        tf.close();
        return traits;
    }

    private Triple<HashSet<String>, HashSet<MetaQTL4MetaTrait>, HashSet<Pair<String, MetaQTL4MetaTrait>>> loadgenotypeMarkerCombinations(String location) throws IOException {
        HashSet<Pair<String, MetaQTL4MetaTrait>> pairs = new HashSet<Pair<String, MetaQTL4MetaTrait>>();
        HashSet<String> snps = new HashSet<String>();
        HashSet<MetaQTL4MetaTrait> probes = new HashSet<MetaQTL4MetaTrait>();
        TextFile tf = new TextFile(location, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[0];
            String probe = elems[1];

            MetaQTL4MetaTrait t = traitAnnotation.getMetaTraitNameToObj().get(probe);
            if (t != null) {
                Pair<String, MetaQTL4MetaTrait> p = new Pair<String, MetaQTL4MetaTrait>(snp, t);
                pairs.add(p);
                snps.add(snp);
                probes.add(t);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        return new Triple<HashSet<String>, HashSet<MetaQTL4MetaTrait>, HashSet<Pair<String, MetaQTL4MetaTrait>>>(snps, probes, pairs);
    }

    private void runProbeFilter() {
        // if we confine to a single chromosome, throw out all the other probes...
        if (m_settings.getConfineToProbesThatMapToChromosome() != null) {
            HashSet<MetaQTL4MetaTrait> finalTraitList = new HashSet<MetaQTL4MetaTrait>();
            if (traitsToInclude != null) {
                for (MetaQTL4MetaTrait metaprobe : traitsToInclude) {
                    String chr = metaprobe.getChr();
                    String chrToInclude = "" + m_settings.getConfineToProbesThatMapToChromosome();
                    if (chr.equals(chrToInclude) && (!m_settings.isExpressionDataLoadOnlyProbesThatMapToChromosome() || (m_settings.isExpressionDataLoadOnlyProbesThatMapToChromosome() && MetaQTL4MetaTrait.mapsToChromosome()))) {
                        finalTraitList.add(metaprobe);
                    }
                }
            } else {
                MetaQTL4MetaTrait[] metaProbes = traitAnnotation.getMetatraits().toArray(new MetaQTL4MetaTrait[0]);
                for (MetaQTL4MetaTrait metaprobe : metaProbes) {
                    String chr = metaprobe.getChr();
                    String chrToInclude = "" + m_settings.getConfineToProbesThatMapToChromosome();
                    if (chr.equals(chrToInclude) && (!m_settings.isExpressionDataLoadOnlyProbesThatMapToChromosome() || (m_settings.isExpressionDataLoadOnlyProbesThatMapToChromosome() && MetaQTL4MetaTrait.mapsToChromosome()))) {
                        finalTraitList.add(metaprobe);
                    }
                }
            }
            traitsToInclude = finalTraitList;

            if (genotypeTraitCombinations != null) {
                // find out whether there are some snp-trait combo's that we cannot test now..
                HashSet<Pair<String, MetaQTL4MetaTrait>> finalgenotypetraitpairs = new HashSet<Pair<String, MetaQTL4MetaTrait>>();
                for (Pair<String, MetaQTL4MetaTrait> p : genotypeTraitCombinations) {
                    MetaQTL4MetaTrait probe = p.getRight();
                    if (traitsToInclude.contains(probe)) {
                        finalgenotypetraitpairs.add(p);
                    }
                }
                genotypeTraitCombinations = finalgenotypetraitpairs;
            }
        }

        if (traitsToInclude != null && traitsToInclude.isEmpty()) {
            throw new IllegalStateException("No traits remain after filtering.");
        }

        if (genotypeTraitCombinations != null && genotypeTraitCombinations.isEmpty()) {
            throw new IllegalStateException("No snp-trait combinations remain after filtering.");
        }
    }

    // create a map from metaProbeId to 
    private void createTraitIndex() {
        // link them together in an index
        HashSet<MetaQTL4MetaTrait> tmpAvailableTraits = new HashSet<MetaQTL4MetaTrait>();
        for (int d = 0; d < datasets.length; d++) {
            String platform = m_settings.getDatasetSettings().get(d).getTraitPlatform();
            Integer platformId = traitAnnotation.getPlatformId(platform);
            String[] availableTraitsInDataset = datasets[d].getTraits();
            for (String trait : availableTraitsInDataset) {
                MetaQTL4MetaTrait metatrait = traitAnnotation.getTraitForPlatformId(platformId, trait);
                if (metatrait != null) {
                    availableTraits.add(metatrait);
                }
            }
        }

        availableTraits = new ArrayList<MetaQTL4MetaTrait>(tmpAvailableTraits);
        availableTraitsHash = new HashMap<MetaQTL4MetaTrait, Integer>();
        traitIndex = new Integer[datasets.length][availableTraits.size()];

        for (int ds = 0; ds < datasets.length; ds++) {
            Integer platformId = traitAnnotation.getPlatformId(m_settings.getDatasetSettings().get(ds).getTraitPlatform());
            for (int i = 0; i < availableTraits.size(); i++) {
                MetaQTL4MetaTrait trait = availableTraits.get(i);
                availableTraitsHash.put(trait, i);
                String probeName = trait.getPlatformIds()[platformId];
                if (probeName != null) {
                    traitIndex[ds][i] = datasets[ds].getTraitId(probeName);
                }
            }
        }

    }

    private void createSNPIndex() {
        // create a list of shared snps
        ArrayList<HashMap<String, GeneticVariant>> tmpVariantsArr = new ArrayList<HashMap<String, GeneticVariant>>();
        HashMap<String, Integer> datasetCounter = new HashMap<String, Integer>();
        HashSet<String> allAvailableVariants = new HashSet<String>();
        // here we'll also filter for snps not mapping to chromosomes etc...
        int maxNrVariants = 0;
        System.out.println("Starting to index.");
        for (int d = 0; d < datasets.length; d++) {
            RandomAccessGenotypeData data = datasets[d].getGenotypeData();
            HashMap<String, GeneticVariant> tmpVariants = new HashMap<String, GeneticVariant>();
            tmpVariantsArr.add(tmpVariants);
            for (GeneticVariant variant : data) {
                String primaryId = variant.getPrimaryVariantId();
                if (variantsToInclude == null || variantsToInclude.contains(primaryId)) {
                    Integer counter = datasetCounter.get(primaryId);
                    if (counter == null) {
                        counter = 1;
                    } else {
                        counter++;
                    }
                    datasetCounter.put(primaryId, counter);
                    tmpVariants.put(primaryId, variant);
                    allAvailableVariants.add(primaryId);
                }
            }
            if (tmpVariants.size() > maxNrVariants) {
                maxNrVariants = tmpVariants.size();
            }
        }

        System.out.println("Max nr variants: " + maxNrVariants);
        // now check whether we actually should exclude variants because they are not in both datasets..
        if (m_settings.getConfineSNPsToSNPsPresentInAllDatasets() != null && m_settings.getConfineSNPsToSNPsPresentInAllDatasets()) {
            HashSet<String> newTmpVariants = new HashSet<String>();
            for (String variant : allAvailableVariants) {
                Integer ctr = datasetCounter.get(variant);
                if (ctr == datasets.length) {

                    newTmpVariants.add(variant);
                }
            }
            allAvailableVariants = newTmpVariants;
            maxNrVariants = newTmpVariants.size();
        }
        System.out.println("Max nr variants: " + maxNrVariants);
        geneticVariantIndex = new GeneticVariant[datasets.length][maxNrVariants];
        int variantCounter = 0;
        for (String variant : allAvailableVariants) {
            for (int d = 0; d < datasets.length; d++) {
                GeneticVariant var = tmpVariantsArr.get(d).get(variant);
                geneticVariantIndex[d][variantCounter] = var;
            }
            variantCounter++;
        }
        System.out.println("Final variants: " + variantCounter);
    }
}
