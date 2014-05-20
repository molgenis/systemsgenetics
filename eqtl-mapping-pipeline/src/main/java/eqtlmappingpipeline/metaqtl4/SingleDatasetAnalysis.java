/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.apache.log4j.Logger;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.math.stats.Correlation;

/**
 *
 * @author Harm-Jan
 */
public class SingleDatasetAnalysis {

    private Set<MetaQTL4MetaTrait> traitsToInclude = null;
    private Set<String> variantsToInclude = null;
    private Set<Pair<String, MetaQTL4MetaTrait>> genotypeTraitCombinations = null;
    private final MetaQTL4Settings m_settings;
    private MetaQTL4TraitAnnotation traitAnnotation;
    private MetaQTL4Dataset dataset;
    private ArrayList<MetaQTL4MetaTrait> availableTraits;
    private TObjectIntHashMap<MetaQTL4MetaTrait> availableTraitsHash; //int as value
    private static final Logger LOG = Logger.getLogger(SingleDatasetAnalysis.class);
    private Integer[] traitIndex;

    public SingleDatasetAnalysis(MetaQTL4Settings settings) throws IOException, Exception {
        LOG.info("WARNING: MetaQTL4 is experimental code!");
        m_settings = settings;
        initialize();
        run();
    }

    public SingleDatasetAnalysis(String settings, String replaceText, String testToReplaceWith) throws IOException, Exception {
        LOG.warn("WARNING: MetaQTL4 is experimental code!");
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

        if (m_settings.getConfineProbeFile() != null) {
            if (!Gpio.existsAndReadable(m_settings.getConfineProbeFile())) {
                throw new IOException("Failed to read file with probes to include: " + m_settings.getConfineProbeFile().getAbsolutePath());
            }
            traitsToInclude = loadMetaTraitList(m_settings.getConfineProbeFile());
        }

        if (m_settings.getConfineSNPFile() != null) {
            if (!Gpio.existsAndReadable(m_settings.getConfineSNPFile())) {
                throw new IOException("Failed to read file with SNPs to include: " + m_settings.getConfineSNPFile().getAbsolutePath());
            }
            TextFile tf = new TextFile(m_settings.getConfineSNPFile(), TextFile.R);
            variantsToInclude = tf.readAsSet(0, TextFile.tab);
            tf.close();
        }

        // load list of specific SNP-trait combinations
        if (m_settings.getSNPProbeConfineFile() != null) {
            if (!Gpio.existsAndReadable(m_settings.getSNPProbeConfineFile())) {
                throw new IOException("Failed to read file with SNP Probe combinations to include: " + m_settings.getSNPProbeConfineFile().getAbsolutePath());
            }
            Triple<HashSet<String>, HashSet<MetaQTL4MetaTrait>, HashSet<Pair<String, MetaQTL4MetaTrait>>> combinations = loadgenotypeMarkerCombinations(m_settings.getSNPProbeConfineFile());
            variantsToInclude = combinations.getLeft();
            traitsToInclude = combinations.getMiddle();
            genotypeTraitCombinations = combinations.getRight();
        }

        // PROBE FILTER...
        LOG.info("Running probe filter");
        runProbeFilter();

        // load the datasets..
        HashSet<String> platformSpecificProbesToInclude = null;
        if (traitsToInclude != null) {
            String platform = m_settings.getDatasetSettings().get(0).getTraitPlatform();
            Integer platformId = traitAnnotation.getPlatformId(platform);
            platformSpecificProbesToInclude = new HashSet<String>();
            for (MetaQTL4MetaTrait metaprobe : traitsToInclude) {
                String platformProbe = metaprobe.getPlatformIds()[platformId];
                if (platformProbe != null) {
                    platformSpecificProbesToInclude.add(platformProbe);
                }
            }
        }

        dataset = new MetaQTL4Dataset(m_settings.getDatasetSettings().get(0), platformSpecificProbesToInclude);
        if (m_settings.isPerformParametricAnalysis()) {
            dataset.rankTraitData(m_settings.isEqualRankForTies());
        }

        createTraitIndex();

        // create lookup table for ZScores 
        // TODO: put this in a nice test encapsulation thing
        int nrSamples = dataset.getGenotypeToTraitCouplingInt().length;

        LOG.info("Meta-analysis will have " + nrSamples + " samples");
        Correlation.correlationToZScore(nrSamples);
    }

    public final void run() {
        // create threadpool
        LOG.info("Running software!");
        int nrPermutations = m_settings.getNrPermutationsFDR();

        // initialize random seed array
        long[] randomizationSeeds = new long[nrPermutations];
        Random rand = new Random();
        for (int permutation = 0; permutation < nrPermutations; permutation++) {
            randomizationSeeds[permutation] = rand.nextLong();
        }

        // create result thread
        ExecutorService resultPool = Executors.newFixedThreadPool(1);
        CompletionService resultPoolService = new ExecutorCompletionService(resultPool);

        // run threads
        int nrThreads = m_settings.getNrThreads();
        ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
        CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);
        MetaQTL4MetaTraitTreeSet traitTreeSet = new MetaQTL4MetaTraitTreeSet();
        for (MetaQTL4MetaTrait trait : availableTraits) {
            traitTreeSet.add(trait);
        }

        int ctr = 0;
        for (GeneticVariant variant : dataset.getGenotypeData()) {
            // submit task

            SingleDatasetAnalysisTask task = new SingleDatasetAnalysisTask(nrThreads, randomizationSeeds, availableTraits, traitTreeSet, traitIndex, dataset, m_settings, traitsToInclude, variant, resultPoolService);
            pool.submit(task);

            ctr++;
        }
//        pool.submit(task);

        int returned = 0;
        ProgressBar pb = new ProgressBar(ctr);
        while (returned < ctr) {
            try {
                Boolean result = pool.take().get();
                if (result != null) {
                    if (result) {
                        returned++;
                        pb.set(returned);
                    }
                }
                if (returned % 10000 == 0) {
                    System.out.println(returned + " have returned..");
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        System.out.println(returned);
        // close threadpool
        threadPool.shutdown();

        // perform FDR estimation
        // shutdown
    }

    // create a map from metaProbeId to dataset probes
    private void createTraitIndex() {
        // link them together in an index
        HashSet<MetaQTL4MetaTrait> tmpAvailableTraits = new HashSet<MetaQTL4MetaTrait>();

        String platform = m_settings.getDatasetSettings().get(0).getTraitPlatform();
        Integer platformId = traitAnnotation.getPlatformId(platform);
        String[] availableTraitsInDataset = dataset.getTraits();
        LOG.info(m_settings.getDatasetSettings().get(0).getName() + " has " + availableTraitsInDataset.length + " traits on platform " + platform + " (" + platformId + ")");
        for (String trait : availableTraitsInDataset) {
            MetaQTL4MetaTrait metatrait = traitAnnotation.getTraitForPlatformId(platformId, trait);
            if (metatrait != null) {
                String chr = metatrait.getChr();
                if (m_settings.isCisAnalysis() && m_settings.isTransAnalysis()) {
                    // include all traits when performing a cis+trans analysis
                    tmpAvailableTraits.add(metatrait);
                } else {
                    // otherwise, only include probes with a proper  annotation.
                    // should we include ChrY traits?
                    byte chrByte = ChrAnnotation.parseChr(chr);
                    if (chr.equals("X") || (chrByte < 23 && chrByte > 0)) {
                        tmpAvailableTraits.add(metatrait);
                    }
                }
            }
        }

        availableTraits = new ArrayList<MetaQTL4MetaTrait>(tmpAvailableTraits);
        availableTraitsHash = new TObjectIntHashMap<MetaQTL4MetaTrait>();
        traitIndex = new Integer[availableTraits.size()];

        for (int i = 0; i < availableTraits.size(); i++) {
            MetaQTL4MetaTrait trait = availableTraits.get(i);
            trait.setMetaTraitId(i);
            availableTraitsHash.put(trait, i);
            String probeName = trait.getPlatformIds()[platformId];
            if (probeName != null) {
                traitIndex[i] = dataset.getTraitId(probeName);
            }
        }

        LOG.info(availableTraits.size() + " traits available to test.");

    }

    // this filters probes on the basis of their annotation
    private void runProbeFilter() {
        // TODO: improve this code..
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

    private Triple<HashSet<String>, HashSet<MetaQTL4MetaTrait>, HashSet<Pair<String, MetaQTL4MetaTrait>>> loadgenotypeMarkerCombinations(File file) throws IOException {
        HashSet<Pair<String, MetaQTL4MetaTrait>> pairs = new HashSet<Pair<String, MetaQTL4MetaTrait>>();
        HashSet<String> snps = new HashSet<String>();
        HashSet<MetaQTL4MetaTrait> probes = new HashSet<MetaQTL4MetaTrait>();
        TextFile tf = new TextFile(file, TextFile.R);
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

    private Set<MetaQTL4MetaTrait> loadMetaTraitList(File location) throws IOException {
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
}
