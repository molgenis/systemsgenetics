/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
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
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.math.stats.Correlation;

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
    private TObjectIntHashMap<MetaQTL4MetaTrait> availableTraitsHash; //int as value
    private GeneticVariant[][] geneticVariantIndex;
    private static final Logger LOG = Logger.getLogger(MetaQTL4.class);

    public MetaQTL4(MetaQTL4Settings settings) throws IOException, Exception {
        LOG.info("WARNING: MetaQTL4 is experimental code!");
        m_settings = settings;
        initialize();
        run();
    }

    public MetaQTL4(String settings, String replaceText, String testToReplaceWith) throws IOException, Exception {
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
        // TODO: there is one exception for the probe filter, 
        // which is when all the datasets are on the same platform. 
        // will implement this at a later stage, or not at all.
        LOG.info("Running probe filter");
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
        LOG.info("Creating trait index");
        createTraitIndex();

        // create SNP index
        LOG.info("Create SNP index");
        createSNPIndex();

        // create lookup table for ZScores 
        // TODO: put this in a nice test encapsulation thing
        int nrSamples = 0;
        for (MetaQTL4Dataset dataset : datasets) {
            nrSamples += dataset.getGenotypeToTraitCouplingInt().length;
        }
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
        int bufferSize = 100;
        int nrThreads = m_settings.getNrThreads();
        ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
        CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);
        int start = 0;
        int stop = geneticVariantIndex.length;
        MetaQTL4ExecutionTask task = new MetaQTL4ExecutionTask(nrThreads, randomizationSeeds, availableTraits, availableTraitsHash, datasets, geneticVariantIndex, m_settings, traitAnnotation, traitIndex, traitsToInclude, variantsToInclude, start, stop, bufferSize, resultPoolService);
        
//        MetaQTL4CorrelationTask task = new MetaQTL4CorrelationTask(nrThreads, distributionSize, randomizationSeeds, availableTraits, availableTraitsHash, datasets, geneticVariantIndex, m_settings, traitAnnotation, traitIndex, traitsToInclude, variantsToInclude, 1);
        pool.submit(task);

        int returned = 0;
        while (returned < nrThreads) {
            try {
                Boolean result = pool.take().get();
                if (result != null) {
                    if(result){
                    returned++;    
                    }
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

    // create a map from metaProbeId to dataset probes
    private void createTraitIndex() {
        // link them together in an index
        HashSet<MetaQTL4MetaTrait> tmpAvailableTraits = new HashSet<MetaQTL4MetaTrait>();
        for (int d = 0; d < datasets.length; d++) {
            String platform = m_settings.getDatasetSettings().get(d).getTraitPlatform();
            Integer platformId = traitAnnotation.getPlatformId(platform);
            String[] availableTraitsInDataset = datasets[d].getTraits();
            LOG.info(m_settings.getDatasetSettings().get(d).getName() + " has " + availableTraitsInDataset.length + " traits on platform " + platform + " (" + platformId + ")");
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
        }

        availableTraits = new ArrayList<MetaQTL4MetaTrait>(tmpAvailableTraits);
        availableTraitsHash = new TObjectIntHashMap<MetaQTL4MetaTrait>();
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

        LOG.info(availableTraits.size() + " traits available to test.");

    }

    private void createSNPIndex() {

        if (traitIndex == null) {
            createTraitIndex();
        }

        // if this is a cis-eQTL analysis, get all the unique variants within the probes vicinity
        LOG.info("Indexing SNPs");
        int distance = m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance();
        Set<Pair<String, Integer>> uniquePositions = new HashSet<Pair<String, Integer>>();

        // at a later stage, we may want to include gene-begin and gene-end as well..
        if (m_settings.isCisAnalysis() && !m_settings.isTransAnalysis()) {
            LOG.info("Performing cis-eQTL analysis, so keeping only SNPs within " + m_settings.getCiseQTLAnalysMaxSNPProbeMidPointDistance());
            LOG.info("");
            int tctr = 0;
            HashMap<Pair<String, Integer>, Integer> allPositions = new HashMap<Pair<String, Integer>, Integer>();
            for (int d = 0; d < datasets.length; d++) {
                RandomAccessGenotypeData data = datasets[d].getGenotypeData();
                ProgressBar pb = new ProgressBar(availableTraits.size(), "Dataset " + d);
                for (MetaQTL4MetaTrait t : availableTraits) {
                    String chr = t.getChr();
                    int midpoint = t.getChrMidpoint();

                    Iterable<GeneticVariant> variants = data.getVariantsByRange(chr, midpoint - distance, midpoint + distance);
                    for (GeneticVariant variant : variants) {
                        String primaryId = variant.getPrimaryVariantId();
                        if (variantsToInclude == null || variantsToInclude.contains(primaryId)) {
                            // we don't need to check the variant chromosome.. This was done already when creating the trait index
                            Pair<String, Integer> position = new Pair<String, Integer>(variant.getSequenceName(), variant.getStartPos());
                            Integer ctr = allPositions.get(position);
                            if (ctr == null) {
                                ctr = 1;
                            } else {
                                ctr++;
                            }
                            allPositions.put(position, ctr);
                        }
                    }

                    pb.iterate();
                    tctr++;
                }
                pb.close();
            }

            LOG.info(allPositions.size() + " variants amongst all datasets.");

            if (m_settings.getConfineSNPsToSNPsPresentInAllDatasets() != null && m_settings.getConfineSNPsToSNPsPresentInAllDatasets()) {
                Set<Pair<String, Integer>> tmpPositions = new HashSet<Pair<String, Integer>>();
                for (Pair<String, Integer> pair : allPositions.keySet()) {
                    Integer ctr = allPositions.get(pair);
                    if (ctr == datasets.length) {
                        // include snp
                        tmpPositions.add(pair);
                    }
                }
                uniquePositions = tmpPositions;
            } else {
                uniquePositions = allPositions.keySet();
            }
            allPositions = null;

        } else { // if cistrans or trans.. include all the SNPs..
            // per chromosome.. (for huge datasets, this saves alot of Integer objects)
            for (int i = 0; i < 23; i++) {
                String chr = "" + i;
                if (i == 23) {
                    chr = "X";
                }

                // get all the unique positions for this chromosome.
                HashMap<Integer, Integer> allPositions = new HashMap<Integer, Integer>();
                for (int d = 0; d < datasets.length; d++) {
                    RandomAccessGenotypeData data = datasets[d].getGenotypeData();
                    Iterable<GeneticVariant> variants = data.getSequenceGeneticVariants(chr);
                    for (GeneticVariant variant : variants) {
                        String primaryId = variant.getPrimaryVariantId();
                        if (variantsToInclude == null || variantsToInclude.contains(primaryId)) {
                            // this time we do need to check the chromosome, etc
                            String variantSequence = variant.getSequenceName();
                            if (variantSequence.equals(chr)) {
                                Integer position = variant.getStartPos();
                                if (position > 0) {
                                    Integer ctr = allPositions.get(position);
                                    if (ctr == null) {
                                        ctr = 1;
                                    } else {
                                        ctr++;
                                    }
                                    allPositions.put(position, ctr);
                                }
                            }
                        }
                    }
                }

                for (Integer position : allPositions.keySet()) {
                    Integer ctr = allPositions.get(position);
                    if (m_settings.getConfineSNPsToSNPsPresentInAllDatasets() != null
                            && (m_settings.getConfineSNPsToSNPsPresentInAllDatasets() && ctr == datasets.length) || (!m_settings.getConfineSNPsToSNPsPresentInAllDatasets())) {
                        uniquePositions.add(new Pair<String, Integer>(chr, position));
                    }
                }
            }
        }
        // ready to index
        int numberFinalPositions = uniquePositions.size();
        if (numberFinalPositions == 0) {
            System.err.println("Error: no SNPs found to test");
        }
        LOG.info("Creating final index");
        geneticVariantIndex = new GeneticVariant[datasets.length][numberFinalPositions];
        ProgressBar pb2 = new ProgressBar(datasets.length * numberFinalPositions);
        for (int dataset = 0; dataset < datasets.length; dataset++) {
            int ctr = 0;
            RandomAccessGenotypeData data = datasets[dataset].getGenotypeData();
            for (Pair<String, Integer> pair : uniquePositions) {
                geneticVariantIndex[dataset][ctr] = data.getSnpVariantByPos(pair.getLeft(), pair.getRight());
                ctr++;
                pb2.iterate();
            }
        }
        pb2.close();
        LOG.info("Done");
    }
}
