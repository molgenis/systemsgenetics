/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.sampleFilter.SampleIncludedFilter;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilterBiAllelic;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.RankArray;

/**
 *
 * @author harmjan
 */
public class MetaQTL4Dataset {

    private final DoubleMatrixDataset<String, String> traitData;
    private DoubleMatrixDataset traitVarianceAndMeanData;
    private final RandomAccessGenotypeData genotypeData;
    private final HashMap<String, String> genotypeToTraitCoupling;
    private int[] intGenotypeToTraitCoupling;
    private final String genotypeDataLocation;

    public MetaQTL4Dataset(MetaQTL4DatasetSettings settings) throws IOException {
        this(settings.getTraitDataLocation(), null, settings.getGenotypeLocation(), settings.getGenotypeFormat(), null, settings.getGenotypeToTraitCoupling());
    }

    public MetaQTL4Dataset(MetaQTL4DatasetSettings settings, Set<String> markersToInclude) throws IOException {
        this(settings.getTraitDataLocation(), markersToInclude, settings.getGenotypeLocation(), settings.getGenotypeFormat(), null, settings.getGenotypeToTraitCoupling());
    }

    public MetaQTL4Dataset(String traitDataFile, Set<String> markersToInclude,
            String genotypeDatalocation, Set<String> snpsToInclude,
            RandomAccessGenotypeDataReaderFormats format) throws IOException {
        this(traitDataFile, markersToInclude, genotypeDatalocation, format, snpsToInclude, null);
    }

    public MetaQTL4Dataset(String traitDataFile, Set<String> markersToInclude, String genotypeDatalocation, RandomAccessGenotypeDataReaderFormats format, Set<String> snpsToInclude, String genotypeToTraitCouplingFile) throws IOException {
        this.genotypeDataLocation = genotypeDatalocation;
        RandomAccessGenotypeData tmpGenoData;
        if (snpsToInclude == null) {
            tmpGenoData = format.createFilteredGenotypeData(genotypeDatalocation, 0, null, new SampleIncludedFilter());
        } else {
            tmpGenoData = format.createFilteredGenotypeData(genotypeDatalocation, 0, new VariantIdIncludeFilter(snpsToInclude), new SampleIncludedFilter());
        }

//        VariantFilter qcFilter = new VariantCombinedFilter(new VariantQcChecker(0.05f, 0.001f, 0.95d), new VariantFilterBiAllelic());
//        for (GeneticVariant variant : tmpGenoData) {
//            if (qcFilter.doesVariantPassFilter(variant)) {
//                //pass
//            }
//        }
        genotypeData = new VariantFilterableGenotypeDataDecorator(tmpGenoData, new VariantFilterBiAllelic());

        // load genotype to expression couplings
        genotypeToTraitCoupling = new HashMap<String, String>();
        Map<String, String> gte = null;
        if (genotypeToTraitCouplingFile != null && Gpio.exists(genotypeToTraitCouplingFile)) {
            TextFile tf = new TextFile(genotypeToTraitCouplingFile, TextFile.R);
            gte = tf.readAsHashMap(0, 1);
            tf.close();

            HashSet<String> traitSamplesToSelect = new HashSet<String>();
            List<Sample> samples = genotypeData.getSamples();
            for (int i = 0; i < samples.size(); i++) {
                String sampleName = samples.get(i).getId();
                String linkedTrait = gte.get(sampleName);
                if (linkedTrait != null) {
                    traitSamplesToSelect.add(linkedTrait);
                }
            }
            traitData = new DoubleMatrixDataset<String, String>(traitDataFile, markersToInclude, traitSamplesToSelect);
            linkSamples(gte);
        } else {
            List<Sample> genotypeSamples = genotypeData.getSamples();
            HashSet<String> traitSamplesToSelect = new HashSet<String>();
            for (int s = 0; s < genotypeSamples.size(); s++) {
                traitSamplesToSelect.add(genotypeSamples.get(s).getId());
            }
            traitData = new DoubleMatrixDataset<String, String>(traitDataFile, null, traitSamplesToSelect);
            linkSamples(null);
        }
        pruneGenotypeSamples(tmpGenoData);
        reorderAndPruneTraitSamples();
        linkSamples(gte);

        determineProbeMeanAndVariance();
    }

    private void linkSamples(Map<String, String> gte) {
        if (gte != null) {
            List<Sample> genotypeSamples = genotypeData.getSamples();
            for (Sample s : genotypeSamples) {
                String genotypeSampleName = s.getId();
                String traitSampleName = gte.get(genotypeSampleName);
                if (traitSampleName != null && traitData.hashCols.containsKey(traitSampleName)) {
                    genotypeToTraitCoupling.put(genotypeSampleName, traitSampleName);
                }
            }
        } else {
            List<Sample> genotypeSamples = genotypeData.getSamples();
            for (Sample s : genotypeSamples) {
                String genotypeSampleName = s.getId();
                if (traitData.hashCols.containsKey(genotypeSampleName)) {
                    genotypeToTraitCoupling.put(genotypeSampleName, genotypeSampleName);
                }
            }
        }

        if (genotypeToTraitCoupling.isEmpty()) {
            throw new IllegalStateException("No links between genotypes and traits found for dataset: " + genotypeDataLocation);
        } else {
            System.out.println("Found "+genotypeToTraitCoupling.size()+ " links between genotypes and traits");
        }
//        intGenotypeToTraitCoupling = getPermutedSampleLinksInt();
    }

    // quite possibly improves performance when no missing data is present..
    private void reorderAndPruneTraitSamples() {

        // transpose: put samples on rows
        traitData.transposeDataset();

        // check the sample variance..
        boolean[] excludeSampleBecauseOfLowVariance = new boolean[traitData.nrRows];

        for (int row = 0; row < traitData.nrRows; row++) {
            double mean = Descriptives.mean(traitData.rawData[row]);
            double variance = Descriptives.variance(traitData.rawData[row], mean);
            if (variance < 1E-200 || Double.isNaN(variance)) {
                System.err.println("WARNING: sample\t" + traitData.rowObjects.get(row) + "\twill be excluded because sample variance " + variance + " < 1E-200");
                excludeSampleBecauseOfLowVariance[row] = true;
            } else {
                excludeSampleBecauseOfLowVariance[row] = false;
            }
        }

        // transpose: put samples back on columns
        traitData.transposeDataset();

        // currently, the trait matrix only contains samples present in both trait and genotype data..
        double[][] traitDataMatrix = traitData.getRawData();

        // make a list of samples to include
        List<Sample> genotypeSamples = genotypeData.getSamples();
        ArrayList<String> samplesToInclude = new ArrayList<String>();
        for (int sampleNr = 0; sampleNr < genotypeSamples.size(); sampleNr++) {
            String traitSample = genotypeToTraitCoupling.get(genotypeSamples.get(sampleNr).getId());
            Integer traitSampleId = traitData.hashCols.get(traitSample);
            if (!excludeSampleBecauseOfLowVariance[traitSampleId]) {
                samplesToInclude.add(traitSample);
            }
        }

        if (samplesToInclude.isEmpty()) {
            throw new IllegalStateException("After removing samples with low trait variance, no samples remain.");
        }

        // create an index of samples..
        int[] sampleIndexes = new int[samplesToInclude.size()];
        for (int sampleNr = 0; sampleNr < samplesToInclude.size(); sampleNr++) {
            int index = traitData.hashCols.get(samplesToInclude.get(sampleNr));
            sampleIndexes[sampleNr] = index;
        }

        // before creating a new trait matrix, determine the mean and variance for each probe using the samples to be included..
        boolean[] includeTrait = new boolean[traitData.nrRows];
        int nrTraitsToInclude = 0;
        for (int row = 0; row < traitData.nrRows; row++) {
            double sum = 0;
            int ctr = 0;
            for (int col = 0; col < traitData.nrCols; col++) {
                if (!excludeSampleBecauseOfLowVariance[col]) {
                    sum += traitData.rawData[row][col];
                    ctr++;
                }
            }
            double mean = sum / ctr;
            double var = 0.0;
            ctr = 0;
            for (int col = 0; col < traitData.nrCols; col++) {
                if (!excludeSampleBecauseOfLowVariance[col]) {
                    var += (traitData.rawData[row][col] - mean) * (traitData.rawData[row][col] - mean);
                    ctr++;
                }
            }
            var = var / (ctr - 1);

            if (var < 1E-200 || Double.isNaN(var)) {
                System.err.println("WARNING: trait\t" + traitData.rowObjects.get(row) + "\twill be excluded because variance " + var + " < 1E-200");
                includeTrait[row] = false;
            } else {
                nrTraitsToInclude++;
                includeTrait[row] = true;
            }
        }

        if (nrTraitsToInclude == 0) {
            throw new IllegalStateException("After removing traits with low variance, no traits remain.");
        }

        // create new trait matrix..
        double[][] newTraitMatrix = new double[nrTraitsToInclude][samplesToInclude.size()];
        int traitCtr = 0;
        ArrayList<String> newTraitIds = new ArrayList<String>();
        for (int row = 0; row < traitDataMatrix.length; row++) {
            if (includeTrait[row]) {
                for (int col = 0; col < sampleIndexes.length; col++) {
                    newTraitMatrix[traitCtr][col] = traitDataMatrix[row][sampleIndexes[col]];
                }
                newTraitIds.add(traitData.rowObjects.get(row));
                traitCtr++;
            }
        }

        traitData.setRawData(newTraitMatrix);
        traitData.rowObjects = newTraitIds;
        traitData.colObjects = samplesToInclude;
        traitData.recalculateHashMaps();
    }

    public HashMap<String, String> getSampleLinks() {
        return genotypeToTraitCoupling;
    }

//    public void permuteSampleLinks() {
//        if (permutedGenotypeToTraitCoupling == null) {
//            permutedGenotypeToTraitCoupling = new HashMap<String, String>();
//            intPermutedGenotypeToTraitCoupling = new int[intGenotypeToTraitCoupling.length];
//        }
//        // permute samples..
//        List<Sample> availableGenotypeSamples = genotypeData.getSamples();
//        ArrayList<Integer> availableTraitSamples = new ArrayList<Integer>();
//        for (int i = 0; i < intGenotypeToTraitCoupling.length; i++) {
//            if (intGenotypeToTraitCoupling[i] != -1) {
//                availableTraitSamples.add(intGenotypeToTraitCoupling[i]);
//            }
//        }
//
//        for (int i = 0; i < intGenotypeToTraitCoupling.length; i++) {
//            if (intGenotypeToTraitCoupling[i] != -1) {
//                intPermutedGenotypeToTraitCoupling[i] = availableTraitSamples.remove((int) (Math.random() * (double) availableTraitSamples.size()));
//                permutedGenotypeToTraitCoupling.put(availableGenotypeSamples.get(i).getId(), traitData.colObjects.get(intPermutedGenotypeToTraitCoupling[i]));
//            } else {
//                intPermutedGenotypeToTraitCoupling[i] = -1;
//            }
//        }
//    }
//    public HashMap<String, String> getPermutedSampleLinks() {
//        return permutedGenotypeToTraitCoupling;
//    }
    public int[] getGenotypeToTraitCouplingInt() {
        if (intGenotypeToTraitCoupling == null) {
            List<Sample> genotypeSamples = genotypeData.getSamples();
            intGenotypeToTraitCoupling = new int[genotypeSamples.size()];
            for (int i = 0; i < genotypeSamples.size(); i++) {
                String traitSample = genotypeToTraitCoupling.get(genotypeSamples.get(i).getId());
                if (traitSample == null) {
                    intGenotypeToTraitCoupling[i] = -1;
                } else {
                    intGenotypeToTraitCoupling[i] = traitData.hashCols.get(traitSample);
                }
            }
        }
        return intGenotypeToTraitCoupling;
    }

//    public int[] getPermutedSampleLinksInt() {
//        return intPermutedGenotypeToTraitCoupling;
//    }
    public void close() throws IOException {
        genotypeData.close();
    }

    public void rankTraitData(boolean rankWithTies) {

        RankArray r = new RankArray();
        for (int p = 0; p < traitData.nrRows; ++p) {
            double[] probeData = traitData.rawData[p];

            // because we have checked low variance traits before, we should never have to enter this loop.
            if (traitVarianceAndMeanData.rawData[0][p] == 0) {
                System.out.println("Trait that has no variance in expression:\t" + traitData.rowObjects.get(p));
            } else {
                if (Double.isNaN(traitVarianceAndMeanData.rawData[0][p]) || Double.isNaN(traitVarianceAndMeanData.rawData[1][p])) {
                    throw new IllegalStateException("Error ranking trait data: mean or variance is NaN!:\t" + p + "\t" + traitData.rowObjects.get(p) + "\tMean: " + traitVarianceAndMeanData.rawData[0][p] + "\tVariance: " + traitVarianceAndMeanData.rawData[1][p]);
                } else {
                    probeData = r.rank(probeData, rankWithTies);
                    traitData.rawData[p] = probeData;
                    traitVarianceAndMeanData.rawData[0][p] = Descriptives.mean(probeData);
                    traitVarianceAndMeanData.rawData[1][p] = Descriptives.variance(probeData, traitVarianceAndMeanData.rawData[0][p]);
                }
            }
        }
    }

    private void determineProbeMeanAndVariance() {
        traitVarianceAndMeanData = new DoubleMatrixDataset();
        traitVarianceAndMeanData.rawData = new double[2][traitData.nrRows];
        for (int i = 0; i < traitData.nrRows; i++) {
            double mean = Descriptives.variance(traitData.rawData[i]);
            traitVarianceAndMeanData.rawData[0][i] = mean;
            traitVarianceAndMeanData.rawData[1][i] = Descriptives.variance(traitData.rawData[i], mean);
        }
    }

    public String[] getTraits() {
        return traitData.rowObjects.toArray(new String[0]);
    }

    public Integer getTraitId(String probeName) {
        return traitData.hashRows.get(probeName);
    }

    public RandomAccessGenotypeData getGenotypeData() {
        return genotypeData;
    }

    // remove genotype samples from the genotypeData that are not linked to trait samples
    private void pruneGenotypeSamples(RandomAccessGenotypeData tmpGenoData) {
//        genotypeToTraitCoupling
        List<Sample> samples = tmpGenoData.getSamples();
        for (Sample sample : samples) {
            String sampleId = sample.getId();
            Boolean includeSample = false;
            if (genotypeToTraitCoupling.containsKey(sampleId)) {
                includeSample = true;
            }
            sample.putAnnotationValues(GenotypeData.BOOL_INCLUDE_SAMPLE, includeSample);
        }
    }

    public synchronized float[] getSampleDosages(GeneticVariant variant) {
        return variant.getSampleDosages();
    }

    public double[] getTraitData(Integer traitId) {
        return traitData.get(traitId);
    }

    public double getTraitMean(Integer traitId) {
        return traitVarianceAndMeanData.get(traitId, 0);
    }

    public double getTraitVariance(Integer traitId) {
        return traitVarianceAndMeanData.get(traitId, 1);
    }
}
