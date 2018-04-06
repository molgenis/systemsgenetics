package org.molgenis.genotype.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.VcfAnnotation;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.tabix.TabixIndex;
import org.molgenis.genotype.tabix.TabixIndex.TabixIterator;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.vcf.VcfInfo;
import org.molgenis.vcf.VcfReader;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaContig;
import org.molgenis.vcf.meta.VcfMetaInfo;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import java.util.Arrays;
import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;

public class VcfGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

    private static final org.apache.log4j.Logger LOG = org.apache.log4j.Logger.getLogger(VcfGenotypeData.class);
    private final File bzipVcfFile;
    private final TabixIndex tabixIndex;
    private final int sampleVariantProviderUniqueId;
    private final SampleVariantsProvider variantProvider;
    private final VcfMeta vcfMeta;
    private transient Map<String, Annotation> cachedSampleAnnotationsMap;
    private transient GeneticVariant cachedGeneticVariant;
    private transient VcfRecord cachedVcfRecord;
    private static int totalRandomAccessRequest = 0;
    private static int currentlyOpenFileHandlers = 0;
    private static int closedFileHandlers = 0;
    private final double minimumPosteriorProbabilityToCall;

    /**
     * VCF genotype reader
     *
     * @param bzipVcfFile
     * @throws IOException
     * @throws FileNotFoundException
     */
    //public VcfGenotypeData(File bzipVcfFile, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
    //	this(bzipVcfFile, 100, minimumPosteriorProbabilityToCall);
    //}
    public VcfGenotypeData(File bzipVcfFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
        this(bzipVcfFile, new File(bzipVcfFile.getAbsolutePath() + ".tbi"), cacheSize, minimumPosteriorProbabilityToCall);
    }

    /**
     * VCF genotype reader with default cache of 100
     *
     * @param bzipVcfFile
     * @param tabixIndexFile
     * @throws IOException
     * @throws FileNotFoundException
     */
    public VcfGenotypeData(File bzipVcfFile, File tabixIndexFile, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
        this(bzipVcfFile, tabixIndexFile, 100, minimumPosteriorProbabilityToCall);
    }

    public VcfGenotypeData(File bzipVcfFile, File tabixIndexFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException,
            IOException {

        if (!bzipVcfFile.isFile()) {
            throw new FileNotFoundException("VCF file not found at " + bzipVcfFile.getAbsolutePath());
        }

        if (!bzipVcfFile.canRead()) {
            throw new IOException("Cannot access VCF file at: " + bzipVcfFile.getAbsolutePath());
        }

        if (!tabixIndexFile.isFile()) {
            throw new TabixFileNotFoundException(tabixIndexFile.getAbsolutePath(), "VCF tabix file not found at " + tabixIndexFile.getAbsolutePath());
        }

        if (!tabixIndexFile.canRead()) {
            throw new IOException("Cannot read tabix file for VCF at: " + tabixIndexFile.getAbsolutePath());
        }

        if (minimumPosteriorProbabilityToCall < 0 || minimumPosteriorProbabilityToCall > 1) {
            throw new GenotypeDataException("Min posterior probability to call must be >0 and <=1 not:" + minimumPosteriorProbabilityToCall);
        }

        this.bzipVcfFile = bzipVcfFile;
        this.tabixIndex = new TabixIndex(tabixIndexFile, bzipVcfFile, null);
        this.minimumPosteriorProbabilityToCall = minimumPosteriorProbabilityToCall;

        VcfReader vcfReader = new VcfReader(new BlockCompressedInputStream(bzipVcfFile));
        try {
            this.vcfMeta = vcfReader.getVcfMeta();
        } finally {
            vcfReader.close();
        }

        if (cacheSize <= 0) {
            variantProvider = this;
        } else {
            variantProvider = new CachedSampleVariantProvider(this, cacheSize);
        }

        sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
    }

    @Override
    public Iterator<GeneticVariant> iterator() {
        final BlockCompressedInputStream inputStream;
        try {
            inputStream = new BlockCompressedInputStream(bzipVcfFile);
        } catch (IOException e) {
            throw new GenotypeDataException(e);
        }

        final VcfReader vcfReader = new VcfReader(inputStream);
        Iterator<GeneticVariant> iterator = new Iterator<GeneticVariant>() {
            private final Iterator<VcfRecord> it = vcfReader.iterator();

            @Override
            public boolean hasNext() {
                boolean hasNext = it.hasNext();
                if (!hasNext) {
                    // close vcf reader
                    try {
                        vcfReader.close();
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }
                return hasNext;
            }

            @Override
            public GeneticVariant next() {
                VcfRecord vcfRecord = it.next();
                return toGeneticVariant(vcfRecord);
            }

            @Override
            public void remove() {
                it.remove();
            }
        };
        iterator.hasNext();//needed to properly initiate
        return iterator;
    }

    @Override
    public List<Alleles> getSampleVariants(final GeneticVariant variant) {
        // get vcf record for variant
        VcfRecord vcfRecord = getVcfRecord(variant);
        int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return Collections.emptyList();
        }

        int idx = vcfRecord.getFormatIndex("GT");
        if (idx != -1) {

            // convert vcf sample alleles to Alleles
            List<Alleles> alleles = new ArrayList<Alleles>(nrSamples);

            try {
                for (VcfSample vcfSample : vcfRecord.getSamples()) {
                    List<Allele> vcfAlleles = vcfSample.getAlleles();
                    alleles.add(Alleles.createAlleles(vcfAlleles));
                }
            } catch (NumberFormatException ex) {
                throw new GenotypeDataException("Error parsing variant: " + variant.getPrimaryVariantId() + " at " + variant.getSequenceName() + ":" + variant.getStartPos(), ex);
            }
            return alleles;

        } else if (vcfRecord.getFormatIndex("GP") != -1) {

            return ProbabilitiesConvertor.convertProbabilitiesToAlleles(getSampleProbilities(variant), variant.getVariantAlleles(), minimumPosteriorProbabilityToCall);

        } else if (vcfRecord.getFormatIndex("DS") != -1) {

            return CalledDosageConvertor.convertDosageToAlleles(getSampleDosage(variant), variant.getVariantAlleles());

        } else {

            ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>(nrSamples);
            for (int i = 0; i < nrSamples; ++i) {
                sampleAlleles.add(Alleles.BI_ALLELIC_MISSING);
            }
            return sampleAlleles;
        }

    }

    @Override
    public Map<String, Annotation> getVariantAnnotationsMap() {
        if (cachedSampleAnnotationsMap == null) {
            cachedSampleAnnotationsMap = new LinkedHashMap<String, Annotation>();

            for (VcfMetaInfo info : vcfMeta.getInfoMeta()) {
                cachedSampleAnnotationsMap.put(info.getId(), VcfAnnotation.fromVcfInfo(info));
            }
        }

        return cachedSampleAnnotationsMap;
    }

    @Override
    public List<Sequence> getSequences() {
        List<String> seqNames = getSeqNames();

        // get sequence length by sequence name
        Map<String, Integer> sequenceLengthById = new HashMap<String, Integer>();
        for (VcfMetaContig contig : vcfMeta.getContigMeta()) {
            sequenceLengthById.put(contig.getId(), contig.getLength());
        }

        List<Sequence> sequences = new ArrayList<Sequence>(seqNames.size());
        for (String seqName : seqNames) {
            sequences.add(new SimpleSequence(seqName, sequenceLengthById.get(seqName), this));
        }

        return sequences;
    }

    @Override
    public List<Sample> getSamples() {
        List<Sample> samples = new ArrayList<Sample>();
        for (String sampleName : vcfMeta.getSampleNames()) {
            samples.add(new Sample(sampleName, null, null));
        }
        return samples;
    }

    @Override
    public int cacheSize() {
        return 0;
    }

    @Override
    public List<Boolean> getSamplePhasing(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return Collections.emptyList();
        }

        List<Boolean> phasing = new ArrayList<Boolean>(nrSamples);
        for (VcfSample vcfSample : vcfRecord.getSamples()) {

            List<Boolean> genotypePhasings = vcfSample.getPhasings();

            if (genotypePhasings == null || genotypePhasings.isEmpty()) {
                phasing.add(Boolean.FALSE);
            } else if (genotypePhasings.size() == 1) {
                phasing.add(genotypePhasings.get(0));
            } else if (genotypePhasings.contains(Boolean.FALSE)) {
                phasing.add(Boolean.FALSE);
            } else {
                phasing.add(Boolean.TRUE);
            }

        }
        return phasing;
    }

    @Override
    public int getSampleVariantProviderUniqueId() {
        return sampleVariantProviderUniqueId;
    }

    @Override
    public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
        return Collections.emptyMap();
    }

    @Override
    public byte[] getSampleCalledDosage(GeneticVariant variant) {
        return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
                variant.getVariantAlleles(), variant.getRefAllele());
    }

    @Override
    public float[] getSampleDosage(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return new float[0];
        }

        float[] dosages;

        int idx = vcfRecord.getFormatIndex("DS");
        if (idx != -1) {
            // retrieve sample dosage from sample info
            dosages = new float[nrSamples];
            int i = 0;
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String dosage = vcfSample.getData(idx);
                if (dosage == null) {
                    //throw new GenotypeDataException("Missing DS format value for sample [" + vcfMeta.getSampleName(i) + "] at variant [" + variant.getPrimaryVariantId() + "]");
                    dosages[i++] = -1;
                } else {
                    try {
                        //Math abs to prevent -0 due to rounding
                        dosages[i++] = Math.abs((Float.parseFloat(dosage) - 2) * -1);
                    } catch (NumberFormatException e) {
                        throw new GenotypeDataException("Error in sample dosage (DS) value for sample [" + vcfMeta.getSampleName(i) + "], found value: " + dosage);
                    }
                }

            }
        } else if (vcfRecord.getFormatIndex("GP") != -1) {
            dosages = ProbabilitiesConvertor.convertProbabilitiesToDosage(getSampleProbilities(variant), minimumPosteriorProbabilityToCall);
        } else if (vcfRecord.getFormatIndex("GT") != -1) {
            // calculate sample dosage from called alleles
            dosages = CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
                    variant.getVariantAlleles(), variant.getRefAllele());
        } else {
            dosages = new float[nrSamples];
            for (int i = 0; i < nrSamples; ++i) {
                dosages[i] = -1;
            }
        }
        return dosages;
    }

    @Override
    public void close() throws IOException {
        // noop
    }

    @Override
    public boolean isOnlyContaingSaveProbabilityGenotypes() {
        return false;
    }

    @Override
    public float[][] getSampleProbilities(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return new float[0][0];
        }

        float[][] probs = null;

        int idx = vcfRecord.getFormatIndex("GP");
        if (idx != -1) {
            // retrieve sample probabilities from sample info
            probs = new float[nrSamples][3];
            int i = 0;
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String probabilitiesStr = vcfSample.getData(idx);
                if (probabilitiesStr == null) {
                    //throw new GenotypeDataException("Missing GP format value for sample [" + vcfMeta.getSampleName(i) + "]");
                    probs[i] = new float[]{0, 0, 0};
                } else {
                    if (probabilitiesStr.matches(".*,+\\.,+.*")) {
//                        System.out.println(probabilitiesStr);
                        probabilitiesStr = probabilitiesStr.replaceAll("\\.", "0");
//                        System.out.println(probabilitiesStr);
                    }
                    String[] probabilities = StringUtils.split(probabilitiesStr, ',');
                    if (probabilities.length != 3) {
                        throw new GenotypeDataException("Error in sample prob (GP) value for sample [" + vcfMeta.getSampleName(i) + "], found value: " + probabilitiesStr);
                    }

                    for (int j = 0; j < 3; ++j) {
                        try {
                            probs[i][j] = Float.parseFloat(probabilities[j]);
                        } catch (NumberFormatException e) {
                            throw new GenotypeDataException("Error in sample prob (GP) value for sample [" + vcfMeta.getSampleName(i) + "], found value: " + probabilitiesStr);
                        }
                    }
                }
                ++i;
            }

        } else if (vcfRecord.getFormatIndex("GT") != -1) {
            probs = ProbabilitiesConvertor.convertCalledAllelesToProbability(getSampleVariants(variant), variant.getVariantAlleles());
        } else if (vcfRecord.getFormatIndex("DS") != -1) {
            // calculate sample probabilities from sample dosage
            probs = ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(getSampleDosage(variant));
        } else {
            probs = new float[nrSamples][3];
        }
        return probs;
    }

    @Override
    public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
        final VcfRecord vcfRecord = getVcfRecord(variant);

        return new FixedSizeIterable<GenotypeRecord>() {
            @Override
            public Iterator<GenotypeRecord> iterator() {
                Iterable<VcfSample> samples = vcfRecord.getSamples();
                return Iterators.transform(samples.iterator(), new Function<VcfSample, GenotypeRecord>() {
                    @Override
                    public GenotypeRecord apply(VcfSample vcfSample) {
                        return toGenotypeRecord(vcfRecord, vcfSample);
                    }
                });
            }

            @Override
            public int size() {
                return vcfRecord.getNrSamples();
            }
        };
    }

    @Override
    public List<String> getSeqNames() {
        return new ArrayList<String>(tabixIndex.getSeqNames());
    }

    @Override
    public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
        return getVariantsByRange(seqName, 0, Integer.MAX_VALUE);
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByPos(final String seqName, final int startPos) {
        return getVariantsByRange(seqName, startPos - 1, startPos);
    }

    @Override
    public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {

        //NOTE this special version will complete the iteration and therefor the file connection will be closed
        //This is kind of an ugly hack but is functional
        Iterable<GeneticVariant> variants = getVariantsByPos(seqName, startPos);
        GeneticVariant snp = null;
        for (GeneticVariant variant : variants) {
            if (snp == null && variant.isSnp()) {
                snp = variant;
            }
        }

        return snp;

    }

    @Override
    public Iterable<GeneticVariant> getVariantsByRange(final String seqName, final int rangeStart, final int rangeEnd) {

        if (rangeStart < 0) {
            throw new GenotypeDataException("Illegal start pos for VCF variant query: " + rangeStart);
        }

        return new Iterable<GeneticVariant>() {
            @Override
            public Iterator<GeneticVariant> iterator() {

                try {

                    ++totalRandomAccessRequest;
                    ++currentlyOpenFileHandlers;

                    return new Iterator<GeneticVariant>() {
                        private final BlockCompressedInputStream stream = new BlockCompressedInputStream(bzipVcfFile);
                        private final TabixIterator it = tabixIndex.queryTabixIndex(seqName, rangeStart, rangeEnd, stream);
                        private String line = readFirst(it);

                        private String readFirst(TabixIterator it) {
                            if (it == null) {
                                try {
                                    stream.close();
                                } catch (IOException e) {
                                    throw new GenotypeDataException(e);
                                }
                                --currentlyOpenFileHandlers;
                                ++closedFileHandlers;
                                return null;
                            } else {
                                try {
                                    String firstLine = it.next();
                                    if (firstLine == null) {
                                        --currentlyOpenFileHandlers;
                                        ++closedFileHandlers;
                                        it.close();
                                    }
                                    return firstLine;
                                } catch (IOException e) {
                                    throw new GenotypeDataException(e);
                                }
                            }
                        }

                        @Override
                        public boolean hasNext() {
                            return line != null;
                        }

                        @Override
                        public GeneticVariant next() {
                            VcfRecord vcfRecord = new VcfRecord(vcfMeta, StringUtils.split(line, '\t'));
                            try {
                                line = it.next();
                                if (line == null) {
                                    --currentlyOpenFileHandlers;
                                    ++closedFileHandlers;
                                    it.close();
                                }
                            } catch (IOException e) {
                                throw new GenotypeDataException(e);
                            }
                            return toGeneticVariant(vcfRecord);
                        }

                        @Override
                        public void remove() {
                            throw new UnsupportedOperationException();
                        }
                    };
                } catch (FileNotFoundException e) {
                    if (e.getMessage().endsWith("(Too many open files)")) {
                        throw new GenotypeDataException("VCF reader trying to open more file connections than allowed by operating system. Currently open connections: " + currentlyOpenFileHandlers + " total opened: " + totalRandomAccessRequest + " total closed: " + closedFileHandlers, e);
                    } else {
                        throw new GenotypeDataException(e);
                    }
                } catch (IOException e) {
                    throw new GenotypeDataException(e);
                }
            }
        };
    }

    private VcfRecord getVcfRecord(GeneticVariant variant) {
        if (!variant.equals(cachedGeneticVariant)) {
            TabixIterator it;
            String line;
            BlockCompressedInputStream stream = null;
            ++totalRandomAccessRequest;
            ++currentlyOpenFileHandlers;

            try {
                stream = new BlockCompressedInputStream(bzipVcfFile);
                it = tabixIndex.queryTabixIndex(variant.getSequenceName(), variant.getStartPos() - 1, variant.getStartPos(), stream);
                while ((line = it.next()) != null) {
                    VcfRecord vcfRecord = new VcfRecord(vcfMeta, StringUtils.split(line, '\t'));
                    if (variant.equals(toGeneticVariant(vcfRecord))) {
                        cachedVcfRecord = vcfRecord;
                        cachedGeneticVariant = variant;
                        break;
                    }
                }
                stream.close();
            } catch (FileNotFoundException e) {
                if (e.getMessage().endsWith("(Too many open files)")) {
                    throw new GenotypeDataException("VCF reader trying to open more file connections than allowed by operating system. Currently open connections: " + currentlyOpenFileHandlers + " total opened: " + totalRandomAccessRequest + " total closed: " + closedFileHandlers, e);
                } else {
                    throw new GenotypeDataException(e);
                }
            } catch (IOException e) {
                throw new GenotypeDataException(e);
            } finally {
                --currentlyOpenFileHandlers;
                ++closedFileHandlers;
                IOUtils.closeQuietly(stream);
            }
        }
        return cachedVcfRecord;
    }

    /**
     * Convert VcfRecord to GeneticVariant
     *
     * @param vcfRecord
     * @return
     */
    private GeneticVariant toGeneticVariant(VcfRecord vcfRecord) {
        List<String> identifiers = vcfRecord.getIdentifiers();
        int pos = vcfRecord.getPosition();
        String sequenceName = vcfRecord.getChromosome();
        Allele refAllele = vcfRecord.getReferenceAllele();
        List<Allele> altAlleles = vcfRecord.getAlternateAlleles();

        Map<String, Object> annotationMap = new HashMap<String, Object>();
        for (VcfInfo vcfInfo : vcfRecord.getInformation()) {
            annotationMap.put(vcfInfo.getKey(), vcfInfo.getVal());
        }
        annotationMap.put("VCF_Filter", vcfRecord.getFilterStatus());
        annotationMap.put("VCF_Qual", vcfRecord.getQuality());

        Alleles alleles;
        if (altAlleles == null || altAlleles.isEmpty()) {
            alleles = Alleles.createAlleles(refAllele);
        } else {
            ArrayList<Allele> allelesList = new ArrayList<Allele>(altAlleles.size() + 1);
            allelesList.add(refAllele);
            allelesList.addAll(altAlleles);
            alleles = Alleles.createAlleles(allelesList);
        }

        GeneticVariantMeta geneticVariantMeta = new VcfGeneticVariantMeta(vcfMeta, Arrays.asList(vcfRecord.getFormat()));
        GeneticVariant geneticVariant = ReadOnlyGeneticVariant.createVariant(geneticVariantMeta, identifiers, pos, sequenceName, annotationMap, variantProvider, alleles, refAllele);

        cachedVcfRecord = vcfRecord;
        cachedGeneticVariant = geneticVariant;
        return geneticVariant;
    }

    /**
     * Convert VcfSample to GenotypeRecord
     *
     * @param vcfRecord
     * @param vcfSample
     * @return
     */
    private GenotypeRecord toGenotypeRecord(VcfRecord vcfRecord, VcfSample vcfSample) {
        return new VcfGenotypeRecord(vcfMeta, vcfRecord, vcfSample);
    }
}
