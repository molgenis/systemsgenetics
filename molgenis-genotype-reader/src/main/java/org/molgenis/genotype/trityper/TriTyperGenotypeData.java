/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.trityper;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import static org.molgenis.genotype.GenotypeData.BOOL_INCLUDE_SAMPLE;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.CaseControlAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.GeneticVariantTreeSet;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariantTriTyper;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.genotype.variantFilter.VariantFilter;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class TriTyperGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

    private ArrayList<Sample> samples;
    private final List<Boolean> samplePhasing;
    private GeneticVariantTreeSet<GeneticVariant> snps = new GeneticVariantTreeSet<GeneticVariant>();
    private final SampleVariantsProvider variantProvider;
    private HashMap<GeneticVariant, Integer> snpToIndex;
    private final File genotypeDataFile;
    private final File imputedDosageDataFile;
    private final File snpFile;
    private final File snpMapFile;
    private final File individualFile;
    private final File phenotypeAnnotationFile;
    private final File baseDir;
    private static final Logger LOG = Logger.getLogger(TriTyperGenotypeData.class);
    private final int cacheSize;
    private final RandomAccessFile dosageHandle;
    private final RandomAccessFile genotypeHandle;
    private final FileChannel genotypeChannel;
    private final FileChannel dosageChannel;
    private final int sampleVariantProviderUniqueId;
    private HashMap<String, SampleAnnotation> sampleAnnotationMap;
    private HashMap<String, Sequence> sequences;
	private final VariantFilter variantFilter;

    public TriTyperGenotypeData(String location) throws IOException {
        this(location, 1024);
    }

	public TriTyperGenotypeData(String location, int cacheSize) throws IOException {
		this(location, cacheSize, null);
	}	

    public TriTyperGenotypeData(String location, int cacheSize, VariantFilter variantFilter) throws IOException {

		this.variantFilter = variantFilter;
		
        if (cacheSize <= 0) {
            variantProvider = this;
        } else {
            variantProvider = new CachedSampleVariantProvider(this, cacheSize);
        }
        this.cacheSize = cacheSize;

        baseDir = new File(location);
        if (!baseDir.exists() || !baseDir.isDirectory()) {
            throw new IOException("Folder not found for TriTyper data: " + location);
        }

        genotypeDataFile = new File(baseDir, "GenotypeMatrix.dat");
        if (!genotypeDataFile.exists()) {
            throw new IOException("GenotypeMatrix.dat not found in TriTyper data folder: " + location);
        }


        File tmpFile = new File(baseDir, "ImputedDosageMatrix.dat");
        if (tmpFile.exists()) {
            imputedDosageDataFile = tmpFile;
        } else {
            imputedDosageDataFile = null;
        }

        tmpFile = new File(baseDir, "SNPs.txt");
        if (!tmpFile.exists()) {
            tmpFile = new File(baseDir, "SNPs.txt.gz");
            if (!tmpFile.exists()) {
                throw new IOException("SNPs.txt or SNPs.txt.gz not found in TriTyper data folder: " + location);
            }
        }
        snpFile = tmpFile;

        tmpFile = new File(baseDir, "SNPMappings.txt");
        if (!tmpFile.exists()) {
            tmpFile = new File(baseDir, "SNPMappings.txt.gz");
            if (!tmpFile.exists()) {
                throw new IOException("SNPMappings.txt or SNPMappings.txt.gz not found in TriTyper data folder: " + location);
            }
        }
        snpMapFile = tmpFile;


        tmpFile = new File(baseDir, "Individuals.txt");
        if (!tmpFile.exists()) {
            tmpFile = new File(baseDir, "Individuals.txt.gz");
            if (!tmpFile.exists()) {
                throw new IOException("Individuals.txt or Individuals.txt.gz not found in TriTyper data folder: " + location);
            }
        }
        individualFile = tmpFile;

        tmpFile = new File(baseDir, "PhenotypeInformation.txt");
        if (!tmpFile.exists()) {
            tmpFile = new File(baseDir, "PhenotypeInformation.txt.gz");
            if (!tmpFile.exists()) {
                throw new IOException("PhenotypeInformation.txt or PhenotypeInformation.txt.gz not found in TriTyper data folder: " + location);
            }
        }
        phenotypeAnnotationFile = tmpFile;

        loadSamples();
        samplePhasing = Collections.nCopies(samples.size(), false);
        loadSNPAnnotation();
        checkFileSize();

        // create file handles //
        if (imputedDosageDataFile != null) {
            dosageHandle = new RandomAccessFile(imputedDosageDataFile, "r");
            dosageChannel = dosageHandle.getChannel();
        } else {
            dosageHandle = null;
            dosageChannel = null;
        }

        genotypeHandle = new RandomAccessFile(genotypeDataFile, "r");
        genotypeChannel = genotypeHandle.getChannel();

        sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

    }

    private void loadSamples() throws IOException {
        // load sample list
        int i = 0;
        TextFile t = new TextFile(individualFile, TextFile.R);
        String[] lineElems = t.readLineElemsReturnReference(TextFile.tab);
        samples = new ArrayList<Sample>();
        HashMap<String, Sample> sampleNameToSampleObj = new HashMap<String, Sample>();
        while (lineElems != null) {
            String individual = new String(lineElems[0].getBytes("UTF-8"));
            Sample sample = new Sample(individual, null, null);

            sampleNameToSampleObj.put(individual, sample);
            samples.add(sample);
            i++;
            lineElems = t.readLineElemsReturnReference(TextFile.tab);
        }
        t.close();

        // load annotation
        t = new TextFile(phenotypeAnnotationFile, TextFile.R);

        int numIncluded = 0;
        int numFemale = 0;
        int numMale = 0;
        int numUnknownSex = 0;
        int numCase = 0;
        int numControl = 0;
        int numUnknown = 0;
        int numAnnotated = 0;

        lineElems = t.readLineElemsReturnReference(TextFile.tab);

        HashSet<Sample> visitedSamples = new HashSet<Sample>();
        while (lineElems != null) {
            String individual = lineElems[0];
            Sample sampleObj = sampleNameToSampleObj.get(individual);
            if (sampleObj != null) {
                if (visitedSamples.contains(sampleObj)) {
                    LOG.warn("Sample " + sampleObj.getId() + " may have duplicate annotation in PhenotypeInformation.txt.");
                } else {
                    Boolean includeSample = parseIncludeExcludeStatus(lineElems[2]);
                    if (!includeSample) {
                        numIncluded++;
                    }
                    CaseControlAnnotation caseControlStatus = CaseControlAnnotation.getCaseAnnotationForTriTyper(lineElems[1]);
                    if (caseControlStatus == CaseControlAnnotation.CASE) {
                        numCase++;
                    } else if (caseControlStatus == CaseControlAnnotation.CONTROL) {
                        numControl++;
                    } else {
                        numUnknown++;
                    }
                    SexAnnotation sex = SexAnnotation.getSexAnnotationForTriTyper(lineElems[3]);
                    if (sex == SexAnnotation.FEMALE) {
                        numFemale++;
                    } else if (sex == SexAnnotation.MALE) {
                        numMale++;
                    } else {
                        numUnknownSex++;
                    }
                    sampleObj.putAnnotationValues(GenotypeData.BOOL_INCLUDE_SAMPLE, includeSample);
                    sampleObj.putAnnotationValues(GenotypeData.CASE_CONTROL_SAMPLE_ANNOTATION_NAME, caseControlStatus);
                    sampleObj.putAnnotationValues(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME, sex);
                    visitedSamples.add(sampleObj);
                }
            }
            lineElems = t.readLineElemsReturnReference(TextFile.tab);
        }
        t.close();
        LOG.info("Loaded " + samples.size() + " samples.\n" + visitedSamples.size() + " samples have annotation: " + numIncluded + " included, " + numFemale + " female, " + numMale + " male, " + numUnknownSex + " with unknown sex");

        sampleAnnotationMap = new HashMap<String, SampleAnnotation>(3);
        sampleAnnotationMap.put(GenotypeData.BOOL_INCLUDE_SAMPLE, new SampleAnnotation(BOOL_INCLUDE_SAMPLE, BOOL_INCLUDE_SAMPLE, null, Annotation.Type.BOOLEAN, SampleAnnotation.SampleAnnotationType.OTHER, false));
        sampleAnnotationMap.put(GenotypeData.CASE_CONTROL_SAMPLE_ANNOTATION_NAME, new SampleAnnotation(CASE_CONTROL_SAMPLE_ANNOTATION_NAME, CASE_CONTROL_SAMPLE_ANNOTATION_NAME, null, Annotation.Type.CASECONTROL, SampleAnnotation.SampleAnnotationType.PHENOTYPE, false));
        sampleAnnotationMap.put(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME, new SampleAnnotation(SEX_SAMPLE_ANNOTATION_NAME, SEX_SAMPLE_ANNOTATION_NAME, null, Annotation.Type.SEX, SampleAnnotation.SampleAnnotationType.COVARIATE, false));
    }

    private Boolean parseIncludeExcludeStatus(String status) {
        if (status == null) {
            return false;
        } else if (status.toLowerCase().equals("exclude")) {
            return false;
        } else if (status.toLowerCase().equals("include")) {
            return true;
        } else {
            return false;
        }
    }

    private void loadSNPAnnotation() throws IOException {
        TextFile tf = new TextFile(snpFile, TextFile.R);
        ArrayList<String> allSNPs = tf.readAsArrayList();
        HashSet<String> allSNPHash = new HashSet<String>();
        allSNPHash.addAll(allSNPs);
        tf.close();

        HashMap<String, PosChr> snpToChr = new HashMap<String, PosChr>();
        TextFile tfSNPMap = new TextFile(snpMapFile, TextFile.R);
        String[] tfStrings = tfSNPMap.readLineElems(TextFile.tab);
        while (tfStrings != null) {
            String snp = tfStrings[2];
            if (allSNPHash.contains(snp)) {
                int pos = 0;
                String chr = tfStrings[0];
                try {
                    pos = Integer.parseInt(tfStrings[1]);
                } catch (NumberFormatException e) {
                    LOG.warn("Position defined for " + snp + " on chromosome " + chr + " is not an integer!");
                }
                snpToChr.put(snp, new PosChr(chr, pos));
            }
            tfStrings = tfSNPMap.readLineElems(TextFile.tab);
        }
        tfSNPMap.close();

        snpToIndex = new HashMap<GeneticVariant, Integer>();


        int index = 0;
        int numberOfSNPsWithAnnotation = 0;
        sequences = new HashMap<String, Sequence>();
        for (String snp : allSNPs) {
            PosChr chrPos = snpToChr.get(snp);
            GeneticVariant variant;
            if (chrPos != null) {

                String chr;
                if (!sequences.containsKey(chrPos.chr)) {
                    chr = chrPos.chr;
                    sequences.put(chrPos.chr, new SimpleSequence(chrPos.chr, 0, this));
                } else {
                    chr = sequences.get(chrPos.chr).getName();
                }

				variant = new ReadOnlyGeneticVariantTriTyper(snp, chrPos.getPos(), chr, variantProvider);



                numberOfSNPsWithAnnotation++;
            } else {
                variant = new ReadOnlyGeneticVariantTriTyper(snp, 0, "0", variantProvider);
            }
			
			if(variantFilter == null || variantFilter.doesVariantPassFilter(variant)){
				snps.add(variant);
				snpToIndex.put(variant, index);
			}
			
			index++;
            
        }
        tf.close();

        LOG.info("Loaded " + allSNPs.size() + " SNPs, " + numberOfSNPsWithAnnotation + " have annotation.");
    }

    private void checkFileSize() {
        long expectedfilesize = (long) (snpToIndex.size() * 2) * (long) samples.size();
        long detectedsize = genotypeDataFile.length();
        if (expectedfilesize != detectedsize) {
            throw new GenotypeDataException("Size of GenotypeMatrix.dat does not match size defined by Indivuals.txt and SNPs.txt. Expected size: " + expectedfilesize + " (" + Gpio.humanizeFileSize(expectedfilesize) + ")\tDetected size: " + detectedsize + " (" + Gpio.humanizeFileSize(detectedsize) + ")\tDiff: " + Math.abs(expectedfilesize - detectedsize));
        }

        if (imputedDosageDataFile != null) {
            expectedfilesize = (long) (snpToIndex.size()) * (long) samples.size();
            detectedsize = imputedDosageDataFile.length();

            if (expectedfilesize != detectedsize) {
                throw new GenotypeDataException("Size of ImputedDosageMatrix.dat does not match size defined by Indivuals.txt and SNPs.txt. Expected size: " + expectedfilesize + " (" + Gpio.humanizeFileSize(expectedfilesize) + ")\tDetected size: " + detectedsize + " (" + Gpio.humanizeFileSize(detectedsize) + ")\tDiff: " + Math.abs(expectedfilesize - detectedsize));
            }
        }
    }

    private class PosChr {

        private final String chr;
        private final int pos;

        public PosChr(String chr, int pos) {
            this.chr = chr;
            this.pos = pos;
        }

        public String getChr() {
            return chr;
        }

        public int getPos() {
            return pos;
        }
    }

    @Override
    public List<Alleles> getSampleVariants(GeneticVariant variant) {
        Integer index = snpToIndex.get(variant);
        if (index == null) {
            throw new GenotypeDataException("Variant " + variant.getPrimaryVariantId() + " does not exist.");
        }

        int numIndividuals = samples.size();
        long indexLong = (long) (index) * (numIndividuals * 2);

        // load bytes using NIO
        ByteBuffer buffer = ByteBuffer.allocate(2 * numIndividuals);
        try {
            genotypeChannel.read(buffer, indexLong);
        } catch (IOException e) {
            throw new GenotypeDataException("Could not read bytes from: " + indexLong + " in genotype file " + genotypeDataFile.getAbsolutePath() + " (size: " + genotypeDataFile.length() + ")");
        }

        List<Alleles> alleles = new ArrayList<Alleles>();
        byte[] bufferArr = buffer.array();
        for (int i = 0; i < numIndividuals; i++) {
            int allele2Pos = numIndividuals + i;
            Alleles a = Alleles.createAlleles(TriTyperAlleleAnnotation.convertByteToAllele(bufferArr[i]), TriTyperAlleleAnnotation.convertByteToAllele(bufferArr[allele2Pos]));
            alleles.add(a);
        }

        return alleles;
    }

    @Override
    public float[] getSampleDosage(GeneticVariant variant) {
        // if there is a dosage file, read from there.. if not, conver genotypes.
        // now transcode into dosage..

        // TODO: optimize this step: no need to get ALL alleles.
        float[] genotypes = CalledDosageConvertor.convertCalledAllelesToDosage(variantProvider.getSampleVariants(variant),
                variant.getVariantAlleles(), variant.getRefAllele());
        if (imputedDosageDataFile != null) {
            Integer index = snpToIndex.get(variant);
            if (index == null) {
                throw new GenotypeDataException("Variant " + variant.getPrimaryVariantId() + " does not exist.");
            }
            int numIndividuals = samples.size();
            long indexLong = (long) index * (long) numIndividuals * 1;
            ByteBuffer buffer = ByteBuffer.allocate(numIndividuals);
            try {
                dosageChannel.read(buffer, indexLong);
            } catch (IOException e) {
                throw new GenotypeDataException("Could not read bytes from: " + indexLong
                        + " in genotype file " + genotypeDataFile.getAbsolutePath()
                        + " (size: " + genotypeDataFile.length() + ")");
            }

            byte[] dosageValues = buffer.array();

            boolean takeComplement = false;
            for (int ind = 0; ind < dosageValues.length; ind++) {
                double dosagevalue = ((double) (-Byte.MIN_VALUE + dosageValues[ind])) / 100;
                if (genotypes[ind] == 0 && dosagevalue > 1) {
                    takeComplement = true;
                    break;
                }
                if (genotypes[ind] == 2 && dosagevalue < 1) {
                    takeComplement = true;
                    break;
                }
            }
            if (takeComplement) {
                for (int ind = 0; ind < dosageValues.length; ind++) {
                    byte dosageValue = (byte) (200 - (-Byte.MIN_VALUE + dosageValues[ind]) + Byte.MIN_VALUE);
                    dosageValues[ind] = dosageValue;
                }
            }
            float[] dosageValuesFloat = new float[dosageValues.length];
            for (int i = 0; i < dosageValues.length; i++) {
                dosageValuesFloat[i] = ((float) (-Byte.MIN_VALUE + dosageValues[ i])) / 100;
            }
            return dosageValuesFloat;
        } else {
            return genotypes;
        }
    }

    @Override
    public byte[] getSampleCalledDosage(GeneticVariant variant) {
        byte[] genotypes = CalledDosageConvertor.convertCalledAllelesToCalledDosage(variantProvider.getSampleVariants(variant),
                variant.getVariantAlleles(), variant.getRefAllele());

        return genotypes;
    }

    @Override
    public List<Sample> getSamples() {
        return samples;
    }

    @Override
    public List<Boolean> getSamplePhasing(GeneticVariant variant) {
        return samplePhasing;
    }

    @Override
    public int cacheSize() {
        return cacheSize;
    }

    public void close() {
        try {
            if (imputedDosageDataFile != null) {
                dosageChannel.close();
                dosageHandle.close();
            }
            genotypeChannel.close();
            genotypeHandle.close();
        } catch (IOException e) {
            throw new GenotypeDataException("Could not close file handle to TriTyper file: " + genotypeDataFile);
        }
    }

    @Override
    public int getSampleVariantProviderUniqueId() {
        return sampleVariantProviderUniqueId;
    }

    @Override
    protected Map<String, SampleAnnotation> getSampleAnnotationsMap() {
        return sampleAnnotationMap;
    }

    @Override
    protected Map<String, ? extends Annotation> getVariantAnnotationsMap() {
        return Collections.emptyMap();
    }

    @Override
    public List<String> getSeqNames() {
        return new ArrayList<String>(sequences.keySet());
    }

    @Override
    public Iterable<Sequence> getSequences() {
        return sequences.values();
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
       return snps.getSequencePosVariants(seqName, startPos);
    }

    @Override
    public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
        return snps.getSequenceVariants(seqName);
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
        return snps.getSequenceRangeVariants(seqName, rangeStart, rangeEnd);
    }
}
