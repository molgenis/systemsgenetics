/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.plink;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.GeneticVariantTreeSet;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class BedBimFamGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	private static final Pattern SEPARATOR_PATTERN = Pattern.compile("[ \\t]+");
	private static final byte MAGIC_NUMBER_1 = 108;
	private static final byte MAGIC_NUMBER_2 = 27;
	private static final byte MODE = 1; //We only write snp major mode
	private static final int READER_MASK = 3;
	private static final int HOMOZYGOTE_FIRST = 0;
	private static final int HOMOZYGOTE_SECOND = 3;
	private static final int HETEROZYGOTE = 2;
	private static final int MISSING = 1;
	private static final Alleles BI_ALLELIC_MISSING = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
	
	private final ArrayList<Sample> samples;
	private final Map<String, SampleAnnotation> sampleAnnotations;
	private final HashMap<String, Sequence> sequences;
	private final GeneticVariantTreeSet<GeneticVariant> snps;
	private final HashMap<GeneticVariant, Integer> snpIndexces;
	private final RandomAccessFile bedFileReader;
	private final SampleVariantsProvider sampleVariantProvider;
	private final int sampleVariantProviderUniqueId;
	private final int cacheSize;
	private final List<Boolean> phasing;
	private final int bytesPerVariant;

	public BedBimFamGenotypeData(String basePath) throws IOException {
		this(basePath, 100);
	}
	
	public BedBimFamGenotypeData(String basePath, int cacheSize) throws IOException {
		this(new File(basePath + ".bed"), new File(basePath + ".bim"), new File(basePath + ".fam"), cacheSize);
	}

	public BedBimFamGenotypeData(File bedFile, File bimFile, File famFile, int cacheSize) throws IOException {

		if (bedFile == null) {
			throw new IllegalArgumentException("BedFile is null");
		}
		if (bimFile == null) {
			throw new IllegalArgumentException("BimFile is null");
		}
		if (famFile == null) {
			throw new IllegalArgumentException("FamFile is null");
		}

		if (!bedFile.isFile()) {
			throw new FileNotFoundException("BED index file not found at "
					+ bedFile.getAbsolutePath());
		}
		if (!bedFile.canRead()) {
			throw new IOException("BED index file not found at " + bedFile.getAbsolutePath());
		}

		if (!bimFile.isFile()) {
			throw new FileNotFoundException("BIM file not found at " + bimFile.getAbsolutePath());
		}
		if (!bimFile.canRead()) {
			throw new IOException("BIM file not found at " + bimFile.getAbsolutePath());
		}

		if (!famFile.isFile()) {
			throw new FileNotFoundException("FAM file not found at " + famFile.getAbsolutePath());
		}
		if (!famFile.canRead()) {
			throw new IOException("FAM file not found at " + famFile.getAbsolutePath());
		}
		
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

		if (cacheSize <= 0) {
			sampleVariantProvider = this;
		} else {
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		}
		this.cacheSize = cacheSize;
		
		sampleAnnotations = PlinkSampleAnnotations.getSampleAnnotations();
		samples = new ArrayList<Sample>();
		readFamFile(famFile);
		
		phasing = Collections.unmodifiableList(Collections.nCopies((int) samples.size(), false));
		
		snpIndexces = new HashMap<GeneticVariant, Integer>();
		snps = new GeneticVariantTreeSet<GeneticVariant>();
		sequences = new HashMap<String, Sequence>();
		readBimFile(bimFile);
		
		bytesPerVariant = samples.size() % 4 == 0 ? samples.size() / 4 : (samples.size() / 4 + 1);
		
		//Check file size of bed file
		if(bedFile.length() != (bytesPerVariant * snpIndexces.size() + 3) ){
			throw new GenotypeDataException("Invalid plink BED file not the expected file size. " + bedFile.getAbsolutePath());
		}
		
		//Check first two bytes for magic number
		bedFileReader = new RandomAccessFile(bedFile, "r");
		if(bedFileReader.read() != MAGIC_NUMBER_1 || bedFileReader.read() != MAGIC_NUMBER_2){
			throw new GenotypeDataException("Error reading plink BED file, magic number not found. " + bedFile.getAbsolutePath());
		}
		
		int bedFileMode = bedFileReader.read();
		if(bedFileMode != MODE){
			if(bedFileMode == 0){
				throw new GenotypeDataException("Error reading BED file, only SNP major mode is supported. " + bedFile.getAbsolutePath());
			} else {
				throw new GenotypeDataException("Error reading BED file, ivalid mode byte detected. " + bedFile.getAbsolutePath());
			}
		}

	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return sampleAnnotations;
	}

	public Map<String, ? extends Annotation> getVariantAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public List<Sample> getSamples() {
		return samples;
	}

	@Override
	public void close() throws IOException {
		bedFileReader.close();
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

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		
		int index = snpIndexces.get(variant);
		
		long startByte = (index * bytesPerVariant) + 3;

		long stopByte = startByte + bytesPerVariant;

		byte[] variantBytes = new byte[(int) (stopByte - startByte)];
		try {
			bedFileReader.seek(startByte);
			if(bedFileReader.read(variantBytes) != variantBytes.length){
				throw new GenotypeDataException("Error reading bed file");
			}
		} catch (IOException ex) {
			throw new GenotypeDataException("Error reading bed file", ex);
		}
		
		ArrayList<Alleles> alleles = new ArrayList<Alleles>(samples.size());
		
		Alleles heterozygote = variant.getVariantAlleles();
		Alleles homozygoteFirst = Alleles.createAlleles(heterozygote.get(0), heterozygote.get(0));
		Alleles homozygoteSecond = Alleles.createAlleles(heterozygote.get(1), heterozygote.get(1));
		
		int sampleCounter = 0;
		
		for(int variantByte : variantBytes){
			
			for(int i = 0 ; i < 4 ; ++i){
				
				if(sampleCounter < samples.size()){
					switch (variantByte & READER_MASK){
						case HOMOZYGOTE_FIRST:
							alleles.add(homozygoteFirst);
							break;
						case HOMOZYGOTE_SECOND:
							alleles.add(homozygoteSecond);
							break;
						case HETEROZYGOTE:
							alleles.add(heterozygote);
							break;
						case MISSING:
							alleles.add(BI_ALLELIC_MISSING);
							break;
						default: 
							throw new GenotypeDataException("Error reading BED, this should not be reachable");
					}
				} else {
					if( (variantByte & READER_MASK) != 0){
						throw new GenotypeDataException("Error reading BED file, found data in padding bits of variant: " + variant.getPrimaryVariantId());
					}
				}
				variantByte = variantByte >>> 2;
				++sampleCounter;
			}
			
		}
		
		
		return Collections.unmodifiableList(alleles);
		
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		return phasing;
	}

	@Override
	public int cacheSize() {
		return cacheSize;
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		return CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	private void readFamFile(File famFile) throws FileNotFoundException, IOException {
		
		BufferedReader famFileReader = new BufferedReader(new FileReader(famFile));
		
		String line;
		while( (line = famFileReader.readLine()) != null ){
			
			String[] elements = SEPARATOR_PATTERN.split(line);
			
			Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
			annotationValues.put(FATHER_SAMPLE_ANNOTATION_NAME, elements[2]);
			annotationValues.put(MOTHER_SAMPLE_ANNOTATION_NAME, elements[3]);
			annotationValues.put(SEX_SAMPLE_ANNOTATION_NAME, SexAnnotation.getSexAnnotationForPlink((byte)elements[4].charAt(0)));
			annotationValues.put(DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME, Double.parseDouble(elements[5]));
			
			samples.add(new Sample(elements[1], elements[0], annotationValues));
			
		}
		
		famFileReader.close();
		
	}

	private void readBimFile(File bimFile) throws FileNotFoundException, IOException {
		
		BufferedReader bimFileReader = new BufferedReader(new FileReader(bimFile));
		
		String line;
		int snpIndex = 0;
		while( (line = bimFileReader.readLine()) != null ){
			
			String[] elements = SEPARATOR_PATTERN.split(line);
			
			String sequenceName = elements[0];
			
			if(!sequences.containsKey(sequenceName)){
				sequences.put(sequenceName, new SimpleSequence(sequenceName, 0, this));
			}
			
			GeneticVariant variant = ReadOnlyGeneticVariant.createVariant(elements[1], Integer.parseInt(elements[3]), sequenceName, sampleVariantProvider, elements[4], elements[5]);
			
			snps.add(variant);
			
			snpIndexces.put(variant, snpIndex);
			
			++snpIndex;
			
		}
		
		bimFileReader.close();
		
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return snps.iterator();
	}
	
	

}