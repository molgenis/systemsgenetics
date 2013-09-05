package org.molgenis.genotype.impute2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SampleAnnotation.SampleAnnotationType;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.GeneticVariantTreeSet;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.io.csv.CsvReader;
import org.molgenis.util.tuple.Tuple;

/**
 * GenotypeData for haps/sample files see http://www.shapeit.fr/
 *
 * First run
 * <code>index-haps.sh yourfile.haps<code> to create the tabix index file
 *
 * The two character string 'NA' is treated as missing when encountered in the sample file
 *
 * @author Patrick Deelen
 *
 */
public class Impute2GenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	private final RandomAccessFile hapsFileReader;
	private Map<String, SampleAnnotation> sampleAnnotations = new LinkedHashMap<String, SampleAnnotation>();
	private final int sampleVariantProviderUniqueId;
	private final SampleVariantsProvider sampleVariantProvider;
	private final GeneticVariantTreeSet<GeneticVariant> variants;
	private final LinkedHashMap<GeneticVariant, Long> variantSampleAllelesIndex;
	private final List<Sample> samples;
	private final HashSet<String> sequenceNames;
	private final int byteToReadForSampleAlleles;
	private static final Logger LOGGER = Logger.getLogger(Impute2GenotypeData.class);
	private AllelesAndPhasing last;

	public Impute2GenotypeData(File hapsFile, File sampleFile) throws IOException {
		this(hapsFile, sampleFile, 100);
	}

	public Impute2GenotypeData(File hapsFile, File sampleFile, int cacheSize) throws IOException{
		this(hapsFile, sampleFile, cacheSize, null);
	}
	
	public Impute2GenotypeData(File hapsFile, File sampleFile, String forceSeqName) throws IOException{
		this(hapsFile, sampleFile, 100, forceSeqName);
	}
	
	public Impute2GenotypeData(File hapsFile, File sampleFile, int cacheSize, String forceSeqName)
			throws IOException {
		if (hapsFile == null) {
			throw new IllegalArgumentException("hapsFile is null");
		}
		if (!hapsFile.isFile()) {
			throw new FileNotFoundException("haps file file not found at "
					+ hapsFile.getAbsolutePath());
		}
		if (!hapsFile.canRead()) {
			throw new IOException("Can not read haps file at "
					+ hapsFile.getAbsolutePath());
		}

		if (sampleFile == null) {
			throw new IllegalArgumentException("sampleFile is null");
		}
		if (!sampleFile.isFile()) {
			throw new FileNotFoundException("sample file file not found at "
					+ sampleFile.getAbsolutePath());
		}
		if (!sampleFile.canRead()) {
			throw new IOException("Can not read sample file " + sampleFile.getAbsolutePath());
		}

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
		if (cacheSize > 0) {
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		} else {
			sampleVariantProvider = this;
		}


		loadAnnotations(sampleFile);
		LOGGER.info("Annotations loaded");

		samples = new ArrayList<Sample>();
		loadSamples(sampleFile);


		variants = new GeneticVariantTreeSet<GeneticVariant>();
		variantSampleAllelesIndex = new LinkedHashMap<GeneticVariant, Long>();
		sequenceNames = new HashSet<String>();
		hapsFileReader = new RandomAccessFile(hapsFile, "r");

		byteToReadForSampleAlleles = loadVariants(forceSeqName);

	}

	@Override
	public Iterable<Sequence> getSequences() {
		List<Sequence> sequences = new ArrayList<Sequence>();
		for (String seqName : getSeqNames()) {
			sequences.add(new SimpleSequence(seqName, null, this));
		}

		return sequences;
	}

	@Override
	public List<Sample> getSamples() {
		return Collections.unmodifiableList(samples);
	}

	private void loadSamples(File sampleFile) {

		CsvReader reader = null;
		try {
			reader = new CsvReader(sampleFile, ' ');
			Iterator<Tuple> it = reader.iterator();

			it.next();// Datatype row

			while (it.hasNext()) {
				Tuple tuple = it.next();

				String familyId = tuple.getString(0);
				String sampleId = tuple.getString(1);

				Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
				annotationValues.put("missing", tuple.getDouble(2));

				for (String colName : sampleAnnotations.keySet()) {
					SampleAnnotation annotation = sampleAnnotations.get(colName);

					Object value = null;
					if (!tuple.getString(colName).equalsIgnoreCase("NA")) {
						switch (annotation.getType()) {
							case STRING:
								value = tuple.getString(colName);
								break;
							case INTEGER:
								value = tuple.getInt(colName);
								break;
							case BOOLEAN:
								value = tuple.getBoolean(colName);
								break;
							case FLOAT:
								value = tuple.getDouble(colName);
								break;
							default:
								LOGGER.warn("Unsupported data type encountered for column [" + colName + "]");
						}
					}

					annotationValues.put(colName, value);
				}

				samples.add(new Sample(sampleId, familyId, annotationValues));
			}

		} catch (FileNotFoundException e) {
			throw new RuntimeException("File [" + sampleFile.getAbsolutePath() + "] does not exists", e);
		} finally {
			IOUtils.closeQuietly(reader);
		}

	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return Collections.emptyMap();
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return sampleAnnotations;
	}

	private void loadAnnotations(File sampleFile) throws IOException {
		sampleAnnotations.clear();

		SampleAnnotation missingAnnotation = new SampleAnnotation("missing", "missing",
				"Missing data proportion of each individual", Annotation.Type.FLOAT, SampleAnnotationType.OTHER, false);
		sampleAnnotations.put(missingAnnotation.getId(), missingAnnotation);

		CsvReader reader = null;
		try {
			reader = new CsvReader(sampleFile, ' ');

			List<String> colNames = Utils.iteratorToList(reader.colNamesIterator());
			Tuple dataTypes = reader.iterator().next();
			for (int i = 3; i < colNames.size(); i++) {
				SampleAnnotation annotation = null;
				if (dataTypes.getString(i).equalsIgnoreCase("D")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.STRING,
							SampleAnnotationType.COVARIATE, false);

				} else if (dataTypes.getString(i).equalsIgnoreCase("C")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotationType.COVARIATE, false);
				} else if (dataTypes.getString(i).equalsIgnoreCase("P")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.FLOAT,
							SampleAnnotationType.PHENOTYPE, false);
				} else if (dataTypes.getString(i).equalsIgnoreCase("B")) {
					annotation = new SampleAnnotation(colNames.get(i), colNames.get(i), "", Annotation.Type.BOOLEAN,
							SampleAnnotationType.PHENOTYPE, false);
				} else {
					LOGGER.warn("Unknown datatype [" + dataTypes.getString(i) + "]");
				}

				if (annotation != null) {
					sampleAnnotations.put(annotation.getId(), annotation);
				}
			}
		} finally {
			IOUtils.closeQuietly(reader);
		}

	}

	@Override
	public synchronized List<Alleles> getSampleVariants(GeneticVariant variant) {
				
		if(last == null || !last.getVariant().equals(variant)){
			last = loadAllelesAndPhasing(variant);
		}

		return last.getAlleles();
		
	}

	@Override
	public synchronized List<Boolean> getSamplePhasing(GeneticVariant variant) {
		if(last == null || !last.getVariant().equals(variant)){
			last = loadAllelesAndPhasing(variant);
		}

		return last.getPhasing();
	}

	@Override
	public int cacheSize() {
		return 0;
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

	@Override
	public void close() throws IOException {
		hapsFileReader.close();
	}

	@Override
	public Iterator<GeneticVariant> iterator() {
		return variants.iterator();
	}

	@Override
	public List<String> getSeqNames() {
		return new ArrayList<String>(sequenceNames);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		return variants.getSequencePosVariants(seqName, startPos);
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		return variants.getSequenceVariants(seqName);
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		return variants.getSequenceRangeVariants(seqName, rangeStart, rangeEnd);
	}

	/**
	 * Very complicated low level reading of the variants. Also stores for each
	 * variant the start pos of the sample alleles in the file
	 *
	 * @return the number of bytes of the longest chunk of sample alleles
	 * @throws IOException
	 */
	private int loadVariants(String forceSeqName) throws IOException {
		StringBuilder stringBuilder = new StringBuilder();
		byte[] buffer = new byte[8192];
		boolean eol = false;

		String seqName = null;
		String variantId = null;
		int position = 0;
		String allele1 = null;
		String allele2;

		int longestedChunk = 0;
		int currentChunk = 0;

		int column = 0;
		long posBeforeBufferRead;

		while (!eol) {

			posBeforeBufferRead = hapsFileReader.getFilePointer();
			int bytesRead = hapsFileReader.read(buffer);

			if (bytesRead == -1) {
				eol = true;
				stringBuilder = new StringBuilder();
			} else {
				for (int i = 0; i < bytesRead; ++i) {
					switch (buffer[i]) {
						case '\n':
						case '\r':
							longestedChunk = longestedChunk < currentChunk ? currentChunk : longestedChunk;
							column = 0;
							break;
						case ' ':
							switch (column) {
								case 0:
									seqName = forceSeqName == null ? stringBuilder.toString().intern() : forceSeqName;
									sequenceNames.add(seqName);
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 1:
									variantId = stringBuilder.toString();
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 2:
									try {
										position = Integer.parseInt(stringBuilder.toString());
									} catch (NumberFormatException ex) {
										throw new GenotypeDataException("Error parsing \"" + stringBuilder.toString() + "\" as position in haps file");
									}
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 3:
									allele1 = stringBuilder.toString();
									++column;
									stringBuilder = new StringBuilder();
									break;
								case 4:
									allele2 = stringBuilder.toString();
									GeneticVariant variant = ReadOnlyGeneticVariant.createVariant(variantId, position, seqName, sampleVariantProvider, allele1, allele2);
									variants.add(variant);
									variantSampleAllelesIndex.put(variant, posBeforeBufferRead + i + 1);
									++column;
									stringBuilder = new StringBuilder();
									currentChunk = -1; // because otherwise we would count this space
									break;
								default:
									//We do not care about the rest of the information. We simply seak the end off line
									break;
							}
							break;
						default:
							if (column < 5) {
								//We do not care about the rest of the information. We simply seak the end off line
								stringBuilder.append((char) buffer[i]);
							}
							break;
					}
					++currentChunk;
				}
			}
		}
		return longestedChunk;

	}

	private AllelesAndPhasing loadAllelesAndPhasing(GeneticVariant geneticVariant) {

		ArrayList<Alleles> alleles = new ArrayList<Alleles>();
		ArrayList<Boolean> phasing = new ArrayList<Boolean>();

		Alleles variantAlleles = geneticVariant.getVariantAlleles();

		long start = variantSampleAllelesIndex.get(geneticVariant);

		byte[] buffer = new byte[byteToReadForSampleAlleles];

		int bytesRead;

		try {
			hapsFileReader.seek(start);
			bytesRead = hapsFileReader.read(buffer);
		} catch (IOException ex) {
			throw new GenotypeDataException("Error loading alleles for variant: " + geneticVariant.getPrimaryVariantId() + " from haps file");
		}

		if (bytesRead == -1) {
			throw new GenotypeDataException("Error loading alleles for variant: " + geneticVariant.getPrimaryVariantId() + " from haps file");
		}

		int i = 0;
		for (int s = 0; s < samples.size(); ++s) {

			Allele allele1;
			Allele allele2;
			Boolean phased;

			switch ((char) buffer[i]) {
				case '0':
					allele1 = variantAlleles.get(0);
					break;
				case '1':
					allele1 = variantAlleles.get(1);
					break;
				case '?':
					allele1 = Allele.ZERO;
					break;
				default:
					throw new GenotypeDataException("Error loading alleles for variant: " + geneticVariant.getPrimaryVariantId() + " and sample " + samples.get(s).getId() + " from haps file found allele: " + (char) buffer[i]);
			}
			++i;
			
			switch ((char) buffer[i]) {
				case ' ':
					phased = Boolean.TRUE;
					break;
				case '*':
					phased = Boolean.FALSE;
					++i;//skip space
					break;
				default:
					throw new GenotypeDataException("Error loading alleles for variant: " + geneticVariant.getPrimaryVariantId() + " and sample " + samples.get(s).getId() + " from haps file expected space or * but found: " + (char) buffer[i]);
			}

			++i;
			switch ((char) buffer[i]) {
				case '0':
					allele2 = variantAlleles.get(0);
					break;
				case '1':
					allele2 = variantAlleles.get(1);
					break;
				case '?':
					allele2 = Allele.ZERO;
					break;
				default:
					throw new GenotypeDataException("Error loading alleles for variant: " + geneticVariant.getPrimaryVariantId() + " and sample " + samples.get(s).getId() + " from haps file found allele: " + (char) buffer[i]);
			}
			
			if(!phased){
				++i;
			}

			++i;
			++i;
			
			alleles.add(Alleles.createAlleles(allele1, allele2));
			phasing.add(phased);
			
		}

		return new AllelesAndPhasing(geneticVariant, alleles, phasing);

	}

	private static class AllelesAndPhasing {

		private final GeneticVariant variant;
		private final List<Alleles> alleles;
		private final List<Boolean> phasing;

		public AllelesAndPhasing(GeneticVariant variant, List<Alleles> alleles, List<Boolean> phasing) {
			this.variant = variant;
			this.alleles = alleles;
			this.phasing = phasing;
		}

		public List<Alleles> getAlleles() {
			return alleles;
		}

		public List<Boolean> getPhasing() {
			return phasing;
		}

		public GeneticVariant getVariant() {
			return variant;
		}
	}
}
