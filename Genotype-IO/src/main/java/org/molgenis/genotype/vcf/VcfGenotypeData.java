package org.molgenis.genotype.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.annotation.Nullable;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeDataIndex;
import org.molgenis.genotype.IndexedGenotypeData;
import org.molgenis.genotype.RawLineQueryResult;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.VcfAnnotation;
import org.molgenis.genotype.tabix.TabixIndex;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

import com.google.common.base.Function;
import com.google.common.collect.Lists;
import java.util.regex.Pattern;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.util.ProbabilitiesConvertor;

public class VcfGenotypeData extends IndexedGenotypeData implements SampleVariantsProvider {

	private final GenotypeDataIndex index;
	private final VcfReader reader;
	private Map<String, Annotation> sampleAnnotationsMap;
	private Map<String, String> altDescriptions;
	private final int sampleVariantProviderUniqueId;
	private static final Pattern COMMA_PATTERNS = Pattern.compile(",");

	/**
	 * VCF genotype reader with default cache of 100
	 *
	 * @param bzipVcfFile
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public VcfGenotypeData(File bzipVcfFile) throws FileNotFoundException, IOException {
		this(bzipVcfFile, 100);
	}

	public VcfGenotypeData(File bzipVcfFile, int cacheSize) throws FileNotFoundException, IOException {
		this(bzipVcfFile, new File(bzipVcfFile.getAbsolutePath() + ".tbi"), cacheSize);
	}

	/**
	 * VCF genotype reader with default cache of 100
	 *
	 * @param bzipVcfFile
	 * @param tabixIndexFile
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public VcfGenotypeData(File bzipVcfFile, File tabixIndexFile) throws FileNotFoundException, IOException {
		this(bzipVcfFile, tabixIndexFile, 100);
	}

	public VcfGenotypeData(File vcfFile, File tabixIndexFile, int cacheSize) throws FileNotFoundException,
			IOException {

		if (!vcfFile.isFile()) {
			throw new FileNotFoundException("VCF file not found at " + vcfFile.getAbsolutePath());
		}

		if (!vcfFile.canRead()) {
			throw new IOException("Cannot read VCF file at: " + vcfFile.getAbsolutePath());
		}

		if (!tabixIndexFile.isFile()) {
			throw new TabixFileNotFoundException("VCF tabix file not found at " + tabixIndexFile.getAbsolutePath(), tabixIndexFile.getAbsolutePath());
		}

		if (!tabixIndexFile.canRead()) {
			throw new IOException("Cannot read tabix file for VCF at: " + tabixIndexFile.getAbsolutePath());
		}

		try {
			reader = new VcfReader(new BlockCompressedInputStream(vcfFile));

			try {
				SampleVariantsProvider sampleVariantProvider = cacheSize <= 0 ? this : new CachedSampleVariantProvider(
						this, cacheSize);

				VariantLineMapper variantLineMapper = new VcfVariantLineMapper(reader.getColNames(),
						getVariantAnnotations(), getAltDescriptions(), sampleVariantProvider);
				index = new TabixIndex(tabixIndexFile, vcfFile, variantLineMapper);
			} finally {
				IOUtils.closeQuietly(reader);
			}

		} catch (IOException e) {
			throw new GenotypeDataException(e);
		}
		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	@Override
	public List<Alleles> getSampleVariants(final GeneticVariant variant) {
		List<VcfSampleGenotype> sampleGenotypes = getSampleGenotypes(variant);
		return Lists.transform(sampleGenotypes, new Function<VcfSampleGenotype, Alleles>() {
			@Override
			@Nullable
			public Alleles apply(@Nullable VcfSampleGenotype input) {
				if (input == null) {
					return null;
				}

				return Alleles.createBasedOnString(input.getSamleVariants(variant.getVariantAlleles()
						.getAllelesAsString()));
			}
		});
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		// TODO: is this the correct annotation map that is returned?????
		if (sampleAnnotationsMap == null) {
			List<VcfInfo> infos = reader.getInfos();
			sampleAnnotationsMap = new LinkedHashMap<String, Annotation>(infos.size());

			for (VcfInfo info : infos) {
				sampleAnnotationsMap.put(info.getId(), VcfAnnotation.fromVcfInfo(info));
			}
		}

		return sampleAnnotationsMap;
	}

	@Override
	public List<Sequence> getSequences() {
		List<String> seqNames = getSeqNames();
		Map<String, Integer> seqLengths = getSequenceLengths();

		List<Sequence> sequences = new ArrayList<Sequence>(seqNames.size());
		for (String seqName : seqNames) {
			sequences.add(new SimpleSequence(seqName, seqLengths.get(seqName), this));
		}

		return sequences;
	}

	@Override
	public List<Sample> getSamples() {
		List<String> sampleNames;
		try {
			sampleNames = reader.getSampleNames();
		} catch (IOException e) {
			throw new GenotypeDataException("IOException getting samplenames from VcfReader", e);
		}

		List<Sample> samples = new ArrayList<Sample>(sampleNames.size());
		for (String sampleName : sampleNames) {
			Sample sample = new Sample(sampleName, null, null);
			samples.add(sample);
		}

		return samples;
	}

	/**
	 * Get sequence length by sequence name
	 */
	private Map<String, Integer> getSequenceLengths() {
		List<VcfContig> contigs = reader.getContigs();
		Map<String, Integer> sequenceLengthById = new HashMap<String, Integer>(contigs.size());

		for (VcfContig contig : contigs) {
			sequenceLengthById.put(contig.getId(), contig.getLength());
		}

		return sequenceLengthById;
	}

	private Map<String, String> getAltDescriptions() {
		if (altDescriptions == null) {
			List<VcfAlt> alts = reader.getAlts();
			altDescriptions = new HashMap<String, String>(alts.size());

			for (VcfAlt alt : alts) {
				altDescriptions.put(alt.getId(), alt.getDescription());
			}
		}

		return altDescriptions;
	}

	@Override
	protected GenotypeDataIndex getIndex() {
		return index;
	}

	@Override
	public int cacheSize() {
		return 0;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		List<VcfSampleGenotype> sampleGenotypes = getSampleGenotypes(variant);
		return Lists.transform(sampleGenotypes, new Function<VcfSampleGenotype, Boolean>() {
			@Override
			@Nullable
			public Boolean apply(
					VcfSampleGenotype input) {
				if (input == null) {
					return false;
				}

				return input.getPhasing().get(0);
			}
		});
	}

	private List<VcfSampleGenotype> getSampleGenotypes(GeneticVariant variant) {
		RawLineQueryResult queryResult = index.createRawLineQuery().executeQuery(variant.getSequenceName(),
				variant.getStartPos());

		List<VcfSampleGenotype> genotypes = new ArrayList<VcfSampleGenotype>();
		try {
			for (String line : queryResult) {
				List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
				VcfRecord record = new VcfRecord(line, reader.getColNames());
				if (record.getChrom().equalsIgnoreCase(variant.getSequenceName())
						&& (record.getPos() == variant.getStartPos()) && record.getAlleles().equals(alleles)) {
					List<String> sampleNames = reader.getSampleNames();
					for (String sampleName : sampleNames) {
						VcfSampleGenotype geno = record.getSampleGenotype(sampleName);
						if (geno == null) {
							throw new GenotypeDataException("Missing GT format value for sample ["
									+ sampleName + "]");
						}
						genotypes.add(geno);
					}

				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			IOUtils.closeQuietly(queryResult);
		}

		return genotypes;
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

		RawLineQueryResult queryResult = index.createRawLineQuery().executeQuery(variant.getSequenceName(),
				variant.getStartPos());

		float[] dosages = null;

		try {
			for (String line : queryResult) {
				List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
				VcfRecord record = new VcfRecord(line, reader.getColNames());
				if (record.getChrom().equalsIgnoreCase(variant.getSequenceName())
						&& (record.getPos() == variant.getStartPos()) && record.getAlleles().equals(alleles)) {


					if (record.getFormat().contains("DS")) {

						List<String> sampleNames = reader.getSampleNames();
						dosages = new float[sampleNames.size()];
						int i = 0;
						for (String sampleName : sampleNames) {
							String dosageString = record.getSampleValue(sampleName, "DS");
							if (dosageString == null) {
								throw new GenotypeDataException("Missing DS format value for sample ["
										+ sampleName + "]");
							}
							float sampleDosage;
							try {
								sampleDosage = (Float.parseFloat(dosageString)-2)*-1;
							} catch (NumberFormatException e) {
								throw new GenotypeDataException("Error in sample dosage (DS) value for sample ["
										+ sampleName + "], found value: " + dosageString);
							}
							dosages[i] = sampleDosage;
							++i;
						}

					}

				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			IOUtils.closeQuietly(queryResult);
		}

		if (dosages == null) {
			dosages = CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
					variant.getVariantAlleles(), variant.getRefAllele());
		}

		return dosages;
	}

	@Override
	public void close() throws IOException {
		reader.close();
	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return false;
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		
		RawLineQueryResult queryResult = index.createRawLineQuery().executeQuery(variant.getSequenceName(),
				variant.getStartPos());

		float[][] probs = null;

		try {
			for (String line : queryResult) {
				List<String> alleles = variant.getVariantAlleles().getAllelesAsString();
				VcfRecord record = new VcfRecord(line, reader.getColNames());
				if (record.getChrom().equalsIgnoreCase(variant.getSequenceName())
						&& (record.getPos() == variant.getStartPos()) && record.getAlleles().equals(alleles)) {


					if (record.getFormat().contains("GP")) {

						List<String> sampleNames = reader.getSampleNames();
						probs = new float[sampleNames.size()][3];
						int i = 0;
						for (String sampleName : sampleNames) {
							String probStrings = record.getSampleValue(sampleName, "GP");
							if (probStrings == null) {
								throw new GenotypeDataException("Missing DS format value for sample ["
										+ sampleName + "]");
							}
							
							String[] probStringsArray = COMMA_PATTERNS.split(probStrings);
							
							if(probStringsArray.length != 3){
								throw new GenotypeDataException("Error in sample prob (GP) value for sample ["
											+ sampleName + "], found value: " + probStrings);
							}
													
							for(int j = 0 ; j < 3 ; ++j){
								float sampleProb;
								try {
									sampleProb = Float.parseFloat(probStringsArray[j]);
								} catch (NumberFormatException e) {
									throw new GenotypeDataException("Error in sample prob (GP) value for sample ["
											+ sampleName + "], found value: " + probStrings);
								}
								probs[i][j] = sampleProb;
							}
							
							
							++i;
						}

					}

				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			IOUtils.closeQuietly(queryResult);
		}

		if (probs == null) {
			probs = ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(variant.getSampleDosages());
		}

		return probs;
		
	}
}
