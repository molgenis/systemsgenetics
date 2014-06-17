package eqtlmappingpipeline.ase;

import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.vcf.VcfGenotypeData;

/**
 *
 * @author Patrick Deelen
 */
public class ReadCountsLoader implements Runnable {

	private final Iterator<File> inputFileIterator;
	private final AseResults aseResults;
	private final AtomicInteger sampleCounter;
	private final AtomicInteger fileCounter;
	private final AseConfiguration configuration;
	private final RandomAccessGenotypeData genotypeReference;
	Map<String, String> refToStudySampleId;
	private static final Logger LOGGER = Logger.getLogger(ReadCountsLoader.class);

	public ReadCountsLoader(Iterator<File> inputFileIterator, AseResults aseResults, AtomicInteger sampleCounter, AtomicInteger fileCounter, AseConfiguration configuration) {
		this(inputFileIterator, aseResults, sampleCounter, fileCounter, configuration, null);
	}

	public ReadCountsLoader(Iterator<File> inputFileIterator, AseResults aseResults, AtomicInteger sampleCounter, AtomicInteger fileCounter, AseConfiguration configuration, RandomAccessGenotypeData genotypeReference) {
		this(inputFileIterator, aseResults, sampleCounter, fileCounter, configuration, genotypeReference, null);
	}

	public ReadCountsLoader(Iterator<File> inputFileIterator, AseResults aseResults, AtomicInteger sampleCounter, AtomicInteger fileCounter, AseConfiguration configuration, RandomAccessGenotypeData genotypeReference, Map<String, String> refToStudySampleId) {
		this.inputFileIterator = inputFileIterator;
		this.aseResults = aseResults;
		this.sampleCounter = sampleCounter;
		this.fileCounter = fileCounter;
		this.configuration = configuration;
		this.genotypeReference = genotypeReference;
		this.refToStudySampleId = refToStudySampleId == null ? (Map<String, String>) Collections.EMPTY_MAP : refToStudySampleId;
	}

	@Override
	public void run() {

		TObjectIntMap<String> referenceSamples = null;
		try {
			if (genotypeReference != null) {
				String[] referenceSampleNames = genotypeReference.getSampleNames();
				referenceSamples = new TObjectIntHashMap<String>(referenceSampleNames.length, 0.2f, -1);

				int i = 0;
				for (String sample : referenceSampleNames) {

					String translatedSampleId = refToStudySampleId.get(sample);
					if(translatedSampleId != null){
						referenceSamples.put(translatedSampleId, i);
					} else {
						referenceSamples.put(sample, i);
					}

					++i;
				}

			}
		} catch (Exception ex) {
			System.err.println("Error reading samples from reference data: " + ex.getMessage());
			LOGGER.fatal("Error reading samples from reference data.", ex);
			System.exit(1);
		}

		File inputFile = null;

		try {

			while (inputFileIterator.hasNext()) {

				synchronized (inputFileIterator) {
					if (inputFileIterator.hasNext()) {
						inputFile = inputFileIterator.next();
					} else {
						break;
					}
				}

				//Loading genotype files:
				GenotypeData genotypeData;
				if (inputFile.isDirectory()) {
					try {
						genotypeData = MultiPartGenotypeData.createFromVcfFolder(inputFile, 100, 0.8);
					} catch (IncompatibleMultiPartGenotypeDataException ex) {
						System.err.println("Error reading folder with VCF files: " + ex.getMessage());
						LOGGER.fatal("Error reading folder with VCF files: ", ex);
						System.exit(1);
						return;
					}
				} else {
					genotypeData = new VcfGenotypeData(inputFile, 100, 0.8);
				}

				//TODO test if VCF contains the read depth field

				List<Sample> samples = genotypeData.getSamples();
				ArrayList<String> sampleIds = new ArrayList<String>(samples.size());
				for(Sample sample : samples){
					sampleIds.add(sample.getId().intern());
				}

				variants:
				for (GeneticVariant variant : genotypeData) {

					

					//Only if variant contains read depth field
					if (variant.getVariantMeta().getRecordType("AD") == GeneticVariantMeta.Type.INTEGER_LIST) {
						//include variant

						GeneticVariantId variantId = variant.getVariantId();
						
						Iterator<String> sampleIdIterator = sampleIds.iterator();

						List<Alleles> referenceVariantAlleles = null;

						for (GenotypeRecord record : variant.getSampleGenotypeRecords()) {

							String sampleId = sampleIdIterator.next();

							try {
								Alleles alleles;
								if (genotypeReference == null) {
									//Use this VCF to check if sample is hetrozygous for this variant
									alleles = record.getSampleAlleles();

//									if (alleles == null) {
//										throw new AseException("When using VCF file with out GT field you must provide a dataset with genotypes");
//									}

								} else {
									//Use reference VCF to check if sample is hetrozygous for this variant
									if (referenceVariantAlleles == null) {

										synchronized (genotypeReference) {

											GeneticVariant referenceVariant = genotypeReference.getSnpVariantByPos(variant.getSequenceName(), variant.getStartPos());
											if (referenceVariant == null) {
												//LOGGER.debug("Variant not found in reference " + variant.getSequenceName() + ":" + variant.getStartPos());
												continue variants;
											}

											if (!referenceVariant.getVariantAlleles().sameAlleles(variant.getVariantAlleles())) {
												continue variants;
											}
											
											if(variantId == null || variantId == GeneticVariantId.BLANK_GENETIC_VARIANT_ID){
												variantId = referenceVariant.getVariantId();
											}

											referenceVariantAlleles = referenceVariant.getSampleVariants();

										}

									}
									int sampleIndexRef = referenceSamples.get(sampleId);
									if (sampleIndexRef == -1) {
										throw new GenotypeDataException("Sample " + sampleId + " not found in data with reference genotypes");
									}
									alleles = referenceVariantAlleles.get(sampleIndexRef);

								}

								if (alleles != null && (alleles.getAlleleCount() != 2 || alleles.get(0) == alleles.get(1) || alleles.contains(Allele.ZERO))) {
									continue;
								}

								List<Integer> counts = (List<Integer>) record.getGenotypeRecordData("AD");							

								int a1Count = counts.get(0);
								int a2Count = counts.get(1);
								int totalReads = a1Count + a2Count;

								//save as double for divide below
								double minReads = Math.min(a1Count, a2Count);

								if (totalReads >= configuration.getMinTotalReads()
										&& minReads >= configuration.getMinAlleleReads()
										&& totalReads <= configuration.getMaxTotalReads()
										&& minReads / totalReads >= configuration.getMinAlleleReadFraction()) {

									if (variant.getVariantMeta().getRecordType("RQ") == GeneticVariantMeta.Type.FLOAT_LIST) {
										List<Double> meanAlleleBaseQualties = (List<Double>) record.getGenotypeRecordData("RQ");
										aseResults.addResult(variant.getSequenceName(), variant.getStartPos(), variantId, variant.getVariantAlleles().get(0), variant.getVariantAlleles().get(1), a1Count, a2Count, sampleId, meanAlleleBaseQualties.get(0), meanAlleleBaseQualties.get(1));
									} else {
										aseResults.addResult(variant.getSequenceName(), variant.getStartPos(), variantId, variant.getVariantAlleles().get(0), variant.getVariantAlleles().get(1), a1Count, a2Count, sampleId);
									}
										
								}

							} catch (GenotypeDataException ex) {
								System.err.println("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sampleId + " " + ex.getMessage());
								LOGGER.fatal("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sampleId, ex);
								System.exit(1);
								return;
//							} catch (AseException ex) {
//								System.err.println("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sample.getId() + " " + ex.getMessage());
//								LOGGER.fatal("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sample.getId(), ex);
//								System.exit(1);
//								return;
							}


						}

					}

				}

				fileCounter.incrementAndGet();

				genotypeData.close();

				sampleCounter.addAndGet(genotypeData.getSampleNames().length);

			}

		} catch (IOException ex) {
			String inputFilePath = inputFile != null ? inputFile.getAbsolutePath() : "?";
			System.err.println("Error reading input data at " + inputFilePath + " error: " + ex.getMessage());
			LOGGER.fatal("Error reading input data at " + inputFilePath, ex);
			System.exit(1);
		} catch (GenotypeDataException ex) {
			String inputFilePath = inputFile != null ? inputFile.getAbsolutePath() : "?";
			System.err.println("Error reading input data at " + inputFilePath + " error: " + ex.getMessage());
			LOGGER.fatal("Error reading input data at " + inputFilePath, ex);
			System.exit(1);
		}



	}
}
