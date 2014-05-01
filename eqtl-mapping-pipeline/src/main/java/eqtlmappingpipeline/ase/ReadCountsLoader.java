/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.ase;

import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.log4j.Logger;
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
	private static final Logger LOGGER = Logger.getLogger(ReadCountsLoader.class);

	public ReadCountsLoader(Iterator<File> inputFileIterator, AseResults aseResults, AtomicInteger sampleCounter, AtomicInteger fileCounter, AseConfiguration configuration) {
		this(inputFileIterator, aseResults, sampleCounter, fileCounter, configuration, null);
	}

	public ReadCountsLoader(Iterator<File> inputFileIterator, AseResults aseResults, AtomicInteger sampleCounter, AtomicInteger fileCounter, AseConfiguration configuration, RandomAccessGenotypeData genotypeReference) {
		this.inputFileIterator = inputFileIterator;
		this.aseResults = aseResults;
		this.sampleCounter = sampleCounter;
		this.fileCounter = fileCounter;
		this.configuration = configuration;
		this.genotypeReference = genotypeReference;
	}

	@Override
	public void run() {

		try {

			while (inputFileIterator.hasNext()) {

				File inputFile;
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
						genotypeData = MultiPartGenotypeData.createFromVcfFolder(inputFile, 100);
					} catch (IncompatibleMultiPartGenotypeDataException ex) {
						System.err.println("Error reading folder with VCF files: " + ex.getMessage());
						LOGGER.fatal("Error reading folder with VCF files: ", ex);
						System.exit(1);
						return;
					}
				} else {
					genotypeData = new VcfGenotypeData(inputFile, 100);
				}

				//TODO test if VCF contains the read depth field

				List<Sample> samples = genotypeData.getSamples();

				TObjectIntMap<String> referenceSamples = null;
				if (genotypeReference != null) {
					String[] referenceSampleNames = genotypeReference.getSampleNames();
					referenceSamples = new TObjectIntHashMap<String>(referenceSampleNames.length, 0.2f, -1);

					int i = 0;
					for (String sample : referenceSampleNames) {
						referenceSamples.put(sample, i);
						++i;
					}

				}

				for (GeneticVariant variant : genotypeData) {


					//Only if variant contains read depth field
					if (variant.getVariantMeta().getRecordType("AD") == GeneticVariantMeta.Type.INTEGER_LIST) {
						//include variant

						Iterator<Sample> sampleIterator = samples.iterator();

						List<Alleles> referenceVariantAlleles = null;

						for (GenotypeRecord record : variant.getSampleGenotypeRecords()) {

							Sample sample = sampleIterator.next();

							if (!record.containsGenotypeRecord("AD")) {
								continue;
							}

							try {
								Alleles alleles;
								if (genotypeReference == null) {
									//Use this VCF to check if sample is hetrozygous for this variant
									alleles = record.getSampleAlleles();

								} else {
									//Use reference VCF to check if sample is hetrozygous for this variant
									if (referenceVariantAlleles == null) {

										synchronized (genotypeReference) {

											GeneticVariant referenceVariant = genotypeReference.getSnpVariantByPos(variant.getSequenceName(), variant.getStartPos());
											if(referenceVariant == null){
												//LOGGER.debug("Variant not found in reference " + variant.getSequenceName() + ":" + variant.getStartPos());
												continue;
											}

											if (!referenceVariant.getVariantAlleles().sameAlleles(variant.getVariantAlleles())) {
												continue;
											}

											referenceVariantAlleles = referenceVariant.getSampleVariants();

										}

									}
									int sampleIndexRef = referenceSamples.get(sample.getId());
									if(sampleIndexRef == -1){
										throw new GenotypeDataException("Sample " + sample.getId() + " not found in data with refernece genotypes");
									}
									alleles = referenceVariantAlleles.get(sampleIndexRef);

								}

								if (alleles.getAlleleCount() != 2 || alleles.get(0) == alleles.get(1)) {
									continue;
								}

								List<Integer> counts = (List<Integer>) record.getGenotypeRecordData("AD");

								int a1Count = counts.get(0);
								int a2Count = counts.get(1);

								if (a1Count + a2Count > configuration.getMinTotalReads()
										&& (a1Count >= configuration.getMinAlleleReads()
										|| a2Count >= configuration.getMinAlleleReads())) {

									aseResults.addResult(variant.getSequenceName(), variant.getStartPos(), variant.getVariantId(), variant.getVariantAlleles().get(0), variant.getVariantAlleles().get(1), a1Count, a2Count);

								}

							} catch (GenotypeDataException ex) {
								System.err.println("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sample.getId() + " " + ex.getMessage());
								LOGGER.fatal("Error parsing " + variant.getSequenceName() + ":" + variant.getStartPos() + " for sample " + sample.getId(), ex);
								System.exit(1);
								return;
							}


						}

					}

				}

				fileCounter.incrementAndGet();

				genotypeData.close();

				sampleCounter.addAndGet(genotypeData.getSampleNames().length);

			}

		} catch (IOException ex) {
			System.err.println("Error reading input data: " + ex.getMessage());
			LOGGER.fatal("Error reading input data", ex);
			System.exit(1);
		} catch (GenotypeDataException ex) {
			System.err.println("Error reading input data: " + ex.getMessage());
			LOGGER.fatal("Error reading input data", ex);
			System.exit(1);
		}



	}
}
