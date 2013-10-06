package org.molgenis.genotype.oxford;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.List;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 * Export a GenotypeData object to an impute2 haps/sample files
 *
 * Missing values in the sample file are written as 'NA'
 *
 * @author erwin
 *
 */
public class HapsGenotypeWriter implements GenotypeWriter {

	public static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	public static final char LINE_ENDING = '\n';
	private static final char SEPARATOR = ' ';
	private static final char UNPHASED_INDICATOR = '*';
	private static final char MISSING_INDICATOR = '?';
	private static final Logger LOG = Logger.getLogger(HapsGenotypeWriter.class);
	private GenotypeData genotypeData;

	public HapsGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String basePath) throws IOException {
		write(new File(basePath + ".haps"), new File(basePath + ".sample"));
	}

	public void write(File hapsFile, File sampleFile) throws IOException {
		LOG.info("Writing haps file [" + hapsFile.getAbsolutePath() + "] and sample file ["
				+ sampleFile.getAbsolutePath() + "]");

		Utils.createEmptyFile(hapsFile, "haps");
		Utils.createEmptyFile(sampleFile, "sample");
		
		HashMap<Sample, Float> sampleMissingness = writeHapsFile(hapsFile);
		OxfordSampleFileWriter.writeSampleFile(sampleFile, genotypeData, sampleMissingness);

	}

	private HashMap<Sample, Float> writeHapsFile(File hapsFile) throws IOException {
		
		BufferedWriter hapsFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(hapsFile), FILE_ENCODING));
		
		float[] sampleMissingCount = new float[genotypeData.getSamples().size()];
		int totalVariants = 0;
		
		for (GeneticVariant variant : genotypeData) {
			
			++totalVariants;

			if (variant.getAlleleCount() > 2) {
				LOG.warn("Skipping variant: " + variant.getPrimaryVariantId() + " at " + variant.getSequenceName() + ":" + variant.getStartPos() + " with more than 2 alleles: " + variant.getVariantAlleles());
			}

			Allele allele0 = variant.getVariantAlleles().get(0);
			Allele allele1 = variant.getAlleleCount() == 1 ? Allele.ZERO : variant.getVariantAlleles().get(1);

			hapsFileWriter.append(variant.getSequenceName());
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(variant.getPrimaryVariantId());
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(String.valueOf(variant.getStartPos()));
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(allele0.getAlleleAsString());
			hapsFileWriter.append(SEPARATOR);
			hapsFileWriter.append(allele1.getAlleleAsString());

			List<Alleles> sampleAlleles = variant.getSampleVariants();
			List<Boolean> phasing = variant.getSamplePhasing();

			if ((sampleAlleles != null) && !sampleAlleles.isEmpty()) {
				for (int i = 0; i < sampleAlleles.size(); i++) {
					hapsFileWriter.append(SEPARATOR);

					Alleles alleles = sampleAlleles.get(i);
					if ((alleles == null) || alleles.getAllelesAsString().isEmpty() || (alleles.get(0) == Allele.ZERO)
							|| (alleles.get(1) == Allele.ZERO)) {
						hapsFileWriter.append(MISSING_INDICATOR);
						hapsFileWriter.append(SEPARATOR);
						hapsFileWriter.append(MISSING_INDICATOR);
						sampleMissingCount[i]++;
					} else {
						if (alleles.get(0).equals(allele0)) {
							hapsFileWriter.append('0');
						} else if (alleles.get(0).equals(allele1)) {
							hapsFileWriter.append('1');
						} else {
							throw new GenotypeDataException("SampleAllele [" + alleles.get(0) + "] for SNP ["
									+ variant.getPrimaryVariantId() + "] does not match one of the variant alleles " + variant.getVariantAlleles());
						}

						if (!phasing.get(i)) {
							hapsFileWriter.append(UNPHASED_INDICATOR);
						}

						hapsFileWriter.append(SEPARATOR);

						if (alleles.get(1).equals(allele0)) {
							hapsFileWriter.append('0');
						} else if (alleles.get(1).equals(allele1)) {
							hapsFileWriter.append('1');
						} else {
							throw new RuntimeException("SampleAllele [" + alleles.get(1) + "] for SNP ["
									+ variant.getPrimaryVariantId() + "] does not match one of the variant alleles " + variant.getVariantAlleles());
						}

						if (!phasing.get(i)) {
							hapsFileWriter.append(UNPHASED_INDICATOR);
						}
					}
				}
			}

			hapsFileWriter.append(LINE_ENDING);

		}
		
		hapsFileWriter.close();
		
		HashMap<Sample, Float> sampleMissingness = new HashMap<Sample, Float>();
		for(int i = 0 ; i < sampleMissingCount.length ; ++i){
			sampleMissingness.put(genotypeData.getSamples().get(i), sampleMissingCount[i] / (float) totalVariants);
		}
		return sampleMissingness;
		
	}
}
