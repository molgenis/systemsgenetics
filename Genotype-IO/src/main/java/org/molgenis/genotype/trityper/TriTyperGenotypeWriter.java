package org.molgenis.genotype.trityper;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.plink.PedMapGenotypeWriter;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;

/**
 *
 * @author Patrick Deelen
 */
public class TriTyperGenotypeWriter implements GenotypeWriter {

	private static Logger LOGGER = Logger.getLogger(PedMapGenotypeWriter.class);
	private final GenotypeData genotypeData;

	public TriTyperGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String path) throws IOException, NotASnpException {
		write(new File(path));
	}

	public void write(File folder) throws IOException {

		if (!folder.isDirectory()) {
			if (!folder.mkdirs()) {
				throw new GenotypeDataException("Failed to create trityper dir at: " + folder.getAbsolutePath());
			}
		}

		File genotypeDataFile = new File(folder, "GenotypeMatrix.dat");
		File imputedDosageDataFile = new File(folder, "ImputedDosageMatrix.dat");
		File snpFile = new File(folder, "SNPs.txt");
		File snpMapFile = new File(folder, "SNPMappings.txt");
		File individualFile = new File(folder, "Individuals.txt");
		File phenotypeAnnotationFile = new File(folder, "PhenotypeInformation.txt");

		writeSnps(snpFile, snpMapFile);
		writeSamples(individualFile, phenotypeAnnotationFile);
		writeGenotypes(genotypeDataFile, imputedDosageDataFile);


	}

	private void writeSnps(File snpFile, File snpMapFile) throws IOException {

		BufferedWriter snpFileWriter = new BufferedWriter(new FileWriter(snpFile));
		BufferedWriter snpMapFileWriter = new BufferedWriter(new FileWriter(snpMapFile));

		for (GeneticVariant variant : genotypeData) {

			if (!variant.isSnp()) {
				LOGGER.warn("Skipping variant: " + variant.getPrimaryVariantId() + ", it is not a SNP");
				continue;
			}

			snpFileWriter.append(variant.getPrimaryVariantId());
			snpFileWriter.append('\n');

			snpMapFileWriter.append(variant.getSequenceName());
			snpMapFileWriter.append('\t');
			snpMapFileWriter.append(String.valueOf(variant.getStartPos()));
			snpMapFileWriter.append('\t');
			snpMapFileWriter.append(variant.getPrimaryVariantId());
			snpMapFileWriter.append('\n');
		}

		snpFileWriter.close();
		snpMapFileWriter.close();

	}

	private void writeSamples(File individualFile, File phenotypeAnnotationFile) throws IOException {

		BufferedWriter individualFileWriter = new BufferedWriter(new FileWriter(individualFile));
		BufferedWriter phenotypeAnnotationFileWriter = new BufferedWriter(new FileWriter(phenotypeAnnotationFile));

		for (Sample sample : genotypeData.getSamples()) {

			individualFileWriter.append(sample.getId());
			individualFileWriter.append('\n');

			phenotypeAnnotationFileWriter.append(sample.getId());
			phenotypeAnnotationFileWriter.append('\t');
			phenotypeAnnotationFileWriter.append(sample.getCaseControlAnnotation().getTriTyperName());
			phenotypeAnnotationFileWriter.append('\t');
			phenotypeAnnotationFileWriter.append(sample.isIncluded() ? "include" : "exclude");
			phenotypeAnnotationFileWriter.append('\t');
			phenotypeAnnotationFileWriter.append(sample.getSex().getGender().toLowerCase());
			phenotypeAnnotationFileWriter.append('\n');

		}

		individualFileWriter.close();
		phenotypeAnnotationFileWriter.close();

	}

	private void writeGenotypes(File genotypeDataFile, File imputedDosageDataFile) throws IOException {

		// no need for buffered stream writer. data we write per SNP.
		FileOutputStream genotypeDataFileWriter = new FileOutputStream(genotypeDataFile);
		FileOutputStream genotypeDosageDataFileWriter = new FileOutputStream(imputedDosageDataFile);

		String[] samples = genotypeData.getSampleNames();
		int sampleCount = samples.length;

		byte[] snpBuffer = new byte[sampleCount * 2];
		byte[] dosageBuffer = new byte[sampleCount];

		for (GeneticVariant variant : genotypeData) {

			if (!variant.isSnp()) {
				continue;
			}

			float[] dosageValues = variant.getSampleDosages();
			int i = 0;
			for (Alleles sampleAlleles : variant.getSampleVariants()) {

				if (sampleAlleles.getAlleleCount() != 2) {
					LOGGER.debug("variant at: " + variant.getSequenceName() + ":" + variant.getStartPos() + " set to missing for " + samples[i]);
					sampleAlleles = Alleles.BI_ALLELIC_MISSING;
				}

				try {
					snpBuffer[i] = sampleAlleles.get(0).isSnpAllele() && sampleAlleles.get(0) != Allele.ZERO ? (byte) sampleAlleles.get(0).getAlleleAsSnp() : 0;
					snpBuffer[i + sampleCount] = sampleAlleles.get(1).isSnpAllele() && sampleAlleles.get(1) != Allele.ZERO ? (byte) sampleAlleles.get(1).getAlleleAsSnp() : 0;
				} catch (Exception e) {
					throw new GenotypeDataException("Error writing TriTyper data: " + e.getMessage(), e);
				}

				float dosage = dosageValues[i];
				if (dosage == -1) {
					dosageBuffer[i] = (byte) 127;
				} else {
					int dosageInt = (int) Math.round(dosage * 100d);
					dosageBuffer[i] = (byte) (Byte.MIN_VALUE + dosageInt);
				}

				++i;
			}

			genotypeDataFileWriter.write(snpBuffer);
			genotypeDosageDataFileWriter.write(dosageBuffer);

		}

		genotypeDataFileWriter.close();
		genotypeDosageDataFileWriter.close();

	}
}
