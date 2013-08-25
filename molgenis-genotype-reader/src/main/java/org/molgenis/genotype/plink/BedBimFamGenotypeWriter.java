/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.plink;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;

/**
 *
 * @author Patrick Deelen
 */
public class BedBimFamGenotypeWriter implements GenotypeWriter {

	private static final byte MAGIC_NUMBER_1 = 108;
	private static final byte MAGIC_NUMBER_2 = 27;
	private static final byte MODE = 1; //We only write snp major mode
	private static final int HOMOZYGOTE_SECOND_BITMASK = 192;
	private static final int HETEROZYGOTE_BITMASK = 128;
	private static final int MISSING_BIT_MASK = 64;
	public static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	public static final char SEPARATOR = ' ';
	private final GenotypeData genotypeData;
	private int writtenSamplesCounter = 0;
	private int writtenVariantsCounter = 0;
	private int excludedVariantsCounter = 0;

	public BedBimFamGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}
	private static final Logger LOGGER = Logger.getLogger(BedBimFamGenotypeWriter.class);

	@Override
	public void write(String path) throws IOException, NotASnpException {
		write(new File(path + ".bed"), new File(path + ".bim"), new File(path + ".fam"));
	}

	public void write(File bedFile, File bimFile, File famFile) throws IOException {

		if (bedFile == null) {
			throw new IllegalArgumentException("No bed file specified to write to");
		}
		if (bimFile == null) {
			throw new IllegalArgumentException("No bim file specified to write to");
		}
		if (famFile == null) {
			throw new IllegalArgumentException("No fam file specified to write to");
		}

		writeBimFile(bimFile);
		writeFamFile(famFile);
		writeBedFile(bedFile);

		LOGGER.info("Binary plink data write completed.\n"
				+ "Number of samples: " + writtenSamplesCounter + "\n"
				+ "Number of SNPs: " + writtenVariantsCounter + "\n"
				+ "Excluded non biallelic SNPs: " + excludedVariantsCounter);


	}

	private void writeBimFile(File bimFile) throws IOException {
		createEmptyFile(bimFile, "bim");

		BufferedWriter bimFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(bimFile), FILE_ENCODING));

		for (GeneticVariant variant : genotypeData) {

			if (!variant.isBiallelic() || !variant.isSnp()) {
				LOGGER.warn("Skipping variant: " + variant.getPrimaryVariantId() + ", it is not a biallelic SNP");
				++excludedVariantsCounter;
				continue;
			}

			bimFileWriter.append(variant.getSequenceName());
			bimFileWriter.append(SEPARATOR);
			bimFileWriter.append(variant.getPrimaryVariantId());
			bimFileWriter.append(SEPARATOR);
			bimFileWriter.append('0');
			bimFileWriter.append(SEPARATOR);
			bimFileWriter.append(String.valueOf(variant.getStartPos()));
			bimFileWriter.append(SEPARATOR);
			bimFileWriter.append(variant.getVariantAlleles().get(0).toString());
			bimFileWriter.append(SEPARATOR);
			bimFileWriter.append(variant.getVariantAlleles().get(1).toString());
			bimFileWriter.append('\n');

			++writtenVariantsCounter;
		}

		bimFileWriter.close();

	}

	private void writeFamFile(File famFile) throws IOException {
		createEmptyFile(famFile, "fam");

		BufferedWriter famFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(famFile), FILE_ENCODING));

		for (Sample sample : genotypeData.getSamples()) {
			famFileWriter.append(sample.getFamilyId() != null ? sample.getFamilyId() : "0");
			famFileWriter.append(SEPARATOR);
			famFileWriter.append(sample.getId());
			famFileWriter.append(SEPARATOR);
			famFileWriter.append(sample.getFatherId());
			famFileWriter.append(SEPARATOR);
			famFileWriter.append(sample.getMotherId());
			famFileWriter.append(SEPARATOR);
			famFileWriter.append(Byte.toString(sample.getSex().getPlinkSex()));
			famFileWriter.append(SEPARATOR);
			famFileWriter.append(Double.toString(getPhenotype(sample)));
			famFileWriter.append('\n');

			++writtenSamplesCounter;
		}

		famFileWriter.close();

	}

	private void writeBedFile(File bedFile) throws IOException {
		createEmptyFile(bedFile, "bed");

		final DataOutputStream bedStreamWriter;

		try {
			bedStreamWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(bedFile)));
		} catch (FileNotFoundException ex) {
			throw new RuntimeException("This should never happen. File will be present at this time", ex);
		}

		bedStreamWriter.writeByte(MAGIC_NUMBER_1);
		bedStreamWriter.writeByte(MAGIC_NUMBER_2);
		bedStreamWriter.writeByte(MODE);

		for (GeneticVariant variant : genotypeData) {

			if (!variant.isBiallelic() || !variant.isSnp()) {
				continue; //Can only write biallelic snps to binary plink. This is logged when writing bim
			}

			Alleles variantAlleles = variant.getVariantAlleles();

			Alleles homozygoteFirst = Alleles.createAlleles(variantAlleles.get(0), variantAlleles.get(0));
			Alleles homozygoteSecond = Alleles.createAlleles(variantAlleles.get(1), variantAlleles.get(1));

			int currentByte = 0; //Bit operations are on int level
			byte counterCurrentByte = 0;

			for (Alleles alleles : variant.getSampleVariants()) {
				if (alleles == homozygoteFirst) {
					//Do nothing, already 00
				} else if (alleles == homozygoteSecond) {
					currentByte = currentByte | HOMOZYGOTE_SECOND_BITMASK;
				} else if (alleles.sameAlleles(variantAlleles)) {
					currentByte = currentByte | HETEROZYGOTE_BITMASK;
				} else if (alleles.contains(Allele.ZERO)) {
					currentByte = currentByte | MISSING_BIT_MASK;
				} else {
					throw new GenotypeDataException("Trying to write alleles " + alleles.getAllelesAsString() + " for " + variantAlleles + " SNP");
				}
				++counterCurrentByte;
				if (counterCurrentByte == 4) {
					bedStreamWriter.writeByte(currentByte);
					currentByte = 0;
					counterCurrentByte = 0;
				} else {
					currentByte = currentByte >>> 2;
				}
			}

			if (counterCurrentByte != 0) {
				while (counterCurrentByte < 3) {
					++counterCurrentByte;
					currentByte = currentByte >>> 2;
				}
				bedStreamWriter.writeByte(currentByte);
			}

		}
		bedStreamWriter.close();

	}

	private void createEmptyFile(File file, String fileName) {

		if (file.exists()) {
			if (file.isDirectory()) {
				throw new GenotypeDataException("Can not overwrite dir with " + fileName + " file:" + file.getAbsolutePath());
			}
			if (file.isFile()) {
				LOGGER.warn("Overriding " + fileName + " file" + file.getAbsolutePath());
				if (!file.delete()) {
					throw new GenotypeDataException("Failed to overwrite " + fileName + " file: " + file.getAbsolutePath());
				}
			}
		}

		if (!file.getParentFile().exists()) {
			if (!file.getParentFile().mkdirs()) {
				throw new GenotypeDataException("Failed to create parent dir for " + fileName + " file: " + file.getAbsolutePath());
			}
		}

		try {
			if (!file.createNewFile()) {
				throw new GenotypeDataException("Error creating " + fileName + " file: " + file.getAbsolutePath());
			}
		} catch (IOException ex) {
			throw new GenotypeDataException("Error creating " + fileName + " file: " + file.getAbsolutePath(), ex);
		}

		if (!file.canWrite()) {
			throw new GenotypeDataException("Created " + fileName + " file but can not write to file: " + file.getAbsolutePath());
		}

	}

	private double getPhenotype(Sample sample) {

		Object value = sample.getAnnotationValues().get(GenotypeData.DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME);
		if (value == null) {
			return -9;
		}

		if (value instanceof Double) {
			return (Double) value;
		}

		return Double.valueOf(value.toString());
	}
}
