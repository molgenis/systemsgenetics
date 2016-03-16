package org.molgenis.genotype.plink;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.Iterator;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;

public class PedMapGenotypeWriter implements GenotypeWriter {

	private static Logger LOGGER = Logger.getLogger(PedMapGenotypeWriter.class);
	private final GenotypeData genotypeData;
	private final char SEPARATOR = ' ';
	private static final DecimalFormat PHENO_FORMATTER = new DecimalFormat("0.#####");
	private static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	private static final Alleles BI_ALLELIC_MISSING = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
	private int writtenSamplesCounter;
	private int writtenVariantsCounter;
	private int excludedVariantsCounter;

	public PedMapGenotypeWriter(GenotypeData genotypeData) {
		this.genotypeData = genotypeData;
	}

	@Override
	public void write(String basePath) throws IOException, NotASnpException {
		write(new File(basePath + ".ped"), new File(basePath + ".map"));
	}

	public void write(File pedFile, File mapFile) throws IOException, NotASnpException {

		if (pedFile == null) {
			throw new IllegalArgumentException("No ped file specified to write to");
		}
		if (mapFile == null) {
			throw new IllegalArgumentException("No map file specified to write to");
		}
		
		writtenSamplesCounter = 0;
		writtenVariantsCounter = 0;
		excludedVariantsCounter = 0;

		writeMapFile(mapFile);
		writePedFile(pedFile);

		LOGGER.info("PED/MAP plink data write completed.\n"
				+ " - Number of samples: " + writtenSamplesCounter + "\n"
				+ " - Number of SNPs: " + writtenVariantsCounter + "\n"
				+ " - Excluded non biallelic SNPs: " + excludedVariantsCounter);

	}

	private void writeMapFile(File mapFile) throws IOException {
		Utils.createEmptyFile(mapFile, "map");

		BufferedWriter mapFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(mapFile), FILE_ENCODING));

		for (GeneticVariant variant : genotypeData) {

			if (variant.getAlleleCount() > 2 || !variant.isSnp()) {
				LOGGER.warn("Skipping variant: " + variant.getPrimaryVariantId() + ", it is not a biallelic SNP");
				++excludedVariantsCounter;
				continue;
			}

			mapFileWriter.append(FormatPlinkChr.formatChr(variant.getSequenceName()));
			mapFileWriter.append(SEPARATOR);
			mapFileWriter.append(variant.getPrimaryVariantId() == null ? variant.getSequenceName() + ":" + variant.getStartPos() : variant.getPrimaryVariantId());
			mapFileWriter.append(SEPARATOR);
			mapFileWriter.append('0');
			mapFileWriter.append(SEPARATOR);
			mapFileWriter.append(String.valueOf(variant.getStartPos()));
			mapFileWriter.append('\n');

			++writtenVariantsCounter;
		}

		mapFileWriter.close();

	}

	private void writePedFile(File pedFile) throws IOException {
		LOGGER.info("Writing genotype data to: " + pedFile);

		Utils.createEmptyFile(pedFile, "bim");

		BufferedWriter pedFilewriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(pedFile), FILE_ENCODING));;

		int sampleInt = 0;
		for (Sample sample : genotypeData.getSamples()) {

			if (sampleInt == 0) {
				++writtenSamplesCounter;
			}

			pedFilewriter.append(sample.getFamilyId() != null ? sample.getFamilyId() : "0");
			pedFilewriter.append(SEPARATOR);
			pedFilewriter.append(sample.getId());
			pedFilewriter.append(SEPARATOR);
			pedFilewriter.append(sample.getFatherId());
			pedFilewriter.append(SEPARATOR);
			pedFilewriter.append(sample.getMotherId());
			pedFilewriter.append(SEPARATOR);
			pedFilewriter.append(Byte.toString(sample.getSex().getPlinkSex()));
			pedFilewriter.append(SEPARATOR);
			pedFilewriter.append(PHENO_FORMATTER.format(getPhenotype(sample)));

			for (GeneticVariant variant : genotypeData) {
				
				

				Alleles variantAlleles = variant.getVariantAlleles();

				if (variant.getAlleleCount() > 2 || !variant.isSnp()) {
					continue;
				}
					
				Iterator<Alleles> sampleAllelesIterator = variant.getSampleVariants().iterator();

				for (int i = 0; i < sampleInt; ++i) {
					sampleAllelesIterator.next();
				}

				Alleles sampleAlleles = sampleAllelesIterator.next();

				if (sampleAlleles.contains(Allele.ZERO) || sampleAlleles.getAlleleCount() < 2) {
					//set both alleles to missing
					sampleAlleles = BI_ALLELIC_MISSING;
				} else if (sampleAlleles.getAlleleCount() > 2 || !variantAlleles.containsAll(sampleAlleles)) {
					throw new GenotypeDataException("Trying to write alleles " + sampleAlleles.getAllelesAsString() + " for " + variantAlleles + " SNP");
				}

				pedFilewriter.append(SEPARATOR);
				pedFilewriter.append(sampleAlleles.getAlleles().get(0).toString());
				pedFilewriter.append(SEPARATOR);
				pedFilewriter.append(sampleAlleles.getAlleles().get(1).toString());

				
			}

			pedFilewriter.append('\n');

			++sampleInt;

			if (sampleInt % 100 == 0) {
				System.out.println(sampleInt + " samples writen to ped file");
			}


		}

		pedFilewriter.close();		
		
		LOGGER.info("All samples and genotypes written to ped file");

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
