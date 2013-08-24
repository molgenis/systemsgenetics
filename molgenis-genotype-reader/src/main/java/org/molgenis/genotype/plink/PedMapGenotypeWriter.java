package org.molgenis.genotype.plink;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.plink.datatypes.MapEntry;
import org.molgenis.genotype.plink.datatypes.PedEntry;
import org.molgenis.genotype.plink.writers.MapFileWriter;
import org.molgenis.genotype.plink.writers.PedFileWriter;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;

public class PedMapGenotypeWriter implements GenotypeWriter
{
	private static Logger LOG = Logger.getLogger(PedMapGenotypeWriter.class);
	private GenotypeData genotypeData;

	public PedMapGenotypeWriter(GenotypeData genotypeData)
	{
		this.genotypeData = genotypeData;
	}

	public void write(String basePath) throws IOException, NotASnpException
	{
		write(new File(basePath + ".ped"), new File(basePath + ".map"));
	}

	public void write(File pedFile, File mapFile) throws IOException, NotASnpException
	{
		writeMapFile(mapFile);
		writePedFile(pedFile);
	}

	@SuppressWarnings("resource")
	private void writeMapFile(File mapFile) throws IOException, NotASnpException
	{
		LOG.info("Going to create [" + mapFile + "]");
		MapFileWriter writer = null;
		try
		{
			writer = new MapFileWriter(mapFile);
			int count = 0;

			for (GeneticVariant variant : genotypeData)
			{
				if (!variant.isSnp())
				{
					throw new NotASnpException(variant);
				}

				MapEntry mapEntry = new MapEntry(variant.getSequenceName(), variant.getPrimaryVariantId(), 0,
						variant.getStartPos());
				writer.write(mapEntry);
				count++;
				if ((count % 100000) == 0)
				{
					LOG.info("Written " + count + " snps");
				}

			}
			LOG.info("Total written " + count + " snps");

		}
		finally
		{
			IOUtils.closeQuietly(writer);
		}
	}

	private void writePedFile(File pedFile) throws IOException
	{
		LOG.info("Going to create [" + pedFile + "]");

		PedFileWriter writer = null;
		try
		{
			writer = new PedFileWriter(pedFile);
			final List<Sample> samples = genotypeData.getSamples();
			int count = samples.size();

			for (int i = 0; i < count; i++)
			{
				Sample sample = samples.get(i);

				PedEntry pedEntry = new PedEntry(getFamilyId(sample), sample.getId(), getFather(sample),
						getMother(sample), getSex(sample), getPhenotype(sample), new BialleleIterator(genotypeData, i));

				writer.write(pedEntry);
				if ((i % 100) == 0)
				{
					LOG.info("Written " + (i + 1) + "/" + count + " samples");
				}
			}

			LOG.info("All samples written");

		}
		catch (IndexOutOfBoundsException e)
		{
			throw new GenotypeDataException(e);
		}
		finally
		{
			IOUtils.closeQuietly(writer);
		}

	}

	private String getFamilyId(Sample sample)
	{
		return sample.getFamilyId() != null ? sample.getFamilyId() : "0";
	}

	private String getFather(Sample sample)
	{

		Object value = sample.getAnnotationValues().get(GenotypeData.FATHER_SAMPLE_ANNOTATION_NAME);
		if (value == null)
		{
			return "0";
		}

		return value.toString();
	}

	private String getMother(Sample sample)
	{

		Object value = sample.getAnnotationValues().get(GenotypeData.MOTHER_SAMPLE_ANNOTATION_NAME);
		if (value == null)
		{
			return "0";
		}

		return value.toString();
	}

	private byte getSex(Sample sample)
	{

		Object value = sample.getAnnotationValues().get(GenotypeData.SEX_SAMPLE_ANNOTATION_NAME);
		if (value == null)
		{
			return 0;
		}

		SexAnnotation sex = (SexAnnotation) value;

		return sex.getPlinkSex();
	}

	private double getPhenotype(Sample sample)
	{

		Object value = sample.getAnnotationValues().get(GenotypeData.DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME);
		if (value == null)
		{
			return -9;
		}

		if (value instanceof Double)
		{
			return (Double) value;
		}

		return Double.valueOf(value.toString());
	}

	private static class BialleleIterator implements Iterator<Alleles>
	{
		private Iterator<GeneticVariant> variantsIterator;
		private int sampleIndex;

		public BialleleIterator(GenotypeData genotypeData, int sampleIndex)
		{
			this.variantsIterator = genotypeData.iterator();
			this.sampleIndex = sampleIndex;
		}

		@Override
		public boolean hasNext()
		{
			return variantsIterator.hasNext();
		}

		@Override
		public Alleles next()
		{
			GeneticVariant variant = variantsIterator.next();
			Alleles variantAlleles = variant.getSampleVariants().get(sampleIndex);
			return variantAlleles;
		}

		@Override
		public void remove()
		{
		}

	}
}
