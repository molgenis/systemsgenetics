package org.molgenis.genotype.tabix;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertTrue;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.VariantQuery;
import org.molgenis.genotype.VariantQueryResult;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.VariantLineMapper;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.genotype.vcf.VcfReader;
import org.molgenis.genotype.vcf.VcfVariantLineMapper;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
public class TabixQueryTest extends ResourceTest implements SampleVariantsProvider
{
	private VariantQuery query;
	private TabixIndex index;
	private VariantQueryResult result;
	private final int sampleVariantProviderUniqueId = 0;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		VcfReader reader = new VcfReader(new BlockCompressedInputStream(getTestVcfGz()));

		try
		{
			VariantLineMapper variantLineMapper = new VcfVariantLineMapper(reader.getColNames(),
					Collections.<Annotation> emptyList(), Collections.<String, String> emptyMap(), this);
			index = new TabixIndex(getTestVcfGzTbi(), getTestVcfGz(), variantLineMapper);
		}
		finally
		{
			reader.close();
		}

	}

	@BeforeMethod
	public void beforeMethod() throws IOException
	{
		query = index.createQuery();
	}

	@AfterMethod
	public void afterMethod() throws IOException
	{
		if (result != null)
		{
			result.close();
		}
	}

	@Test
	public void queryPos() throws IOException
	{
		result = query.executeQuery("1", 565286, 6097450);
		Iterator<GeneticVariant> it = result.iterator();

		int i = 0;
		while (it.hasNext())
		{
			i++;
			GeneticVariant variant = it.next();
			assertNotNull(variant);
			assertEquals(variant.getSequenceName(), "1");
		}

		// 565286 is excluded, 6097450 is included
		assertEquals(i, 5);
	}

	@Test
	public void querySeq1() throws IOException
	{
		result = query.executeQuery("1");
		Iterator<GeneticVariant> it = result.iterator();

		int i = 0;
		while (it.hasNext())
		{
			i++;
			GeneticVariant variant = it.next();
			assertNotNull(variant);
			assertEquals(variant.getSequenceName(), "1");
		}
		assertEquals(i, 6);
	}

	@Test
	public void querySeq2() throws IOException
	{
		result = query.executeQuery("2");
		Iterator<GeneticVariant> it = result.iterator();

		int i = 0;
		while (it.hasNext())
		{
			i++;
			GeneticVariant variant = it.next();
			assertNotNull(variant);
			assertEquals(variant.getSequenceName(), "2");
		}
		assertEquals(i, 1);
	}

	@Test
	public void queryStartPos()
	{
		result = query.executeQuery("1", 3172273);
		Iterator<GeneticVariant> variants = result.iterator();

		assertNotNull(variants);
		assertTrue(variants.hasNext());
		assertEquals(variants.next().getPrimaryVariantId(), "rs2455100");
		assertFalse(variants.hasNext());

		assertNotNull(query.executeQuery("x", 3172273));
		assertFalse(query.executeQuery("x", 3172273).iterator().hasNext());
		assertNotNull(query.executeQuery("1", 31722730));
		assertFalse(query.executeQuery("1", 31722730).iterator().hasNext());
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)
	{
		return Collections.emptyList();
	}

	@Override
	public int cacheSize()
	{
		return 0;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant)
	{
		return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant)
	{
		return CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
				variant.getVariantAlleles(), variant.getRefAllele());
	}
}
