package org.molgenis.genotype.impute2;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;

import java.io.IOException;
import java.lang.reflect.Array;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.ResourceTest;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class Impute2GenotypeDataTest extends ResourceTest
{
	private Impute2GenotypeData genotypeData;

	@BeforeClass
	public void beforeClass() throws IOException, URISyntaxException
	{
		genotypeData = new Impute2GenotypeData(getTestImpute2Haps(), getTestImpute2Sample());
	}

	@Test
	public void getSeqNames()
	{
		List<String> seqNames = genotypeData.getSeqNames();
		assertNotNull(seqNames);
		assertEquals(seqNames, Arrays.asList("7", "8"));
	}

	@Test
	public void iterator()
	{
		Iterator<GeneticVariant> it = genotypeData.iterator();
		assertNotNull(it);

		GeneticVariant var = it.next();
		assertNotNull(var);
		assertEquals(var.getPrimaryVariantId(), "SNP1");
		assertEquals(var.getSequenceName(), "7");
		assertEquals(var.getStartPos(), 123);
		assertEquals(var.getAlleleCount(), 2);
		assertEquals(var.getAllIds(), Arrays.asList("SNP1"));
		assertEquals(var.getVariantAlleles(), Alleles.createBasedOnChars('A', 'G'));
		assertEquals(
				var.getSampleVariants(),
				Arrays.asList(Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'A'),
						Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'G')));
		assertEquals(var.getMinorAllele(), Allele.G);
		assertEquals(var.getMinorAlleleFrequency(), 0.375, 0.001);
		assertEquals(Utils.iteratorToList(genotypeData.iterator()).size(), 4);
	}

	@Test
	public void getSamples()
	{
		List<Sample> samples = genotypeData.getSamples();
		assertEquals(samples.size(), 4);

		Sample sample = samples.get(0);
		assertEquals(sample.getId(), "1");
		assertEquals(sample.getFamilyId(), "1");

		Map<String, ?> annotations = sample.getAnnotationValues();
		assertEquals(annotations.size(), 7);
		assertEquals((Double) annotations.get("missing"), 0.007d, 0.0001);
		assertEquals(annotations.get("cov_1"), "1");
		assertEquals(annotations.get("cov_2"), "2");
		assertEquals((Double) annotations.get("cov_3"), 0.0019d, 0.00001);
		assertNull(annotations.get("cov_4"));
		assertEquals((Double) annotations.get("pheno1"), 1.233d, 0.0001);
		assertEquals(annotations.get("bin1"), true);
	}

	@Test
	public void getSequenceByName()
	{
		assertNotNull(genotypeData.getSequenceByName("7"));
		assertNull(genotypeData.getSequenceByName("bogus"));
	}

	@Test
	public void getSequenceGeneticVariants()
	{
		Iterator<GeneticVariant> it = genotypeData.getSequenceGeneticVariants("7").iterator();
		List<GeneticVariant> variants = Utils.iteratorToList(it);
		assertEquals(variants.size(), 3);

		GeneticVariant var = variants.get(0);
		assertEquals(var.getPrimaryVariantId(), "SNP1");
		assertEquals(var.getSequenceName(), "7");
		assertEquals(var.getStartPos(), 123);
		assertEquals(var.getVariantAlleles(), Alleles.createBasedOnChars('A', 'G'));
		assertEquals(
				var.getSampleVariants(),
				Arrays.asList(Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'A'),
						Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'G')));
		
		List<Boolean> expectedPhasing = Arrays.asList(true,false,true,true);
		assertEquals(var.getSamplePhasing(), expectedPhasing);
		
	}

	@Test
	public void getSnpVariantByPos()
	{
		GeneticVariant var = genotypeData.getSnpVariantByPos("7", 789);
		assertNotNull(var);
		assertEquals(var.getPrimaryVariantId(), "SNP3");
		assertEquals(var.getSequenceName(), "7");
		assertEquals(var.getStartPos(), 789);
		assertEquals(var.getVariantAlleles(), Alleles.createBasedOnChars('A', 'T'));
		assertEquals(
				var.getSampleVariants(),
				Arrays.asList(Alleles.createBasedOnChars('A', 'T'), Alleles.createBasedOnChars('T', 'A'),
						Alleles.createBasedOnChars('T', 'T'), Alleles.createBasedOnChars('T', 'T')));
		
		List<Boolean> expectedPhasing = Arrays.asList(true,true,true,true);
		assertEquals(var.getSamplePhasing(), expectedPhasing);
	}
	
	@Test
	public void testForceSeqName() throws URISyntaxException, IOException{
		
		String testHapsPath = getTestImpute2Haps().toString();
		String hapsBasePath = testHapsPath.substring(0, testHapsPath.length() - 5);
		
		RandomAccessGenotypeData genotypeDataForceSeq = RandomAccessGenotypeDataReaderFormats.SHAPEIT2.createGenotypeData(hapsBasePath, 100, "1");
		
		List<String> seqNames = genotypeDataForceSeq.getSeqNames();
		assertNotNull(seqNames);
		assertEquals(seqNames, Arrays.asList("1"));
		
		

		
		Iterator<GeneticVariant> it = genotypeDataForceSeq.iterator();
		assertNotNull(it);

		GeneticVariant var = it.next();
		assertNotNull(var);
		assertEquals(var.getPrimaryVariantId(), "SNP1");
		assertEquals(var.getSequenceName(), "1");
		assertEquals(var.getStartPos(), 123);
		assertEquals(var.getAlleleCount(), 2);
		assertEquals(var.getAllIds(), Arrays.asList("SNP1"));
		assertEquals(var.getVariantAlleles(), Alleles.createBasedOnChars('A', 'G'));
		assertEquals(
				var.getSampleVariants(),
				Arrays.asList(Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'A'),
						Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'G')));
		assertEquals(Utils.iteratorToList(genotypeDataForceSeq.iterator()).size(), 4);
		
		it = genotypeDataForceSeq.getSequenceGeneticVariants("1").iterator();
		List<GeneticVariant> variants = Utils.iteratorToList(it);
		assertEquals(variants.size(), 4);
		
		var = variants.get(0);
		assertEquals(var.getPrimaryVariantId(), "SNP1");
		assertEquals(var.getSequenceName(), "1");
		assertEquals(var.getStartPos(), 123);
		assertEquals(var.getVariantAlleles(), Alleles.createBasedOnChars('A', 'G'));
		assertEquals(
				var.getSampleVariants(),
				Arrays.asList(Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'A'),
						Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'G')));
		
		List<Boolean> expectedPhasing = Arrays.asList(true,false,true,true);
		assertEquals(var.getSamplePhasing(), expectedPhasing);
		
	}
}
