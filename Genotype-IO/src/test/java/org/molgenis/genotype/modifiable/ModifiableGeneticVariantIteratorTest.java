package org.molgenis.genotype.modifiable;

import static org.mockito.Mockito.mock;
import static org.testng.Assert.assertEquals;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

public class ModifiableGeneticVariantIteratorTest
{
	public static final ModifiableGenotypeData dummyModifiableGenotypeData = new ModifiableGenotypeDataInMemory(null);

	private GeneticVariantMeta variantMeta;
	
	@BeforeMethod
	public void setUp() {
		this.variantMeta = mock(GeneticVariantMeta.class);
	}
	
	@Test
	public void createModifiableGeneticVariantIterable()
	{
		GeneticVariant variant1 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs1", 1, "chr1", null, 'A', 'T');
		GeneticVariant variant2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs2", 20, "chr1", null, 'G', 'C');
		ArrayList<GeneticVariant> variants = new ArrayList<GeneticVariant>(2);
		variants.add(variant1);
		variants.add(variant2);

		HashSet<ModifiableGeneticVariant> excludeList = new HashSet<ModifiableGeneticVariant>();

		Iterable<ModifiableGeneticVariant> modifiableVariants = ModifiableGeneticVariantIterator
				.createModifiableGeneticVariantIterable(variants.iterator(), dummyModifiableGenotypeData, excludeList);

		Iterator<ModifiableGeneticVariant> modifiableVariantsIterator = modifiableVariants.iterator();

		assertEquals(modifiableVariantsIterator.next().getOriginalVariant(), variant1);
		assertEquals(modifiableVariantsIterator.next().getOriginalVariant(), variant2);

		assertEquals(modifiableVariantsIterator.hasNext(), false);

	}

	@Test
	public void createModifiableGeneticVariantIterableFiltered()
	{
		GeneticVariant variant1 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs1", 1, "chr1", null, 'A', 'T');
		GeneticVariant variant2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs2", 20, "chr1", null, 'G', 'C');
		ArrayList<GeneticVariant> variants = new ArrayList<GeneticVariant>(2);
		variants.add(variant1);
		variants.add(variant2);

		HashSet<ModifiableGeneticVariant> excludeList = new HashSet<ModifiableGeneticVariant>();
		excludeList.add(new ModifiableGeneticVariant(variant1, dummyModifiableGenotypeData));

		Iterable<ModifiableGeneticVariant> modifiableVariants = ModifiableGeneticVariantIterator
				.createModifiableGeneticVariantIterable(variants.iterator(), dummyModifiableGenotypeData, excludeList);

		Iterator<ModifiableGeneticVariant> modifiableVariantsIterator = modifiableVariants.iterator();

		assertEquals(modifiableVariantsIterator.next().getOriginalVariant(), variant2);

		assertEquals(modifiableVariantsIterator.hasNext(), false);

	}

	@Test
	public void createModifiableGeneticVariantIterableFiltered2()
	{
		GeneticVariant variant1 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs1", 1, "chr1", null, 'A', 'T');
		GeneticVariant variant2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs2", 20, "chr1", null, 'G', 'C');
		ArrayList<GeneticVariant> variants = new ArrayList<GeneticVariant>(2);
		variants.add(variant1);
		variants.add(variant2);

		HashSet<ModifiableGeneticVariant> excludeList = new HashSet<ModifiableGeneticVariant>();
		excludeList.add(new ModifiableGeneticVariant(variant2, dummyModifiableGenotypeData));

		Iterable<ModifiableGeneticVariant> modifiableVariants = ModifiableGeneticVariantIterator
				.createModifiableGeneticVariantIterable(variants.iterator(), dummyModifiableGenotypeData, excludeList);

		Iterator<ModifiableGeneticVariant> modifiableVariantsIterator = modifiableVariants.iterator();

		assertEquals(modifiableVariantsIterator.next().getOriginalVariant(), variant1);

		assertEquals(modifiableVariantsIterator.hasNext(), false);

	}

	@Test
	public void createModifiableGeneticVariantIterableFiltered3()
	{
		GeneticVariant variant1 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs1", 1, "chr1", null, 'A', 'T');
		GeneticVariant variant2 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs2", 20, "chr1", null, 'G', 'C');
		GeneticVariant variant3 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs3", 30, "chr1", null, 'G', 'C');
		GeneticVariant variant4 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs4", 40, "chr1", null, 'G', 'C');
		GeneticVariant variant5 = ReadOnlyGeneticVariant.createSnp(variantMeta, "Rs5", 50, "chr1", null, 'G', 'C');

		ArrayList<GeneticVariant> variants = new ArrayList<GeneticVariant>(2);
		variants.add(variant1);
		variants.add(variant2);
		variants.add(variant3);
		variants.add(variant4);
		variants.add(variant5);

		HashSet<ModifiableGeneticVariant> excludeList = new HashSet<ModifiableGeneticVariant>();
		excludeList.add(new ModifiableGeneticVariant(variant3, dummyModifiableGenotypeData));
		excludeList.add(new ModifiableGeneticVariant(variant2, dummyModifiableGenotypeData));

		Iterable<ModifiableGeneticVariant> modifiableVariants = ModifiableGeneticVariantIterator
				.createModifiableGeneticVariantIterable(variants.iterator(), dummyModifiableGenotypeData, excludeList);

		Iterator<ModifiableGeneticVariant> modifiableVariantsIterator = modifiableVariants.iterator();

		assertEquals(modifiableVariantsIterator.next().getOriginalVariant(), variant1);
		assertEquals(modifiableVariantsIterator.next().getOriginalVariant(), variant4);
		assertEquals(modifiableVariantsIterator.next().getOriginalVariant(), variant5);

		assertEquals(modifiableVariantsIterator.hasNext(), false);

	}

}
