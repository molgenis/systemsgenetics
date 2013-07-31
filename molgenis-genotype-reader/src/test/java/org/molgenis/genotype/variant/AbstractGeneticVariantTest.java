package org.molgenis.genotype.variant;

import java.util.ArrayList;
import org.molgenis.genotype.Alleles;
import static org.testng.AssertJUnit.assertEquals;

import org.molgenis.genotype.DummySampleVariantsProvider;
import org.molgenis.genotype.util.GeneticVariantTreeSet;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

public class AbstractGeneticVariantTest extends GeneticVariantTreeSet<GeneticVariant>
{

	GeneticVariant variant1;
	GeneticVariant variant2;
	GeneticVariant variant3;
	GeneticVariant variant4;
	GeneticVariant variant5;
	GeneticVariant variant6;
	GeneticVariant variant7;
	GeneticVariant variant8;
	GeneticVariant variant9;
	GeneticVariant variant10;
	GeneticVariant variant11;
	GeneticVariant variant12;

	@BeforeMethod
	public void setup()
	{

		SampleVariantsProvider provider1 = new DummySampleVariantsProvider(new ArrayList<Alleles>());
		SampleVariantsProvider provider2 = new DummySampleVariantsProvider(new ArrayList<Alleles>());

		variant1 = ReadOnlyGeneticVariant.createSnp("rs1", 1, "1", provider1, 'A', 'T');
		variant2 = ReadOnlyGeneticVariant.createSnp("rs2", 1, "1", provider1, 'A', 'T');
		variant3 = ReadOnlyGeneticVariant.createSnp("rs3", 2, "1", provider1, 'A', 'T');
		variant4 = ReadOnlyGeneticVariant.createSnp("rs4", 1, "2", provider1, 'A', 'T');
		variant5 = ReadOnlyGeneticVariant.createSnp("rs5", 1, "1", provider1, 'G', 'T');
		variant6 = ReadOnlyGeneticVariant.createVariant("rs6", 1, "1", provider1, "G", "AA");
		variant7 = ReadOnlyGeneticVariant.createVariant("rs7", 1, "1", provider1, "GG", "T");
		variant8 = ReadOnlyGeneticVariant.createVariant("rs8", 1, "X", provider1, "GG", "T");
		variant9 = ReadOnlyGeneticVariant.createSnp("rs1", 1, "1", provider2, 'A', 'T');
		variant10 = ReadOnlyGeneticVariant.createSnp("rs1", 3, "1", provider2, 'T', 'T');
		variant11 = new GeneticVariantTreeSet.DummyGenticVariant("22", Integer.MIN_VALUE);
		variant12 = new GeneticVariantTreeSet.DummyGenticVariant("22", Integer.MAX_VALUE);

	}

	@Test
	public void compareTo()
	{
		assertEquals(variant1.compareTo(variant3) < 0, true);
		assertEquals(variant3.compareTo(variant1) > 0, true);
		assertEquals(variant1.compareTo(variant4) < 0, true);
		assertEquals(variant4.compareTo(variant1) > 0, true);
		assertEquals(variant3.compareTo(variant4) < 0, true);
		assertEquals(variant4.compareTo(variant3) > 0, true);
		assertEquals(variant3.compareTo(variant8) < 0, true);
		assertEquals(variant8.compareTo(variant3) > 0, true);
		assertEquals(variant1.compareTo(variant10) < 0, true);
		assertEquals(variant10.compareTo(variant1) > 0, true);
		assertEquals(variant11.compareTo(variant12) < 0, true);
//		assertEquals(variant12.compareTo(variant11) > 0, true);
	}

	@Test
	public void equals()
	{

		assertEquals(variant1.equals(variant1), true);
		assertEquals(variant1.equals(variant2), true);
		assertEquals(variant1.equals(variant3), false);
		assertEquals(variant1.equals(variant4), false);
		assertEquals(variant1.equals(variant5), false);
		assertEquals(variant1.equals(variant6), false);
		assertEquals(variant1.equals(variant7), false);
		assertEquals(variant1.equals(variant8), false);
		assertEquals(variant1.equals(variant9), false);
		assertEquals(variant2.equals(variant3), false);
		assertEquals(variant5.equals(variant6), false);
		assertEquals(variant7.equals(variant8), false);
		assertEquals(variant11.equals(variant12), false);
		assertEquals(variant12.equals(variant11), false);

	}
	
	@Test
	public void getHwePvalueTest() {
		
		//Expected values calculated using R package hwde
		
		assertEquals(Double.isNaN(variant1.getHwePvalue()), true);
		
		GeneticVariant variant;
		ArrayList<Alleles> sampleAlleles;
		SampleVariantsProvider sampleAllelesProvider;
				
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(Double.isNaN(variant.getHwePvalue()), true);
		
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		assertEquals(variant.getHwePvalue(), 1d);
		
		
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		assertEquals(variant.getHwePvalue(), 0.04761905d, 0.000001d);
		
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		assertEquals(variant.getHwePvalue(), 0.04761905d, 0.000001d);
		
		
		
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		assertEquals(variant.getHwePvalue(), 0.3650794d, 0.000001d);
		
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		assertEquals(variant.getHwePvalue(), 0.3650794d, 0.000001d);
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', '\0'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C'}));
		assertEquals(variant.getHwePvalue(), 0.3650794d, 0.000001d);
		
	}
	
	@Test
	public void getCallRateTest(){
		
		assertEquals(Double.isNaN(variant1.getCallRate()), true);
		
		GeneticVariant variant;
		ArrayList<Alleles> sampleAlleles;
		SampleVariantsProvider sampleAllelesProvider;
				
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('C', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(variant.getCallRate(), 1d, 0.000001d);
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('\0', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(variant.getCallRate(), 0.8d, 0.000001d);
		
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('\0', '\0'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(variant.getCallRate(), 0.8d, 0.000001d);
		
			sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', '0'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(variant.getCallRate(), 0.8d, 0.000001d);
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('\0', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', '\0'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(variant.getCallRate(), 0.6d, 0.000001d);
		
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('\0', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('\0', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', '\0'));
		sampleAlleles.add(Alleles.createBasedOnChars('\0', '\0'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', 'A'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(variant.getCallRate(), 0.2d, 0.000001d);
		
		
		sampleAlleles = new ArrayList<Alleles>();
		sampleAlleles.add(Alleles.createBasedOnChars('\0', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('\0', 'C'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', '\0'));
		sampleAlleles.add(Alleles.createBasedOnChars('\0', '\0'));
		sampleAlleles.add(Alleles.createBasedOnChars('A', '\0'));
		sampleAllelesProvider = new DummySampleVariantsProvider(sampleAlleles);	
		variant = ReadOnlyGeneticVariant.createVariant("rs", 1, "chr", sampleAllelesProvider, Alleles.createBasedOnChars(new char[]{'A', 'C', 'G'}));
		
		assertEquals(variant.getCallRate(), 0.0d, 0.000001d);
		
	}
}
