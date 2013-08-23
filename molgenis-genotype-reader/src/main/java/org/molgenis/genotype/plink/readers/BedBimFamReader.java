package org.molgenis.genotype.plink.readers;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.plink.datatypes.BimEntry;
import org.molgenis.genotype.plink.datatypes.FamEntry;
import org.molgenis.genotype.plink.drivers.BedFileDriver;
import org.molgenis.genotype.plink.drivers.BimFileDriver;
import org.molgenis.genotype.plink.drivers.FamFileDriver;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 * Plink binary reader/converter. See:
 * https://github.com/molgenis/molgenis/blob/standalone_tools
 * /src/plinkbintocsv/PlinkbinToCsv.java
 * 
 * @author joeri
 * 
 */
public class BedBimFamReader implements SampleVariantsProvider
{

	private BedFileDriver bedfd;
	private BimFileDriver bimfd;
	private FamFileDriver famfd;

	private long nrOfIndividuals;
	private long nrOfSnps;
	private long nrOfGenotypes;
	private int paddingPerSnp;
	private List<String> individualNames;
	private List<String> snpNames;
	private List<BimEntry> bimEntries;
	private List<FamEntry> famEntries;
	private List<String> sequences; // usually the chromosome numbers, unique
									// within this list
	private HashMap<String, Alleles> snpCoding;

	// helper variables
	private Map<String, Integer> snpIndexById = new HashMap<String, Integer>();
	private final int sampleVariantProviderUniqueId;

	// first lookup is on sequence (usually chromosome, ie. "chr1" ), second on
	// position (basepair, ie. 9345352)
	private Map<String, Map<Long, Integer>> snpIndexByPosition = new HashMap<String, Map<Long, Integer>>();

	// sample phasing
	private final List<Boolean> phasing;

	public BedBimFamReader(File bed, File bim, File fam) throws IOException
	{
		bedfd = new BedFileDriver(bed);
		bimfd = new BimFileDriver(bim);
		famfd = new FamFileDriver(fam);

		nrOfIndividuals = famfd.getNrOfElements();
		nrOfSnps = bimfd.getNrOfElements();
		nrOfGenotypes = nrOfIndividuals * nrOfSnps;
		paddingPerSnp = (int) ((bedfd.getNrOfElements() - nrOfGenotypes) / nrOfSnps);

		phasing = Collections.nCopies((int) nrOfIndividuals, false);

		sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
	}

	public void setIndividuals() throws IOException
	{
		List<String> individualNames = new ArrayList<String>();
		List<FamEntry> famEntries = famfd.getAllEntries();
		if (famEntries.size() != nrOfIndividuals)
		{
			throw new GenotypeDataException("Problem with FAM file: scanned number of elements (" + nrOfIndividuals
					+ ") does not match number of parsed elements (" + famEntries.size() + ")");
		}
		for (FamEntry fe : famEntries)
		{
			if (individualNames.contains(fe.getIndividual()))
			{
				throw new GenotypeDataException("Problem with FAM file: Individual '" + fe.getIndividual() + "' is not unique!");
			}
			individualNames.add(fe.getIndividual());
		}
		this.individualNames = individualNames;
		this.famEntries = famEntries;
	}

	public void setSnps() throws IOException
	{
		nrOfSnps = bimfd.getNrOfElements();
		List<String> snpNames = new ArrayList<String>();
		List<String> uniqueChromosomes = new ArrayList<String>();
		List<BimEntry> bimEntries = bimfd.getAllEntries();
		if (bimEntries.size() != nrOfSnps)
		{
			throw new GenotypeDataException(
					"Problem with BIM file: scanned number of elements does not match number of parsed elements");
		}
		this.bimEntries = bimEntries;
		snpCoding = new HashMap<String, Alleles>();
		int index = 0;
		for (BimEntry be : bimEntries)
		{
			if (snpCoding.containsKey(be.getSNP()))
			{
				throw new GenotypeDataException("Problem with BIM file: SNP '" + be.getSNP() + "' is not unique!");
			}
			snpCoding.put(be.getSNP(), be.getBiallele());
			snpNames.add(be.getSNP());
			if (!uniqueChromosomes.contains(be.getChromosome()))
			{
				uniqueChromosomes.add(be.getChromosome());
			}
			snpIndexById.put(be.getSNP(), index);

			// add to sequence -> position map, create one if missing for this
			// sequence (or 'chromosome')
			if (snpIndexByPosition.get(be.getChromosome()) == null)
			{
				Map<Long, Integer> mapForThisSequence = new HashMap<Long, Integer>();
				snpIndexByPosition.put(be.getChromosome(), mapForThisSequence);
			}
			snpIndexByPosition.get(be.getChromosome()).put(be.getBpPos(), index);

			index++;
		}
		this.snpNames = snpNames;
		this.sequences = uniqueChromosomes;
	}

	public List<GeneticVariant> loadVariantsForSequence(String seq)
	{
		List<GeneticVariant> variants = new ArrayList<GeneticVariant>();

		int index = 0;
		for (BimEntry entry : this.bimEntries)
		{
			String sequenceName = entry.getChromosome();

			if (!sequenceName.equals(seq))
			{
				continue;
			}

			String id = entry.getSNP();
			int startPos = (int) entry.getBpPos();
			Alleles alleles = bimEntries.get(index).getBiallele();
			GeneticVariant snp = ReadOnlyGeneticVariant.createVariant(id, startPos, sequenceName, this, alleles);
			variants.add(snp);
			index++;
		}

		return variants;
	}

	public void extractGenotypes(File writeTo) throws Exception
	{
		setIndividuals();
		setSnps();

		Writer genotypesOut = new OutputStreamWriter(new FileOutputStream(writeTo), "UTF-8");

		// /header: all individual names
		for (String indvName : individualNames)
		{
			genotypesOut.write("\t" + indvName);
		}
		genotypesOut.write("\n");

		// elements: snp name + genotypes
		int snpCounter = 0;
		for (int genotypeStart = 0; genotypeStart < nrOfGenotypes; genotypeStart += nrOfIndividuals)
		{
			System.out.println((int) ((genotypeStart / (double) nrOfGenotypes) * 100) + "% of genotypes done");

			String snpName = snpNames.get(snpCounter);
			String[] allIndividualsForThisSNP = bedfd.getElements(genotypeStart, genotypeStart + nrOfIndividuals,
					paddingPerSnp, snpCounter);
			String a1 = snpCoding.get(snpName).get(0).toString();
			String a2 = snpCoding.get(snpName).get(1).toString();
			String hom1 = a1 + a1;
			String hom2 = a2 + a2;
			String hetr = a1 + a2;

			StringBuilder lineOfGenotypesBuilder = new StringBuilder(snpName);
			for (String s : allIndividualsForThisSNP)
			{
				lineOfGenotypesBuilder.append('\t').append(bedfd.convertGenoCoding(s, hom1, hom2, hetr, ""));
			}
			lineOfGenotypesBuilder.append('\n');
			genotypesOut.write(lineOfGenotypesBuilder.toString());

			snpCounter++;
		}
		genotypesOut.close();
	}

	public List<FamEntry> getFamEntries()
	{
		return famEntries;
	}

	public List<BimEntry> getBimEntries()
	{
		return bimEntries;
	}

	public List<String> getSequences()
	{
		return sequences;
	}

	public Integer getSnpIndexByPosition(String seq, long pos)
	{
		if(snpIndexByPosition.containsKey(seq)){
			return snpIndexByPosition.get(seq).get(pos);
		} else {
			return null;
		}
		
	}

	public static void main(String[] args) throws Exception
	{
		File bed = new File(Alleles.class.getResource("../testfiles/test.bed").getFile());
		File bim = new File(Alleles.class.getResource("../testfiles/test.bim").getFile());
		File fam = new File(Alleles.class.getResource("../testfiles/test.fam").getFile());
		BedBimFamReader bbfr = new BedBimFamReader(bed, bim, fam);
		File out = new File("geno_tmp.txt");
		System.out.println("going to write to: " + out.getAbsolutePath());
		bbfr.extractGenotypes(out);
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant)
	{
		if (variant.getPrimaryVariantId() == null)
		{
			throw new IllegalArgumentException("Not a snp, missing primaryVariantId");
		}

		Integer index = snpIndexById.get(variant.getPrimaryVariantId());

		if (index == null)
		{
			throw new IllegalArgumentException("Unknown primaryVariantId [" + variant.getPrimaryVariantId() + "]");
		}

		List<Alleles> bialleles = new ArrayList<Alleles>();

		try
		{
			String[] allIndividualsForThisSNP = bedfd.getSNPs(index.longValue(), (int) nrOfIndividuals);

			Allele a1 = snpCoding.get(snpNames.get(index)).get(0);
			Allele a2 = snpCoding.get(snpNames.get(index)).get(1);

			for (int i = 0; i < this.nrOfIndividuals; i++)
			{
				Alleles b;
				if (allIndividualsForThisSNP[i].equals("00"))
				{
					b = Alleles.createAlleles(a1, a1);
				}
				else if (allIndividualsForThisSNP[i].equals("01"))
				{
					b = Alleles.createAlleles(a1, a2);
				}
				else if (allIndividualsForThisSNP[i].equals("11"))
				{
					b = Alleles.createAlleles(a2, a2);
				}
				else
				{
					b = Alleles.createAlleles(Allele.ZERO, Allele.ZERO);
				}
				bialleles.add(b);
			}

		}
		catch (Exception e)
		{
			throw new GenotypeDataException(e);
		}

		return bialleles;
	}

	@Override
	public int cacheSize()
	{
		return 0;
	}

	@Override
	public int getSampleVariantProviderUniqueId()
	{
		return sampleVariantProviderUniqueId;
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant)
	{
		return phasing;
	}

	public List<GeneticVariant> loadVariantsForIndex(Integer index)
	{
		if(index == null){
			throw new NullPointerException("Error in index accessing binary plink data");
		}
		BimEntry be = bimEntries.get(index);
		List<GeneticVariant> variants = new ArrayList<GeneticVariant>();
		Alleles alleles = be.getBiallele();
		GeneticVariant snp = ReadOnlyGeneticVariant.createVariant(be.getSNP(), (int) be.getBpPos(), be.getChromosome(),
				this, alleles);
		variants.add(snp);
		return variants;
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
