package org.molgenis.genotype.multipart;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;

public class MultiPartGenotypeData implements RandomAccessGenotypeData
{

	private List<Sample> samples = null;
	private static final Pattern VCF_PATTERN = Pattern.compile("vcf.gz$", Pattern.CASE_INSENSITIVE);

	/**
	 * Can map multiple times to same genotype dataset if a genotype dataset
	 * contains multiple sequences
	 * 
	 * seqName, GenotypeData
	 */
	private LinkedHashMap<String, RandomAccessGenotypeData> genotypeDatasets = new LinkedHashMap<String, RandomAccessGenotypeData>();
	private Set<RandomAccessGenotypeData> genotypeDataCollection;

	public MultiPartGenotypeData(Collection<RandomAccessGenotypeData> genotypeDataCollection)
			throws IncompatibleMultiPartGenotypeDataException
	{
		this(new HashSet<RandomAccessGenotypeData>(genotypeDataCollection));
	}

	public MultiPartGenotypeData(HashSet<RandomAccessGenotypeData> genotypeDataCollection)
			throws IncompatibleMultiPartGenotypeDataException
	{
		for (RandomAccessGenotypeData genotypeData : genotypeDataCollection)
		{
			if (samples != null)
			{
				if (!genotypeData.getSamples().equals(samples))
				{
					throw new IncompatibleMultiPartGenotypeDataException(
							"Incompatible multi part genotype data. All files should contain identical samples in same order.");
				}
			}
			else
			{
				samples = genotypeData.getSamples();
			}

			for (String seqName : genotypeData.getSeqNames())
			{
				if (genotypeDatasets.containsKey(seqName))
				{
					throw new IncompatibleMultiPartGenotypeDataException(
							"Incompatible multi part genotype data. A seq/chr can not be present in multiple files.");
				}
				genotypeDatasets.put(seqName, genotypeData);
			}

		}

		this.genotypeDataCollection = genotypeDataCollection;
	}

	/**
	 * Folder with VCF files. Matches all vcf.gz (case insensitive). Can only
	 * handle one file per chr. vcf.gz.tbi should be present. All files must
	 * have the same samples in the same order.
	 * 
	 * @param vcfFolder
	 *            folder with vcf files
	 * @param cacheSize
	 *            size of the cache per vcf file.
	 * @throws IOException
	 * @throws IncompetibleMultiPartGenotypeDataException
	 *             if the datasets are not compatible
	 * @throws Exception
	 *             If multiple files for one chr found
	 */
	public static MultiPartGenotypeData createFromVcfFolder(File vcfFolder, int cacheSize) throws IOException,
			IncompatibleMultiPartGenotypeDataException
	{

		Set<RandomAccessGenotypeData> genotypeDataSets = new HashSet<RandomAccessGenotypeData>();

		if (!vcfFolder.isDirectory())
		{
			throw new IOException("This is not a directory: " + vcfFolder.getAbsolutePath());
		}

		for (File file : vcfFolder.listFiles())
		{
			if (file.isDirectory())
			{
				continue;
			}
			Matcher matcher = VCF_PATTERN.matcher(file.getName());

			if (matcher.matches())
			{
				genotypeDataSets.add(new VcfGenotypeData(file, cacheSize));
			}
		}

		return new MultiPartGenotypeData(genotypeDataSets);

	}

	@Override
	public List<String> getSeqNames()
	{
		return new ArrayList<String>(genotypeDatasets.keySet());
	}

	@Override
	public Iterable<Sequence> getSequences()
	{
		return new MultiPartSequencesIterable(genotypeDatasets.values());
	}

	@Override
	public Sequence getSequenceByName(String name)
	{
		if (genotypeDatasets.containsKey(name))
		{
			return genotypeDatasets.get(name).getSequenceByName(name);
		}
		else
		{
			return null;
		}

	}

	@Override
	public List<Annotation> getVariantAnnotations()
	{
		throw new UnsupportedOperationException(
				"Not yet implemented for multipart genotype data. Feel free to contact the developers");
	}

	@Override
	public Annotation getVariantAnnotation(String annotationId)
	{
		throw new UnsupportedOperationException(
				"Not yet implemented for multipart genotype data. Feel free to contact the developers");
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos)
	{
		if (genotypeDatasets.containsKey(seqName))
		{
			return genotypeDatasets.get(seqName).getVariantsByPos(seqName, startPos);
		}
		else
		{
			return Collections.emptyList();
		}
	}

	@Override
	public GeneticVariant getSnpVariantByPos(String seqName, int startPos)
	{
		if (genotypeDatasets.containsKey(seqName))
		{
			return genotypeDatasets.get(seqName).getSnpVariantByPos(seqName, startPos);
		}
		else
		{
			return null;
		}
	}

	@Override
	public List<Sample> getSamples()
	{
		return Collections.unmodifiableList(samples);
	}

	@Override
	public Iterator<GeneticVariant> iterator()
	{
		return new MultiPartVariantsIterable(genotypeDataCollection).iterator();
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName)
	{
		if (genotypeDatasets.containsKey(seqName))
		{
			return genotypeDatasets.get(seqName).getSequenceGeneticVariants(seqName);
		}
		else
		{
			return Collections.emptyList();
		}
	}

	@Override
	public List<SampleAnnotation> getSampleAnnotations()
	{
		throw new UnsupportedOperationException(
				"Not yet implemented for multipart genotype data. Feel free to contact the developers");
	}

	@Override
	public Annotation getSampleAnnotation(String annotationId)
	{
		throw new UnsupportedOperationException(
				"Not yet implemented for multipart genotype data. Feel free to contact the developers");
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd)
	{
		if (genotypeDatasets.containsKey(seqName))
		{
			return genotypeDatasets.get(seqName).getVariantsByRange(seqName, rangeStart, rangeEnd);
		}
		else
		{
			return null;
		}
	}
}
