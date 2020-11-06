package org.molgenis.genotype.multipart;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.oxford.GenGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;

public class MultiPartGenotypeData extends AbstractRandomAccessGenotypeData
{

	private List<Sample> samples = null;
	private static final Pattern VCF_PATTERN = Pattern.compile(".*vcf\\.gz$", Pattern.CASE_INSENSITIVE);
	private final Map<String, Annotation> variantAnnotationsMap;
	private final Map<String, SampleAnnotation> sampleAnnotationsMap;
	private static final Logger LOGGER = Logger.getLogger(MultiPartGenotypeData.class);

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
		this(new LinkedHashSet<RandomAccessGenotypeData>(genotypeDataCollection));
	}
	
	public MultiPartGenotypeData(RandomAccessGenotypeData... genotypeDataCollection){
		this(new LinkedHashSet<RandomAccessGenotypeData>(Arrays.asList(genotypeDataCollection)));
	}

	public MultiPartGenotypeData(Set<RandomAccessGenotypeData> genotypeDataCollection)
			throws IncompatibleMultiPartGenotypeDataException
	{
		
		variantAnnotationsMap = new LinkedHashMap<>();
		sampleAnnotationsMap = new LinkedHashMap<>();
		
		for (RandomAccessGenotypeData genotypeData : genotypeDataCollection)
		{
			
			String sequenceName = genotypeData.getSeqNames().get(0);
			
			//LOGGER.debug("Started loading chr " + sequenceName + " to multipart data");
			
			if (samples != null)
			{
				if (!genotypeData.getSamples().equals(samples))
				{



					if(genotypeData.getSamples().size() != samples.size()){
						throw new IncompatibleMultiPartGenotypeDataException(
							"Incompatible multi part genotype data. All files should contain identical samples in same order. Number of samples is not identical for chr: "+ sequenceName);
					}

					Iterator<Sample> newSampleIterator = genotypeData.getSamples().iterator();
					Iterator<Sample> totalSampleIterator = samples.iterator();

					while(newSampleIterator.hasNext()){
						Sample newSample = newSampleIterator.next();
						Sample existingSample = totalSampleIterator.next();

						if(!newSample.equals(existingSample)){
							throw new IncompatibleMultiPartGenotypeDataException(
							"Incompatible multi part genotype data. All files should contain identical samples in same order. Found sample: " + newSample + " expected: " + existingSample + " for chr: " + sequenceName);
						}

					}

					throw new IncompatibleMultiPartGenotypeDataException(
							"Incompatible multi part genotype data. All files should contain identical samples in same order. Cause of difference unkown for chr: " + sequenceName);
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
			
			variantAnnotationsMap.putAll(genotypeData.getVariantAnnotationsMap());
			sampleAnnotationsMap.putAll(genotypeData.getSampleAnnotationsMap());

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
	 * @throws IncompatibleMultiPartGenotypeDataException
	 *             if the datasets are not compatible
	 * @throws Exception
	 *             If multiple files for one chr found
	 */
	public static MultiPartGenotypeData createFromVcfFolder(File vcfFolder, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException,
			IncompatibleMultiPartGenotypeDataException
	{

		Set<RandomAccessGenotypeData> genotypeDataSets = new LinkedHashSet<RandomAccessGenotypeData>();

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
				//LOGGER.debug("Adding to multipart data: " + file.getAbsolutePath());
				genotypeDataSets.add(new VcfGenotypeData(file, cacheSize, minimumPosteriorProbabilityToCall));
			} 
		}
		
		if(genotypeDataSets.isEmpty()){
			throw new GenotypeDataException("Did not detect any vcf.gz files at: " + vcfFolder.getAbsolutePath());
		}

		return new MultiPartGenotypeData(genotypeDataSets);

	}
	
	public static MultiPartGenotypeData createFromGenFolder(File genFolder, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException,
			IncompatibleMultiPartGenotypeDataException
	{

		Set<RandomAccessGenotypeData> genotypeDataSets = new LinkedHashSet<RandomAccessGenotypeData>();

		if (!genFolder.isDirectory())
		{
			throw new IOException("This is not a directory: " + genFolder.getAbsolutePath());
		}

		for (File file : genFolder.listFiles())
		{
			if (file.isDirectory())
			{
				continue;
			}
			File sampleFile = new File(file.getAbsolutePath() + ".sample");
			if (sampleFile.exists())
			{
				//LOGGER.debug("Adding to multipart data: " + file.getAbsolutePath());
				genotypeDataSets.add(new GenGenotypeData(file, sampleFile, cacheSize, minimumPosteriorProbabilityToCall));
			} 
		}
		
		if(genotypeDataSets.isEmpty()){
			throw new GenotypeDataException("Did not detect any gen files at: " + genFolder.getAbsolutePath());
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
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd)
	{
		if (genotypeDatasets.containsKey(seqName))
		{
			return genotypeDatasets.get(seqName).getVariantsByRange(seqName, rangeStart, rangeEnd);
		}
		else
		{
			return Collections.emptyList();
		}
	}

	@Override
	public void close() throws IOException {
		for(RandomAccessGenotypeData g : genotypeDataCollection){
			g.close();
		}
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		return variantAnnotationsMap;
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		return sampleAnnotationsMap;
	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		for(RandomAccessGenotypeData g : genotypeDataCollection){
			if(g.isOnlyContaingSaveProbabilityGenotypes() == false){
				return false;
			}
		}
		return true;
	}
}
