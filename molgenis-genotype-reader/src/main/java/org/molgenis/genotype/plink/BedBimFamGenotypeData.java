package org.molgenis.genotype.plink;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.SexAnnotation;
import org.molgenis.genotype.plink.datatypes.FamEntry;
import org.molgenis.genotype.plink.readers.BedBimFamReader;
import org.molgenis.genotype.variant.GeneticVariant;

public class BedBimFamGenotypeData extends AbstractRandomAccessGenotypeData
{

	private final BedBimFamReader reader;
	private Map<String, SampleAnnotation> sampleAnnotations;

	public BedBimFamGenotypeData(File bedFile, File bimFile, File famFile) throws IOException
	{

		if (bedFile == null) {
			throw new IllegalArgumentException("BedFile is null");
		}
		if (bimFile == null) {
			throw new IllegalArgumentException("BimFile is null");
		}
		if (famFile == null) {
			throw new IllegalArgumentException("FamFile is null");
		}

		if (!bedFile.isFile()) {
			throw new FileNotFoundException("BED index file not found at "
+ bedFile.getAbsolutePath());
		}
		if (!bedFile.canRead()) {
			throw new IOException("BED index file not found at " + bedFile.getAbsolutePath());
		}

		if (!bimFile.isFile()) {
			throw new FileNotFoundException("BIM file not found at " + bimFile.getAbsolutePath());
		}
		if (!bimFile.canRead()) throw new IOException("BIM file not found at " + bimFile.getAbsolutePath());

		if (!famFile.isFile()) throw new FileNotFoundException("FAM file not found at " + famFile.getAbsolutePath());
		if (!famFile.canRead()) throw new IOException("FAM file not found at " + famFile.getAbsolutePath());

		this.reader = new BedBimFamReader(bedFile, bimFile, famFile);

		reader.setIndividuals();
		reader.setSnps();

		sampleAnnotations = PlinkSampleAnnotations.getSampleAnnotations();
	}

	@Override
	public List<Sequence> getSequences()
	{
		List<String> seqNames = getSeqNames();

		List<Sequence> sequences = new ArrayList<Sequence>(seqNames.size());
		for (String seqName : seqNames)
		{
			sequences.add(new SimpleSequence(seqName, null, this));
		}

		return sequences;
	}

	@Override
	public List<Sample> getSamples()
	{
		List<Sample> samples = new ArrayList<Sample>();
		for (FamEntry famEntry : reader.getFamEntries())
		{
			Map<String, Object> annotationValues = new LinkedHashMap<String, Object>();
			annotationValues.put(FATHER_SAMPLE_ANNOTATION_NAME, famEntry.getFather());
			annotationValues.put(MOTHER_SAMPLE_ANNOTATION_NAME, famEntry.getMother());
			annotationValues.put(SEX_SAMPLE_ANNOTATION_NAME, SexAnnotation.getSexAnnotationForPlink(famEntry.getSex()));
			annotationValues.put(DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME, famEntry.getPhenotype());

			samples.add(new Sample(famEntry.getIndividual(), famEntry.getFamily(), annotationValues));
		}
		return samples;
	}

	@Override
	protected Map<String, Annotation> getVariantAnnotationsMap()
	{
		return Collections.emptyMap();
	}

	@Override
	public List<String> getSeqNames()
	{
		return this.reader.getSequences();
	}

	@Override
	public List<GeneticVariant> getVariantsByPos(String seqName, int startPos)
	{
		int index = this.reader.getSnpIndexByPosition(seqName, startPos);
		return this.reader.loadVariantsForIndex(index);
	}

	@Override
	public Iterator<GeneticVariant> iterator()
	{
		Iterator<GeneticVariant> it = new Iterator<GeneticVariant>()
		{
			private int currentIndex = 0;

			@Override
			public boolean hasNext()
			{
				return currentIndex < reader.getBimEntries().size();
			}

			@Override
			public GeneticVariant next()
			{
				return reader.loadVariantsForIndex(currentIndex++).get(0);
			}

			@Override
			public void remove()
			{
				throw new GenotypeDataException("Operation 'remove' not implemented for binary Plink iterator");
			}
		};
		return it;
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName)
	{
		List<GeneticVariant> variants = reader.loadVariantsForSequence(seqName);
		if (variants == null)
		{
			throw new IllegalArgumentException("Unknown sequence [" + seqName + "]");
		}
		return variants;
	}

	@Override
	protected Map<String, org.molgenis.genotype.annotation.SampleAnnotation> getSampleAnnotationsMap()
	{
		return sampleAnnotations;
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd)
	{
		throw new UnsupportedOperationException();
	}

}
