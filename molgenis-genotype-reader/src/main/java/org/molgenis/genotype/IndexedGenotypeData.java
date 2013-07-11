package org.molgenis.genotype;

import java.util.List;

import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;

public abstract class IndexedGenotypeData extends AbstractRandomAccessGenotypeData
{

	@Override
	public List<String> getSeqNames()
	{
		return getIndex().getSeqNames();
	}

	@Override
	public List<GeneticVariant> getVariantsByPos(String seqName, int startPos)
	{
		VariantQueryResult result = getIndex().createQuery().executeQuery(seqName, startPos);
		try
		{
			return Utils.iteratorToList(result.iterator());
		}
		finally
		{
			IOUtils.closeQuietly(result);
		}
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName)
	{
		return getIndex().createQuery().executeQuery(seqName);
	}

	protected abstract GenotypeDataIndex getIndex();

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd)
	{
		// query is start exclusive and end inclusive. So -1
		VariantQueryResult result = getIndex().createQuery().executeQuery(seqName, rangeStart - 1, rangeEnd - 1);
		try
		{
			return Utils.iteratorToList(result.iterator());
		}
		finally
		{
			IOUtils.closeQuietly(result);
		}
	}
}
