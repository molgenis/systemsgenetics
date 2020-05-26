package org.molgenis.genotype;

import java.util.List;

public interface GenotypeDataIndex
{
	List<String> getSeqNames();

	VariantQuery createQuery();

	RawLineQuery createRawLineQuery();
}
