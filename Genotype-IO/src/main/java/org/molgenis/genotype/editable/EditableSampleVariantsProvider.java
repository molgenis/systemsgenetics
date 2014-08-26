package org.molgenis.genotype.editable;

import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public interface EditableSampleVariantsProvider extends SampleVariantsProvider{

	GeneticVariantMeta getGeneticVariantMeta();
	
}
