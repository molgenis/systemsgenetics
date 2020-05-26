/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.examples;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleFilterableGenotypeDataDecorator;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterBiAllelic;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantQcChecker;
import org.molgenis.genotype.vcf.VcfGenotypeData;

/**
 * Untested examples
 *
 *
 * @author Patrick Deelen
 */
@edu.umd.cs.findbugs.annotations.SuppressWarnings(value={"DLS_DEAD_LOCAL_STORE", "DMI_HARDCODED_ABSOLUTE_FILENAME"}, justification="It is just an example")
public class FilteringSamplesAndVariants {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {

		/*
		 * Note: Filtering a dataset on variants or samples will create a new
		 * RandomAccessGenotypeData that when accessed will only reviel the included
		 * variants and samples as if the others do not even exist. Data is not copied,
		 * instead the filtered datasets are backed by the orignal
		 */

		String datasetPath = "/some/path/file"; //Note omition of extentions
		
		RandomAccessGenotypeData randomAccessGenotypeData = new VcfGenotypeData(new File(datasetPath), 100, 0.8);

		HashSet<String> includedSnps = new HashSet<String>();
		includedSnps.add("rs1");
		includedSnps.add("rs2");

		//Create filter
		VariantFilter variantFilter = new VariantIdIncludeFilter(includedSnps);
		
		//Filter to only include selected variants
		RandomAccessGenotypeData randomAccessGenotypeDataVariantFilter1 = new VariantFilterableGenotypeDataDecorator(randomAccessGenotypeData, variantFilter);

		HashSet<String> includedSsamples = new HashSet<String>();
		includedSsamples.add("pietje");
		includedSsamples.add("jansje");
		
		SampleFilter sampleFilter = new SampleIdIncludeFilter(includedSsamples);

		//Filter the data with the selected variants to only include the selected samples
		RandomAccessGenotypeData randomAccessGenotypeDataSampleFilter = new SampleFilterableGenotypeDataDecorator(randomAccessGenotypeDataVariantFilter1, sampleFilter);

		//Here we define a combined variant filter consting of 2 filters. This can contain as many filters as requered
		VariantFilter combinedFilter = new VariantCombinedFilter(new VariantQcChecker(0.05f, 0.95f, 0.001d), new VariantFilterBiAllelic());

		//Now we also do some QC of the variants after the samples are filtered. maf 0.05, call rate 0.95 and hwe-p 0.001
		RandomAccessGenotypeData randomAccessGenotypeDataVariantFilter2 = new VariantFilterableGenotypeDataDecorator(randomAccessGenotypeDataSampleFilter, combinedFilter);

		//Note: although it posible to stack variant filters as done above it is more efficient to this at once.


		/*
		 * It is also possible to use the auto loader from the basicUsage examples in combination with filters.
		 * This is faster and more memory efficient for some filetypes, like TriTyper and in the futere 
		 * (or now if this comment is not updated) also for ped_map and binary plink
		 */

		//Note: first samples are filtered and then the variants are filted. 
		RandomAccessGenotypeDataReaderFormats.VCF.createFilteredGenotypeData(datasetPath, 1000, combinedFilter, sampleFilter);

	}
}
