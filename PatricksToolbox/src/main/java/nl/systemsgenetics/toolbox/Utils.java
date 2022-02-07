/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.toolbox;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterBiAllelic;
import org.molgenis.genotype.variantFilter.VariantFilterMaf;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;

/**
 *
 * @author patri
 */
public class Utils {

	private static final Logger LOGGER = Logger.getLogger(Utils.class);
	
	public static RandomAccessGenotypeData loadGenotypes(Options options, Collection<String> variantsToInclude) throws IOException {
		final RandomAccessGenotypeData referenceGenotypeData;

		final SampleFilter sampleFilter;
		if (options.getGenotypeSamplesFile() != null) {
			sampleFilter = readSampleFile(options.getGenotypeSamplesFile());
		} else {
			sampleFilter = null;
		}

		ArrayList<VariantFilter> variantFilters = new ArrayList<>();
		variantFilters.add(new VariantFilterBiAllelic());

		if (variantsToInclude != null) {
			variantFilters.add(new VariantIdIncludeFilter(new HashSet<>(variantsToInclude)));
		}

		if (options.getMafFilter() != 0) {
			variantFilters.add(new VariantFilterMaf(options.getMafFilter()));
		}

		final VariantFilter variantFilter;
		if (variantFilters.size() == 1) {
			variantFilter = variantFilters.get(0);
		} else {
			variantFilter = new VariantCombinedFilter(variantFilters);
		}

		referenceGenotypeData = options.getGenotypeType().createFilteredGenotypeData(options.getGenotypeBasePath(), 10000, variantFilter, sampleFilter, null, 0.34f);

		return referenceGenotypeData;
	}

	public static SampleIdIncludeFilter readSampleFile(File sampleFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withCSVParser(parser).withSkipLines(0).build();

		final HashSet<String> samples = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			samples.add(nextLine[0]);

		}
		LOGGER.debug("Loaded " + samples.size() + " IDs from: " + sampleFile.getAbsolutePath());

		return new SampleIdIncludeFilter(samples);

	}

}
