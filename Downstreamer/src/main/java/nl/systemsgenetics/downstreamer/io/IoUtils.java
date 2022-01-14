package nl.systemsgenetics.downstreamer.io;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.gene.Gene;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterMaf;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;

import java.io.*;
import java.util.*;
import nl.systemsgenetics.downstreamer.containers.LeadVariant;
import org.molgenis.genotype.variantFilter.VariantFilterBiAllelic;
import umcg.genetica.collections.ChrPosTreeMap;

public class IoUtils {

	private static final Logger LOGGER = Logger.getLogger(IoUtils.class);

	private static Map<String, ChrPosTreeMap<LeadVariant>> leadVariantsPerTraitCache = null;

	/**
	 * Load genotype data matching GWAS matrix and MAF filter.
	 *
	 * @param options Depict options object
	 * @return RandomAccesGenotypeData for all SNPs in GWAS matrix and MAF
	 * @throws IOException
	 */
	public static RandomAccessGenotypeData readReferenceGenotypeDataMatchingGwasSnps(DownstreamerOptions options) throws IOException {
		return readReferenceGenotypeDataMatchingGwasSnps(options, null);
	}

	/**
	 * Load genotype data matching GWAS matrix and MAF filter.
	 *
	 * @param options Depict options object
	 * @return RandomAccesGenotypeData for all SNPs in GWAS matrix and MAF
	 * @throws IOException
	 */
	public static RandomAccessGenotypeData readReferenceGenotypeDataMatchingGwasSnps(DownstreamerOptions options, Set<String> variantSubset) throws IOException {

		final List<String> variantsInZscoreMatrix;
		if (variantSubset == null) {
			variantsInZscoreMatrix = IoUtils.readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".rows.txt"));
		} else {
			variantsInZscoreMatrix = new ArrayList<>(variantSubset);
		}

		final List<String> phenotypesInZscoreMatrix = IoUtils.readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".cols.txt"));

		LOGGER.info("Number of phenotypes in GWAS matrix: " + Downstreamer.LARGE_INT_FORMAT.format(phenotypesInZscoreMatrix.size()));
		LOGGER.info("Number of variants in GWAS matrix: " + Downstreamer.LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));

		if (options.getVariantFilterFile() != null) {
			HashSet<String> variantsToInclude = IoUtils.readVariantFilterFile(options.getVariantFilterFile());
			Iterator<String> variantsInZscoreMatrixIt = variantsInZscoreMatrix.iterator();
			while (variantsInZscoreMatrixIt.hasNext()) {
				String variant = variantsInZscoreMatrixIt.next();
				if (!variantsToInclude.contains(variant)) {
					variantsInZscoreMatrixIt.remove();
				}
			}
			LOGGER.info("Number of variants after filtering on selected variants: " + Downstreamer.LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));
		}

		return IoUtils.loadGenotypes(options, variantsInZscoreMatrix);
	}

	public static final List<String> readMatrixAnnotations(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).build();

		ArrayList<String> identifiers = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			identifiers.add(nextLine[0]);
		}

		return identifiers;
	}

	public static List<Gene> readGenes(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

		final ArrayList<Gene> genes = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]));

		}

		return genes;

	}

	public static IntervalTreeMap<Gene> readGenesAsIntervalTree(File geneFile) throws Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

		IntervalTreeMap<Gene> intervalTree = new IntervalTreeMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			Gene gene = new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]);
			intervalTree.put(gene, gene);
		}
//		
//		LOGGER.info("Loaded genes in interval tree:" + intervalTree.size());
//		
//		for(Gene gene : intervalTree.getOverlapping(new Interval("1", 67457690, 67957690))){
//			LOGGER.info(gene.getGene());
//		}

		return intervalTree;
	}

	public static Map<String, ChrPosTreeMap<LeadVariant>> loadLeadVariantsPerTrait(DownstreamerOptions options) throws IOException {

		if (leadVariantsPerTraitCache != null) {
			//this is save because options should not change during a run
			return leadVariantsPerTraitCache;
		} else {

			Map<String, File> alternativeLeadVariantFiles = options.getAlternativeTopHitFiles();

			Map<String, ChrPosTreeMap<LeadVariant>> leadVariantsPerTrait = new HashMap<>();

			for (Map.Entry<String, File> alternativeEntry : alternativeLeadVariantFiles.entrySet()) {

				String trait = alternativeEntry.getKey();
				File file = alternativeEntry.getValue();

				ChrPosTreeMap<LeadVariant> leadVariants = readLeadVariantFile(file);

				LOGGER.info("Loaded " + leadVariants.size() + " for " + trait + " lead variants from: " + file.getAbsolutePath());

				leadVariantsPerTrait.put(trait, leadVariants);

			}

			//Now read all files in lead variants folder.
			File leadVariantFolder = options.getLeadSnpsFolder();

			File[] leadVariantFiles = leadVariantFolder.listFiles(new LeadVariantFileNameFilter());

			if (leadVariantFiles == null) {
				LOGGER.info("Lead variant folder not found");
			} else {

				for (File leadVariantFile : leadVariantFiles) {

					String trait = leadVariantFile.getName().substring(0, leadVariantFile.getName().length() - 10);

					if (!leadVariantsPerTrait.containsKey(trait)) {
						ChrPosTreeMap<LeadVariant> leadVariants = readLeadVariantFile(leadVariantFile);
						leadVariantsPerTrait.put(trait, leadVariants);
						LOGGER.info("Loaded " + leadVariants.size() + " for " + trait + " lead variants");
					}

				}
			}
			leadVariantsPerTraitCache = leadVariantsPerTrait;

			return leadVariantsPerTrait;
		}
	}

	public static ChrPosTreeMap<LeadVariant> readLeadVariantFile(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).withSkipLines(1).build();

		ChrPosTreeMap<LeadVariant> leadVariants = new ChrPosTreeMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			String id = nextLine[0];
			String chr = nextLine[1];
			int pos = Integer.parseInt(nextLine[2]);
			double pvalue = Double.parseDouble(nextLine[3]);
			leadVariants.put(chr, pos, new LeadVariant(id, chr, pos, pvalue));
		}

		return leadVariants;

	}

	public static LinkedHashMap<String, Gene> readGenesMap(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

		final LinkedHashMap<String, Gene> genes = new LinkedHashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.put(nextLine[0], new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]));

		}

		return genes;

	}

	public static LinkedHashMap<String, List<Gene>> readGenesAsChrMap(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

		final LinkedHashMap<String, List<Gene>> genes = new LinkedHashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (!genes.keySet().contains(nextLine[1])) {
				genes.put(nextLine[1], new ArrayList<>());
			}
			genes.get(nextLine[1]).add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]));

		}

		return genes;

	}

	public static SampleIdIncludeFilter readSampleFile(File sampleFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withCSVParser(parser).withSkipLines(0).build();

		final HashSet<String> samples = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			samples.add(nextLine[0]);

		}

		return new SampleIdIncludeFilter(samples);

	}

	public static HashSet<String> readVariantFilterFile(File variantFilterFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(variantFilterFile))).withCSVParser(parser).withSkipLines(0).build();

		final HashSet<String> variants = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			variants.add(nextLine[0]);

		}

		return variants;

	}

	public static RandomAccessGenotypeData loadGenotypes(DownstreamerOptions options, Collection<String> variantsToInclude) throws IOException {
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

	public static RandomAccessGenotypeData loadGenotypes(DownstreamerOptions options) throws IOException {
		final RandomAccessGenotypeData referenceGenotypeData;

		referenceGenotypeData = options.getGenotypeType().createGenotypeData(options.getGenotypeBasePath());

		return referenceGenotypeData;

	}

	private static class LeadVariantFileNameFilter implements FilenameFilter {

		public LeadVariantFileNameFilter() {
		}

		@Override
		public boolean accept(File dir, String name) {
			return name.toLowerCase().endsWith("_leads.txt");
		}

	}

}
