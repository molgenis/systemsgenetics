package nl.systemsgenetics.downstreamer.io;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import htsjdk.samtools.util.IntervalTreeMap;
import nl.systemsgenetics.downstreamer.DownstreamerDeprecated;
import nl.systemsgenetics.downstreamer.DownstreamerStep2Results;
import nl.systemsgenetics.downstreamer.pathway.PathwayDatabase;
import nl.systemsgenetics.downstreamer.pathway.PathwayEnrichments;
import nl.systemsgenetics.downstreamer.runners.DownstreamerUtilities;
import nl.systemsgenetics.downstreamer.runners.options.DownstreamerOptionsDeprecated;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.runners.options.GenotypeFileProvider;
import nl.systemsgenetics.downstreamer.summarystatistic.LdScore;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterMaf;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

import nl.systemsgenetics.downstreamer.containers.LeadVariant;
import org.molgenis.genotype.variantFilter.VariantFilterBiAllelic;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

public class IoUtils {

	private static final Logger LOGGER = LogManager.getLogger(IoUtils.class);
	private static Map<String, ChrPosTreeMap<LeadVariant>> leadVariantsPerTraitCache = null;
	private static class LeadVariantFileNameFilter implements FilenameFilter {

		public LeadVariantFileNameFilter() {
		}

		@Override
		public boolean accept(File dir, String name) {
			return (name.toLowerCase().endsWith("_leads.txt") || name.toLowerCase().endsWith("_leads.txt.gz"));
		}

	}

	public static List<String> readMatrixAnnotations(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (file.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))))).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).build();
		}

		ArrayList<String> identifiers = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			identifiers.add(nextLine[0]);
		}

		return identifiers;
	}

	public static List<Gene> readGenes(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(geneFile)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();
		}
		final ArrayList<Gene> genes = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]));

		}

		return genes;

	}

	public static IntervalTreeMap<Gene> readGenesAsIntervalTree(File geneFile) throws Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(geneFile)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();
		}

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

	public static ChrPosTreeMap<LeadVariant> readLeadVariantFile(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

		CSVReader reader = null;
		if (file.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).withSkipLines(1).build();
		}


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

		CSVReader reader = null;
		if (geneFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(geneFile)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();
		}

		final LinkedHashMap<String, Gene> genes = new LinkedHashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.put(nextLine[0], new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]));

		}

		return genes;

	}

	public static LinkedHashMap<String, List<Gene>> readGenesAsChrMap(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(geneFile)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();
		}

		final LinkedHashMap<String, List<Gene>> genes = new LinkedHashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (!genes.containsKey(nextLine[1])) {
				genes.put(nextLine[1], new ArrayList<>());
			}
			genes.get(nextLine[1]).add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]));

		}

		return genes;

	}

	public static LinkedHashMap<String, List<Gene>> readGenesAsChrArmMap(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(geneFile)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();
		}

		final LinkedHashMap<String, List<Gene>> genes = new LinkedHashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			Gene curGene = new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5], nextLine[6]);

			if (!genes.containsKey(curGene.getChrAndArm())) {
				genes.put(curGene.getChrAndArm(), new ArrayList<>());
			}
			genes.get(curGene.getChrAndArm()).add(curGene);

		}

		return genes;

	}

	public static SampleIdIncludeFilter readSampleFile(File sampleFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

		CSVReader reader = null;
		if (sampleFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(sampleFile)))))).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withCSVParser(parser).withSkipLines(0).build();
		}


		final HashSet<String> samples = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			samples.add(nextLine[0]);

		}

		return new SampleIdIncludeFilter(samples);

	}

	public static HashSet<String> readVariantFilterFile(File variantFilterFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (variantFilterFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(variantFilterFile)))))).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(variantFilterFile))).withCSVParser(parser).withSkipLines(0).build();
		}


		final HashSet<String> variants = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			variants.add(nextLine[0]);

		}

		return variants;

	}

	public static RandomAccessGenotypeData loadGenotypes(GenotypeFileProvider options, Collection<String> variantsToInclude) throws IOException {
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

	public static RandomAccessGenotypeData loadGenotypes(GenotypeFileProvider options) throws IOException {
		final RandomAccessGenotypeData referenceGenotypeData;

		referenceGenotypeData = options.getGenotypeType().createGenotypeData(options.getGenotypeBasePath());

		return referenceGenotypeData;

	}

	/**
	 * Read LD score files into an IntervalTreeMap in the format provided by
	 * https://github.com/bulik/ldsc
	 *
	 * @param options
	 * @return
	 * @throws IOException
	 */
	public static IntervalTreeMap<LdScore> readLdScores(String path) throws IOException {

		IntervalTreeMap<LdScore> output = new IntervalTreeMap<>();
		for (int i = 1; i < 23; i++) {
			BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path + "/" + i + ".l2.ldscore.gz"))));
			String[] line = reader.readLine().split("\t");

			if (line[0].equals("CHR")) {
				continue;
			}

			LdScore curLdscore = new LdScore(line[0],
					Integer.parseInt(line[2]),
					line[1],
					Double.parseDouble(line[5]));

			output.put(curLdscore, curLdscore);
		}

		return output;
	}

	// Specalized readers

	/**
	 * Load existing results from step 2 from storage
	 *
	 * @param options
	 * @return
	 * @throws Exception
	 */
	public static DownstreamerStep2Results loadExistingStep2Results(DownstreamerOptionsDeprecated options) throws Exception {
		return loadExistingStep2Results(options, false);
	}

	/**
	 * Load existing results from step 2 from storage. If
	 * matchToNormalizedPvalues = true, the genePvalues and
	 * normalizedGenePvalues are matched.
	 *
	 * @param options
	 * @return
	 * @throws Exception
	 */
	public static DownstreamerStep2Results loadExistingStep2Results(DownstreamerOptionsDeprecated options, boolean matchToNormalizedPvalues) throws Exception {

		DoubleMatrixDataset<String, String> genePvalues = DoubleMatrixDataset.loadDoubleBinaryData(options.getRun1BasePath() + "_genePvalues");
		DoubleMatrixDataset<String, String> normalizedGenePvalues;
		if (options.isForceNormalGenePvalues()) {
			normalizedGenePvalues = DownstreamerUtilities.getNormalizedGwasGenePvaluesReturn(options);
		} else {
			normalizedGenePvalues = null;
		}

		if (matchToNormalizedPvalues) {
			genePvalues = genePvalues.viewRowSelection(normalizedGenePvalues.getRowObjects());
		}

		final List<PathwayDatabase> pathwayDatabases = options.getPathwayDatabases();
		ArrayList<PathwayEnrichments> pathwayEnrichments = new ArrayList<>(pathwayDatabases.size());
		for (PathwayDatabase pathwayDatabase : pathwayDatabases) {
			pathwayEnrichments.add(new PathwayEnrichments(pathwayDatabase, options.getIntermediateFolder(), options.isExcludeHla()));
		}
		return new DownstreamerStep2Results(pathwayEnrichments, genePvalues, normalizedGenePvalues);

	}

	/**
	 * Load genotype data matching GWAS matrix and MAF filter.
	 *
	 * @param options Depict options object
	 * @return RandomAccesGenotypeData for all SNPs in GWAS matrix and MAF
	 * @throws IOException
	 */
	public static RandomAccessGenotypeData readReferenceGenotypeDataMatchingGwasSnps(DownstreamerOptionsDeprecated options) throws IOException {
		return readReferenceGenotypeDataMatchingGwasSnps(options, null);
	}

	/**
	 * Load genotype data matching GWAS matrix and MAF filter.
	 *
	 * @param options Depict options object
	 * @return RandomAccesGenotypeData for all SNPs in GWAS matrix and MAF
	 * @throws IOException
	 */
	public static RandomAccessGenotypeData readReferenceGenotypeDataMatchingGwasSnps(DownstreamerOptionsDeprecated options, Set<String> variantSubset) throws IOException {

		final List<String> variantsInZscoreMatrix;
		if (variantSubset == null) {
			String varSubSetFile = options.getGwasZscoreMatrixPath() + ".rows.txt";
			if (!new File(varSubSetFile).canRead()) {
				varSubSetFile += ".gz";
				if (!new File(varSubSetFile).canRead()) {
					throw new FileNotFoundException("Cannot find file: " + options.getGwasZscoreMatrixPath() + ".rows.txt or " + options.getGwasZscoreMatrixPath() + ".rows.txt.gz");
				}
			}
			variantsInZscoreMatrix = IoUtils.readMatrixAnnotations(new File(varSubSetFile));
		} else {
			variantsInZscoreMatrix = new ArrayList<>(variantSubset);
		}

		String gwasColFile = options.getGwasZscoreMatrixPath() + ".cols.txt";
		if (!new File(gwasColFile).canRead()) {
			gwasColFile += ".gz";
			if (!new File(gwasColFile).canRead()) {
				throw new FileNotFoundException("Cannot find file: " + options.getGwasZscoreMatrixPath() + ".cols.txt or " + options.getGwasZscoreMatrixPath() + ".cols.txt.gz");
			}
		}

		final List<String> phenotypesInZscoreMatrix = IoUtils.readMatrixAnnotations(new File(gwasColFile));

		LOGGER.info("Number of phenotypes in GWAS matrix: " + DownstreamerDeprecated.LARGE_INT_FORMAT.format(phenotypesInZscoreMatrix.size()));
		LOGGER.info("Number of variants in GWAS matrix: " + DownstreamerDeprecated.LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));

		if (options.getVariantFilterFile() != null) {
			HashSet<String> variantsToInclude = IoUtils.readVariantFilterFile(options.getVariantFilterFile());
			Iterator<String> variantsInZscoreMatrixIt = variantsInZscoreMatrix.iterator();
			while (variantsInZscoreMatrixIt.hasNext()) {
				String variant = variantsInZscoreMatrixIt.next();
				if (!variantsToInclude.contains(variant)) {
					variantsInZscoreMatrixIt.remove();
				}
			}
			LOGGER.info("Number of variants after filtering on selected variants: " + DownstreamerDeprecated.LARGE_INT_FORMAT.format(variantsInZscoreMatrix.size()));
		}

		return IoUtils.loadGenotypes(options, variantsInZscoreMatrix);
	}

	/**
	 * Load lead variants for all traits in the run.
	 * @param options
	 * @return
	 * @throws IOException
	 */
	public static Map<String, ChrPosTreeMap<LeadVariant>> loadLeadVariantsPerTrait(DownstreamerOptionsDeprecated options) throws IOException {

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


	public static void writeBlockDiagonalIndices(List<int[]> blockDiagonalIndices, List<String> colObjects, String basePath) throws IOException {
		File blockdiagonalFile =  new File(basePath);

		BufferedWriter writer = new BufferedWriter(new FileWriter(blockdiagonalFile));
		int curBlock = 0;
		writer.write("block\tindex\tgene");
		writer.newLine();
		for (int[] cur :blockDiagonalIndices) {

			for (int i : cur) {
				writer.write(curBlock + "\t" + i + "\t" + colObjects.get(i));
				writer.newLine();
			}
		}
		writer.flush();
		writer.close();

	}





}
