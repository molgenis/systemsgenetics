package nl.systemsgenetics.downstreamer.io;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.gene.Gene;
import nl.systemsgenetics.downstreamer.summarystatistic.SummaryStatisticRecord;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterMaf;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;

import java.io.*;
import java.util.*;

public class IoUtils {

    private static final Logger LOGGER = Logger.getLogger(IoUtils.class);


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


    public static PerChrIntervalTree<Gene> readGenesAsIntervalTree(File geneFile) throws Exception {
        return readGenesAsIntervalTreeWindowed(geneFile, 0);
    }

    public static PerChrIntervalTree<Gene> readGenesAsIntervalTreeWindowed(File geneFile, int window) throws Exception {

        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

        PerChrIntervalTree<Gene> intervalTree = new PerChrIntervalTree<>(Gene.class);

        String[] nextLine;
        while ((nextLine = reader.readNext()) != null) {
            final ArrayList<Gene> genes = new ArrayList<>();
            genes.add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]) - window, Integer.parseInt(nextLine[3]) + window, nextLine[5]));

            intervalTree.addChrElements(nextLine[1], genes);
        }

        return intervalTree;
    }


    public static Set<String> readAlternativeIndependentVariants(String path) throws IOException {
        return readAlternativeIndependentVariants(new File(path));
    }

    public static Set<String> readAlternativeIndependentVariants(File path) throws IOException {
        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(path))).withCSVParser(parser).withSkipLines(0).build();

        HashSet<String> variantIds = new HashSet<>();

        String[] nextLine;
        while ((nextLine = reader.readNext()) != null) {
            //SummaryStatisticRecord curRec = new SummaryStatisticRecord(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Double.parseDouble(nextLine[3]));
            variantIds.add(nextLine[0]);
        }

        return variantIds;
    }

    public static Map<String, Set<String>> readAlternativeIndependentVariants(Map<String, File> files) throws IOException {
        Map<String, Set<String>> output = new HashMap<>();

        for (String file: files.keySet()) {
            output.put(file, readAlternativeIndependentVariants(files.get(file)));
        }

        return output;
    }



    public static List<SummaryStatisticRecord> readAlternativeIndependentVariantsAsRecords(String path) throws IOException {
        return readAlternativeIndependentVariantsAsRecords(new File(path));
    }

    public static List<SummaryStatisticRecord> readAlternativeIndependentVariantsAsRecords(File path) throws IOException {
        final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(path))).withCSVParser(parser).withSkipLines(0).build();

        List<SummaryStatisticRecord> variants = new ArrayList<>();

        String[] nextLine;
        while ((nextLine = reader.readNext()) != null) {
            SummaryStatisticRecord curRec = new SummaryStatisticRecord(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Double.parseDouble(nextLine[3]));
            variants.add(curRec);
        }

        return variants;
    }


    public static Map<String, Set<String>> readIndependentVariants(String path) throws IOException {
        return readIndependentVariants(new File(path));
    }

    public static Map<String, Set<String>> readIndependentVariants(File path) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(path));
        Map<String, Set<String>> output = new HashMap<>();

        String line;
        while ((line = reader.readLine()) != null) {
            List<String> data = Arrays.asList(line.split("\t"));
            Set<String> rsids = new HashSet<>();
            if(!data.isEmpty()){
                rsids.addAll(data.subList(1, data.size() - 1));
            }

            output.put(data.get(0), rsids);
        }

        reader.close();
        return output;
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

        VariantFilter variantFilter;
        if (variantsToInclude == null) {
            variantFilter = null;
        } else {
            variantFilter = new VariantIdIncludeFilter(new HashSet<>(variantsToInclude));
        }

        if (options.getMafFilter() != 0) {
            VariantFilter mafFilter = new VariantFilterMaf(options.getMafFilter());
            if (variantFilter == null) {
                variantFilter = mafFilter;
            } else {
                variantFilter = new VariantCombinedFilter(variantFilter, mafFilter);
            }
        }

        referenceGenotypeData = options.getGenotypeType().createFilteredGenotypeData(options.getGenotypeBasePath(), 10000, variantFilter, sampleFilter, null, 0.34f);

        return referenceGenotypeData;
    }
}
