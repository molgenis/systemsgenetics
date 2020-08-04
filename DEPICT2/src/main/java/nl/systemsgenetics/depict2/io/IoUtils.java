package nl.systemsgenetics.depict2.io;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import nl.systemsgenetics.depict2.Depict2Options;
import nl.systemsgenetics.depict2.gene.Gene;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterMaf;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;

import java.io.*;
import java.util.*;

public class IoUtils {

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

            genes.add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5]));

        }

        return genes;

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
            rsids.addAll(data.subList(1, data.size() - 1));

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

            genes.put(nextLine[0], new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5]));

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
            genes.get(nextLine[1]).add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3]), nextLine[5]));

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

    public static RandomAccessGenotypeData loadGenotypes(Depict2Options options, Collection<String> variantsToInclude) throws IOException {
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
