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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

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

    public static RandomAccessGenotypeData loadGenotypes(Depict2Options options, List<String> variantsToInclude) throws IOException {
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
