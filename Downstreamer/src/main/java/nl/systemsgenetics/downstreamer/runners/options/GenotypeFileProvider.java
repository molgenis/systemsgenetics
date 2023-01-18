package nl.systemsgenetics.downstreamer.runners.options;

import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;

import java.io.File;

public interface GenotypeFileProvider {

    String[] getGenotypeBasePath();
    RandomAccessGenotypeDataReaderFormats getGenotypeType();
    File getGenotypeSamplesFile();
    double getMafFilter();

}
