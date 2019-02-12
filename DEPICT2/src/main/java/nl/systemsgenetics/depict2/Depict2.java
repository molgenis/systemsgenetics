package nl.systemsgenetics.depict2;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ResourceBundle;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author Patrick Deelen
 */
public class Depict2 {

	private static final String VERSION = ResourceBundle.getBundle("verion").getString("application.version");
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Logger LOGGER = Logger.getLogger(Depict2.class);
	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |                DEPICT2                |\n"
			+ "  |                                       |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	/**
	 * @param args the command line arguments
	 * @throws java.lang.InterruptedException
	 * @throws java.io.IOException
	 */
	public static void main(String[] args) throws InterruptedException, IOException {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		Date currentDataTime = new Date();
		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		Thread.sleep(25); //Allows flush to complete

		Depict2Options options;

		if (args.length == 0) {
			Depict2Options.printHelp();
			return;
		}

		try {
			options = new Depict2Options(args);
		} catch (ParseException ex) {
			System.err.println("Error parsing commandline: " + ex.getMessage());
			Depict2Options.printHelp();
			return;
		}


		if (options.getLogFile().getParentFile() != null && !options.getLogFile().getParentFile().isDirectory()) {
            if (!options.getLogFile().getParentFile().mkdirs()) {
                System.err.println("Failed to create output folder: " + options.getLogFile().getParent());
                System.exit(1);
            }
        }

        try {
            FileAppender logAppender = new FileAppender(new SimpleLayout(), options.getLogFile().getCanonicalPath(), false);
            LOGGER.getRootLogger().removeAllAppenders();
            LOGGER.getRootLogger().addAppender(logAppender);
            if (options.isDebugMode()) {
                LOGGER.setLevel(Level.DEBUG);
            } else {
                LOGGER.setLevel(Level.INFO);
            }
        } catch (IOException e) {
            System.err.println("Failed to create logger: " + e.getMessage());
            System.exit(1);
        }
		
		options.printOptions();

		switch (options.getMode()) {
			case CONVERT_EQTL:
				convertEqtlToBin(options);
				break;
			case CONVERT_TXT:
				convertTxtToBin(options);
				break;
			case RUN:
				run(options);
				break;
		}

		System.out.println("Completed mode: " + options.getMode());

		currentDataTime = new Date();
		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));

	}

	public static void run(Depict2Options options) throws IOException {
		
		final List<String> variantsInZscoreMatrix = readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".rows.txt"));
		final List<String> phenotypesInZscoreMatrix = readMatrixAnnotations(new File(options.getGwasZscoreMatrixPath() + ".cols.txt"));

		System.out.println("Number of phenotypes in GWAS matrix: " + phenotypesInZscoreMatrix.size());
		System.out.println("Number of variants in GWAS matrix: " + variantsInZscoreMatrix.size());
		
		RandomAccessGenotypeData referenceGenotypeData = loadGenotypes(options, variantsInZscoreMatrix);
		
		
//		Iterator<GeneticVariant> it = referenceGenotypeData.iterator();
//		for(int i = 0 ; i < 5 ; i++){
//			if(!it.hasNext()){
//				System.out.println("No more variants");
//			}
//			GeneticVariant v = it.next();
//			System.out.println(v.getSequenceName() + ":" + v.getStartPos() + " " + v.getVariantId().getPrimairyId());
//		}
//		System.out.println("---");
//		GeneticVariant v = referenceGenotypeData.getVariantsByPos("1", 10177).iterator().next();
//		System.out.println(v.getSequenceName() + ":" + v.getStartPos() + " " + v.getVariantId().getPrimairyId());
//		System.out.println("---");
//		for (GeneticVariant v2 : referenceGenotypeData.getVariantsByRange("1", 10100, 10200)) {
//			System.out.println(v2.getSequenceName() + ":" + v2.getStartPos() + " " + v2.getVariantId().getPrimairyId());
//		}
//		System.out.println("---");
//		for (GeneticVariant v2 : referenceGenotypeData.getVariantsByRange("1", 11919466, 12132045)) {
//			System.out.println(v2.getSequenceName() + ":" + v2.getStartPos() + " " + v2.getVariantId().getPrimairyId());
//		}
//		

		System.out.println("Done loading genotype data");

		List<Gene> genes = readGenes(options.getGeneInfoFile());

		System.out.println("Loaded " + genes.size() + " genes");

		DoubleMatrixDataset<String, String> genePvalues = CalculateGenePvalues.calculatorGenePvalues(options.getGwasZscoreMatrixPath(), new GenotypeCorrelationGenotypes(referenceGenotypeData), genes, options.getWindowExtend(), options.getMaxRBetweenVariants(), options.getNumberOfPermutations());

		System.out.println("Finished calculating gene p-values");

		genePvalues.save(options.getOutputBasePath() + "-genePvalues.txt");
	}

	private static RandomAccessGenotypeData loadGenotypes(Depict2Options options, List<String> variantsToInclude) throws IOException {
		final RandomAccessGenotypeData referenceGenotypeData;

		final SampleFilter sampleFilter;
		if (options.getGenotypeSamplesFile() != null) {
			sampleFilter = readSampleFile(options.getGenotypeSamplesFile());
		} else {
			sampleFilter = null;
		}
		
		try {
			referenceGenotypeData = options.getGenotypeType().createFilteredGenotypeData(options.getGenotypeBasePath(), 10000, new VariantIdIncludeFilter(new HashSet<>(variantsToInclude)), sampleFilter, null, 0.34f);
		} catch (TabixFileNotFoundException e) {
			System.err.println("Tabix file not found for input data at: " + e.getPath() + "\n"
					+ "Please see README on how to create a tabix file");
			System.exit(1);
			return null;
		} catch (IOException e) {
			System.err.println("Error accessing input data: " + e.getMessage());
			System.exit(1);
			return null;
		} catch (IncompatibleMultiPartGenotypeDataException e) {
			System.err.println("Error combining the impute genotype data files: " + e.getMessage());
			System.exit(1);
			return null;
		} catch (GenotypeDataException e) {
			System.err.println("Error reading input data: " + e.getMessage());
			System.exit(1);
			return null;
		}
		return referenceGenotypeData;
	}

	protected static final List<String> readMatrixAnnotations(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).build();

		ArrayList<String> identifiers = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			identifiers.add(nextLine[0]);

		}

		return identifiers;

	}

	private static List<Gene> readGenes(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

		final ArrayList<Gene> genes = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.add(new Gene(nextLine[0], nextLine[1], Integer.parseInt(nextLine[2]), Integer.parseInt(nextLine[3])));

		}

		return genes;

	}

	private static SampleIdIncludeFilter readSampleFile(File sampleFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withCSVParser(parser).withSkipLines(0).build();

		final HashSet<String> samples = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			samples.add(nextLine[0]);

		}

		return new SampleIdIncludeFilter(samples);

	}

	private static void convertTxtToBin(Depict2Options options) throws IOException {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadDoubleTextData(options.getGwasZscoreMatrixPath(), '\t');
		matrix.saveBinary(options.getOutputBasePath());

	}

	private static void convertEqtlToBin(Depict2Options options) throws IOException {

		DoubleMatrixDataset<String, String> matrix = DoubleMatrixDataset.loadTransEqtlExpressionMatrix(options.getGwasZscoreMatrixPath());
		matrix.saveBinary(options.getOutputBasePath());

	}

}
