package eqtlmappingpipeline.ase;

import eqtlmappingpipeline.Main;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.lang.Thread.UncaughtExceptionHandler;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;
import umcg.genetica.io.gtf.GffElement;
import umcg.genetica.io.gtf.GtfReader;

/**
 *
 * @author Patrick Deelen
 */
public class Ase {

	private static final String VERSION = Main.VERSION;
	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |  Allele Specific Expression Mapper    |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  | Dasha Zhernakova, Marijke v/d Sijde,  |\n"
			+ "  |   Marc Jan Bonder, Harm-Jan Westra,   |\n"
			+ "  |      Lude Franke, Morris Swertz       |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Logger LOGGER = Logger.getLogger(Ase.class);
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();

	public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();

	public static void main(String[] args) {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://molgenis.org/systemsgenetics");
		System.out.println();

		System.out.println("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		System.out.println();

		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		if (args.length == 0) {
			AseConfiguration.printHelp();
			System.exit(1);
		}

		final AseConfiguration configuration;
		try {
			configuration = new AseConfiguration(args);
		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			AseConfiguration.printHelp();
			System.exit(1);
			return;
		}

		if (!configuration.getOutputFolder().isDirectory() && !configuration.getOutputFolder().mkdirs()) {
			System.err.println("Failed to create output folder: " + configuration.getOutputFolder().getAbsolutePath());
			System.exit(1);
		}

		startLogging(configuration.getLogFile(), configuration.isDebugMode());
		configuration.printOptions();

		final AseResults aseResults = new AseResults();
		final AtomicInteger sampleCounter = new AtomicInteger(0);
		final AtomicInteger fileCounter = new AtomicInteger(0);
		final RandomAccessGenotypeData referenceGenotypes;

		if (configuration.isRefSet()) {
			try {
				referenceGenotypes = configuration.getRefDataType().createGenotypeData(configuration.getRefBasePaths(), configuration.getRefDataCacheSize());
				System.out.println("Loading reference data complete");
				LOGGER.info("Loading reference data complete");
			} catch (IOException ex) {
				System.err.println("Unable to load reference genotypes file.");
				LOGGER.fatal("Unable to load reference genotypes file.", ex);
				System.exit(1);
				return;
			} catch (IncompatibleMultiPartGenotypeDataException ex) {
				System.err.println("Unable to load reference genotypes file.");
				LOGGER.fatal("Unable to load reference genotypes file.", ex);
				System.exit(1);
				return;
			} catch (GenotypeDataException ex) {
				System.err.println("Unable to load reference genotypes file.");
				LOGGER.fatal("Unable to load reference genotypes file.", ex);
				System.exit(1);
				return;
			}
		} else {
			referenceGenotypes = null;
		}

		if (configuration.isGtfSet()) {
			if (!configuration.getGtf().canRead()) {
				System.err.println("Cannot read GENCODE gft file.");
				LOGGER.fatal("Cannot read GENCODE gft file");
				System.exit(1);
				return;
			}
		}

		final Iterator<File> inputFileIterator = configuration.getInputFiles().iterator();
		
		int threadCount = configuration.getInputFiles().size() < configuration.getThreads() ? configuration.getInputFiles().size() : configuration.getThreads();
		List<Thread> threads = new ArrayList<Thread>(threadCount);
		final ThreadErrorHandler threadErrorHandler = new ThreadErrorHandler(threads);
		for (int i = 0; i < threadCount; ++i) {

			Thread worker = new Thread(new ReadCountsLoader(inputFileIterator, aseResults, sampleCounter, fileCounter, configuration, referenceGenotypes));
			worker.setUncaughtExceptionHandler(threadErrorHandler);
			worker.start();
			threads.add(worker);

		}

		int running;
		int nextReport = 100;
		do {
			running = 0;
			for (Thread thread : threads) {
				if (thread.isAlive()) {
					running++;
				}
			}
			try {
				Thread.sleep(500);
			} catch (InterruptedException ex) {
			}
			int currentCount = fileCounter.get();

			if (currentCount > nextReport) {
				//sometimes we skiped over report because of timing. This solved this
				System.out.println("Loaded " + DEFAULT_NUMBER_FORMATTER.format(nextReport) + " out of " + DEFAULT_NUMBER_FORMATTER.format(configuration.getInputFiles().size()) + " files");
				nextReport += 100;
			}

		} while (running > 0);


		LOGGER.info("Loading files complete. Detected " + DEFAULT_NUMBER_FORMATTER.format(sampleCounter) + " samples.");
		System.out.println("Loading files complete. Detected " + DEFAULT_NUMBER_FORMATTER.format(sampleCounter) + " samples.");

		Iterator<AseVariant> aseIterator = aseResults.iterator();
		while (aseIterator.hasNext()) {
			if (aseIterator.next().getSampleCount() < configuration.getMinSamples()) {
				aseIterator.remove();
			}
		}

		AseVariant[] aseVariants = new AseVariant[aseResults.getCount()];
		{
			int i = 0;
			for (AseVariant aseVariant : aseResults) {

				//This can be made multithreaded if needed
				aseVariant.calculateStatistics();

				aseVariants[i] = aseVariant;
				++i;

			}
		}

		int numberOfTests = aseVariants.length;
		double bonferroniCutoff = 0.05 / numberOfTests;
		System.out.println("Performed " + DEFAULT_NUMBER_FORMATTER.format(numberOfTests) + " tests. Bonferroni FWER 0.05 cut-off: " + bonferroniCutoff);
		LOGGER.info("Performed " + DEFAULT_NUMBER_FORMATTER.format(numberOfTests) + " tests. Bonferroni FWER 0.05 cut-off: " + bonferroniCutoff);


		Arrays.sort(aseVariants);

		final PerChrIntervalTree<GffElement> gtfAnnotations;
		if (configuration.isGtfSet()) {
			try {
				System.out.println("Started loading GTF file.");
				gtfAnnotations = new GtfReader(configuration.getGtf()).createIntervalTree();
				System.out.println("Loaded " + DEFAULT_NUMBER_FORMATTER.format(gtfAnnotations.size()) + " annotations from GTF file.");
				LOGGER.info("Loaded " + DEFAULT_NUMBER_FORMATTER.format(gtfAnnotations.size()) + " annotations from GTF file.");
			} catch (FileNotFoundException ex) {
				System.err.println("Cannot read GENCODE gft file.");
				LOGGER.fatal("Cannot read GENCODE gft file", ex);
				System.exit(1);
				return;
			} catch (Exception ex) {
				System.err.println("Cannot read GENCODE gft file. Error: " + ex.getMessage());
				LOGGER.fatal("Cannot read GENCODE gft file", ex);
				System.exit(1);
				return;
			}
		} else {
			gtfAnnotations = null;
		}


		File outputFileAll = new File(configuration.getOutputFolder(), "ase.txt");
		try {

			//print all restuls
			printAseResults(outputFileAll, aseVariants, gtfAnnotations);

		} catch (UnsupportedEncodingException ex) {
			throw new RuntimeException(ex);
		} catch (FileNotFoundException ex) {
			System.err.println("Unable to create output file at " + outputFileAll.getAbsolutePath());
			LOGGER.fatal("Unable to create output file at " + outputFileAll.getAbsolutePath(), ex);
			System.exit(1);
			return;
		} catch (IOException ex) {
			System.err.println("Unable to create output file at " + outputFileAll.getAbsolutePath());
			LOGGER.fatal("Unable to create output file at " + outputFileAll.getAbsolutePath(), ex);
			System.exit(1);
			return;
		}

		File outputFileBonferroni = new File(configuration.getOutputFolder(), "ase_bonferroni.txt");
		try {

			//print bonferroni significant results
			int bonferroniSignificantCount = printAseResults(outputFileBonferroni, aseVariants, gtfAnnotations, bonferroniCutoff, false);

			System.err.println("Completed writing " + DEFAULT_NUMBER_FORMATTER.format(bonferroniSignificantCount) + " bonferroni significant ASE variants");
			LOGGER.fatal("Completed writing " + DEFAULT_NUMBER_FORMATTER.format(bonferroniSignificantCount) + " bonferroni significant ASE variants");
			
		} catch (UnsupportedEncodingException ex) {
			throw new RuntimeException(ex);
		} catch (FileNotFoundException ex) {
			System.err.println("Unable to create output file at " + outputFileBonferroni.getAbsolutePath());
			LOGGER.fatal("Unable to create output file at " + outputFileBonferroni.getAbsolutePath(), ex);
			System.exit(1);
			return;
		} catch (IOException ex) {
			System.err.println("Unable to create output file at " + outputFileBonferroni.getAbsolutePath());
			LOGGER.fatal("Unable to create output file at " + outputFileBonferroni.getAbsolutePath(), ex);
			System.exit(1);
			return;
		}


		System.out.println("Program completed");
		LOGGER.info("Program completed");


	}

	private static void startLogging(File logFile, boolean debugMode) {

		try {
			FileAppender logAppender = new FileAppender(new SimpleLayout(), logFile.getCanonicalPath(), false);
			Logger.getRootLogger().removeAllAppenders();
			Logger.getRootLogger().addAppender(logAppender);
			if (debugMode) {
				LOGGER.setLevel(Level.DEBUG);
			} else {
				LOGGER.setLevel(Level.INFO);
			}
		} catch (IOException e) {
			System.err.println("Failed to create logger: " + e.getMessage());
			System.exit(1);
		}

		LOGGER.info(
				"\n" + HEADER);
		LOGGER.info("Version: " + VERSION);
		LOGGER.info("Current date and time: " + DATE_TIME_FORMAT.format(currentDataTime));
		LOGGER.info("Log level: " + LOGGER.getLevel());

		System.out.println("Started logging");
		System.out.println();
	}

	/**
	 *
	 * @param outputFile
	 * @param aseVariants
	 * @throws UnsupportedEncodingException
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private static int printAseResults(File outputFile, AseVariant[] aseVariants, final PerChrIntervalTree<GffElement> gtfAnnotations) throws UnsupportedEncodingException, FileNotFoundException, IOException {
		return printAseResults(outputFile, aseVariants, gtfAnnotations, 1);
	}

	/**
	 *
	 * @param outputFile
	 * @throws UnsupportedEncodingException
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private static int printAseResults(File outputFile, AseVariant[] aseVariants, final PerChrIntervalTree<GffElement> gtfAnnotations, double maxPvalue) throws UnsupportedEncodingException, FileNotFoundException, IOException {

		return printAseResults(outputFile, aseVariants, gtfAnnotations, maxPvalue, false);

	}

	private static int printAseResults(File outputFile, AseVariant[] aseVariants, final PerChrIntervalTree<GffElement> gtfAnnotations, double maxPvalue, boolean excludeNegativeCountR) throws UnsupportedEncodingException, FileNotFoundException, IOException {

		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), AseConfiguration.ENCODING));

		outputWriter.append("Meta_P\tMeta_Z\tChr\tPos\tSnpId\tSample_Count\tRef_Allele\tAlt_Allele\tCount_Pearson_R\tGenes\tRef_Counts\tAlt_Counts\n");
		
		int counter = 0;

		HashSet<String> genesPrinted = new HashSet<String>();
		for (AseVariant aseVariant : aseVariants) {

			if (aseVariant.getMetaPvalue() > maxPvalue) {
				continue;
			}

			if (excludeNegativeCountR && aseVariant.getCountPearsonR() < 0) {
				continue;
			}
			
			++counter;

			outputWriter.append(String.valueOf(aseVariant.getMetaPvalue()));
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getMetaZscore()));
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getChr());
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getPos()));
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getId().getPrimairyId() == null ? "." : aseVariant.getId().getPrimairyId());
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getSampleCount()));
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getA1().getAlleleAsString());
			outputWriter.append('\t');
			outputWriter.append(aseVariant.getA2().getAlleleAsString());
			outputWriter.append('\t');

			outputWriter.append(String.valueOf(aseVariant.getCountPearsonR()));
			outputWriter.append('\t');

			if (gtfAnnotations != null) {

				genesPrinted.clear();

				List<GffElement> elements = gtfAnnotations.searchPosition(aseVariant.getChr(), aseVariant.getPos());


				boolean first = true;
				for (GffElement element : elements) {

					String geneId = element.getAttributeValue("gene_id");

					if(genesPrinted.contains(geneId)){
						continue;
					}

					if (first) {
						first = false;
					} else {
						outputWriter.append(',');
					}
					outputWriter.append(geneId);
					genesPrinted.add(geneId);
				}

			}

			outputWriter.append('\t');

			for (int i = 0; i < aseVariant.getA1Counts().size(); ++i) {
				if (i > 0) {
					outputWriter.append(',');
				}
				outputWriter.append(String.valueOf(aseVariant.getA1Counts().getQuick(i)));
			}
			outputWriter.append('\t');
			for (int i = 0; i < aseVariant.getA2Counts().size(); ++i) {
				if (i > 0) {
					outputWriter.append(',');
				}
				outputWriter.append(String.valueOf(aseVariant.getA2Counts().getQuick(i)));
			}
			outputWriter.append('\n');

		}


		outputWriter.close();
		return counter;

	}
	
	private static class ThreadErrorHandler implements UncaughtExceptionHandler {

		private final List<Thread> threads;

		public ThreadErrorHandler(List<Thread> threads) {
			this.threads = threads;
		}
		
		@Override
		public void uncaughtException(Thread t, Throwable e) {
			
			System.err.println("Fatal error: " + e.getMessage());
			LOGGER.fatal("Fatal error: ", e);
			System.exit(1);
		}
		
	}
	
}