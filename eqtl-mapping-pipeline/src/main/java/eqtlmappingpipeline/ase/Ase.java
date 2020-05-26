package eqtlmappingpipeline.ase;

import eqtlmappingpipeline.Main;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.lang.Thread.UncaughtExceptionHandler;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;
import umcg.genetica.io.bedgraph.BedGraphEntry;
import umcg.genetica.io.bedgraph.BedGraphFile;
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
			+ "  |     Genomics Coordination Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Logger LOGGER = Logger.getLogger(Ase.class);
	private static final DateFormat DATE_TIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	private static final Date currentDataTime = new Date();
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();

	public static void main(String[] args) {

		System.out.println(HEADER);
		System.out.println();
		System.out.println("          --- Version: " + VERSION + " ---");
		System.out.println();
		System.out.println("More information: http://www.molgenis.org/systemsgenetics/Allele-Specific-Expression");
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
			System.err.println("Invalid command line arguments: " + ex.getMessage());
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

		LOGGER.debug("Java version: " + System.getProperty("java.version"));

		try {

			final AseResults aseResults = new AseResults();
			final Set<String> detectedSampleSet = Collections.synchronizedSet(new HashSet<String>());
			final RandomAccessGenotypeData referenceGenotypes;
			final Map<String, String> refToStudySampleId;
			final int minimumNumberSamples = configuration.getMinSamples();

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

				if (configuration.isSampleToRefSampleFileSet()) {
					try {
						refToStudySampleId = readSampleMapping(configuration.getSampleToRefSampleFile());
						System.out.println("Found " + refToStudySampleId.size() + " sample mappings");
						LOGGER.info("Found " + refToStudySampleId.size() + " sample mappings");
					} catch (FileNotFoundException ex) {
						System.err.println("Cannot find samples mapping file at: " + configuration.getSampleToRefSampleFile().getAbsolutePath());
						LOGGER.fatal("Cannot find samples mapping file at: " + configuration.getSampleToRefSampleFile().getAbsolutePath(), ex);
						System.exit(1);
						return;
					} catch (Exception ex) {
						System.err.println("Error reading sample mapping file: " + ex.getMessage());
						LOGGER.fatal("Error reading sample mapping file", ex);
						System.exit(1);
						return;
					}
				} else {
					refToStudySampleId = null;
				}

			} else {
				referenceGenotypes = null;
				refToStudySampleId = null;
			}

			if (configuration.isGtfSet()) {
				if (!configuration.getGtf().canRead()) {
					System.err.println("Cannot read GENCODE gft file.");
					LOGGER.fatal("Cannot read GENCODE gft file");
					System.exit(1);
				}
			}

			final List<File> inputFiles = configuration.getInputFiles();

			System.out.println("Loading sample allele counts");
			if (referenceGenotypes == null) {
				loadAseData(inputFiles, aseResults, detectedSampleSet, configuration, null, refToStudySampleId, configuration.getChrFilter(), true);

				Iterator<AseVariantAppendable> aseIterator = aseResults.iterator();
				while (aseIterator.hasNext()) {
					if (aseIterator.next().getSampleCount() < minimumNumberSamples) {
						aseIterator.remove();
					}
				}

			} else {

				//We have reference genotypes. This means a lot of random access to this single genotype data which slowed down the analyis a lot. 
				//Chucking the analysis will be used so that we can efficiently cache the reference genotypes.

				final int chunkSize = configuration.getChunkSize();

				for (String chr : referenceGenotypes.getSeqNames()) {

					if (configuration.getChrFilter() != null && !configuration.getChrFilter().equals(chr)) {
						continue;
					}

					for (int startChunk = 0; referenceGenotypes.getVariantsByRange(chr, startChunk, Integer.MAX_VALUE).iterator().hasNext(); startChunk += chunkSize) {
						if (chunkSize == Integer.MAX_VALUE) {
							System.out.println("Chr: " + chr);
						} else {
							System.out.println("Chr: " + chr + " chunk: " + DEFAULT_NUMBER_FORMATTER.format(startChunk) + "-" + DEFAULT_NUMBER_FORMATTER.format(startChunk + chunkSize));
						}

						loadAseData(inputFiles, aseResults, detectedSampleSet, configuration, referenceGenotypes, refToStudySampleId, chr, startChunk, (startChunk + chunkSize - 1), false);
						//System.out.println("Current number of ASE targets: " + aseResults.getCount());
					}
					//Clean up ASE that do not meet minimum number of samples	
					if (aseResults.chrIterator(chr) != null) {
						for (Iterator<AseVariantAppendable> aseChrIterator = aseResults.chrIterator(chr); aseChrIterator.hasNext();) {
							if (aseChrIterator.next().getSampleCount() < minimumNumberSamples) {
								aseChrIterator.remove();
							}
						}
					}
					//System.out.println("Current number of ASE targets after checking sample count: " + aseResults.getCount());

				}


			}

			final boolean encounteredBaseQuality = aseResults.isEncounteredBaseQuality();

			LOGGER.info("Loading files complete. Detected " + DEFAULT_NUMBER_FORMATTER.format(detectedSampleSet.size()) + " samples.");
			System.out.println("Loading files complete. Detected " + DEFAULT_NUMBER_FORMATTER.format(detectedSampleSet.size()) + " samples.");

			if (configuration.isMappabilityTrackSet()) {

				System.out.println("Start loading mappability track");
				
				int removeNoMappabilityInfo = 0;
				int removeLowMappability = 0;

				double minimumMappability = configuration.getMappabilityMinimum();

				PerChrIntervalTree<BedGraphEntry> mappabilities;
				try {
					BedGraphFile bedGraphFile = new BedGraphFile(configuration.getMappabilityTrackFile());
					mappabilities = bedGraphFile.createIntervalTree();
				} catch (FileNotFoundException ex) {
					System.err.println("Could not find mappability file at: " + configuration.getMappabilityTrackFile().getAbsolutePath());
					LOGGER.fatal("Could not find mappability file at: " + configuration.getMappabilityTrackFile().getAbsolutePath());
					System.exit(1);
					return;
				} catch (IOException ex) {
					System.err.println("Could not read mappability file at: " + configuration.getMappabilityTrackFile().getAbsolutePath());
					LOGGER.fatal("Could not read mappability file at: " + configuration.getMappabilityTrackFile().getAbsolutePath());
					System.exit(1);
					return;
				} catch (Exception ex) {
					System.err.println("Error reading mappability file " + ex.getMessage());
					LOGGER.fatal("Error reading mappability file " + ex.getMessage());
					System.exit(1);
					return;
				}


				for (Iterator<AseVariantAppendable> aseIterator = aseResults.iterator(); aseIterator.hasNext();) {
					AseVariantAppendable ase = aseIterator.next();
					List<BedGraphEntry> aseMappabilities = mappabilities.searchPosition(ase.getChr(), ase.getPos());
					if (aseMappabilities.isEmpty()) {
						aseIterator.remove();
						++removeNoMappabilityInfo;
						continue;
					} else if (aseMappabilities.size() > 1) {
						System.err.println("Error reading mappability file " + ase.getChr() + ":" + ase.getPos() + " contains multiple mappability scores");
						LOGGER.fatal("Error reading mappability file " + ase.getChr() + ":" + ase.getPos() + " contains multiple mappability scores");
						System.exit(1);
						return;
					} else if (aseMappabilities.get(0).getValue() < minimumMappability) {
						++removeLowMappability;
						aseIterator.remove();
					}

				}

				System.out.println("ASE removed due to low mappability: " + DEFAULT_NUMBER_FORMATTER.format(removeLowMappability));
				LOGGER.info("ASE removed due to low mappability: " + DEFAULT_NUMBER_FORMATTER.format(removeLowMappability));

				if (removeNoMappabilityInfo != 0) {
					System.out.println("ASE removed due to no mappability info present: " + DEFAULT_NUMBER_FORMATTER.format(removeNoMappabilityInfo));
					LOGGER.info("ASE removed due to no mappability info present: " + DEFAULT_NUMBER_FORMATTER.format(removeNoMappabilityInfo));
				}

			}


			AseVariantAppendable[] aseVariants = new AseVariantAppendable[aseResults.getCount()];
			{
				int i = 0;
				for (AseVariantAppendable aseVariant : aseResults) {
					aseVariants[i] = aseVariant;
					++i;
				}
			}
			
			AseCalculator.startAseCalculators(aseVariants, configuration.getThreads());

			System.out.println("Completed ASE calculations");
			
			//int numberOfTests = aseVariants.length;
			//double bonferroniCutoff = 0.05 / numberOfTests;
			//System.out.println("Performed " + DEFAULT_NUMBER_FORMATTER.format(numberOfTests) + " tests. Bonferroni FWER 0.05 cut-off: " + bonferroniCutoff);
			//LOGGER.info("Performed " + DEFAULT_NUMBER_FORMATTER.format(numberOfTests) + " tests. Bonferroni FWER 0.05 cut-off: " + bonferroniCutoff);


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


			for (MultipleTestingCorrectionMethod correctionMethod : EnumSet.of(MultipleTestingCorrectionMethod.NONE, MultipleTestingCorrectionMethod.BONFERRONI, MultipleTestingCorrectionMethod.HOLM, MultipleTestingCorrectionMethod.BH)) {

				File outputFileBonferroni = new File(configuration.getOutputFolder(), correctionMethod == MultipleTestingCorrectionMethod.NONE ? "ase.txt" : "ase_" + correctionMethod.toString().toLowerCase() + ".txt");
				try {

					int writtenResults = printAseResults(outputFileBonferroni, aseVariants, gtfAnnotations, correctionMethod, encounteredBaseQuality);

					if (correctionMethod == MultipleTestingCorrectionMethod.NONE) {
						System.out.println("Completed writing all " + DEFAULT_NUMBER_FORMATTER.format(writtenResults) + " ASE variants");
						LOGGER.info("Completed writing all " + DEFAULT_NUMBER_FORMATTER.format(writtenResults) + " ASE variants");
					} else {
						System.out.println("Completed writing " + DEFAULT_NUMBER_FORMATTER.format(writtenResults) + " " + correctionMethod.toString().toLowerCase() + " significant ASE variants");
						LOGGER.info("Completed writing " + DEFAULT_NUMBER_FORMATTER.format(writtenResults) + " " + correctionMethod.toString().toLowerCase() + " significant ASE variants");
					}


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
				} catch (AseException ex) {
					System.err.println("Error creating output file: " + ex.getMessage());
					LOGGER.fatal("Error creating output file.", ex);
					System.exit(1);
					return;
				}
			}

			System.out.println("Program completed");
			LOGGER.info("Program completed");
		} catch (OutOfMemoryError ex) {
			System.err.println("Error: out of memory, use -Xmx##g -Xms##g to reserve more memory, see manual for more information.");
			LOGGER.fatal("Out of memory, use -Xmx##g -Xms##g to reserve more memory, see manual for more information.", ex);
			System.exit(1);
		}

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
	 * @param aseVariants must be sorted on significance
	 * @throws UnsupportedEncodingException
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private static int printAseResults(File outputFile, AseVariantAppendable[] aseVariants, final PerChrIntervalTree<GffElement> gtfAnnotations, final boolean encounteredBaseQuality) throws UnsupportedEncodingException, FileNotFoundException, IOException, AseException {
		return printAseResults(outputFile, aseVariants, gtfAnnotations, MultipleTestingCorrectionMethod.NONE, encounteredBaseQuality);
	}

	/**
	 *
	 * @param outputFile
	 * @param aseVariants must be sorted on significance
	 * @param gtfAnnotations
	 * @param onlyHolmSignificant
	 * @return
	 * @throws UnsupportedEncodingException
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private static int printAseResults(final File outputFile, final AseVariantAppendable[] aseVariants, final PerChrIntervalTree<GffElement> gtfAnnotations, final MultipleTestingCorrectionMethod multipleTestingCorrectionMethod, final boolean encounteredBaseQuality) throws UnsupportedEncodingException, FileNotFoundException, IOException, AseException {

		final BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), AseConfiguration.ENCODING));

		outputWriter.append("LikelihoodRatioP\tLikelihoodRatioD\tLikelihoodRatioP\tLikelihoodRatioTheta\tMeta_P\tMeta_Z\tChr\tPos\tSnpId\tSample_Count\tRef_Allele\tAlt_Allele\tCount_Pearson_R\tGenes\tRef_Counts\tAlt_Counts\tBinom_P\tSampleIds");

//		if (encounteredBaseQuality) {
//			outputWriter.append("\tRef_MeanBaseQuality\tAlt_MeanBaseQuality\tRef_MeanBaseQualities\tAlt_MeanBaseQualities");
//		}
		outputWriter.append('\n');

		final double significance = 0.05;

		int counter = 0;
		double lastRatioD = Double.POSITIVE_INFINITY;

		final int totalNumberOfTests = aseVariants.length;
		final double bonferroniCutoff = significance / totalNumberOfTests;

		HashSet<String> genesPrinted = new HashSet<String>();
		aseVariants:
		for (AseVariantAppendable aseVariant : aseVariants) {

			double ratioD = aseVariant.getMle().getRatioD();
			if (ratioD > lastRatioD) {
				throw new AseException("ASE results not sorted");
			}
			lastRatioD = ratioD;

			switch (multipleTestingCorrectionMethod) {
				case NONE:
					break;
				case NOMINAL:
					if (aseVariant.getMle().getRatioP() > significance) {
						break aseVariants;
					}
					break;
				case BONFERRONI:
					if (aseVariant.getMle().getRatioP() > bonferroniCutoff) {
						break aseVariants;
					}
					break;
				case HOLM:
					final double holmCutoff = significance / (totalNumberOfTests - counter);
					if (aseVariant.getMle().getRatioP() > holmCutoff) {
						break aseVariants;
					}
					break;
				case BH:
					final double qvalue = ((counter + 1d) / totalNumberOfTests) * significance;
					if (aseVariant.getMle().getRatioP() > qvalue) {
						break aseVariants;
					}
					break;
				default:
					throw new AseException("Multiple testing method: " + multipleTestingCorrectionMethod + " is not supported");

			}

			++counter;


			outputWriter.append(String.valueOf(aseVariant.getMle().getRatioP()));
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getMle().getRatioD()));
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getMle().getMaxLikelihoodP()));
			outputWriter.append('\t');
			outputWriter.append(String.valueOf(aseVariant.getMle().getMaxLogLikelihoodTheta()));
			outputWriter.append('\t');
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

					if (genesPrinted.contains(geneId)) {
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

			outputWriter.append('\t');
			for (int i = 0; i < aseVariant.getPValues().size(); ++i) {
				if (i > 0) {
					outputWriter.append(',');
				}
				outputWriter.append(String.valueOf(aseVariant.getPValues().getQuick(i)));
			}

			outputWriter.append('\t');
			for (int i = 0; i < aseVariant.getSampleIds().size(); ++i) {
				if (i > 0) {
					outputWriter.append(',');
				}
				outputWriter.append(aseVariant.getSampleIds().get(i));
			}

//			if (encounteredBaseQuality) {
//
//				StringBuilder refMeanBaseQualities = new StringBuilder();
//				double sumRefMeanBaseQualities = 0;
//				for (int i = 0; i < aseVariant.getA1MeanBaseQualities().size(); ++i) {
//					sumRefMeanBaseQualities += aseVariant.getA1MeanBaseQualities().getQuick(i);
//					if (i > 0) {
//						refMeanBaseQualities.append(',');
//					}
//					refMeanBaseQualities.append(String.valueOf(aseVariant.getA1MeanBaseQualities().getQuick(i)));
//				}
//				outputWriter.append('\t');
//				outputWriter.append(String.valueOf( sumRefMeanBaseQualities / aseVariant.getA1MeanBaseQualities().size() ));
//				
//				StringBuilder altMeanBaseQualities = new StringBuilder();
//				double sumAtMeanBaseQualities = 0;
//				for (int i = 0; i < aseVariant.getA2MeanBaseQualities().size(); ++i) {
//					sumAtMeanBaseQualities += aseVariant.getA2MeanBaseQualities().getQuick(i);
//					if (i > 0) {
//						altMeanBaseQualities.append(',');
//					}
//					altMeanBaseQualities.append(String.valueOf(aseVariant.getA2MeanBaseQualities().getQuick(i)));
//				}
//				outputWriter.append('\t');
//				outputWriter.append(String.valueOf( sumAtMeanBaseQualities / aseVariant.getA2MeanBaseQualities().size() ));
//				outputWriter.append('\t');
//				outputWriter.append(refMeanBaseQualities);
//				outputWriter.append('\t');
//				outputWriter.append(altMeanBaseQualities);
//
//			}



			outputWriter.append('\n');

		}


		outputWriter.close();
		return counter;

	}

	/**
	 *
	 *
	 * @param sampleToRefSampleFile reference and value sample ID of study
	 * @return unmodifiable map with key sample ID in
	 */
	private static Map<String, String> readSampleMapping(File sampleToRefSampleFile) throws FileNotFoundException, UnsupportedEncodingException, IOException, Exception {

		HashMap<String, String> sampleMap = new HashMap<String, String>();

		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(sampleToRefSampleFile), "UTF-8"));

		String line;
		String[] elements;

		while ((line = reader.readLine()) != null) {
			elements = TAB_PATTERN.split(line);
			if (elements.length != 2) {
				throw new Exception("Detected " + elements.length + " columns instead of 2 for this line: " + line);
			}
			sampleMap.put(elements[1], elements[0]);
		}

		return Collections.unmodifiableMap(sampleMap);
	}

	protected static class ThreadErrorHandler implements UncaughtExceptionHandler {

		@Override
		public void uncaughtException(Thread t, Throwable e) {

			if(e instanceof OutOfMemoryError){
				System.err.println("Error: out of memory, use -Xmx##g -Xms##g to reserve more memory, see manual for more information.");
				LOGGER.fatal("Out of memory, use -Xmx##g -Xms##g to reserve more memory, see manual for more information.", e);
			} else {
				System.err.println("Fatal error: " + e.getMessage());
				LOGGER.fatal("Fatal error: ", e);
			}
			
			
			System.exit(1);
		}
	}

	private static void loadAseData(List<File> inputFiles, AseResults aseResults, Set<String> detectedSampleSet, AseConfiguration configuration, RandomAccessGenotypeData referenceGenotypes, Map<String, String> refToStudySampleId, String chr, boolean showFileProgress) {
		loadAseData(inputFiles, aseResults, detectedSampleSet, configuration, referenceGenotypes, refToStudySampleId, chr, 0, Integer.MAX_VALUE, showFileProgress);
	}

	private static void loadAseData(List<File> inputFiles, AseResults aseResults, Set<String> detectedSampleSet, AseConfiguration configuration, RandomAccessGenotypeData referenceGenotypes, Map<String, String> refToStudySampleId, String chr, int start, int stop, boolean showFileProgress) {

		final AtomicInteger fileCounter = new AtomicInteger(0);
		int threadCount = configuration.getInputFiles().size() < configuration.getThreads() ? configuration.getInputFiles().size() : configuration.getThreads();
		List<Thread> threads = new ArrayList<Thread>(threadCount);
		final ThreadErrorHandler threadErrorHandler = new ThreadErrorHandler();

		Iterator<File> inputFileIterator = inputFiles.iterator();

		for (int i = 0; i < threadCount; ++i) {

			Thread worker = new Thread(new ReadCountsLoader(inputFileIterator, aseResults, detectedSampleSet, fileCounter, configuration, referenceGenotypes, refToStudySampleId, chr, start, stop));
			worker.setUncaughtExceptionHandler(threadErrorHandler);
			worker.start();
			threads.add(worker);

		}

		boolean running;
		int nextReport = 100;
		do {
			running = false;
			for (Thread thread : threads) {
				if (thread.isAlive()) {
					running = true;
				}
			}
			try {
				Thread.sleep(500);
			} catch (InterruptedException ex) {
			}

			if (showFileProgress) {
				
				int currentCount = fileCounter.get();

				if (currentCount > nextReport) {
					//sometimes we skiped over report because of timing. This solved this
					System.out.println("Loaded " + DEFAULT_NUMBER_FORMATTER.format(nextReport) + " out of " + DEFAULT_NUMBER_FORMATTER.format(configuration.getInputFiles().size()) + " files");
					nextReport += 100;
				}
			}

		} while (running);

	}
	
	
	
}