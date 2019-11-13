package org.molgenis.genotype;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.util.GenotypeCountCalculator;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilterSeqPos;

/**
 *
 * @author Patrick Deelen
 */
public class GenotypeInfo {

	private static final Logger LOGGER;
	private static final Options OPTIONS;
	private static final Pattern CHR_POS_SPLITTER = Pattern.compile("\\s+|:");
	public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();

	static {

		LOGGER = Logger.getLogger(GenotypeInfo.class);

		OPTIONS = new Options();

		Option option;

		option = OptionBuilder.withArgName("basePath")
				.hasArgs()
				.withDescription("The base path of the data to align. The extensions are determined based on the input data type.")
				.withLongOpt("input")
				.isRequired()
				.create("i");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("type")
				.hasArg()
				.withDescription("The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* BGEN - Oxford .bgen & .sampleg\n"
				+ "* TRITYPER - TriTyper format folder")
				.withLongOpt("inputType")
				.create("I");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("basePath")
				.hasArg()
				.withDescription("The output bash path")
				.withLongOpt("output")
				.isRequired()
				.create("o");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("string")
				.hasArg()
				.withDescription("Path to file with samples IDs to include from input data. For plink data only the sample id (column 2) is used")
				.withLongOpt("sampleFilterList")
				.create("sf");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("string")
				.hasArg()
				.withDescription("Path to file with variant CHR\tPOS or CHR:POS to include from input data.")
				.withLongOpt("variantPosFilterList")
				.create("pf");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("string")
				.hasArg()
				.withDescription("Shapeit2 does not output the sequence name in the first column of the haplotype file. Use this option to force the chromosome for all variants. This option is only valid in combination with --inputType SHAPEIT2")
				.withLongOpt("forceChr")
				.create("f");
		OPTIONS.addOption(option);

		option = OptionBuilder.withArgName("double")
				.hasArg()
				.withDescription("The minimum posterior probability to call genotypes in the input data " + 0.8)
				.withLongOpt("inputProb")
				.create("ip");
		OPTIONS.addOption(option);

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		try {

			CommandLineParser parser = new PosixParser();
			final CommandLine commandLine = parser.parse(OPTIONS, args, true);

			final String[] inputBasePaths = commandLine.getOptionValues('i');
			final RandomAccessGenotypeDataReaderFormats inputType;

			try {
				if (commandLine.hasOption('I')) {
					inputType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue('I').toUpperCase());
				} else {
					if (inputBasePaths[0].endsWith(".vcf")) {
						throw new ParseException("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					}
					try {
						inputType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(inputBasePaths[0]);
					} catch (GenotypeDataException e) {
						throw new ParseException("Unable to determine input type based on specified path. Please specify --inputType");
					}
				}
			} catch (IllegalArgumentException e) {
				throw new ParseException("Error parsing --inputType \"" + commandLine.getOptionValue('I') + "\" is not a valid input data format");
			}

			final String outputBasePath = commandLine.getOptionValue('o');
			final File sampleFilterListFile = commandLine.hasOption("sf") ? new File(commandLine.getOptionValue("sf")) : null;
			final File variantPosFilterListFile = commandLine.hasOption("pf") ? new File(commandLine.getOptionValue("pf")) : null;

			final double minimumPosteriorProbability;
			try {
				minimumPosteriorProbability = commandLine.hasOption("ip") ? Double.parseDouble(commandLine.getOptionValue("ip")) : 0.8;
			} catch (NumberFormatException e) {
				throw new ParseException("Error parsing --inputProb \"" + commandLine.getOptionValue("ip") + "\" is not an double");
			}

			StringBuilder inputPaths = new StringBuilder();
			for (String path : inputBasePaths) {
				inputPaths.append(path);
				inputPaths.append(' ');
			}

			String forceSeqName = commandLine.hasOption('f') ? commandLine.getOptionValue('f') : null;

			LOGGER.info("Input base path: " + inputPaths);
			LOGGER.info("Input data type: " + inputType.getName());
			LOGGER.info("Output base path: " + outputBasePath);
			LOGGER.info("Force input sequence name: " + (forceSeqName == null ? "not forcing" : "forcing to: " + forceSeqName));
			if (sampleFilterListFile != null) {
				LOGGER.info("Filter input data to samples present in: " + sampleFilterListFile);
			}
			if (variantPosFilterListFile != null) {
				LOGGER.info("Filter input data to variants present in: " + variantPosFilterListFile);
			}
			LOGGER.info("Minimum posterior probability for input data: " + minimumPosteriorProbability);

			final VariantFilterSeqPos varFilter;
			if (variantPosFilterListFile != null) {
				varFilter = new VariantFilterSeqPos();
				int includeCountVar = 0;
				try {
					BufferedReader variantChrPosFilterReader = new BufferedReader(new FileReader(variantPosFilterListFile));
					String line;
					while ((line = variantChrPosFilterReader.readLine()) != null) {

						String[] elements = CHR_POS_SPLITTER.split(line);

						if (elements.length != 2) {
							LOGGER.error("Error parsing chr pos for line: " + line + " skipping line");
							continue;
						}
						++includeCountVar;
						varFilter.addSeqPos(elements[0], Integer.parseInt(elements[1]));
					}

					LOGGER.info("Included " + DEFAULT_NUMBER_FORMATTER.format(includeCountVar) + " variants from file with chr and pos");

				} catch (FileNotFoundException ex) {
					LOGGER.fatal("Unable to find file with variants to filter on at: " + variantPosFilterListFile.getAbsolutePath());
					System.exit(1);
				} catch (IOException e) {
					LOGGER.fatal("Error reading file with variants to filter on at: " + variantPosFilterListFile, e);
					System.exit(1);
				}
			} else {
				varFilter = null;
			}

			final SampleIdIncludeFilter sampleFilter;

			if (sampleFilterListFile != null) {
				HashSet<String> samples = new HashSet<String>();
				try {
					BufferedReader sampleIdFilterReader = new BufferedReader(new FileReader(sampleFilterListFile));
					String line;
					while ((line = sampleIdFilterReader.readLine()) != null) {
						samples.add(line);
					}
				} catch (FileNotFoundException ex) {
					LOGGER.fatal("Unable to find file with samples to filter on at: " + sampleFilterListFile.getAbsolutePath());
					System.exit(1);
				} catch (IOException e) {
					LOGGER.fatal("Error reading file with samples to filter on at: " + sampleFilterListFile.getAbsolutePath(), e);
					System.exit(1);
				}
				sampleFilter = new SampleIdIncludeFilter(samples);
			} else {
				sampleFilter = null;
			}

			final RandomAccessGenotypeData inputData;

			try {
				inputData = inputType.createFilteredGenotypeData(inputBasePaths, 0, varFilter, sampleFilter, forceSeqName, minimumPosteriorProbability);
			} catch (TabixFileNotFoundException e) {
				LOGGER.fatal("Tabix file not found for input data at: " + e.getPath() + "\n"
						+ "Please see README on how to create a tabix file");
				System.exit(1);
				return;
			} catch (IOException e) {
				LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
				System.exit(1);
				return;
			} catch (IncompatibleMultiPartGenotypeDataException e) {
				LOGGER.fatal("Error combining the impute genotype data files: " + e.getMessage(), e);
				System.exit(1);
				return;
			} catch (GenotypeDataException e) {
				LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
				System.exit(1);
				return;
			}

			LOGGER.info("Data loaded");

			final int[] sampleHetCounts = new int[inputData.getSamples().size()];
			final int[] sampleCallCounts = new int[inputData.getSamples().size()];

			double sumVariantCallRate = 0;

			File variantInfoFile = new File(outputBasePath + ".vars");
			if (variantInfoFile.getParentFile() != null) {
				variantInfoFile.getParentFile().mkdirs();
			}

			BufferedWriter variantInfoWriter = new BufferedWriter(new FileWriter(variantInfoFile));

			variantInfoWriter.append("ID\tCHR\tPOS\tAlleles\tMA\tMAF\tCALL\tHWE\tMACH_R2\tGenotype_Counts\n");

			int variantCount = 0;
			for (GeneticVariant variant : inputData) {
				variantInfoWriter.append(variant.getPrimaryVariantId());
				variantInfoWriter.append('\t');
				variantInfoWriter.append(variant.getSequenceName());
				variantInfoWriter.append('\t');
				variantInfoWriter.append(String.valueOf(variant.getStartPos()));
				variantInfoWriter.append('\t');
				boolean notFirst = false;
				for (Allele a : variant.getVariantAlleles()) {
					if (notFirst) {
						variantInfoWriter.append('/');
					}
					variantInfoWriter.append(a.getAlleleAsString());
					notFirst = true;
				}
				variantInfoWriter.append('\t');
				variantInfoWriter.append(String.valueOf(variant.getMinorAllele()));
				variantInfoWriter.append('\t');
				variantInfoWriter.append(String.valueOf(variant.getMinorAlleleFrequency()));
				variantInfoWriter.append('\t');
				double callrate = variant.getCallRate();
				sumVariantCallRate += callrate;
				variantInfoWriter.append(String.valueOf(callrate));
				variantInfoWriter.append('\t');
				variantInfoWriter.append(String.valueOf(variant.getHwePvalue()));
				variantInfoWriter.append('\t');
				variantInfoWriter.append(String.valueOf(variant.getMachR2()));
				variantInfoWriter.append('\t');

				ArrayList<GenotypeCountCalculator.GenotypeCount> genotypeCounts = GenotypeCountCalculator.countGenotypes(variant);
				for (GenotypeCountCalculator.GenotypeCount genotypeCount : genotypeCounts) {
					variantInfoWriter.append(genotypeCount.getGenotype().toString());
					variantInfoWriter.append(": ");
					variantInfoWriter.append(String.valueOf(genotypeCount.getCount()));
					variantInfoWriter.append(", ");
				}

				variantInfoWriter.append('\n');

				List<Alleles> sampleAlleles = variant.getSampleVariants();
				float[][] probs = variant.getSampleGenotypeProbilities();
				
				sampleAllelesLoop:
				for (int j = 0 ; j < sampleAlleles.size() ; ++j) {
					
					Alleles sampleAllele = sampleAlleles.get(j);
					
					if (sampleAllele.getAlleleCount() == 0 || sampleAllele.contains(Allele.ZERO)) {
//						System.out.println(inputData.getSampleNames()[j]);
//						System.out.println(variant.getPrimaryVariantId());
//						System.out.println(sampleAllele);
//						System.out.println(probs[j][0] + " " + probs[j][1] + " " + probs[j][2]);
//						System.out.println("------");
						continue sampleAllelesLoop;
					}
					sampleCallCounts[j]++;

					if (sampleAllele.getAlleleCount() > 1) {
						Allele a1 = sampleAllele.getAlleles().get(0);
						for (Allele a2 : sampleAllele.getAlleles()) {
							if (a1 != a2) {
								sampleHetCounts[j]++;
							}
						}
					}

				}
				++variantCount;
				if ((variantCount % 100000) == 0) {
					System.out.println("Processed " + DEFAULT_NUMBER_FORMATTER.format(variantCount) + " variants");
				}

			}

			variantInfoWriter.close();


			File sampleInfoFile = new File(outputBasePath + ".samples");
			if (sampleInfoFile.getParentFile() != null) {
				sampleInfoFile.getParentFile().mkdirs();
			}

			double sumSampleCallRate = 0;
			
			BufferedWriter sampleInfoWriter = new BufferedWriter(new FileWriter(sampleInfoFile));

			sampleInfoWriter.append("ID\tCallRate\tHetRate\n");
			int j = 0;
			sampleAlleles:
			for (String sample : inputData.getSampleNames()) {
				sampleInfoWriter.append(sample);
				sampleInfoWriter.append('\t');
				
				double sampleCallRate = sampleCallCounts[j] / (double) variantCount;
				sumSampleCallRate += sampleCallRate;
				sampleInfoWriter.append(String.valueOf(sampleCallRate));
				sampleInfoWriter.append('\t');
				
				double sampleHetRate = sampleHetCounts[j] / (double) variantCount;
				sampleInfoWriter.append(String.valueOf(sampleHetRate));
				sampleInfoWriter.append('\n');
				++j;
			}

			sampleInfoWriter.close();

			double averageVariantCallRate = sumVariantCallRate / variantCount;
			double averageSampleCallRate = sumSampleCallRate / inputData.getSamples().size();
			
			System.out.println("Samples: " + inputData.getSampleNames().length);
			System.out.println("Variants: " + variantCount);
			
			System.out.println("Average variant call rate: " + averageVariantCallRate);
			System.out.println("Average sample call rate: " + averageSampleCallRate);
			
			LOGGER.info("Done writing genotype info");



		} catch (ParseException ex) {
			LOGGER.fatal("Invalid command line arguments: ");
			LOGGER.fatal(ex.getMessage());
			System.err.println();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(" ", OPTIONS);
		}

	}
}
