/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.NavigableMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeInfo;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author patri
 */
public class ModuleEqtlGeuvadisReplication {

	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |                                       |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Options OPTIONS;
	private static Logger LOGGER;

	static {

		LOGGER = Logger.getLogger(GenotypeInfo.class);

		OPTIONS = new Options();

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Genotypes for LD calculations");
		OptionBuilder.withLongOpt("genotypes");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("genotypesFormat");
		OPTIONS.addOption(OptionBuilder.create("G"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Path to output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('o'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Path to replication eQTL file");
		OptionBuilder.withLongOpt("eqtls");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('e'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("File with interaction QTLs");
		OptionBuilder.withLongOpt("interactions");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('i'));
		
		OptionBuilder.withArgName("double");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("LD cutoff");
		OptionBuilder.withLongOpt("ld");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("ld"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("window");
		OptionBuilder.withLongOpt("window");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("w"));

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, LdCalculatorException {

		System.out.println(HEADER);
		System.out.println();
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		CommandLineParser parser = new PosixParser();
		final CommandLine commandLine;
		try {
			commandLine = parser.parse(OPTIONS, args, true);
		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: " + ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		final String[] genotypesBasePaths = commandLine.getOptionValues("g");
		final RandomAccessGenotypeDataReaderFormats genotypeDataType;
		final String replicationQtlFilePath = commandLine.getOptionValue("e");
		final String interactionQtlFilePath = commandLine.getOptionValue("i");
		final String outputFilePath = commandLine.getOptionValue("o");
		final double ldCutoff = Double.parseDouble(commandLine.getOptionValue("ld"));
		final int window = Integer.parseInt(commandLine.getOptionValue("w"));
		
		System.out.println("Genotype: " + Arrays.toString(genotypesBasePaths));
		System.out.println("Interaction file: " + interactionQtlFilePath);
		System.out.println("Replication file: " + replicationQtlFilePath);
		System.out.println("Output: " + outputFilePath);
		System.out.println("LD: " + ldCutoff);
		System.out.println("Window: " + window);

		try {
			if (commandLine.hasOption("G")) {
				genotypeDataType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue("G").toUpperCase());
			} else {
				if (genotypesBasePaths[0].endsWith(".vcf")) {
					System.err.println("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					System.exit(1);
					return;
				}
				try {
					genotypeDataType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(genotypesBasePaths[0]);
				} catch (GenotypeDataException e) {
					System.err.println("Unable to determine input 1 type based on specified path. Please specify -G");
					System.exit(1);
					return;
				}
			}
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --genotypesFormat \"" + commandLine.getOptionValue("G") + "\" is not a valid input data format");
			System.exit(1);
			return;
		}

		final RandomAccessGenotypeData genotypeData;

		try {
			genotypeData = genotypeDataType.createFilteredGenotypeData(genotypesBasePaths, 100, null, null, null, 0.8);
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

		ChrPosTreeMap<EQTL> replicationQtls = new QTLTextFile(replicationQtlFilePath, false).readQtlsAsTreeMap();

		int interactionSnpNotInGenotypeData = 0;
		int noReplicationQtlsInWindow = 0;
		int noReplicationQtlsInLd = 0;
		int multipleReplicationQtlsInLd = 0;
		int replicationTopSnpNotInGenotypeData = 0;

		final CSVWriter outputWriter = new CSVWriter(new FileWriter(new File(outputFilePath)), '\t', '\0');
		final String[] outputLine = new String[11];
		int c = 0;
		outputLine[c++] = "Chr";
		outputLine[c++] = "Pos";
		outputLine[c++] = "SNP";
		outputLine[c++] = "Gene";
		outputLine[c++] = "Module";
		outputLine[c++] = "ReplicationZ";
		outputLine[c++] = "DiscoveryAlleleAssessed";
		outputLine[c++] = "ReplicationAlleleAssessed";
		outputLine[c++] = "bestLd";
		outputLine[c++] = "bestLd_dist";
		outputLine[c++] = "nextLd";
		outputWriter.writeNext(outputLine);

		CSVReader interactionQtlReader = new CSVReader(new FileReader(interactionQtlFilePath),'\t');
		interactionQtlReader.readNext();//skip header
		String[] interactionQtlLine;
		while ((interactionQtlLine = interactionQtlReader.readNext()) != null) {

			String snp = interactionQtlLine[1];
			String chr = interactionQtlLine[2];
			int pos = Integer.parseInt(interactionQtlLine[3]);
			String gene = interactionQtlLine[4];
			String alleleAssessed = interactionQtlLine[9];
			String module = interactionQtlLine[12];

			GeneticVariant interactionQtlVariant = genotypeData.getSnpVariantByPos(chr, pos);

			if (interactionQtlVariant == null) {
				System.err.println("Interaction QTL SNP not found in genotype data: " + chr + ":" + pos);
				++interactionSnpNotInGenotypeData;
				continue;
			}

			NavigableMap<Integer, EQTL> potentionalReplicationQtls = replicationQtls.getChrRange(chr, pos - window, true, pos + window, true);

			EQTL bestMatch = null;
			double bestMatchLd = Double.NaN;
			double nextBestLd = Double.NaN;
			
			

			for (EQTL potentialReplicationQtl : potentionalReplicationQtls.values()) {
				GeneticVariant potentialReplicationQtlVariant = genotypeData.getSnpVariantByPos(potentialReplicationQtl.getRsChr().toString(), potentialReplicationQtl.getRsChrPos());

				if (potentialReplicationQtlVariant == null) {
					System.err.println("Not found: " + potentialReplicationQtl.getRsChr().toString() + ":" + potentialReplicationQtl.getRsChrPos());
					++replicationTopSnpNotInGenotypeData;
					continue;
				}
				
				double ld = interactionQtlVariant.calculateLd(potentialReplicationQtlVariant).getR2();
				
				if(bestMatch == null){
					bestMatch = potentialReplicationQtl;
					bestMatchLd = ld;
				} else if (ld > bestMatchLd){
					bestMatch = potentialReplicationQtl;
					nextBestLd = bestMatchLd;
					bestMatchLd = ld;
				}				
				
			}
			
			c = 0;
			outputLine[c++] = chr;
			outputLine[c++] = String.valueOf(pos);
			outputLine[c++] = snp;
			outputLine[c++] = gene;
			outputLine[c++] = module;
			outputLine[c++] = bestMatch == null ? "NA" : String.valueOf(bestMatch.getZscore());
			outputLine[c++] = alleleAssessed;
			outputLine[c++] = bestMatch == null ? "NA" : String.valueOf(bestMatch.getAlleleAssessed());
			outputLine[c++] = String.valueOf(bestMatchLd);
			outputLine[c++] = bestMatch == null ? "NA" : String.valueOf(Math.abs(pos - bestMatch.getRsChrPos()));
			outputLine[c++] = String.valueOf(nextBestLd);
			outputWriter.writeNext(outputLine);

		}
		
		outputWriter.close();

		System.out.println("interactionSnpNotInGenotypeData: " + interactionSnpNotInGenotypeData);
		System.out.println("noReplicationQtlsInWindow: " + noReplicationQtlsInWindow);
		System.out.println("noReplicationQtlsInLd: " + noReplicationQtlsInLd);
		System.out.println("multipleReplicationQtlsInLd: " + multipleReplicationQtlsInLd);
		System.out.println("replicationTopSnpNotInGenotypeData: " + replicationTopSnpNotInGenotypeData);

	}

}
