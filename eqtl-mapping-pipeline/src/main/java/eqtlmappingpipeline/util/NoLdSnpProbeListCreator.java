package eqtlmappingpipeline.util;

import com.google.common.collect.Lists;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilterBiAllelic;
import org.molgenis.genotype.variantFilter.VariantFilterSeq;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantQcChecker;

/**
 * @author Patrick Deelen
 */
public class NoLdSnpProbeListCreator {
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws UnsupportedEncodingException, FileNotFoundException, IOException, Exception {
		CommandLineParser parser = new GnuParser();
		
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}
		
		Options options = new Options();
		
		Option FileOut = OptionBuilder.withArgName("path").hasArg().withDescription("Location and name of the output file.").withLongOpt("OutputFile").create("o");
		Option ReferenceTypeIn = OptionBuilder.withArgName("type").hasArg().withDescription("Type of reference data.").withLongOpt("ReferenceType").create("rt");
		Option ReferenceIn = OptionBuilder.withArgName("path").hasArg().withDescription("Location for the reference data").withLongOpt("ReferenceLocation").create("ri");
		Option InFile = OptionBuilder.withArgName("path").hasArg().withDescription("Location of the input data").withLongOpt("input").create("i");
		Option WindowSize = OptionBuilder.withArgName("int").hasArg().withDescription("Half window size (default 250000).").withLongOpt("Window_Size").create("ws");
		Option ProbeMargin = OptionBuilder.withArgName("int").hasArg().withDescription("Additional probe margin (default 0).").withLongOpt("Probe_margin").create("pm");
		Option MaxDprime = OptionBuilder.withArgName("int").hasArg().withDescription("Cut-off point for max d'(default 0.2).").withLongOpt("max_dPrime").create("dp");
		Option MaxRsquare = OptionBuilder.withArgName("int").hasArg().withDescription("Cut-off point for max r2 (default 0.2).").withLongOpt("max_rSquare").create("r");
		Option variantFilter = OptionBuilder.withArgName("String").hasArg().withDescription("List with variant ids to select.").withLongOpt("variant_filter").create("vf");
		Option callRate = OptionBuilder.withArgName("int").hasArg().withDescription("Call-rate cut-off.").withLongOpt("min_callRate").create("vc");
		Option mafFilter = OptionBuilder.withArgName("int").hasArg().withDescription("Minor allel cut-off filter.").withLongOpt("min_maf").create("mf");
		Option chrFilter = OptionBuilder.withArgName("string").hasArg().withDescription("Filter input data on chromosome").withLongOpt("chrFilter").create("ch");
		options.addOption(FileOut).addOption(ReferenceTypeIn).addOption(ReferenceIn).addOption(WindowSize).addOption(InFile).addOption(ProbeMargin).addOption(MaxDprime).addOption(MaxRsquare).addOption(variantFilter).addOption(callRate).addOption(mafFilter).addOption(chrFilter);
		
		File probeFile = null;
		String genotypePath = null;
		String genotypeType = null;
		int windowHalfSize = 250000;
		int probeMargin = 1;
		float cRate = 0.0f;
		float maf = 0.0f;
		double HWE = 0.0d;
		String variantFilterList = null;
		double maxDprime = 0.2;
		double maxR2 = 0.2;
		String chrF = null;
		File outputFile = null;
		
		String logOutput = "";
		
		CommandLine cmd;
		try {
			cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			
			if (cmd.hasOption("OutputFile") || cmd.hasOption("o")) {
				// initialise the member variable
				outputFile = new File(cmd.getOptionValue("OutputFile"));
				logOutput += "OutputFile:\t" + outputFile + "\n";
			} else {
				System.out.println("Missing necesarray information");
				formatter.printHelp("ant", options);
				System.exit(0);
			}
			if (cmd.hasOption("ReferenceLocation") || cmd.hasOption("ri")) {
				// initialise the member variable
				genotypePath = cmd.getOptionValue("ReferenceLocation");
				logOutput += "ReferenceLocation:\t" + genotypePath + "\n";
			} else {
				System.out.println("Missing necesarray information");
				formatter.printHelp("ant", options);
				System.exit(0);
			}
			if (cmd.hasOption("ReferenceType") || cmd.hasOption("rt")) {
				// initialise the member variable
				genotypeType = cmd.getOptionValue("ReferenceType");
				logOutput += "ReferenceType:\t" + genotypeType + "\n";
			} else {
				System.out.println("Missing necesarray information");
				formatter.printHelp("ant", options);
				System.exit(0);
			}
			if (cmd.hasOption("input") || cmd.hasOption("i")) {
				// initialise the member variable
				probeFile = new File(cmd.getOptionValue("input"));
				logOutput += "Input probemapping:\t" + probeFile + "\n";
			} else {
				System.out.println("Missing necesarray information");
				formatter.printHelp("ant", options);
				System.exit(0);
			}
			if (cmd.hasOption("Window_Size") || cmd.hasOption("ws")) {
				// initialise the member variable
				windowHalfSize = Integer.parseInt(cmd.getOptionValue("Window_Size"));
				logOutput += "Window Size:\t" + windowHalfSize + "\n";
			}
			if (cmd.hasOption("Probe_margin") || cmd.hasOption("pm")) {
				// initialise the member variable
				probeMargin = Integer.parseInt(cmd.getOptionValue("Probe_margin"));
				logOutput += "Probe margin:\t" + probeMargin + "\n";
			}
			if (cmd.hasOption("variant_filter") || cmd.hasOption("vf")) {
				// initialise the member variable
				variantFilterList = cmd.getOptionValue("variant_filter");
				logOutput += "Variant filter:\t" + variantFilterList + "\n";
			}
			if (cmd.hasOption("min_callRate") || cmd.hasOption("vc")) {
				// initialise the member variable
				cRate = Float.parseFloat(cmd.getOptionValue("min_callRate"));
				logOutput += "Min callRate:\t" + cRate + "\n";
			}
			if (cmd.hasOption("min_maf") || cmd.hasOption("mf")) {
				// initialise the member variable
				maf = Float.parseFloat(cmd.getOptionValue("min_maf"));
				logOutput += "Min maf:\t" + maf + "\n";
			}
			if (cmd.hasOption("MaxDprime") || cmd.hasOption("md")) {
				// initialise the member variable
				maxDprime = Double.parseDouble(cmd.getOptionValue("MaxDprime"));
				logOutput += "Max D':\t" + maxDprime + "\n";
			}
			if (cmd.hasOption("max_rSquare") || cmd.hasOption("mr")) {
				// initialise the member variable
				maxR2 = Double.parseDouble(cmd.getOptionValue("max_rSquare"));
				logOutput += "Max rSquare:\t" + maxR2 + "\n";
			}
			if (cmd.hasOption("chrFilter") || cmd.hasOption("ch")) {
				chrF = cmd.getOptionValue("chrFilter");
				logOutput += "Chr filter:\t" + chrF + "\n";
			}
		} catch (ParseException ex) {
			Logger.getLogger(NoLdSnpProbeListCreator.class.getName()).log(Level.SEVERE, null, ex);
		}
		
		System.out.println(logOutput);
		SampleFilter sf = null;
		
		VariantCombinedFilter varFilter = new VariantCombinedFilter();
		varFilter.add(new VariantFilterBiAllelic());
		
		if (cRate != 0.0f || maf != 0.0f || HWE != 0.0d) {
			VariantQcChecker qc = new VariantQcChecker(maf, cRate, HWE);
			varFilter.add(qc);
		}
		if (variantFilterList != null) {
			HashSet<String> snps = new HashSet<String>();
			
			BufferedReader variantIdFilterReader = new BufferedReader(new FileReader(variantFilterList));
			String line;
			while ((line = variantIdFilterReader.readLine()) != null) {
				snps.add(line);
			}
			VariantIdIncludeFilter snpIdFilter = new VariantIdIncludeFilter(snps);
			varFilter.add(snpIdFilter);
			
		}
		if (chrF != null) {
			VariantFilterSeq seqFilter = new VariantFilterSeq(chrF);
			varFilter.add(seqFilter);
		}
		RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOf(genotypeType).createFilteredGenotypeData(genotypePath, 10000, varFilter, sf);
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(probeFile), "UTF-8"));
		final BufferedWriter snpProbeToTestWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), "UTF-8"));
		
		int rSquareExclusions = 0;
		int dPrimeExclusions = 0;
		int possibleCombinations = 0;
		int locatedInProbe = 0;
		String line;
		String[] elements;
		
		int lineNr = 0;
		while ((line = reader.readLine()) != null) {
			lineNr++;
			elements = StringUtils.splitPreserveAllTokens(line, '\t');
			
			if ((probeFile.getName().endsWith(".bed") && elements.length < 4) || (probeFile.getName().endsWith(".txt"))) {
				throw new Exception("Incorrect file format.");
			}
			
			String chr;
			int probeStartPos;
			int probeStopPos;
			String probeName;
			
			if (probeFile.getName().endsWith(".bed")) {
				chr = elements[0].replace("chr", "");
				probeStartPos = Integer.parseInt(elements[1]);
				probeStopPos = Integer.parseInt(elements[2]);
				probeName = elements[3];
				
			} else if (probeFile.getName().endsWith(".txt") || probeFile.getName().endsWith(".tsv")) {
				if (lineNr == 1) {
					line = reader.readLine();
					elements = StringUtils.splitPreserveAllTokens(line, '\t');
				}
				
				if (elements.length < 6) {
					System.err.println("Expected at least 6 columns in the annotation file.");
					System.exit(-1);
				}
				
				chr = elements[3];
				probeStartPos = Integer.parseInt(elements[4]);
				probeStopPos = Integer.parseInt(elements[5]);
				probeName = elements[2];
			} else {
				throw new Exception("Incorrect file format.");
			}
			
			int probeCenter = probeStartPos + (int) (Math.floor((probeStopPos - probeStartPos) / 2.0d));
			
			int windowsStart = probeCenter - windowHalfSize;
			windowsStart = windowsStart < 0 ? 0 : windowsStart;
			int windowStop = probeCenter + windowHalfSize;
			
			ArrayList<GeneticVariant> probeVariants = Lists.newArrayList(genotypeData.getVariantsByRange(chr, (probeStartPos - probeMargin), (probeStopPos + probeMargin)));
			
			variants:
			for (GeneticVariant variant : genotypeData.getVariantsByRange(chr, windowsStart, windowStop)) {
				possibleCombinations++;
				
				//skip over probe variants
				if (variant.getStartPos() >= (probeStartPos - probeMargin) && variant.getStartPos() <= (probeStopPos + probeMargin)) {
					locatedInProbe++;
					continue variants;
				}
				
				if (!probeVariants.isEmpty()) {
					
					for (GeneticVariant probeVariant : probeVariants) {
						
						Ld ld = variant.calculateLd(probeVariant);
						
						if (ld.getDPrime() >= maxDprime) {
							//Exclude
							dPrimeExclusions++;
							continue variants;
						}
						if (ld.getR2() >= maxR2) {
							//Exclude
							rSquareExclusions++;
							continue variants;
						}
						
					}
				}
				//Probe SNP combination is okay
				String variantPrimaryId = variant.getPrimaryVariantId();
				if (variantPrimaryId == null) {
					snpProbeToTestWriter.append(variant.getSequenceName()).append('-').append(String.valueOf(variant.getStartPos()));
				} else {
					snpProbeToTestWriter.append(variantPrimaryId);
				}
				
				snpProbeToTestWriter.append('-');
				snpProbeToTestWriter.append(probeName);
				snpProbeToTestWriter.append('\n');
				
			}
			
		}
		
		snpProbeToTestWriter.close();
		System.out.println("Number of valid combination: " + possibleCombinations);
		System.out.println("Excluded because located in probe: " + locatedInProbe);
		System.out.println("Excluded based on D': " + dPrimeExclusions);
		System.out.println("Excluded based on r2: " + rSquareExclusions);
	}
}
