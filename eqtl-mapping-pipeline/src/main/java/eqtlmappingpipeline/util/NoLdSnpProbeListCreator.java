package eqtlmappingpipeline.util;

import com.google.common.collect.Lists;
import eqtlmappingpipeline.ase.AseConfiguration;
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
import java.util.ArrayList;
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
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
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
        Option RefferenceTypeIn = OptionBuilder.withArgName("type").hasArg().withDescription("Type of refference data.").withLongOpt("RefferenceType").create("rt");
        Option RefferenceIn = OptionBuilder.withArgName("path").hasArg().withDescription("Location for the reference data").withLongOpt("RefferenceLocation").create("ri");
        Option BedIn = OptionBuilder.withArgName("path").hasArg().withDescription("Location of the input data").withLongOpt("input").create("i");
        Option WindowSize = OptionBuilder.withArgName("int").hasArg().withDescription("Half window size (default 250000).").withLongOpt("Window_Size").create("ws");
        Option ProbeMargin = OptionBuilder.withArgName("int").hasArg().withDescription("Additional probe margin (default 0).").withLongOpt("Probe_margin").create("pm");
        Option MaxDprime = OptionBuilder.withArgName("int").hasArg().withDescription("Cut-off point for max d'(default 0.2).").withLongOpt("max_dPrime").create("dp");
        Option MaxRsquare = OptionBuilder.withArgName("int").hasArg().withDescription("Cut-off point for max r2 (default 0.2).").withLongOpt("max_rSquare").create("r");
        options.addOption(FileOut).addOption(RefferenceTypeIn).addOption(RefferenceIn).addOption(WindowSize).addOption(BedIn).addOption(ProbeMargin).addOption(MaxDprime).addOption(MaxRsquare);
        
		File probeFile = null;
		String genotypePath = null;
		String genotypeType = null;
		int windowHalfSize = 250000;
		int probeMargin = 0;
		double maxDprime = 0.2;
		double maxR2 = 0.2;
		File outputFile = null;
        
        
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
            HelpFormatter formatter = new HelpFormatter();
            
            if (cmd.hasOption("OutputFile") || cmd.hasOption("o")) {
            // initialise the member variable
                outputFile = new File(cmd.getOptionValue("OutputFile"));
            } else {
                System.out.println("Missing necesarray information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("RefferenceLocation")|| cmd.hasOption("ri")) {
                // initialise the member variable
                genotypePath = cmd.getOptionValue("RefferenceLocation");
            } else {
                System.out.println("Missing necesarray information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("RefferenceType")|| cmd.hasOption("rt")) {
                // initialise the member variable
                genotypeType = cmd.getOptionValue("RefferenceType");
            } else {
                System.out.println("Missing necesarray information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("input") || cmd.hasOption("i")) {
            // initialise the member variable
                probeFile = new File(cmd.getOptionValue("input"));
            } else {
                System.out.println("Missing necesarray information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("Window_Size")|| cmd.hasOption("ws")) {
                // initialise the member variable
                windowHalfSize = Integer.parseInt(cmd.getOptionValue("Window_Size"));
            }
            if (cmd.hasOption("Probe_margin")|| cmd.hasOption("pm")) {
                // initialise the member variable
                probeMargin = Integer.parseInt(cmd.getOptionValue("Probe_margin"));
            }
            if (cmd.hasOption("MaxDprime")|| cmd.hasOption("md")) {
                // initialise the member variable
                maxDprime = Double.parseDouble(cmd.getOptionValue("MaxDprime"));
            }
            if (cmd.hasOption("max_rSquare")|| cmd.hasOption("mr")) {
                // initialise the member variable
                maxR2 = Double.parseDouble(cmd.getOptionValue("max_rSquare"));
            }
        } catch (ParseException ex) {
            Logger.getLogger(NoLdSnpProbeListCreator.class.getName()).log(Level.SEVERE, null, ex);
        }
       
		RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOf(genotypeType).createGenotypeData(genotypePath, 10000);

		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(probeFile), "UTF-8"));
		final BufferedWriter snpProbeToTestWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), "UTF-8"));


		String line;
		String[] elements;
		while ((line = reader.readLine()) != null) {
			elements = StringUtils.splitPreserveAllTokens(line, '\t');

			if (elements.length < 4) {
				throw new Exception();
			}

			String chr = elements[0].replace("chr", "");//TODO remove chr
			int probeStartPos = Integer.parseInt(elements[1]);
			int probeStopPos = Integer.parseInt(elements[2]);
			String probeName = elements[3];


			int windowsStart = probeStartPos - windowHalfSize;
            windowsStart = windowsStart < 0 ? 0 : windowsStart;
			int windowStop = probeStopPos + windowHalfSize;


			ArrayList<GeneticVariant> probeVariants = Lists.newArrayList(genotypeData.getVariantsByRange(chr, (probeStartPos - probeMargin), (probeStopPos + probeMargin)));

			variants:
			for (GeneticVariant variant : genotypeData.getVariantsByRange(chr, windowsStart, windowStop)) {

				//skip over probe variants
				if (variant.getStartPos() >= (probeStartPos - probeMargin) && variant.getStartPos() <= (probeStopPos + probeMargin)) {
					continue variants;
				}

				for (GeneticVariant probeVariant : probeVariants) {

					Ld ld = variant.calculateLd(probeVariant);

					if (ld.getDPrime() >= maxDprime || ld.getR2() >= maxR2) {
						//Exclude
						continue variants;
					}

				}

				//Probe SNP combination is okay
                String variantPrimaryId = variant.getPrimaryVariantId();
                if(variantPrimaryId==null){
                    snpProbeToTestWriter.append(variant.getSequenceName()).append('-').append(String.valueOf(variant.getStartPos()));
                } else {
                    snpProbeToTestWriter.append(variantPrimaryId);
                }
                
				snpProbeToTestWriter.append('\t');
				snpProbeToTestWriter.append(probeName);
				snpProbeToTestWriter.append('\n');

			}

		}

        snpProbeToTestWriter.close();
	}

    private static void printHelp() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
