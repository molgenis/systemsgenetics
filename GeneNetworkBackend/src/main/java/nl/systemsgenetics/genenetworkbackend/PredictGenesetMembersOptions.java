/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend;

import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author patri
 */
public class PredictGenesetMembersOptions {

	private static final Options OPTIONS;

	private final File eigenVectorFile;
	private final File pathwayMatrixFile;
	private final File predictionsFile;
	private final File backgroudGenesFile;

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The eigenvectors to use for prediction");
		OptionBuilder.withLongOpt("eigen");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("e"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The annototed pathway matrix");
		OptionBuilder.withLongOpt("pathways");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The pathway predictions per gene");
		OptionBuilder.withLongOpt("predictions");
		OPTIONS.addOption(OptionBuilder.create("o"));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("A file with genes to use as background for predictions");
		OptionBuilder.withLongOpt("background");
		OPTIONS.addOption(OptionBuilder.create("b"));

	}

	public PredictGenesetMembersOptions(String... args) throws ParseException {

		CommandLineParser parser = new PosixParser();
		final CommandLine commandLine = parser.parse(OPTIONS, args, false);
		
		eigenVectorFile = new File(commandLine.getOptionValue("e"));
		pathwayMatrixFile = new File(commandLine.getOptionValue("p"));
		
		predictionsFile = new File(commandLine.getOptionValue("o", eigenVectorFile.getPath() + ".GenesetZScores.txt"));
		
		backgroudGenesFile = commandLine.hasOption("b") ? new File(commandLine.getOptionValue("b")) : null;

	}
	
	public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
    }
	
	public void printOptions() {
		
		System.out.println(" - Eigenvector matrix: " + eigenVectorFile.getPath());
		System.out.println(" - Pathway matrix: " + pathwayMatrixFile.getPath());
		System.out.println(" - Pathway predictions: " + predictionsFile.getPath());
		
		if(backgroudGenesFile == null){
			System.out.println(" - Not using background genes");
		} else {
			System.out.println(" - Background gene list: " + backgroudGenesFile.getPath());
		}
		
	}


	public File getEigenVectorFile() {
		return eigenVectorFile;
	}

	public File getPathwayMatrixFile() {
		return pathwayMatrixFile;
	}

	public File getPredictionsFile() {
		return predictionsFile;
	}

	public File getBackgroudGenesFile() {
		return backgroudGenesFile;
	}
	
	

}
