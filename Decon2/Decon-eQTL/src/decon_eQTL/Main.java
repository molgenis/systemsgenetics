package decon_eQTL;

import java.io.IOException;

import org.apache.commons.cli.ParseException;

public class Main {


	/**
	 * Deconvolute a set of QTLs given the expression levels, genotypes,
	 * and cell counts. Calculates the p-values for the deconvoluted QTLs
	 * and writes them to an out file
	 * 
	 * @param args List of command line arguments
	 * 
	 * @throws ParseException	If cell count file is not in right format to be parsed correctly
	 * @throws IllegalAccessException	If out folder can not be retrieved from commandLineOptions
	 * @throws IOException	If cell counts file can not be found or read
	 */
	public static void main(String[] args) throws ParseException, IllegalAccessException, IOException {
		CommandLineOptions commandLineOptions = new CommandLineOptions(); 
		Deconvolution deconvolution = new Deconvolution();
		commandLineOptions.parseCommandLine(args);
		deconvolution.runDeconPerGeneSnpPair(commandLineOptions);
	}

}
