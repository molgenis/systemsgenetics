package main.java.decon_eQTL_simple;

import java.io.IOException;
import java.util.List;

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
		commandLineOptions.parseCommandLine(args);
		Deconvolution deconvolution = new Deconvolution(commandLineOptions);
		deconvolution.readInputData();
		List<DeconvolutionResult> deconvolutionResults = deconvolution.runDeconPerGeneSnpPair();
		deconvolution.writeDeconvolutionResults(deconvolutionResults);

	}

}
