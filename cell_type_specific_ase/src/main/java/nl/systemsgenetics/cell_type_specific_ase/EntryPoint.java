/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

import java.text.NumberFormat;
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
import org.molgenis.genotype.GenotypeInfo;

/**
 *
 * @author Adriaan van der Graaf
 */
public class EntryPoint {
    /**
     * This is the entrypoint for cell type specific ASE,
     * In the future it will do the correct stuff based command line arguments, 
     * but currently it will just be used for testing purposes.
     * @param args
    */
    private static final Logger LOGGER;
    private static final Options OPTIONS;
    private static final Pattern CHR_POS_SPLITTER = Pattern.compile("\\s+|:");
    public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();

    static {

		LOGGER = Logger.getLogger(GenotypeInfo.class);

		OPTIONS = new Options();

		Option option;

		option = OptionBuilder.withArgName("action")
				.hasArgs()
				.withDescription("Determine what to do in this program, currently the following options are available:\n"+
                                                 "\t1\tDetermine Allele specific reads per SNP:     \tASREADS\n" +  
                                                 "\t2\tPerform a binomial test for ASE:             \tBINOMTEST\n" + 
                                                 "\t3\tEstimate per sample beta binomoverdispersion:\tBETABINOMDISP\n" +
                                                 "Please Run an option based on the number in the first column, or through the name in the third column."
                                                )
				.withLongOpt("action")
				.isRequired()
				.create("A");
		OPTIONS.addOption(option);


	}
    
    public static void main(String[] args){
        
        try {
            CommandLineParser parser = new PosixParser();
            final CommandLine commandLine = parser.parse(OPTIONS, args, true);
            
        
        } catch (ParseException ex) {
            LOGGER.fatal("Invalid command line arguments: ");
            LOGGER.fatal(ex.getMessage());
            System.err.println();
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", OPTIONS);
        }
    }
}
