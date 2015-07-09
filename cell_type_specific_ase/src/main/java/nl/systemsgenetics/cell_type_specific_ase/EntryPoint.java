/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cell_type_specific_ase;

import java.text.NumberFormat;
import java.util.regex.Pattern;
import static nl.systemsgenetics.cell_type_specific_ase.BinomialEntry.BinomialEntry;
import static nl.systemsgenetics.cell_type_specific_ase.readGenoAndAsFromIndividual.readGenoAndAsFromIndividual;
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
     * This will provide the logic, based on what to do based on command line statements.
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

		option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Determine what to do in this program, currently the following options are available:\n"+
                                                 "\t1\tDetermine Allele specific reads per SNP:     \tASREADS\n" +  
                                                 "\t2\tPerform a binomial test for ASE:             \tBINOMTEST\n" + 
                                                 "\t3\tEstimate per sample beta binomoverdispersion:\tBETABINOMTEST\n" +
                                                 "Please Run an option based on the number in the first column, or through the name in the third column."
                                                )
				.withLongOpt("action")
				.isRequired()
				.create('A');
		OPTIONS.addOption(option);
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Path to bamfile from which to load data.\n "
                                               + "Required when action is: ASREADS")                                                
				.withLongOpt("bam_file")
				.create('B');
		OPTIONS.addOption(option);
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Path to a file containing paths from which to load data for testing.\n "
                                               + "Required when action is: BINOMTEST, BETABINOMTEST.")                                                
				.withLongOpt("as_locations")
				.create('L');
		OPTIONS.addOption(option);                
                

	}
    
    public static void main(String[] args) throws Exception{
        
        String bamFile = new String();
        String asLocations  = new String();
        String asFile = new String();
        
        try {
            CommandLineParser parser = new PosixParser();
            final CommandLine commandLine = parser.parse(OPTIONS, args, true);
            
            try{
                if(commandLine.hasOption('A')){
                    String programAction = commandLine.getOptionValue('A').toUpperCase();
                    
                    if(programAction.equals("ASREADS") || programAction.equals("1")){
        
                        //Do the AS part of the program
                        if(commandLine.hasOption('B')){
                            bamFile = commandLine.getOptionValue('B');
                        } else{
                            throw new ParseException("Required input --bamFile when --action is ASREADS");
                        }
                        //TODO, add the genotype folder location, and add the coupling file location.
                        
                        //START reading for AS reads.
                        readGenoAndAsFromIndividual(bamFile);
                        
                    
                    }else if(programAction.equals("BINOMTEST") || programAction.equals("2")){
                        //Do a binomial test
                        
                        if(commandLine.hasOption('L')){
                            asLocations = commandLine.getOptionValue('L');
                        } else{
                            throw new ParseException("Required input --as_location when --action is BINOMTEST");
                        }
                        
                        //START BINOMIAL
                        BinomialEntry(asLocations);
                        
                    }else if(programAction.equals("BETABINOMTEST") || programAction.equals("3")){
                        //Determine Allele specific reads per individual
                        
                        if(commandLine.hasOption('L')){
                            asFile = commandLine.getOptionValue('L');
                        } else{
                            throw new ParseException("Required input --as_location when --action is BETABINOMTEST");
                        }
                        
                        //TODO maybe fix this
                        
                        //START BETABINOMIAL
                        BetaBinomEntry  a = new BetaBinomEntry(asFile);
                    
                    }else{
                        throw new ParseException("Unable to determine what to do. Please specify --Action");
                    }
                }
                
                
                
            }catch (ParseException ex) {
                LOGGER.fatal("Invalid command line arguments");
                LOGGER.fatal(ex.getMessage());
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp(" ", OPTIONS);
            
            }
            
        
        } catch (ParseException ex) {
            LOGGER.fatal("Invalid command line arguments: ");
            LOGGER.fatal(ex.getMessage());
            System.err.println();
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(" ", OPTIONS);
        }
    }
}
