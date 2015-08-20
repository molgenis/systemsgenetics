/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.text.NumberFormat;
import static nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression.BinomialEntry.BinomialEntry;
import static nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression.readGenoAndAsFromIndividual.readGenoAndAsFromIndividual;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;
import org.jdom.IllegalDataException;
import org.molgenis.genotype.GenotypeInfo;


/**
 *
 * @author Adriaan van der Graaf
 */
public class mainEntryPoint {
    /**
     * This is the entrypoint for cell type specific ASE,
     * This will provide the logic, based on what to do based on command line statements.
     * @param args
    */
    private static final Logger LOGGER;
    private static final Options OPTIONS;
    
    public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();

    static {

		LOGGER = Logger.getLogger(GenotypeInfo.class);

		OPTIONS = new Options();

		Option option;

                /*
                    Required Arguments (not optionally)
                */
                
		option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Determine what to do in this program, currently the following options are available:\n"+
                                                 "\t1\tDetermine Allele specific reads per SNP:     \tASREADS\n\n" +  
                                                 "\t2\tPerform a binomial test for ASE:             \tBINOMTEST\n" + 
                                                 "\t3\tEstimate per sample beta binomoverdispersion:\tBETABINOMTEST\n\n" +
                                                 "\t4\tCell type specific Binomial test:            \tCTSBINOMTEST\n" +
                                                 "\t5\tCell type specific Beta Binomial test        \tCTSBETABINOMTEST\n " +
                                                 "\t6\tPhased entry, currently used for testing     \tPHASED"+
                                                 "Please Run an option based on the number in the first column, or through the name in the third column."
                                                )
				.withLongOpt("action")
				.isRequired()
				.create('A');
		OPTIONS.addOption(option);
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("location of file in which to write output.\n "
                                               + "When action is: ASREADS, this will be the output of as_file\n"
                                               + "When action other, it will provide the test statistics per SNP.\n"
                                               + "When test contains dispersion estimatino, dispersion will be saved as: "
                                               + "<output-no-extension>_dispersion.txt")                                                
				.withLongOpt("output")
                                .isRequired()
				.create('O');
		OPTIONS.addOption(option);
                
                /*
                    Arguments that are conditionally required 
                        (based on the specific action that is being done)
                */
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("location of file where couplindg data is stored.\n "
                                               + "Required when action is: ASREADS\n"
                                               + "Coupling data will have the individual name in the first column and \n"
                                               + "the associated sample bam file name (without path and extenstion)in the second column")                                                
				.withLongOpt("coupling_file")
				.create('C');
		OPTIONS.addOption(option);
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Path to directory from which to load genotype data.\n"
                                                + "Currently only trityper data is available"
                                                + "Required when action is: ASREADS")                                                
				.withLongOpt("genotype_path")
				.create('G');
		OPTIONS.addOption(option);
                
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Location of bamfile from which to load data.\n "
                                               + "Required when action is: ASREADS")                                                
				.withLongOpt("bam_file")
				.create('B');
		OPTIONS.addOption(option);
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Location of  a file containing paths from which to load as data created in .\n "
                                               + "Required when action is: BINOMTEST, BETABINOMTEST, CTSBINOMTEST, CTSBETABINOMTEST.")                                                
				.withLongOpt("as_locations")
				.create('L');
		OPTIONS.addOption(option);                
                
                
                option = OptionBuilder.withArgName("string")
				.hasArgs()
				.withDescription("Path to a file phenotype information from which cell proportions are loaded.\n "
                                               + "Required when action is: CTSBINOMTEST,CTSBETABINOMTEST.")                                                
				.withLongOpt("pheno_file")
				.create('P');
		OPTIONS.addOption(option);
                
                
                /*
                    fully optionally arguments
                */
                
                option = OptionBuilder.withArgName("String")
				.hasArgs()
				.withDescription("Location of a file containing a list of snp names (rs123456789 format)\n"
                                               + "Standard setting is none.\n"
                                               + "Used when action is: ASREADS.")                                                
				.withLongOpt("snp_list")
				.create("snp_list");
		OPTIONS.addOption(option);
                
                
                
                option = OptionBuilder.withArgName("String")
				.hasArgs()
				.withDescription("Integer specifying how many heterozygotes should present before to run an AS test.\n"
                                               + "Standard setting is 1, minimum should be 0 but is not checked.\n"
                                               + "Used when action is: BINOMTEST, BETABINOMTEST, CTSBINOMTEST, CTSBETABINOMTEST.")                                                
				.withLongOpt("minimum_heterozygotes")
				.create("minimum_hets");
		OPTIONS.addOption(option);                
                
                
                option = OptionBuilder.withArgName("String")
				.hasArgs()
				.withDescription("Integer specifying how many reads should be overlapping before to run an AS test.\n"
                                               + "Standard setting is 10, minimum should be 0 but is not checked.\n"
                                               + "Used when action is: BINOMTEST, BETABINOMTEST, CTSBINOMTEST, CTSBETABINOMTEST.")                                                
				.withLongOpt("minimum_reads")
				.create("minimum_reads");
		OPTIONS.addOption(option);                
                
                
                

	}
    
    public static void main(String... args) throws Exception{
        
        //Required Arguments
        String outputLocation = new String();
        
        //ASREADS specific arguments
        String bamFile = new String();
        String couplingLocation = new String();
        String genotypeLocation = new String();
        String snpsLocation = new String();
        
                
        //BINOMTEST and BETABINOMTEST specific arguments
        String asFile = new String();

        //Cell type specific locations
        String phenoTypeLocation = new String();
        
        try {
            CommandLineParser parser = new PosixParser();
            final CommandLine commandLine = parser.parse(OPTIONS, args, true);
            
            try{
                //Read outputLocation
                
                if(commandLine.hasOption('O')){
                    outputLocation = commandLine.getOptionValue('O');
                } else{
                    throw new ParseException("Required command line input: --output ");
                }
                
                // Optional arguments that are not passed to the Entry constructors
                // But are saved in the GlobalVariables class.
                
                if(commandLine.hasOption("minimum_hets")){
                    GlobalVariables.minHets = Integer.parseInt(commandLine.getOptionValue("minimum_hets"));
                    //Check if this is bigger than 0, otherwise exit 
                    if(GlobalVariables.minHets <= 0){
                        throw new IllegalDataException("Minimum Number of hets cannot be smaller than one for AS testing\n"
                                         + "Exitting");
                    }
                    
                }
                
                if(commandLine.hasOption("minimum_reads")){
                    GlobalVariables.minReads = Integer.parseInt(commandLine.getOptionValue("minimum_reads"));
                    //Check if this is bigger than 0, otherwise exit
                    if(GlobalVariables.minReads <= 0){
                        throw new IllegalDataException("Minimum Number of reads cannot be smaller than one for AS testing\n"
                                         + "Exitting");
                    }
                }
                
                if(commandLine.hasOption('A')){
                    String programAction = commandLine.getOptionValue('A').toUpperCase();
                    
                    if(programAction.equals("ASREADS") || programAction.equals("1")){
        
                        //Do the AS determination part of the program
                        
                        //read genotype from options 
                        if(commandLine.hasOption('G')){
                             genotypeLocation = commandLine.getOptionValue('G');
                        } else {
                            throw new ParseException("Required command line input --genotype_path when --action is ASREADS");
                        }
                        
                        //Read binomTest arguments
                        if(commandLine.hasOption('C')){
                            couplingLocation = commandLine.getOptionValue('C');
                        } else{
                            throw new ParseException("Required command line input --coupling_file when --action is ASREADS");
                        }                        
                         
                        if(commandLine.hasOption('B')){
                            bamFile = commandLine.getOptionValue('B');
                        } else{
                            throw new ParseException("Required command line input --bam_file when --action is ASREADS");
                        }         
                        if(commandLine.hasOption("snp_list")){
                            snpsLocation = commandLine.getOptionValue("snp_list");
                        }else{
                            snpsLocation = "";
                        
                        }      
                        
                        
                        
                        /*
                            START reading for AS reads.
                        */
                        readGenoAndAsFromIndividual(bamFile, genotypeLocation, couplingLocation, outputLocation, snpsLocation);
                        
                    
                    }else if(programAction.equals("BINOMTEST") || programAction.equals("2")){
                        //Do a binomial test
                        
                        //Read het location
                        if(commandLine.hasOption('L')){
                            asFile = commandLine.getOptionValue('L');
                        } else{
                            throw new ParseException("Required command line input --as_location when --action is BINOMTEST");
                        }
                        
                        
                        /*
                            START BINOMIAL
                        */
                                
                        BinomialEntry(asFile, outputLocation);
                        
                    }else if(programAction.equals("BETABINOMTEST") || programAction.equals("3")){
                        //Determine Allele specific reads per individual
                        
                        if(commandLine.hasOption('L')){
                            asFile = commandLine.getOptionValue('L');
                        } else{
                            throw new ParseException("Required command line input --as_location when --action is BETABINOMTEST");
                        }
                        
                    
                        /*
                            START BETABINOMIAL
                        */
                        BetaBinomEntry  a = new BetaBinomEntry(asFile, outputLocation);
                    
                        
                        
                    }else if(programAction.equals("CTSBINOMTEST") || programAction.equals("4")){
                        

                        if(commandLine.hasOption('L')){
                            asFile = commandLine.getOptionValue('L');
                        } else{
                            throw new ParseException("Required command line input --as_location when --action is CTSBINOMTEST");
                        }
                        if(commandLine.hasOption('P')){
                            phenoTypeLocation = commandLine.getOptionValue('P');
                        } else{
                            throw new ParseException("Required command line input --pheno_file when --action is CTSBINOMTEST");
                        }
                        
                       
                        /*
                            START BINOMIAL CELL TYPE SPECIFIC TEST
                        */
                        
                        CTSbinomialEntry a =  new CTSbinomialEntry(asFile, phenoTypeLocation,  outputLocation);
                        
                    }else if(programAction.equals("CTSBETABINOMTEST") || programAction.equals("5")){

                        if(commandLine.hasOption('L')){
                            asFile = commandLine.getOptionValue('L');
                        } else{
                            throw new ParseException("Required command line input --as_location when --action is CTSBETABINOMTEST");
                        }
                        if(commandLine.hasOption('P')){
                            phenoTypeLocation = commandLine.getOptionValue('P');
                        } else{
                            throw new ParseException("Required command line input --pheno_file when --action is CTSBETABINOMTEST");
                        }
                  
                        /*
                            START BETA BINOMIAL CELL TYPE SPECIFIC TEST
                        */
                        
                        CTSbetaBinomialEntry a =  new CTSbetaBinomialEntry(asFile, phenoTypeLocation,  outputLocation);
                    
                    }else if(programAction.equals("PHASED") || programAction.equals("6")){
                        
                        if(commandLine.hasOption('L')){
                            asFile = commandLine.getOptionValue('L');
                        } else{
                            throw new ParseException("Required command line input --as_location when --action is PHASEDANALYSIS");
                        }
                        if(commandLine.hasOption('C')){
                            couplingLocation = commandLine.getOptionValue('C');
                        } else{
                            throw new ParseException("Required command line input --coupling_file when --action is PHASEDANALYSIS");
                        }           
                        
                       PhasedEntry a = new PhasedEntry(asFile, couplingLocation);
                        
                        
                    }else{
                        throw new ParseException("Unable to determine what to do. Please specify a correct value to --Action");
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
