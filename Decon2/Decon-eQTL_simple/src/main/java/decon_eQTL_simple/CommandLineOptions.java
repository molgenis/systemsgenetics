package main.java.decon_eQTL_simple;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class CommandLineOptions {
	private String expressionFile;
	private String genotypeFile;
	private String cellCountFile;
	private String snpsToTestFile;
	private String outfile = "deconvolutionResults.csv";
	private String outfolder;
	private String programVersion;
	
	/**
	 * Standard command line parsing.
	 * 
	 * @param args A string vector of all arguments given to the command line, e.g. for `java -jar Deconvolution.jar -o` args = ["-o"]
	 * 
	 * @throws ParseException	If command line options can not be parsed
	 * @throws IOException 
	 */
	public void parseCommandLine(String[] args) throws ParseException, IOException {
		// load the properties file so that version can be printed
		Properties properties = new Properties();
		properties.load(this.getClass(). getClassLoader().getResourceAsStream("project.properties"));
		programVersion = properties.getProperty("artifactId")+"_"+properties.getProperty("version");
		Options options = new Options();
		Option help = new Option("help", "print this message");
		Option cellcount = Option.builder("c").required(true).hasArg().longOpt("cellcount").desc("Cellcount file name")
				.argName("file").build();
		Option expression = Option.builder("e").required(true).hasArg().longOpt("expression")
				.desc("Expression file name").argName("file").build();
		Option genotype = Option.builder("g").required(true).hasArg().longOpt("genotype").desc("Genotype file name")
				.argName("file").build();
		Option outfolder = Option.builder("o").required(true).hasArg().longOpt("outfolder").desc("Path to folder to write output to")
				.argName("path").build();
		Option outfile = Option.builder("of").required(false).hasArg().longOpt("outfile").desc("Outfile name of deconvolution results (will be written in outfolder)")
				.argName("file").build();
		Option snpsToTestOption = Option.builder("sn").required(true).hasArg().longOpt("snpsToTest").argName("file")
				.desc("Tab delimited file with first column gene name, second column SNP name. Need to match with names from genotype and expression files.").build();
		Option version = Option.builder("v").required(false).longOpt("version")
				.desc("Print the version of the program").build();
		options.addOption(help);
		options.addOption(outfile);
		options.addOption(expression);
		options.addOption(genotype);
		options.addOption(cellcount);
		options.addOption(outfolder);
		options.addOption(snpsToTestOption);
		options.addOption(version);
		CommandLineParser cmdLineParser = new DefaultParser();
		try{
			CommandLine cmdLine = cmdLineParser.parse(options, args);
			if (cmdLine.hasOption("help")) {
				// automatically generate the help statement
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp("deconvolution", options, true);
			}
			parseOptions (cmdLine);
			printArgumentValues(cmdLine);
		}
			catch(MissingOptionException e){
				HelpFormatter formatter = new HelpFormatter();
				DeconvolutionLogger.log.info(e.toString());
				formatter.printHelp("deconvolution", options, true);
				System.exit(0);
		}
	}
	
	private void parseOptions(CommandLine cmdLine) throws IOException{
		expressionFile = cmdLine.getOptionValue("expression");
		genotypeFile = cmdLine.getOptionValue("genotype");
		cellCountFile = cmdLine.getOptionValue("cellcount");
		snpsToTestFile = cmdLine.getOptionValue("snpsToTest");
		// check if all input files exist before starting the program to return error as early as possible
		if(!new File(expressionFile).exists() || new File(expressionFile).isDirectory()) { 
		    throw new FileNotFoundException(expressionFile+" does not exist");
		}
		if(!new File(genotypeFile).exists() || new File(genotypeFile).isDirectory()) { 
		    throw new FileNotFoundException(genotypeFile+" does not exist");
		}
		if(!new File(cellCountFile).exists() || new File(cellCountFile).isDirectory()) { 
		    throw new FileNotFoundException(cellCountFile+" does not exist");
		}
		if(!new File(snpsToTestFile).exists() || new File(snpsToTestFile).isDirectory()) { 
		    throw new FileNotFoundException(snpsToTestFile+" does not exist");
		}
				
		if (cmdLine.hasOption("outfile")) {
			outfile = cmdLine.getOptionValue("outfile");
		}
		
		outfolder = cmdLine.getOptionValue("outfolder")+"/";
				
		if (cmdLine.hasOption("version")) {
			DeconvolutionLogger.log.info("Version: "+programVersion);
		}
	}
	

	private void printArgumentValues(CommandLine cmdLine){
	    try {
	    	File outfolderDir = new File(outfolder);
	    	// if the directory does not exist, create it
	    	Boolean dirDidNotExist = false;
	    	if (!outfolderDir.exists()) {
	    		outfolderDir.mkdir();
	    		dirDidNotExist = true;
	    	}
	    	DeconvolutionLogger.setup(outfolder, false);
	    	if(dirDidNotExist){
	    		DeconvolutionLogger.log.info("Created directory "+outfolder);
	    	}
	    	DeconvolutionLogger.log.info("Writing output and logfile to "+outfolder);
	    } catch (IOException e) {
	      e.printStackTrace();
	      throw new RuntimeException("Problems with creating the log files");
	    }
		
	    File outfolderDir = new File(outfolder);
		// if the directory does not exist, create it
		Boolean dirDidNotExist = false;
		if (!outfolderDir.exists()) {
			outfolderDir.mkdir();
			dirDidNotExist = true;
		}
		if(dirDidNotExist){
			DeconvolutionLogger.log.info("Created directory "+outfolder);
		}
	    DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	    Date date = new Date();
	    DeconvolutionLogger.log.info("Starting deconvolution");
	    DeconvolutionLogger.log.info(dateFormat.format(date));
	    DeconvolutionLogger.log.info(String.format("Running deconvolution version: %s", programVersion));
	    DeconvolutionLogger.log.info("======= DECONVOLUTION paramater settings =======");
		DeconvolutionLogger.log.info(String.format("Expression file (-e): %s", expressionFile));
		DeconvolutionLogger.log.info(String.format("Genotype file (-g): %s", genotypeFile));
		DeconvolutionLogger.log.info(String.format("Cellcount file (-c): %s", cellCountFile));
		DeconvolutionLogger.log.info(String.format("SNPs to test file (-sn): %s", snpsToTestFile));
		DeconvolutionLogger.log.info(String.format("Outfolder (-o): %s", outfolder));
		DeconvolutionLogger.log.info(String.format("Outfile (-of): %s", outfile));
		DeconvolutionLogger.log.info("=================================================");
	}
	public String getExpressionFile(){
		return (expressionFile);
	}
	public String getSnpsToTestFile(){
		return (snpsToTestFile);
	}	
	public String getGenotypeFile(){
		return genotypeFile;
	}
	public String getCellcountFile(){
		return cellCountFile;
	}
	public String getOutfile(){
		return outfile;
	}

	public String getOutfolder() throws IllegalAccessException{
		if(this.outfolder == null){
			throw new IllegalAccessException("Outfolder has not been set");
		}
		return outfolder;
	}

	public void setOutfolder(String newOutfolder){
		outfolder = newOutfolder;
	}	
}





