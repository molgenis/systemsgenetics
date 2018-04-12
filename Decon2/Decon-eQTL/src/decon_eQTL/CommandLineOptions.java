package decon_eQTL;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

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
	private int minimumSamplesPerGenotype = 0;
	private Boolean roundDosage = false;
	private Boolean allDosages = false;
	private Boolean filterSamples = false;
	private Boolean removeConstraintViolatingSamples = false;
	private Boolean testRun = false;
	private Boolean skipGenotypes = false;
	private Boolean wholeBloodQTL = false;
	private Boolean noConsole = false;
	private Boolean outputPredictedExpression = false;
	private String genotypeConfigurationType = "all";
	
	/**
	 * Standard command line parsing.
	 * 
	 * @param args A string vector of all arguments given to the command line, e.g. for `java -jar Deconvolution.jar -o` args = ["-o"]
	 * 
	 * @throws ParseException	If command line options can not be parsed
	 * 
	 * @throws FileNotFoundException	If deconvolution log file can't be written
	 */
	public void parseCommandLine(String[] args) throws ParseException, FileNotFoundException {
		Options options = new Options();
		Option help = new Option("help", "print this message");
		Option allDosages = Option.builder("ad").required(false).longOpt("all_dosages")
				.desc("Filter out QTLs where not all dosages are present in at least 1 sample").build();
		Option cellcount = Option.builder("c").required(true).hasArg().longOpt("cellcount").desc("Cellcount file name")
				.argName("file").build();
		Option expression = Option.builder("e").required(true).hasArg().longOpt("expression")
				.desc("Expression file name").argName("file").build();
		Option filterSamplesOption =  Option.builder("f").required(false).longOpt("filter_samples")
				.desc("If set, remove samples that are filtered out because of -m, -nn or -ad. By default p-values of these are set to 333.0").build();
		Option genotype = Option.builder("g").required(true).hasArg().longOpt("genotype").desc("Genotype file name")
				.argName("file").build();
		Option genotypeConfigurationTypeOption = Option.builder("gc").required(false).hasArg().longOpt("genotypeConfigurationType")
				.desc("Which genotype configuration type to use (either all or two)").build();
		Option minimumSamplesPerGenotype = Option.builder("m").required(false).hasArg().longOpt("minimum_samples_per_genotype")
				.desc("The minimum amount of samples need for each genotype of a QTL for the QTL to be included in the results")
				.argName("int").build();
		Option noConsoleOption = Option.builder("no").required(false).longOpt("no_console")
				.desc("Do not output logging info to the console").build();
		Option outfolder = Option.builder("o").required(true).hasArg().longOpt("outfolder").desc("Path to folder to write output to")
				.argName("path").build();
		Option outputPredictedExpressionOption = Option.builder("oe").required(false).longOpt("outputPredictedExpression").desc("Write output file with predicted expression")
				.build();
		Option outfile = Option.builder("of").required(false).hasArg().longOpt("outfile").desc("Outfile name of deconvolution results (will be written in outfolder)")
				.argName("file").build();
		Option roundDosage = Option.builder("r").required(false).longOpt("round_dosage")
				.desc("Round the dosage to the closest int").build();
		Option skipGenotypes = Option.builder("sg").required(false).longOpt("skip_genotypes")
				.desc("Skip genotypes that are in the GeneSNP pair file but not in the genotype file.").build();
		Option snpsToTestOption = Option.builder("sn").required(true).hasArg().longOpt("snpsToTest").argName("file")
				.desc("Tab delimited file with first column gene name, second column SNP name. Need to match with names from genotype and expression files.").build();
		Option doTestRun = Option.builder("t").required(false).longOpt("test_run")
				.desc("Only run deconvolution for 100 QTLs for quick test run").build();
		Option wholeBloodQTL = Option.builder("w").required(false).longOpt("whole_blood_qtl")
				.desc("Add whole blood eQTL (pearson correlation genotypes and expression)").build();
		options.addOption(filterSamplesOption);
		options.addOption(help);
		options.addOption(outfile);
		options.addOption(expression);
		options.addOption(genotype);
		options.addOption(cellcount);
		options.addOption(roundDosage);
		options.addOption(minimumSamplesPerGenotype);
		options.addOption(allDosages);
		options.addOption(outfolder);
		options.addOption(doTestRun);
		options.addOption(snpsToTestOption);
		options.addOption(skipGenotypes);
		options.addOption(wholeBloodQTL);
		options.addOption(noConsoleOption);
		options.addOption(outputPredictedExpressionOption);
		options.addOption(genotypeConfigurationTypeOption);
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
	
	private void parseOptions(CommandLine cmdLine) throws FileNotFoundException{

		if(cmdLine.hasOption("remove_constraint_violating_samples")){
			removeConstraintViolatingSamples = !removeConstraintViolatingSamples;
		}
		
		if (cmdLine.hasOption("round_dosage")) {
			roundDosage = !roundDosage;
		}

		if (cmdLine.hasOption("all_dosages")){
			allDosages = !allDosages;
		}
		if (cmdLine.hasOption("skip_genotypes")){
			skipGenotypes = !skipGenotypes;
		}
		
		if (cmdLine.hasOption("minimum_samples_per_genotype")) {
			minimumSamplesPerGenotype = Integer.parseInt(cmdLine.getOptionValue("minimum_samples_per_genotype"));
			if(minimumSamplesPerGenotype < 0){
				minimumSamplesPerGenotype = 0;
			}
		}
		
		if(cmdLine.hasOption("genotypeConfigurationType")){
			genotypeConfigurationType = cmdLine.getOptionValue("genotypeConfigurationType");
			
			if(!(genotypeConfigurationType.equals("all") || genotypeConfigurationType.equals("two") || genotypeConfigurationType.equals("one"))){
				throw new IllegalArgumentException("genotypeConfigurationType should be all or two, not "+genotypeConfigurationType);
			}
		}

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
		
		outfolder = cmdLine.getOptionValue("outfolder");
		
		if (cmdLine.hasOption("no_console")) {
			noConsole = !noConsole;
		}
		if (cmdLine.hasOption("filter_samples")){
			filterSamples = !filterSamples;
		}

		if (cmdLine.hasOption("test_run")) {
			testRun = !testRun;
		}
		
		if (cmdLine.hasOption("whole_blood_qtl")){
			wholeBloodQTL = !wholeBloodQTL;
		}
		
		if (cmdLine.hasOption("outputPredictedExpression")){
			outputPredictedExpression = !outputPredictedExpression;
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
	    	DeconvolutionLogger.setup(outfolder, noConsole);
	    	if(dirDidNotExist){
	    		DeconvolutionLogger.log.info("Created directory "+outfolder);
	    	}
	    	DeconvolutionLogger.log.info("Writing output and logfile to "+outfolder);
	    } catch (IOException e) {
	      e.printStackTrace();
	      throw new RuntimeException("Problems with creating the log files");
	    }
	    DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	    Date date = new Date();
	    DeconvolutionLogger.log.info("Starting deconvolution");
	    DeconvolutionLogger.log.info(dateFormat.format(date));
	    DeconvolutionLogger.log.info("Running deconvolution version 1.0.1, src last changed at 06-APR-2017");
	    DeconvolutionLogger.log.info("======= DECONVOLUTION paramater settings =======");
		DeconvolutionLogger.log.info(String.format("Expression file (-e): %s", expressionFile));
		DeconvolutionLogger.log.info(String.format("Genotype file (-g): %s", genotypeFile));
		DeconvolutionLogger.log.info(String.format("Cellcount file (-c): %s", cellCountFile));
		DeconvolutionLogger.log.info(String.format("SNPs to test file (-sn): %s", snpsToTestFile));
		DeconvolutionLogger.log.info(String.format("Outfolder (-o): %s", outfolder));
		DeconvolutionLogger.log.info(String.format("Outfile (-of): %s", outfile));
		DeconvolutionLogger.log.info(String.format("Round dosage (-r): %s", roundDosage));
		DeconvolutionLogger.log.info(String.format("Filter out QTLs where not all dosages are present in at least 1 sample (-ad): %s", allDosages));
		DeconvolutionLogger.log.info(String.format("Minimum samples per genotype (-m): %s", minimumSamplesPerGenotype));
		DeconvolutionLogger.log.info(String.format("Filter samples from output (-f): %s", filterSamples));
		DeconvolutionLogger.log.info(String.format("Remove constraint violating samples (-rc): %s", removeConstraintViolatingSamples));
		DeconvolutionLogger.log.info(String.format("test run doing only 100 QTL (-t): %s", testRun));
		DeconvolutionLogger.log.info(String.format("Skipping genotypes that are in SNP-gene pair file but not in genotype file (-sg): %s", skipGenotypes));
		DeconvolutionLogger.log.info(String.format("Add whole blood eQTL (pearson correlation genotypes and expression) (-w): %s",wholeBloodQTL));
		DeconvolutionLogger.log.info(String.format("Do not ouput logging info to console (-no): %s", noConsole));
		DeconvolutionLogger.log.info(String.format("Write predicted expression to output file (-oe): %s", outputPredictedExpression));
		DeconvolutionLogger.log.info(String.format("Genotype configuration to use (-gc): %s", genotypeConfigurationType));
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

	public int getMinimumSamplesPerGenotype(){
		return minimumSamplesPerGenotype;
	}
	public Boolean getRoundDosage(){
		return roundDosage;
	}

	public Boolean getAllDosages(){
		return allDosages;
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
	public Boolean getFilterSamples(){
		return filterSamples;
	}
	public Boolean getRemoveConstraintViolatingSamples(){
		return removeConstraintViolatingSamples;
	}

	public Boolean getTestRun(){
		return testRun;
	}
	public Boolean getSkipGenotypes(){
		return skipGenotypes;
	}
	public Boolean getWholeBloodQTL(){
		return wholeBloodQTL;
	}

	public Boolean getOutputPredictedExpression(){
		return outputPredictedExpression;
	}

	public String getGenotypeConfigurationType() {
		return genotypeConfigurationType;
	}
	
}





