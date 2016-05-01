package deconvolution;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class CommandLineOptions {
	private String expressionFile;
	private String genotypeFile;
	private String cellcountFile;
	private String outfile = "deconvolutionResults.csv";
	private String outfolder;
	private int numberOfPermutations = 0;
	private String permutationType = "genotype";
	private int numberOfThreads = 1;
	private Boolean forceNormalExpression = false;
	private Boolean forceNormalCellcount = false;
	private int minimumSamplesPerGenotype = 0;
	private Boolean roundDosage = false;
	private Boolean plot_beta_times_variables = false;
	private Boolean allDosages = false;
	private Boolean plotBetas = false;
	private String normalizationType = "normalizeAddMean";
	public void parseCommandLine(String[] args) throws ParseException {
		/*
		 * Standard command line parsing.
		 * 
		 * @param args A string vector of all arguments given to the command
		 * line, e.g. for `java -jar Deconvolution.jar -o` args = ["-o"]
		 * 
		 * @return A CommandLine object that includes all the options given to
		 * the command line
		 */
		Options options = new Options();
		Option help = new Option("help", "print this message");
		Option roundDosage = Option.builder("r").required(false).longOpt("round_dosage")
				.desc("Round the dosage to the closest int").build();
		Option permute = Option.builder("p").required(false).longOpt("permute").hasArg()
				.desc("Do permutations. Uses more than usual memory.").build();
		Option numberOfThreads = Option.builder("t").required(false).longOpt("threads").hasArg()
				.desc("Number of threads to use.").build();
		Option outfile = Option.builder("of").required(false).hasArg().longOpt("outfile").desc("Outfile name of deconvolution results (will be written in outfolder)")
				.argName("file").build();
		Option outfolder = Option.builder("o").required(true).hasArg().longOpt("outfolder").desc("Path to folder to write output to")
				.argName("path").build();
		Option expression = Option.builder("e").required(true).hasArg().longOpt("expression")
				.desc("Expression file name").argName("file").build();
		Option genotype = Option.builder("g").required(true).hasArg().longOpt("genotype").desc("Genotype file name")
				.argName("file").build();
		Option cellcount = Option.builder("c").required(true).hasArg().longOpt("cellcount").desc("Cellcount file name")
				.argName("file").build();
		Option minimum_samples_per_genotype = Option.builder("m").required(false).hasArg().longOpt("minimum_samples_per_genotype")
				.desc("The minimum amount of samples need for each genotype of a QTL for the QTL to be included in the results")
				.build();
		Option plotBetaTimesVariables = Option.builder("b").required(false).hasArg().longOpt("plot_betas")
				.desc("Plot the B1*X1, B2*X2 etc values if the sum of Bx*CELLTYPEz+By*CELLTYPEz:GT < 0").build();
		Option forceNormalExpression = Option.builder("ne").required(false).hasArg().longOpt("force_normal_expression")
				.desc("Force normal on the expression data").build();
		Option normalizationType = Option.builder("n").required(false).hasArg().longOpt("normalization_type")
				.desc("Type to normalization to use when normalizing expression data (Default: normalizeAddMean)").build();
		Option forceNormalCellcount = Option.builder("nc").required(false).hasArg().longOpt("force_normal_cellcount")
				.desc("Force normal on the expression data").build();
		Option allDosages = Option.builder("ad").required(false).hasArg().longOpt("all_dosages")
				.desc("Filter out QTLs where not all dosages are present in at least 1 sample").build();
		Option permuteType = Option.builder("pt").required(false).hasArg().longOpt("permutation_type")
				.desc("Type to permute on, either genotype or expression (Default: genotype)").build();
		
		options.addOption(normalizationType);
		options.addOption(plotBetaTimesVariables);
		options.addOption(help);
		options.addOption(permute);
		options.addOption(numberOfThreads);
		options.addOption(outfile);
		options.addOption(expression);
		options.addOption(genotype);
		options.addOption(cellcount);
		options.addOption(roundDosage);
		options.addOption(minimum_samples_per_genotype);
		options.addOption(forceNormalExpression);
		options.addOption(forceNormalCellcount);
		options.addOption(permuteType);
		options.addOption(allDosages);
		options.addOption(outfolder);
		CommandLineParser cmdLineParser = new DefaultParser();
		CommandLine cmdLine = cmdLineParser.parse(options, args);
		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		if (cmdLine.hasOption("help")) {
			formatter.printHelp("deconvolution", options, true);
		}
		parseOptions (cmdLine);
		printArgumentValues(cmdLine);
	}
	
	private void parseOptions(CommandLine cmdLine){
		permutationType = "genotype";
		if(cmdLine.hasOption("permutation_type")){
			permutationType = cmdLine.getOptionValue("permutation_type");
			if(!(permutationType.equals("genotye") || permutationType.equals("expression"))){
				throw new IllegalArgumentException("permutation_type should be genotype or expression, not "+cmdLine.getOptionValue("permutation_type"));
			}
		}
		if (cmdLine.hasOption("permute")) {
			numberOfPermutations = Integer.parseInt(cmdLine.getOptionValue("permute"));
			if(numberOfPermutations < 1){
				numberOfPermutations = 1;
			}
		}
		
		if (cmdLine.hasOption("threads")) {
			numberOfThreads = Integer.parseInt(cmdLine.getOptionValue("threads"));
		}

		if (cmdLine.hasOption("round_dosage")) {
			roundDosage = true;
		}
		if (cmdLine.hasOption("force_normal_expression")){
			forceNormalExpression = true;
		}

		if (cmdLine.hasOption("force_normal_cellcount")){
			forceNormalCellcount = true;
		}

		if (cmdLine.hasOption("all_dosages")){
			allDosages = true;
		}

		if (cmdLine.hasOption("minimum_samples_per_genotype")) {
			minimumSamplesPerGenotype = Integer.parseInt(cmdLine.getOptionValue("minimum_samples_per_genotype"));
			if(minimumSamplesPerGenotype < 0){
				minimumSamplesPerGenotype = 0;
			}
		}
		if (cmdLine.hasOption("plot_betas")) {
			plotBetas = true;
		}

		expressionFile = cmdLine.getOptionValue("expression");
		genotypeFile = cmdLine.getOptionValue("genotype");
		cellcountFile = cmdLine.getOptionValue("cellcount");
		if (cmdLine.hasOption("outfile")) {
			outfile = cmdLine.getOptionValue("outfile");
		}
		outfolder = cmdLine.getOptionValue("outfolder");
	}
	

	private void printArgumentValues(CommandLine cmdLine){
		System.out.println("======= DECONVOLUTION paramater settings =======");
		System.out.printf("Expression file (-e): %s\n", expressionFile);
		System.out.printf("Genotype file (-g): %s\n", genotypeFile);
		System.out.printf("Cellcount file (-c): %s\n", cellcountFile);
		System.out.printf("Outfolder (-o): %s\n", outfolder);
		System.out.printf("Outfile (-of): %s\n", outfile);
		System.out.printf("Number of permutations (-p): %s\n", numberOfPermutations);
		if(numberOfPermutations > 0){
			System.out.printf("Permutation type (-pt): %s\n", permutationType);
		}
		System.out.printf("Threads (-t): %d\n", numberOfThreads);
		System.out.printf("Plot beta (-b): %s\n", plot_beta_times_variables);
		System.out.printf("Plot beta x variables (-b): %s\n", plotBetas);
		System.out.printf("Normalize expression (-ne): %s\n", forceNormalExpression);
		System.out.printf("Normalization type (-nt): %s\n", normalizationType);
		System.out.printf("Normalize cellcounts (-ne): %s\n", forceNormalCellcount);
		System.out.printf("Round dosage (-r): %s\n", roundDosage);
		System.out.printf("All dosages (-ad): %s\n", allDosages);
		System.out.printf("Minimum samples per genotype (-m): %s\n", minimumSamplesPerGenotype);
		
		System.out.println("=================================================");
	}
	public String getExpressionFile(){
		return (expressionFile);
	}
	public String getGenotypeFile(){
		return(genotypeFile);
	}
	public String getCellcountFile(){
		return(cellcountFile);
	}
	public String getOutfile(){
		return(outfile);
	}
	public int getNumberOfPermutations(){
		return(numberOfPermutations);
	}
	public String getPermutationType(){
		return(permutationType);
	}
	public int getNumberOfThreads(){
		return(numberOfThreads);
	}
	public Boolean getForceNormalExpression(){
		return(forceNormalExpression);
	}
	public Boolean getForceNormalCellcount(){
		return(forceNormalCellcount);
	}
	public int getMinimumSamplesPerGenotype(){
		return(minimumSamplesPerGenotype);
	}
	public Boolean getRoundDosage(){
		return(roundDosage);
	}
	public Boolean getPlotBetas(){
		return(plotBetas);
	}
	public Boolean getAllDosages(){
		return(allDosages);
	}
	public String getOutfolder(){
		return(outfolder);
	}
	public String getNormalizationType(){
		return(normalizationType);
	}
	public void setOutfolder(String newOutfolder){
		outfolder = newOutfolder;
	}
}





