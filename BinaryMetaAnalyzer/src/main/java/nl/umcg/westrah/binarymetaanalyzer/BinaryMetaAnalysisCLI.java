package nl.umcg.westrah.binarymetaanalyzer;

import org.apache.commons.cli.*;

import java.io.IOException;

public class BinaryMetaAnalysisCLI {
	
	private static final Options OPTIONS;
	
	static {
		OPTIONS = new Options();

//		Option option = Option.builder().longOpt("gwas").build();
//		OPTIONS.addOption(option);

//		option = Option.builder()
//				.longOpt("meta")
//				.hasArg()
//				.desc("Run second iteration conditional on a set of predefined variants. Specify file in format: region iter variant pval")
//				.build();
//		OPTIONS.addOption(option);
		
		
		Option option = Option.builder()
				.longOpt("meta")
				
				.desc("Run meta-analysis")
				.build();
		OPTIONS.addOption(option);
		option = Option.builder()
				.longOpt("qc")
				
				.desc("Run QC")
				.build();
		
		OPTIONS.addOption(option);
		
		
		option = Option.builder("s")
				.longOpt("settings")
				.hasArg()
				.desc("XML settings file")
				.build();
		OPTIONS.addOption(option);
		option = Option.builder("rt")
				.longOpt("replacetext")
				.hasArg()
				.desc("Replace a pattern")
				.build();
		OPTIONS.addOption(option);
		option = Option.builder("rtw")
				.longOpt("replacetextwith")
				.hasArg()
				.desc("Replace the pattern with this text")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("e")
				.longOpt("eqtls")
				.hasArg()
				.desc("eQTL Files (comma separated)")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("g")
				.longOpt("groups")
				.hasArg()
				.desc("File with group definition, tab-separated, format: cohort\tgroup")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("o")
				.longOpt("out")
				.hasArg()
				.desc("Output file (for QC)")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("a")
				.longOpt("annot")
				.hasArg()
				.desc("SNP annotation file (SNP ID should be on first column, tab separated)")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("t")
				.longOpt("threshold")
				.hasArg()
				.desc("Threshold for p-val")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder()
				.longOpt("nonan")
				.desc("Force presence of QTL in all datasets")
				.build();
		OPTIONS.addOption(option);
		
		
	}
	
	public static void main(String[] args) {
		BinaryMetaAnalysisCLI cli = new BinaryMetaAnalysisCLI(args);
	}
	
	public BinaryMetaAnalysisCLI(String[] args) {
		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);
			
			if (cmd.hasOption("meta")) {
				String settings = null;
				String texttoreplace = null;
				String replacewith = null;
				if (cmd.hasOption("settings")) {
					settings = cmd.getOptionValue("settings");
				}
				if (cmd.hasOption("rt")) {
					texttoreplace = cmd.getOptionValue("rt");
				}
				if (cmd.hasOption("rtw")) {
					replacewith = cmd.getOptionValue("rtw");
				}
				if (settings == null) {
					System.err.println("Error: Please provide settings with --meta");
					System.out.println();
					printHelp();
				} else {
					BinaryMetaAnalysis bm = new BinaryMetaAnalysis(settings, texttoreplace, replacewith);
				}
			} else if (cmd.hasOption("qc")) {
				MetaAnalysisQC q = new MetaAnalysisQC();
				String eqtlfile = null;
				String groupdefinition = null;
				String outfile = null;
				String snpannotation = null;
				
				if (cmd.hasOption("e")) {
					eqtlfile = cmd.getOptionValue("e");
				}
				if (cmd.hasOption("o")) {
					outfile = cmd.getOptionValue("o");
				}
				if (cmd.hasOption("g")) {
					groupdefinition = cmd.getOptionValue("g");
				}
				
				double threshold = 0.05;
				if (cmd.hasOption("threshold")) {
					threshold = Double.parseDouble(cmd.getOptionValue("threshold"));
				}
				
				boolean nonan = false;
				if (cmd.hasOption("nonan")) {
					nonan = true;
				}
				
				if (cmd.hasOption("annot")) {
					snpannotation = cmd.getOptionValue("annot");
				}
				
				if (eqtlfile == null || groupdefinition == null || outfile == null) {
					System.err.println("Specify -e, -g and -o with --qc");
					System.out.println();
					printHelp();
				} else {
					q.ComparePBMCWithWholeBloodCohorts(eqtlfile,
							groupdefinition,
							outfile,
							snpannotation,
							threshold,
							nonan);
				}
				
			} else {
				System.out.println("Choose --meta or --qc");
				printHelp();
			}
			
		} catch (ParseException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		
		System.exit(-1);
	}
}
