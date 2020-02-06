/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.conditionalanalysis;

import umcg.genetica.console.ConsoleGUIElems;

/**
 * @author harmjan
 */
public class ConditionalAnalysisConsoleGUI {

	public ConditionalAnalysisConsoleGUI(String[] args) {

		String settingsfile = null;
		String settingstexttoreplace = null;
		String settingstexttoreplacewith = null;
		String in = null;
		String out = null;
		boolean cis = false;
		boolean trans = false;
		int perm = 1;
		String outtype = "text";
		String inexp = null;
		String inexpplatform = null;
		String inexpannot = null;
		String gte = null;
		String snpfile = null;
		Integer threads = null;
		boolean textout = false;
		boolean binout = false;

		boolean iterativeConditional = false;
		boolean skipinitialsnpmap = false;
		boolean skipalleqtlmap = false;
		boolean runSNPCentric = false;
		Integer startiter = 1;

		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			String val = null;

			if (i + 1 < args.length) {
				val = args[i + 1];
			}

			if (arg.equals("--settings")) {
				settingsfile = val;
			} else if (arg.equals("--replacetext")) {
				settingstexttoreplace = val;
			} else if (arg.equals("--replacetextwith")) {
				settingstexttoreplacewith = val;
			} else if (arg.equals("--in")) {
				in = val;
			} else if (arg.equals("--out")) {
				out = val;
			} else if (arg.equals("--text")) {
				textout = true;
			} else if (arg.equals("--binary")) {
				binout = true;
			} else if (arg.equals("--inexp")) {
				inexp = val;
			} else if (arg.equals("--inexpplatform")) {
				inexpplatform = val;
			} else if (arg.equals("--inexpannot")) {
				inexpannot = val;
			} else if (arg.equals("--gte")) {
				gte = val;
			} else if (arg.equals("--cis")) {
				cis = true;
			} else if (arg.equals("--trans")) {
				trans = true;
			} else if (arg.equals("--snps")) {
				snpfile = val;
			} else if (arg.equals("--perm")) {
				try {
					perm = Integer.parseInt(val);
				} catch (NumberFormatException e) {
					System.out.println("Please supply an integer for --perm");
				}
			} else if (arg.equals("--threads")) {
				try {
					threads = Integer.parseInt(val);
				} catch (NumberFormatException e) {
					System.err.println("Error --threads should be an integer");
				}
			} else if (arg.equals("--skip-initial-snpmapping")) {
				skipinitialsnpmap = true;
			} else if (arg.equals("--skip-all-eqtlmapping")) {
				skipalleqtlmap = true;
			} else if (arg.equals("--snpcentric")) {
				runSNPCentric = true;
			} else if (arg.equals("--iterative")) {
				iterativeConditional = true;
			} else if (arg.equals("--startiter")) {
				try {
					startiter = Integer.parseInt(val);
				} catch (NumberFormatException e) {
					System.err.println("Error --threads should be an integer");
				}

			}
		}

		try {
			if (settingsfile == null && in == null) {
				System.out.println("ERROR: Please supply settings file (--settings settings.xml) or --in and --out");
				printUsage();
			} else {
				if (iterativeConditional) {
					IterativeConditionalAnalysis m = new IterativeConditionalAnalysis();
					m.setStartIter(startiter);
					m.run(settingsfile, settingstexttoreplace, settingstexttoreplacewith, in, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads);
				} else {
					ConditionalAnalysis m = new ConditionalAnalysis();
					if (skipinitialsnpmap) {
						m.setSkipInitialSNPMapping();
					}
					if (skipalleqtlmap) {
						m.setSkipAllEQTLMapping();
					}
					if (!binout && !textout) {
						textout = true;
					}
					if (runSNPCentric) {
						m.runGetResultsForAllSNPs(settingsfile, settingstexttoreplace, settingstexttoreplacewith, in, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads);
					} else {
						m.run(settingsfile, settingstexttoreplace, settingstexttoreplacewith, in, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads);
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void printUsage() {
		System.out.print("\nMetaQTL3 - Conditional analysis\n" + ConsoleGUIElems.LINE);
		System.out.println("Perform conditional analysis on a set of SNPs");
		System.out.print("Settings file options:\n" + ConsoleGUIElems.LINE);
		System.out.println("--settings\t\tsettings.xml\tLocation of settings file\n"
				+ "--skip-initial-snpmapping\t\tskips initial SNP mapping, but expects SNP-Initial folder in --out. \n"
				+ "--skip-all-eqtlmapping\t\tskips all eQTL mapping (just summarizes and calculates LD between SNPs), but expects SNP-Initial folder in --out. \n"
				+ "--replacetext\t\ttext\t\tText to replace in settings file\n"
				+ "--replacetextwith\ttext\t\tReplace the text in the settings file, defined by --replacetext with the following text (can be empty)\n"
				+ "--iterative\t\tPerform conditional analysis in iterations.\n"
				+ "--startiter\t\tint\t\tStart iterative analysis at this iteration");

		System.out.println("");
	}
}
