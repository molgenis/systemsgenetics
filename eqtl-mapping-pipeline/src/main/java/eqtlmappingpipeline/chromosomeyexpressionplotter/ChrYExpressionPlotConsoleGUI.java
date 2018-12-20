/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.chromosomeyexpressionplotter;

import eqtlmappingpipeline.metaqtl3.containers.Settings;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;

import java.util.ArrayList;

/**
 * @author harmjan
 */
public class ChrYExpressionPlotConsoleGUI {
	
	public ChrYExpressionPlotConsoleGUI(String[] args) {
		
		
		String xmlSettingsFile = null;
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
		
		
		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			String val = null;
			
			if (i + 1 < args.length) {
				val = args[i + 1];
			}
			
			if (arg.equals("--settings")) {
				xmlSettingsFile = val;
			} else if (arg.equals("--replacetext")) {
				settingstexttoreplace = val;
			} else if (arg.equals("--replacetextwith")) {
				settingstexttoreplacewith = val;
			} else if (arg.equals("--in")) {
				in = val;
			} else if (arg.equals("--out")) {
				out = val;
			} else if (arg.equals("--inexp")) {
				inexp = val;
			} else if (arg.equals("--inexpplatform")) {
				inexpplatform = val;
			} else if (arg.equals("--inexpannot")) {
				inexpannot = val;
			} else if (arg.equals("--gte")) {
				gte = val;
			}
		}
		
		try {
			if (xmlSettingsFile == null && in == null) {
				System.out.println("ERROR: Please supply settings file (--settings settings.xml) or --in and --out");
				printUsage();
			} else {
				
				Settings s = null;
				if (xmlSettingsFile != null) {
					s = new Settings();
					
					if (settingstexttoreplace != null && settingstexttoreplace.length() > 0) {
						s.settingsTextReplaceWith = settingstexttoreplace;
						s.settingsTextToReplace = settingstexttoreplacewith;
					}
					s.load(xmlSettingsFile);
					for (TriTyperGeneticalGenomicsDatasetSettings gs : s.datasetSettings) {
						gs.quantilenormalize = false;
						gs.confineProbesToProbesMappingToAnyChromosome = true;
					}
				} else {
					if (out == null) {
						System.out.println("ERROR: Please supply --out when specifiying parameters from command line");
						printUsage();
					} else {
						// prepare settings object...
						s = new Settings();
						s.expressionLocation = inexp;
						s.expressionplatform = inexpplatform;
						s.probeannotation = inexpannot;
						s.genotypeLocation = in;
						s.genotypeToExpressionCoupling = gte;
						s.datasetSettings = new ArrayList<TriTyperGeneticalGenomicsDatasetSettings>();
						s.outputReportsDir = out;
						TriTyperGeneticalGenomicsDatasetSettings gs = new TriTyperGeneticalGenomicsDatasetSettings();
						gs.expressionLocation = inexp;
						gs.expressionplatform = inexpplatform;
						gs.probeannotation = inexpannot;
						gs.genotypeLocation = in;
						gs.genotypeToExpressionCoupling = gte;
						gs.quantilenormalize = false;
						gs.logtransform = false;
						gs.confineProbesToProbesMappingToAnyChromosome = true;
						gs.name = "Dataset";
						s.datasetSettings.add(gs);
					}
				}
				
				if (s != null) {
					ChromosomeYExpressionPlotter c = new ChromosomeYExpressionPlotter();
					c.run(s);
				} else {
					System.out.println("Error determining settings");
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private void printUsage() {
		
		System.out.print("\nChromosome Y Expression plotter\n" + ConsoleGUIElems.LINE);
		System.out.println("Chromosome Y Expression plotter can be used to get an overview of possible gender mismatches, on the basis of comparing chromosome Y gene expression data to predicted genotypes.");
		System.out.print("\nExamples\n" + ConsoleGUIElems.LINE);
		System.out.println("Example using settingsfile:\tjava -jar eQTLMappingPipeline.jar --mode chryplot --settings settings.xml");
		System.out.println("Example using commandline:\tjava -jar eQTLMappingPipeline.jar --mode chryplot --in /path/to/GenotypeMatrix.dat --out /path/to/output/ --inexp /path/to/expressiondata.txt --inexpannot /path/to/annotation.txt --gte /path/to/genotypetoexpressioncoupling.txt");
		System.out.println("");
		System.out.print("Settings file options:\n" + ConsoleGUIElems.LINE);
		System.out.println("--settings\t\tsettings.xml\tLocation of settings file\n"
				+ "--replacetext\t\ttext\t\tText to replace in settings file\n"
				+ "--replacetextwith\ttext\t\tReplace the text in the settings file, defined by --replacetext with the following text (can be empty)");
		System.out.println("");
		System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
		System.out.println("--in\t\t\tdir\t\tLocation of the genotype data\n"
				+ "--out\t\t\tdir\t\tLocation where the output should be stored\n"
				+ "--inexp\t\t\tstring\t\tLocation of expression data\n"
				+ "--inexpplatform\t\tstring\t\tGene expression platform\n"
				+ "--inexpannot\t\tstring\t\tLocation of annotation file for gene expression data\n"
				+ "--gte\t\t\tstring\t\tLocation of genotype to expression coupling file\n");
		System.out.println("");
	}
	
}
