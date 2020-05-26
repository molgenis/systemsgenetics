/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.mixupmapper;

import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.containers.Settings;
import org.apache.commons.configuration.ConfigurationException;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;

import java.io.IOException;

/**
 * @author harmjan
 */
public class MixupMapperConsoleGUI {

	public MixupMapperConsoleGUI(String[] args) {


		String xmlSettingsFile = null;
		String settingstexttoreplace = null;
		String settingstexttoreplacewith = null;
		String in = null;
		String out = null;
		boolean cis = true;
		boolean trans = false;
		int perm = 10;
		String outtype = "text";
		String inexp = null;
		String inexpplatform = null;
		String inexpannot = null;
		String gte = null;
		Integer threads = null;
		String inputeQTLs = null;
		String snps = null;
		boolean allCombos = false;

		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			String val = null;

			if (i + 1 < args.length) {
				val = args[i + 1];
			}

			if (arg.equals("--settings")) {
				xmlSettingsFile = val;
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
			} else if (arg.equals("--snps")) {
				snps = val;
			} else if (arg.equals("--gte")) {
				gte = val;
			} else if (arg.equals("--eqtls")) {
				inputeQTLs = val;
			} else if (arg.equals("--testall")) {
				allCombos = true;
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

			}
		}

		try {
			if ((in == null || out == null) && xmlSettingsFile == null) {
				System.out.println("ERROR: Please supply settings --in and --out or --settings settings.xml");
				printUsage();
			} else {


				if (xmlSettingsFile != null) {

					eqtlmappingpipeline.metaqtl3.containers.Settings settings = new Settings();
					try {
						settings.load(xmlSettingsFile);
						for (TriTyperGeneticalGenomicsDatasetSettings dataset : settings.datasetSettings) {
							MixupMapper m = new MixupMapper();
							String outdir = settings.outputReportsDir + "/" + dataset.name + "/";
							settings.fdrType = FDR.FDRMethod.FULL;
							System.out.println(settings.numberOfVariantsToBuffer);
							m.run(null, null, null, dataset.genotypeLocation, dataset.expressionLocation, dataset.expressionplatform,
									dataset.probeannotation, dataset.genotypeToExpressionCoupling, outdir, true, false, 10, true, false, null,
									settings.nrThreads, settings.maxNrMostSignificantEQTLs, null, null, inputeQTLs, allCombos, settings.numberOfVariantsToBuffer);
							System.gc();
							System.gc();
						}

					} catch (IOException e) {
						e.printStackTrace();
					} catch (ConfigurationException e) {
						e.printStackTrace();
					} catch (Exception e) {
						e.printStackTrace();
					}


				} else {
					MixupMapper m = new MixupMapper();
					m.run(xmlSettingsFile, settingstexttoreplace, settingstexttoreplacewith, in, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, true, false, snps, threads, 500000, null, null, inputeQTLs, allCombos);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void printUsage() {
		System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
		System.out.println("--in\t\t\tdir\t\tLocation of the genotype data\n"
				+ "--out\t\t\tdir\t\tLocation where the output should be stored\n"
				+ "--inexp\t\t\tstring\t\tLocation of expression data\n"
				+ "--inexpplatform\t\tstring\t\tGene expression platform\n"
				+ "--inexpannot\t\tstring\t\tLocation of annotation file for gene expression data\n"
				+ "--gte\t\t\tstring\t\tLocation of genotype to expression coupling file\n"
				+ "--perm\t\t\tstring\t\tNumber of permutations to perform in order to determine FDR\n"
				+ "--eqtls\t\t\tstring\t\tPath to eQTL file to use for MixupMapper\n"
				+ "--testall\t\t\t\tTest all possible combinations of genotype and gene expression samples\n"
				+ "--threads\t\tinteger\t\tNumber of threads to calculate with. Default is number of processors.\n"
				+ "--snps\t\tstring\t\tList of SNPs to test.\n"
				+ "--settings\t\t\tfile\t\txml settings file\n");
		System.out.println("");
	}
}
