/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import eqtlmappingpipeline.metaqtl3.EQTLRegression;
import eqtlmappingpipeline.metaqtl3.MetaQTL3;
import eqtlmappingpipeline.metaqtl3.containers.Settings;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.ConfigurationException;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 * @author harmjan
 */
public class RegressCisEffectsFromGeneExpressionData extends MetaQTL3 {
	
	public RegressCisEffectsFromGeneExpressionData(String[] args) {
		
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
		String eqtleffectstoregressout = null;
		Double maf = 0.05;
		Double hwe = 0.0001;
		
		Integer nrEQTLsToOutput = null;
		
		
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
			} else if (arg.equals("--regressouteqtls")) {
				eqtleffectstoregressout = val;
			} else if (arg.equals("--perm")) {
				try {
					perm = Integer.parseInt(val);
				} catch (NumberFormatException e) {
					System.out.println("Please supply an integer for --perm");
				}
			} else if (arg.equals(
					"--maf")) {
				try {
					maf = Double.parseDouble(val);
				} catch (NumberFormatException e) {
					System.out.println("Please supply an integer for --maf");
				}
			} else if (arg.equals(
					"--hwe")) {
				try {
					hwe = Double.parseDouble(val);
				} catch (NumberFormatException e) {
					System.out.println("Please supply an integer for --hwe");
				}
			} else if (arg.equals("--threads")) {
				try {
					threads = Integer.parseInt(val);
				} catch (NumberFormatException e) {
					System.err.println("Error --threads should be an integer");
				}
				
			} else if (arg.equals("--maxresults")) {
				try {
					nrEQTLsToOutput = Integer.parseInt(val);
				} catch (NumberFormatException e) {
					System.err.println("Error --maxresults should be an integer");
				}
				
			}
		}
		
		try {
			if (settingsfile == null && in == null) {
				System.out.println("ERROR: Please supply settings file (--settings settings.xml) or --in and --out");
				printUsage();
			} else {
				
				if (!binout && !textout) {
					textout = true;
				}
				this.initialize(settingsfile, settingstexttoreplace, settingstexttoreplacewith, null, null, in, inexp, inexpplatform, inexpannot, gte, out, cis, trans, perm, textout, binout, snpfile, threads, nrEQTLsToOutput, eqtleffectstoregressout, null, false, false, null, maf, hwe);
				
				// now save all the expressiondata to a new file..
				for (int d = 0; d < m_gg.length; d++) {
					TriTyperExpressionData data = m_gg[d].getExpressionData();
					DoubleMatrixDataset<String, String> dmd = new DoubleMatrixDataset<String, String>();
					dmd.rowObjects = Arrays.asList(data.getProbes());
					dmd.colObjects = Arrays.asList(data.getIndividuals());
					
					dmd.rawData = data.getMatrix();
					String outputfile = m_settings.datasetSettings.get(d).expressionLocation + "-CisEffectsRegressedOut.txt";
					System.out.println("Saving to: " + outputfile);
					dmd.save(outputfile);
				}
				System.exit(0);
				
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public RegressCisEffectsFromGeneExpressionData(String settingsFile, String eQTLFile, String textToReplace, String textToReplaceWith) {
		m_settings = new Settings();
		m_settings.settingsTextToReplace = textToReplace;
		m_settings.settingsTextReplaceWith = textToReplaceWith;
		
		try {
			m_settings.load(settingsFile);
		} catch (IOException ex) {
			Logger.getLogger(RegressCisEffectsFromGeneExpressionData.class.getName()).log(Level.SEVERE, null, ex);
		} catch (ConfigurationException ex) {
			Logger.getLogger(RegressCisEffectsFromGeneExpressionData.class.getName()).log(Level.SEVERE, null, ex);
		}
		
		int numDatasets = m_settings.datasetSettings.size();
		m_gg = new TriTyperGeneticalGenomicsDataset[numDatasets];
		numAvailableInds = 0;
		int nrOfDatasetsWithGeneExpressionData = 0;
		
		for (int i = 0; i < numDatasets; i++) {
			System.out.println("- Loading dataset: " + m_settings.datasetSettings.get(i).name + "");
			m_settings.datasetSettings.get(i).confineProbesToProbesMappingToAnyChromosome = m_settings.confineProbesToProbesMappingToAnyChromosome;
			System.out.println(ConsoleGUIElems.LINE);
			try {
				m_gg[i] = new TriTyperGeneticalGenomicsDataset(m_settings.datasetSettings.get(i), null);
			} catch (IOException ex) {
				Logger.getLogger(RegressCisEffectsFromGeneExpressionData.class.getName()).log(Level.SEVERE, null, ex);
			} catch (Exception ex) {
				Logger.getLogger(RegressCisEffectsFromGeneExpressionData.class.getName()).log(Level.SEVERE, null, ex);
			}
			
			if (m_gg[i].isExpressionDataLoadedCorrectly()) {
				nrOfDatasetsWithGeneExpressionData++;
			}
			
		}
		
		if (nrOfDatasetsWithGeneExpressionData == 0 || nrOfDatasetsWithGeneExpressionData != m_gg.length) {
			System.out.println("Error: Something is wrong with the supplied expression data sets,\n please check if all your datasets contain any gene expression data for the settings you have specified");
			System.exit(0);
		}
		
		for (int i = 0; i < numDatasets; i++) {
			m_gg[i].getExpressionData().calcAndSubtractMean();
			m_gg[i].getExpressionData().calcMeanAndVariance();
			numAvailableInds += m_gg[i].getExpressionToGenotypeIdArray().length;
		}
		
		EQTLRegression eqr = new EQTLRegression();
		try {
			eqr.regressOutEQTLEffects(eQTLFile, true, m_gg);
		} catch (IOException ex) {
			Logger.getLogger(RegressCisEffectsFromGeneExpressionData.class.getName()).log(Level.SEVERE, null, ex);
		}
	}
	
	private void printUsage() {
		System.out.print("\nRegressCisEffectsFromGeneExpressionData\n" + ConsoleGUIElems.LINE);
		System.out.println("RegressCisEffectsFromGeneExpressionData removes cis effects from your expression data.");
		System.out.print("\nExamples\n" + ConsoleGUIElems.LINE);
		System.out.println("Example using settingsfile:\tjava -jar eQTLMappingPipeline.jar --mode metaqtl --settings settings.xml");
		System.out.println("Example using commandline:\tjava -jar eQTLMappingPipeline.jar --mode metaqtl --in /path/to/GenotypeMatrix.dat --out /path/to/output/ --cis --perm 10 --text --inexp /path/to/expressiondata.txt --inexpannot /path/to/annotation.txt --gte /path/to/genotypetoexpressioncoupling.txt");
		System.out.println("");
		System.out.print("Settings file options:\n" + ConsoleGUIElems.LINE);
		System.out.println("--settings\t\tsettings.xml\tLocation of settings file\n"
				+ "--replacetext\t\ttext\t\tText to replace in settings file\n"
				+ "--replacetextwith\ttext\t\tReplace the text in the settings file, defined by --replacetext with the following text (can be empty)");
		
		System.out.println("");
		System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
		System.out.println("--in\t\t\tdir\t\tLocation of the genotype data\n"
				+ "--out\t\t\tdir\t\tLocation where the output should be stored\n"
				+ "--cis\t\t\t\t\tPerform cis-eQTL analysis\n"
				+ "--trans\t\t\t\t\tPerform trans-eQTL analysis\n"
				+ "--perm\t\t\tint\t\tNumber of permutations to perform\n"
				+ "--text\t\t\t\t\tOutput results in text format\n"
				+ "--binary\t\t\t\tOutput results in binary format\n"
				+ "--inexp\t\t\tstring\t\tLocation of expression data\n"
				+ "--inexpplatform\t\tstring\t\tGene expression platform\n"
				+ "--inexpannot\t\tstring\t\tLocation of annotation file for gene expression data\n"
				+ "--gte\t\t\tstring\t\tLocation of genotype to expression coupling file\n"
				+ "--snps\t\t\tstring\t\tLocation of file containing SNPs to confine to\n"
				+ "--threads\t\tinteger\t\tNumber of threads to calculate with. Default is number of processors.\n"
				+ "--regressouteqtls\t\tstring\t\tRegress out these eQTL effects before starting the analysis.");
		System.out.println("");
	}
}
