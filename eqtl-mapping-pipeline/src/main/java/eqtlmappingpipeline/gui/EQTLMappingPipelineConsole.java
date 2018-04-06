/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.gui;

import eqtlmappingpipeline.Main;
import eqtlmappingpipeline.ase.AnnotateAseWithSnpEffVcf;
import eqtlmappingpipeline.ase.Ase;
import eqtlmappingpipeline.ase.AseRecalculate;
import eqtlmappingpipeline.binaryInteraction.BinaryInteractionMetaAnalysis;
import eqtlmappingpipeline.binaryInteraction.ConvertTextOutputToBinary;
import eqtlmappingpipeline.binaryInteraction.CovariateImportance;
import eqtlmappingpipeline.binaryInteraction.InvestigateCovariate;
import eqtlmappingpipeline.binaryInteraction.QueryBinaryInteraction;
import eqtlmappingpipeline.binaryInteraction.ReplicateInteractions;
import eqtlmappingpipeline.interactionanalysis.InteractionAnalysisConsoleGUI;
import eqtlmappingpipeline.chromosomeyexpressionplotter.ChrYExpressionPlotConsoleGUI;
import eqtlmappingpipeline.conditionalanalysis.ConditionalAnalysisConsoleGUI;
import eqtlmappingpipeline.eQTLFoldChangeCalculator.QTLFoldChangeCalculatorGUI;
import eqtlmappingpipeline.causalinference.IVConsoleGUI;
import eqtlmappingpipeline.metaqtl3.MetaQTL3ConsoleGUI;
import eqtlmappingpipeline.metaqtl4.MetaQTL4ConsoleUI;
import eqtlmappingpipeline.mixupmapper.MixupMapperConsoleGUI;
import eqtlmappingpipeline.normalization.NormalizationConsoleGUI;
import eqtlmappingpipeline.pcaoptimum.PCAOptimumConsoleGUI;
import eqtlmappingpipeline.qcpca.QCPCAConsoleGui;
import eqtlmappingpipeline.util.ModuleEqtWestraReplication;
import eqtlmappingpipeline.util.ModuleEqtlGeuvadisReplication;
import eqtlmappingpipeline.util.ModuleEqtlNeutrophilReplication;
import eqtlmappingpipeline.util.UtilConsoleGUI;
import eqtlmappingpipeline.util.QTLFileCompare;
import eqtlmappingpipeline.util.QTLmeQTLCompare;

import java.util.Arrays;

import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.pileup.PileupToVcf;

/**
 * @author harmjan
 */
public class EQTLMappingPipelineConsole {
	
	public void main(String[] args) throws Exception {
		
		if (args == null || args.length == 0) {
			
			printHeader();
			printUsage();
			System.out.println("\nERROR: Please supply --mode\n");
			return;
		}
		
		String mode = null;
		
		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			String val = null;
			
			if (i + 1 < args.length) {
				val = args[i + 1];
			}
			
			if (arg.equals("--imputationtool")) {
				
				imputationtool.ImputationTool.main(args);
			}
			
			if (arg.equals("--metamode")) {
				eqtlmappingpipeline.binarymeta.Main.main(args);
				System.exit(0);
			}
			
			if (arg.equals("--utilmode")) {
				UtilConsoleGUI g = new UtilConsoleGUI(args);
				System.exit(0);
			}
			
			if (arg.equals("--mode")) {
				mode = val;
			}
		}
		
		if (mode == null) {
			System.out.println("ERROR: Please supply --mode");
			printUsage();
		} else {
			if (mode.equals("metaqtl")) {
				MetaQTL3ConsoleGUI metaQTL = new MetaQTL3ConsoleGUI(args);
			} else if (mode.equals("metaqtl4")) {
				MetaQTL4ConsoleUI mui = new MetaQTL4ConsoleUI(args);
			} else if (mode.equals("mixupmapper")) {
				MixupMapperConsoleGUI m = new MixupMapperConsoleGUI(args);
			} else if (mode.equals("normalize")) {
				NormalizationConsoleGUI p = new NormalizationConsoleGUI(args);
			} else if (mode.equals("compare")) {
				QTLFileCompare r = new QTLFileCompare(args);
			} else if (mode.equals("compareM")) {
				QTLmeQTLCompare r = new QTLmeQTLCompare(args);
			} else if (mode.equals("chryplot")) {
				ChrYExpressionPlotConsoleGUI r = new ChrYExpressionPlotConsoleGUI(args);
			} else if (mode.equals("pcaoptimum")) {
				PCAOptimumConsoleGUI g = new PCAOptimumConsoleGUI(args);
			} else if (mode.equals("foldchange")) {
				QTLFoldChangeCalculatorGUI g = new QTLFoldChangeCalculatorGUI(args);
			} else if (mode.equals("util")) {
				UtilConsoleGUI g = new UtilConsoleGUI(args);
			} else if (mode.equals("qcpca")) {
				QCPCAConsoleGui q = new QCPCAConsoleGui(args);
			} else if (mode.equals("conditional")) {
				ConditionalAnalysisConsoleGUI q = new ConditionalAnalysisConsoleGUI(args);
			} else if (mode.equals("iv") || mode.equals("causalinference")) {
				IVConsoleGUI q = new IVConsoleGUI(args);
			} else if (mode.equals("celltypespecific") || mode.equals("interaction")) {
				InteractionAnalysisConsoleGUI q = new InteractionAnalysisConsoleGUI(args);
			} else if (mode.equals("ase")) {
				Ase.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("aseR")) {
				AseRecalculate.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("itb")) {
				ConvertTextOutputToBinary.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("queryInteraction") || mode.equals("qi")) {
				QueryBinaryInteraction.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("binaryMeta") || mode.equals("bm")) {
				BinaryInteractionMetaAnalysis.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("replicateInteraction") || mode.equals("ri")) {
				ReplicateInteractions.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("interactionChi2") || mode.equals("ic")) {
				CovariateImportance.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("covariate")) {
				InvestigateCovariate.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("pileupToVcf")) {
				PileupToVcf.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("moduleReplication")) {
				ModuleEqtlGeuvadisReplication.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("moduleNeutrophilReplication")) {
				ModuleEqtlNeutrophilReplication.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("moduleestralReplication")) {
				ModuleEqtWestraReplication.main(Arrays.copyOfRange(args, 2, args.length));
				return;
			} else if (mode.equals("aseSnpEff")) {
				AnnotateAseWithSnpEffVcf.annotateAseWithSnpEffVcf(args[2], args[3], args[4]);
				return;
			} else {
				printUsage();
			}
		}
		
		System.out.println(ConsoleGUIElems.DOUBLELINE);
	}
	
	private void printHeader() {
		
		//Note: Version is null when running from netbeans but will be set when buidling a jar
		System.out.println("\n"
				+ ConsoleGUIElems.DOUBLELINE
				+ "Version: " + Main.VERSION + "\n"
				+ "Department of Genetics, University Medical Center Groningen\n"
				+ "www.molgenis.org/systemsgenetics\n"
				+ "Harm-Jan Westra, Patrick Deelen, Marc Jan Bonder, Dasha Zhernakova and Lude Franke\n"
				+ ConsoleGUIElems.DOUBLELINE + "\n");
	}
	
	/*
	 * prints usage
	 */
	private void printUsage() {
		System.out.println("");
		System.out.print("\nCommand line options\n" + ConsoleGUIElems.LINE);
		System.out.println("       metaqtl\t\tPerform QTL mapping on the dataset\n"
				+ "       mixupmapper\tPerform MixUp mapping on the dataset\n"
				+ "       normalize\tPerform and Remove principal components from expression data\n"
				+ "       pcaoptimum\tDetermine optimum PCs to subtract from the data\n"
				+ "       conditional\tPerform conditional eQTL analysis given a set of SNPs\n"
				+ "       causalinference\t\tPerform instrumental variable analysis\n"
				+ "       celltypespecific\tCell type specific eQTL mapping\n"
				+ "       interaction\teQTL mapping using interaction model\n"
				+ "       util\tSmall utilities\n"
				+ "       ase\t\tAllele Specific Expression mapping\n"
				+ "       binaryMeta\t\tDo interaction meta analysis on binary files\n"
				+ "       queryInteraction\t\tQuery binary interaction file\n"
				+ "       pileupToVcf\tConvert a pileup file to vcf for usage in ASE mapping\n"
				+ "       compare\tCompare QTL files\n"
				+ "       compareM\tCompare eQTL and eQTM files\n");
		System.out.println("");
		System.out.println("More information: www.molgenis.org/systemsgenetics/QTL-mapping-pipeline");
		
	}
}
