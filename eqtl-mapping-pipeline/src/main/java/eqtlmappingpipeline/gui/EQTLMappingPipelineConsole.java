/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.gui;

import eqtlmappingpipeline.Main;
import eqtlmappingpipeline.celltypespecific.CelltypeSpecificeQTLMappingConsoleGUI;
import eqtlmappingpipeline.chromosomeyexpressionplotter.ChrYExpressionPlotConsoleGUI;
import eqtlmappingpipeline.conditionalanalysis.ConditionalAnalysisConsoleGUI;
import eqtlmappingpipeline.eQTLFoldChangeCalculator.eQTLFoldChangeCalculatorGUI;
import eqtlmappingpipeline.causalinference.IVConsoleGUI;
import eqtlmappingpipeline.metaqtl3.MetaQTL3ConsoleGUI;
import eqtlmappingpipeline.mixupmapper.MixupMapperConsoleGUI;
import eqtlmappingpipeline.normalization.NormalizationConsoleGUI;
import eqtlmappingpipeline.pcaoptimum.PCAOptimumConsoleGUI;
import eqtlmappingpipeline.qcpca.QCPCAConsoleGui;
import eqtlmappingpipeline.util.UtilConsoleGUI;
import eqtlmappingpipeline.util.eQTLFileCompare;
import umcg.genetica.console.ConsoleGUIElems;

/**
 *
 * @author harmjan
 *
 */
public class EQTLMappingPipelineConsole {

    public void main(String[] args) {
        printHeader();

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
                //MetaQTL4ConsoleGUI metaQTL = new MetaQTL4ConsoleGUI(args);
            } else if (mode.equals("mixupmapper")) {
                MixupMapperConsoleGUI m = new MixupMapperConsoleGUI(args);
            } else if (mode.equals("normalize")) {
                NormalizationConsoleGUI p = new NormalizationConsoleGUI(args);
            } else if (mode.equals("compare")) {
                eQTLFileCompare r = new eQTLFileCompare(args);
            } else if (mode.equals("chryplot")) {
                ChrYExpressionPlotConsoleGUI r = new ChrYExpressionPlotConsoleGUI(args);
            } else if (mode.equals("pcaoptimum")) {
                PCAOptimumConsoleGUI g = new PCAOptimumConsoleGUI(args);
            } else if (mode.equals("foldchange")) {
                eQTLFoldChangeCalculatorGUI g = new eQTLFoldChangeCalculatorGUI(args);
            } else if (mode.equals("util")) {
                UtilConsoleGUI g = new UtilConsoleGUI(args);
            } else if (mode.equals("qcpca")) {
                QCPCAConsoleGui q = new QCPCAConsoleGui(args);
            } else if (mode.equals("conditional")) {
                ConditionalAnalysisConsoleGUI q = new ConditionalAnalysisConsoleGUI(args);
            } else if (mode.equals("iv")) {
                IVConsoleGUI q = new IVConsoleGUI(args);
            } else if (mode.equals("celltypespecific")) {
                CelltypeSpecificeQTLMappingConsoleGUI q = new CelltypeSpecificeQTLMappingConsoleGUI(args);
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
                + "Department of Genetics, University Medical Centre Groningen\n"
                + "Harm-Jan Westra, Lude Franke, Marjolein Peters, Tonu Esko, Claudia Schurmann\n"
                + ConsoleGUIElems.LINE
                + "Please note that this software comes with no warranty and should not be\n"
                + "distributed to third parties without proper consent of the developers.\n"
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
                + "       qcpca\t\tPerform principal component analysis on both genotype and gene expression data\n"
                + "       conditional\tPerform conditional eQTL analysis given a set of SNPs\n"
                + "       iv\t\tPerform instrumental variable analysis\n"
                + "       celltypespecific\tCell type specific eQTL mapping");
        System.out.println("");

    }
}
