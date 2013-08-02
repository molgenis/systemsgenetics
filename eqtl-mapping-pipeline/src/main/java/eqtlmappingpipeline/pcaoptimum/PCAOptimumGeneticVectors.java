/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.pcaoptimum;

import eqtlmappingpipeline.util.eQTLFileCompare;
import eqtlmappingpipeline.metaqtl3.MetaQTL3Settings;
import eqtlmappingpipeline.normalization.Normalizer;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Properties;
import java.util.Set;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class PCAOptimumGeneticVectors extends PCAOptimum {
    
    public void run(String xmlSettingsFile, String texttoreplace, String texttoreplacewith,
            String ingt, String inexp, String inexpplatform, String inexpannot, String gte,
            String out, boolean cis, boolean trans, int perm, boolean textout, boolean binout, String snpfile, Integer threads, boolean runonlypcqtl, Integer runpcsremoved) throws IOException, Exception {
        if (!out.endsWith("/")) {
            out += "/";
        }
        if (!Gpio.exists(out)) {
            Gpio.createDir(out);
        }

        permutations = perm;

        String origInExp = inexp;

        m_settings = new MetaQTL3Settings();

        int nrProcs = Runtime.getRuntime().availableProcessors();
        if (threads != null && threads > 0 && threads <= nrProcs) {
            //
        } else {
            if (threads != null && threads > nrProcs) {
                System.out.println("The number of threads you set using the command line is not correct for your system. You set " + threads + " threads, while your machine has " + nrProcs + " processors");
            }
            threads = nrProcs;
        }

        m_threads = threads;

        int round = 0;

        EQTL[] originalCisEQTLs = null;
        EQTL[] originalTransEQTLs = null;

        // set standard cis-settings
        this.inexpannot = inexpannot;
        this.inexpplatform = inexpplatform;
        this.ingt = ingt;
        this.gte = gte;


        Properties properties = new Properties();
        InputStream inputStream = getClass().getResourceAsStream("CisSNPs.properties");
        properties.load(inputStream);
        Set keys = properties.keySet();
        HashSet<String> cisSnpsToTest = new HashSet<String>();
        for (Object k : keys) {
            cisSnpsToTest.add("rs" + (String) k);
        }
        inputStream.close();

        properties = new Properties();
        InputStream inputStream2 = getClass().getResourceAsStream("TransSNPs.properties");
        properties.load(inputStream2);
        keys = properties.keySet();
        HashSet<String> transSnpsToTest = new HashSet<String>();
        for (Object k : keys) {
            transSnpsToTest.add("rs" + (String) k);
        }
        inputStream2.close();

//	if (runonlypcqtl || runpcsremoved == null) {
        performeQTLMappingOverEigenvectorMatrixAndReNormalize(origInExp, out);
//	}

        String outputDir = out + "CisTrans-PCAEigenVectors/";
        eQTLFileCompare e = new eQTLFileCompare();
        round = 1;

        // run cis / trans on 0 pcs removed
//	if ((!runonlypcqtl && runpcsremoved == null) || (runpcsremoved != null && runpcsremoved == 0)) {
        String nextInExp = origInExp + ".QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt";
        boolean zeropcexpfileexists = true;
        if (covariatesremoved) {
            nextInExp = origInExp + ".QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt";
            if (!Gpio.exists(nextInExp)) {
                nextInExp += ".gz";
                if (!Gpio.exists(nextInExp)) {
                    zeropcexpfileexists = false;
                }
            }
        } else {
            if (!Gpio.exists(nextInExp)) {
                nextInExp += ".gz";
                if (!Gpio.exists(nextInExp)) {
                    zeropcexpfileexists = false;
                }
            }
        }

        Gpio.createDir(out + "Cis-0PCAsRemoved/");
        Gpio.createDir(out + "Trans-0PCAsRemoved/");

        if (Gpio.isDir(out + "Cis-0PCAsRemoved/")) {
            if (!zeropcexpfileexists) {
                System.out.println("Could not find expression file: " + nextInExp);
                System.exit(0);
            }
            outputDir = out + "Cis-0PCAsRemoved/";
            performeQTLMapping(true, false, nextInExp, outputDir, cisSnpsToTest, null, threads);
            cleanup();
        }

        if (Gpio.isDir(out + "Trans-0PCAsRemoved/")) {
            if (!zeropcexpfileexists) {
                System.out.println("Could not find expression file: " + nextInExp);
                System.exit(0);
            }
            outputDir = out + "Trans-0PCAsRemoved/";
            performeQTLMapping(false, true, nextInExp, outputDir, transSnpsToTest, null, threads);
            cleanup();
        }
//	}

        // run the remaining cis/trans eQTL steps

//	if ((!runonlypcqtl && runpcsremoved == null)) {
        for (int pca = 5; pca < 101; pca += 5) {
//		if (runpcsremoved == null || runpcsremoved == pca) {
            nextInExp = origInExp + "." + pca + "PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz";
            if (Gpio.exists(nextInExp)) {
                outputDir = out + "Cis-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/";
                performeQTLMapping(true, false, nextInExp, outputDir, cisSnpsToTest, null, threads);
                cleanup();

                outputDir = out + "Trans-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/";
                performeQTLMapping(false, true, nextInExp, outputDir, transSnpsToTest, null, threads);
                cleanup();
            } else {
                System.out.println("Error could not find: " + nextInExp);
            }

        }

//            String cisNull = out + "Cis-0PCAsRemoved/eQTLsFDR0.05.txt";
//            String transNull = out + "Trans-0PCAsRemoved/eQTLsFDR0.05.txt";
//
//            String cisOut  = out + "Cis-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/eQTLsFDR0.05.txt";
//            String transOut  = out + "Trans-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/eQTLsFDR0.05.txt";
//            
//            e.compareOverlapAndZScoreDirectionTwoEQTLFiles(cisOut,   cisNull  ,     out+"Cis-" + pca +"PCAsRemoved-GeneticVectorsNotRemoved");
//            e.compareOverlapAndZScoreDirectionTwoEQTLFiles(transOut ,transNull   ,   out+"Trans-" + pca +"PCAsRemoved-GeneticVectorsNotRemoved");
//            
//            eQTLTextFile etf2 = new eQTLTextFile(out + "Cis-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
//            nrCISEQTLsPerRound[round] = etf2.read().length;
//            etf2.close();
//
//            etf2 = new eQTLTextFile(out + "Trans-" + pca + "PCAsRemoved-GeneticVectorsNotRemoved/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
//            nrTransEQTLsPerRound[round] = etf2.read().length;
//            round++;
//	}

//	    if (runpcsremoved == null) {
        PCAOptimumInventorize pi = new PCAOptimumInventorize();

        //TODO: MJ: dit is een hack hier moet ook gekeken worden naar hoeveel PCs ed
        pi.inventorypcqtl(out, true, true, 100, 5);
//	    }
    }
//        System.out.println("PCs\tCis\tShared\tDifferentAllelicDirection\tTrans\tShared\tDifferentAllelicDirection");
//        for (int i = 0; i < 21; i++) {
//            //    System.out.println((i*5)+"\t"+nrCISEQTLsPerRound[i]+"\t"+nrCISEQTLsPerRoundShared[i]+"\t"+nrCISEQTLsPerRoundWithDifferentDirection[i]+"\t"+nrTransEQTLsPerRound[i]+"\t"+nrTranEQTLsPerRoundShared[i]+"\t"+nrTransEQTLsPerRoundWithDifferentDirection[i]);
//            System.out.println((i * 5) + "\t" + nrCISEQTLsPerRound[i] + "\t" + nrTransEQTLsPerRound[i]);
//        }

    private void performeQTLMappingOverEigenvectorMatrixAndReNormalize(String origInExp, String out) throws IOException, Exception {
        // Eigenvector mapping
        TextFile tf = new TextFile(origInExp, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        int nrCols = header.length;
        tf.close();

        int nrToRemove = 101;
        if (nrToRemove > nrCols) {
            nrToRemove = nrCols;
        }

        HashSet<String> probesToTest = new HashSet<String>();
        for (int i = 1; i < nrToRemove; i++) {
            probesToTest.add("Comp" + i);
        }

        // ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.PCAOverSamplesEigenvectorsTransposed
        File f = new File(origInExp + ".PCAOverSamplesEigenvectorsTransposed.txt.gz");
        String startfile = origInExp + ".PCAOverSamplesEigenvectorsTransposed.txt.gz";
        if (!f.exists()) {
            f = new File(origInExp + ".PCAOverSamplesEigenvectorsTransposed.txt");
            if (f.exists()) {
                startfile = origInExp + ".PCAOverSamplesEigenvectorsTransposed.txt";
            } else {
                System.out.println("Error! Could not find " + startfile);
                System.exit(0);
            }

        }

        String nextInExp = startfile;
        performeQTLMapping(true, true, nextInExp, out + "CisTrans-PCAEigenVectors/", null, probesToTest, m_threads);
        cleanup();

        eQTLTextFile etf = new eQTLTextFile(out + "CisTrans-PCAEigenVectors/eQTLProbesFDR0.05.txt", eQTLTextFile.R);
        EQTL[] eigenvectorEQTLs = etf.read();
        etf.close();

        HashSet<Integer> geneticEigenVectors = new HashSet<Integer>();
        for (EQTL e : eigenvectorEQTLs) {
            Double fdr = e.getFDR();
            if (fdr == null) {
                System.out.println("Error with eQTL file!: FDR == null: " + out + "CisTrans-EigenVectors/eQTLProbesFDR0.05.txt");
                System.exit(0);
            }
            if (fdr > 0) {
                break;
            } else {
                String probe = e.getProbe();
                Integer compId = Integer.parseInt(probe.replace("Comp", ""));
                // quick hack: component 1 captures population stratification information...
                if (compId > 1) {
                    geneticEigenVectors.add(compId);
                }
            }
        }

        System.out.println("Repeating PCA analysis, without removal of genetically controlled PCs");
        System.out.println("Components under genetic control: " + Strings.concat(geneticEigenVectors.toArray(new Integer[0]), Strings.comma));
        System.out.println();
        if (geneticEigenVectors.size() > 0) {
            Normalizer n = new Normalizer();
            n.repeatPCAOmitCertainPCAs(geneticEigenVectors, origInExp, 100, 5, covariatesremoved);
        } else {
            System.out.println("No PCA vectors seem to be genetically associated.\n"
                    + "To find the optimum number of PCs to use, rerun this command without --pcqtl, if you have not done so already.");
            System.exit(0);
        }

    }

    
}
