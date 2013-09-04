/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.normalization;

import java.io.IOException;
import umcg.genetica.io.Gpio;

/**
 *
 * @author harmjan
 */
public class NormalizationConsoleGUI {

    public NormalizationConsoleGUI(String[] args) {

        String in = null;
        String out = null;
        String cov = null;

        boolean fullNorm = true;
        boolean runLogTransform = false;
        boolean runMTransform = false;
        boolean runQQNorm = false;
        boolean runCenterScale = false;
        boolean runCovariateAdjustment = false;
        boolean runPCAdjustment = false;
        boolean orthogonalizecovariates = false;

        int maxPcaToRemove = 100;
        int stepSizePcaRemoval = 5;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }
            if (arg.equals("--in")) {
                in = val;
            }
            if (arg.equals("--out")) {
                out = val;
            }
            if (arg.equals("--logtransform")) {
                runLogTransform = true;
                fullNorm = false;
            }
            if (arg.equals("--Mtransform")) {
                runMTransform = true;
                fullNorm = false;
            }
            if (arg.equals("--qqnorm")) {
                runQQNorm = true;
                fullNorm = false;
            }
            if (arg.equals("--centerscale")) {
                runCenterScale = true;
                fullNorm = false;
            }
            if (arg.equals("--adjustcovariates")) {
                runCovariateAdjustment = true;
                fullNorm = false;
            }
            if (arg.equals("--adjustPCA")) {
                runPCAdjustment = true;
                fullNorm = false;
            }

            if (arg.equals("--cov")) {
                cov = val;
            }
            if (arg.equals("--covpca")) {
                orthogonalizecovariates = true;
            }

            if (arg.equals("--maxnrpcaremoved")) {
                maxPcaToRemove = Integer.parseInt(val);
            }
            if (arg.equals("--stepsizepcaremoval")) {
                stepSizePcaRemoval = Integer.parseInt(val);
            }
        }

        if (in == null) {
            System.out.println("Please supply your data file.\n");
            printUsage();
            System.exit(0);
        }

        if (!Gpio.exists(in)) {
            System.out.println("Error: the file you specified does not exist.\n");
            System.out.println("Could not find file: " + in);
            System.exit(0);
        }
        if(runLogTransform && runMTransform){
            System.out.println("Error: cant perform both log and M-value transformation");
            System.exit(0);
        }


        try {
            Normalizer p = new Normalizer();
            
            if (!fullNorm) {
                p.normalize(in, maxPcaToRemove, stepSizePcaRemoval, cov, orthogonalizecovariates, out,
                        runQQNorm, runLogTransform, runMTransform, runCenterScale, runPCAdjustment, runCovariateAdjustment);
            } else {
                // run full normalization
                p.normalize(in, maxPcaToRemove, stepSizePcaRemoval, cov, orthogonalizecovariates, out,
                        true, true, false, true, true, true);
            }

            System.out.println("Done.");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void printUsage() {
        System.out.println("Normalization tool\n"
                + "Parameters for full normalization (which includes the following procedures in order: qqnorm, log2transform, centerscale, [covariate adjustment], PCA adjustment):\n"
                + "--in\t\t\tstring\t\tProvide the location of the expression data matrix\n"
                + "--out\t\t\tdir\t\tProvide the directory to output files [optional; default is --in directory]\n"
                + "\n"
                + "Specific normalization methodologies (if multiple methodologies are specified, they will follow the order of full normalization; see above):\n"
                + "--qqnorm\t\t\t\tQuantile normalization\n"
                + "--logtransform\t\t\t\tRun log2 transformation\n"
                + "--Mtransform\t\t\t\tRun M-val (log) transformation for methylation Beta values\n"
                + "--adjustcovariates\t\t\tRun covariate adjustment\n"
                + "--centerscale\t\t\t\tCenter the mean to 0, linearly scale using standard deviation\n"
                + "--adjustPCA\t\t\t\tRun PCA adjustment \n"
                + "\n"
                + "Covariate adjustment parameters:\n"
                + "--cov\t\t\tstring\t\tCovariates to remove\n"
                + "--covpca\t\t\t\tOrthogonalize covariates using PCA before regression\n"
                + "\n"
                + "PCA parameters\n"
                + "--maxnrpcaremoved\tinteger\t\tMaximum number of PCs to remove\n"
                + "--stepsizepcaremoval\tinteger\t\tStep size for PC removal\n");
    }
}
