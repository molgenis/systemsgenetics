/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.interactionanalysis;

import java.io.IOException;
import umcg.genetica.console.ConsoleGUIElems;

/**
 *
 * @author harm-jan
 */
public class InteractionAnalysisConsoleGUI {

    enum RUNMODE {

        NORMALIZE, CELLTYPESPECIFICEQTLMAPPING, PLOT
    };

    /**
     * @param args the command line arguments
     */
    public InteractionAnalysisConsoleGUI(String[] args) {
        String inexpraw = null;
        String out = null;
        String celltypespecificprobefile = null;
        String mdscomponents = null;
        String cellcountfile = null;
        String in = null;
        String gte = null;
        String snpprobecombofile = null;
        String covariates = null;
        String inexp = null;
        RUNMODE step = null;
        boolean binaryoutput = false;

        boolean matchCovariateNamesToExpressionProbeNames = false;
        Integer nrThreads = null;
        String covariateList = null;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--step")) {
                if (val == null) {

                } else if (val.equals("normalize")) {
                    step = RUNMODE.NORMALIZE;
                } else if (val.equals("mapeqtls")) {
                    step = RUNMODE.CELLTYPESPECIFICEQTLMAPPING;
                } else if (val.equals("plot")) {
                    step = RUNMODE.PLOT;
                }
            } else if (arg.equals("--inexpraw")) {
                inexpraw = val;
            } else if (arg.equals("--covariatelist")){
                covariateList = val;
            } else if (arg.equals("--binary")) {
                binaryoutput = true;
            } else if (arg.equals("--covariates")) {
                covariates = val;
            } else if (arg.equals("--inexp")) {
                inexp = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--celltypespecificprobes")) {
                celltypespecificprobefile = val;
            } else if (arg.equals("--mdscomponents")) {
                mdscomponents = val;
            } else if (arg.equals("--cellcounts")) {
                cellcountfile = val;
            } else if (arg.equals("--gte")) {
                gte = val;
            } else if (arg.equals("--snpprobe")) {
                snpprobecombofile = val;
            } else if (arg.equals("--testMatchingCovariates")) {
                matchCovariateNamesToExpressionProbeNames = true;
            } else if (arg.equals("--threads")) {
                try {
                    nrThreads = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.err.println("ERROR: value supplied for --threads is not a numerical value.");
                    System.exit(-1);
                }
                if (nrThreads != null && nrThreads < 1) {
                    System.err.println("ERROR: value supplied for --threads is smaller than 1.");
                    System.exit(-1);
                }

            }
        }

        if (step == null) {
            System.err.println("ERROR: please select the step to run.");
            printUsage();
        }

        try {
            if (step == RUNMODE.PLOT) {
                System.out.println("Interaction plotter");
                boolean kill = false;
                if (covariates == null) {
                    System.err.println("Error: please supply --covariates");
                    kill = true;
                }
                if (in == null) {
                    System.err.println("Error: please supply --in");
                    kill = true;
                }
                if (inexp == null) {
                    System.err.println("Error: please supply --inexp");
                    kill = true;
                }
                if (out == null) {
                    System.err.println("Error: please supply --out");
                    kill = true;
                }
                if (kill) {
                    System.err.println("");
                    printUsage();
                } else {
                    InteractionPlotter plotter = new InteractionPlotter(snpprobecombofile, in, inexp, covariates, gte, out);
                }
            } else {
                InteractionAnalysisMultiThreaded qmt = new InteractionAnalysisMultiThreaded();
                if (step == RUNMODE.NORMALIZE) {
                    System.out.println("Cell type specific cis-eQTL normalization");
                    boolean kill = false;
                    if (inexpraw == null) {
                        System.err.println("Error: please supply --inexpraw");
                        kill = true;
                    }
                    if (out == null) {
                        System.err.println("Error: please supply --out");
                        kill = true;
                    }
                    if (celltypespecificprobefile == null) {
                        System.err.println("Error: please supply --celltypespecificprobes");
                        kill = true;
                    }
                    if (kill) {
                        System.err.println("");
                        printUsage();
                    } else {
                        qmt.prepareDataForCelltypeSpecificEQTLMapping(inexpraw, out, null, celltypespecificprobefile, mdscomponents, cellcountfile, gte, nrThreads);
                    }
                } else if (step == RUNMODE.CELLTYPESPECIFICEQTLMAPPING) {
                    System.out.println("Cell type specific cis-eQTL mapping");
                    boolean kill = false;
                    if (covariates == null) {
                        System.err.println("Error: please supply --covariates");
                        kill = true;
                    }
                    if (inexp == null) {
                        System.err.println("Error: please supply --inexp");
                        kill = true;
                    }
                    if (out == null) {
                        System.err.println("Error: please supply --out");
                        kill = true;
                    }
                    if (cellcountfile == null) {
//                    System.err.println("Warning: yo please supply --cellcounts");
                        //kill = true;
                    }

                    if (kill) {
                        System.err.println("");
                        printUsage();
                    } else {
                        qmt.runInteractionAnalysis(inexp, 
                                covariates, 
                                in, 
                                gte, 
                                snpprobecombofile, 
                                nrThreads, 
                                out, 
                                covariateList);
//                    qmt.runCelltypeSpecificEQTLMapping(inexppccorrected, inexpraw, in, gte, snpprobecombofile, cellcountfile, nrThreads, out, testAllCovariatesInCovariateData);
                    }

                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void printUsage() {
        System.out.print("\nCell type specific eQTL Mapping\n" + ConsoleGUIElems.LINE);
        System.out.println("This program uses an OLS model to test eQTLs for cell type specificity.");

        System.out.println("");
        System.out.print("Step 1: Normalization\n" + ConsoleGUIElems.LINE);
        System.out.println("--step normalize\t\t\t\tTell the program to run normalization.\n"
                + "--inexpraw\t\t\tdir\t\tLocation of the gene expression data\n"
                + "--out\t\t\t\tdir\t\tLocation where the output should be stored\n"
                + "--celltypespecificprobes\tString\t\tLocation of the file containing list of cell-type specific probes\n"
                + "--mdscomponents\t\t\tString\t\tLocation of the file containing MDS components (optional)\n"
                + "--gte\t\t\t\tString\t\tLocation of the genotype to expression coupling file (optional)\n"
                + "--cellcounts\t\t\tString\t\tLocation of the cell count file (optional)\n");

        System.out.println("");

        System.out.print("Step 2: Mapping eQTLs with interaction model\n" + ConsoleGUIElems.LINE);
        System.out.println("--step mapeqtls\t\t\t\tTell the program to map eQTLs.\n"
                + "--inexp\tdir\t\tLocation of the dependent dataset\n"
                + "--covariates\t\tdir\t\tLocation of covariate file (the raw gene expression data or the matrix containing the covariates to analyze)\n"
                + "--gte\t\t\tString\t\tLocation of the genotype to expression coupling file\n"
                + "--in\t\t\tdir\t\tLocation of the genotype data\n"
                + "--out\t\t\tdir\t\tLocation where the output should be stored\n"
                + "--snpprobe\t\tString\t\tLocation of the SNP-Probe combination file\n"
                + "--threads\t\tInteger\t\tThe number of threads to use for calculations.\n"
                + "--covariatelist\t\tList of covariates to test\n");

        System.out.println("");

        System.out.print("Step 3: Plot effects\n" + ConsoleGUIElems.LINE);
        System.out.println("--step plot\t\t\t\tTell the program to plot interaction effects.\n"
                + "--inexp\tdir\t\tLocation of the dependent dataset\n"
                + "--covariates\t\tdir\t\tLocation of covariate file (the raw gene expression data or the matrix containing the covariates to analyze)\n"
                + "--gte\t\t\tString\t\tLocation of the genotype to expression coupling file\n"
                + "--in\t\t\tdir\t\tLocation of the genotype data\n"
                + "--out\t\t\tdir\t\tLocation where the output should be stored\n"
                + "--snpprobe\t\tString\t\tLocation of the SNP-Covariate-Probe combination file\n"
        );

    }
}
