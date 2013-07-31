/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import eqtlmappingpipeline.textmeta.FixedEffectMetaAnalysis;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.util.eqtlfilesorter.EQTLFileSorter;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class UtilConsoleGUI {

    public static enum MODE {

        GETSNPSFROMREGION, GETPROBESFROMREGION, GETSNPSINPROBEREGION, FDR, GETMAF, MERGE, REGRESS, GETSNPSTATS, PROXYSEARCH, DOTPLOT, META, SORTFILE, CONVERTBINARYMATRIX
    };
    MODE run;

    public UtilConsoleGUI(String[] args) {

        String settingsfile = null;
        String settingstexttoreplace = null;
        String settingstexttoreplacewith = null;
        String in = null;
        String in2 = null;
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
        String probefile = null;
        String region = "";

        Double threshold = null;
        Integer nreqtls = null;

        Double r2 = null;
        Double maf = 0.05;
        Double cr = 0.95;
        Double hwep = 0.001;
        Integer dist = 1000000;

        Integer minnrdatasets = null;
        Integer minnrsamples = null;

        boolean createQQPlot = true;

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
            } else if (arg.equals("--in2")) {
                in2 = val;
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
            } else if (arg.equals("--convertbinarymatrix")) {
                region = val;
                run = MODE.CONVERTBINARYMATRIX;
            } else if (arg.equals("--getsnpsinregion")) {
                region = val;
                run = MODE.GETSNPSFROMREGION;
            } else if (arg.equals("--sortfile")) {
                region = val;
                run = MODE.SORTFILE;
            } else if (arg.equals("--findproxy")) {
                region = val;
                run = MODE.PROXYSEARCH;
            } else if (arg.equals("--getmaf")) {
                region = val;
                run = MODE.GETMAF;
            } else if (arg.equals("--getsnpsinproberegion")) {
                region = val;
                run = MODE.GETSNPSINPROBEREGION;
            } else if (arg.equals("--merge")) {
                run = MODE.MERGE;
            } else if (arg.equals("--fdr")) {
                region = val;
                run = MODE.FDR;
            } else if (arg.equals("--dotplot")) {
                region = val;
                run = MODE.DOTPLOT;
            } else if (arg.equals("--regress")) {
                run = MODE.REGRESS;
            } else if (arg.equals("--snpstats")) {
                run = MODE.GETSNPSTATS;
            } else if (arg.equals("--meta")) {
                run = MODE.META;
            } else if (arg.equals("--snps")) {
                snpfile = val;
            } else if (arg.equals("--probes")) {
                probefile = val;
            } else if (arg.equals("--perm")) {
                perm = Integer.parseInt(val);
            } else if (arg.equals("--nreqtls")) {
                nreqtls = Integer.parseInt(val);
            } else if (arg.equals("--threshold")) {
                threshold = Double.parseDouble(val);
            } else if (arg.equals("--r2")) {
                r2 = Double.parseDouble(val);
            } else if (arg.equals("--maf")) {
                maf = Double.parseDouble(val);
            } else if (arg.equals("--hwep")) {
                hwep = Double.parseDouble(val);
            } else if (arg.equals("--dist")) {
                dist = Integer.parseInt(val);
            } else if (arg.equals("--skipqqplot")) {
                createQQPlot = false;
            }
        }
        if (run == null) {
            System.err.println("Please specify an util.");
            printUsage();
        } else {
            try {
                switch (run) {
                    case CONVERTBINARYMATRIX:
                        if (in == null || out == null) {
                            System.out.println("Usage: --util --convertbinarymatrix --in /path/to/matrix.binary --out /path/to/textoutput.txt");
                        } else {
                            if (in.endsWith(".txt")) {
                                System.out.println("The file provided with --in is already a text file: " + in);
                            } else {
                                if (in.endsWith(".dat")) {
                                    in = in.substring(0, in.length() - 4);
                                }
                                System.out.println("Converting: " + in);
                                DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(in);
                                ds.save(out);
                            }
                        }
                        break;


                    case REGRESS:

                        RegressCisEffectsFromGeneExpressionData r = new RegressCisEffectsFromGeneExpressionData(args);
                        break;
                    case PROXYSEARCH:

                        if (in == null || snpfile == null || out == null || r2 == null) {
                            System.out.println("Usage: --mode util --findproxy --r2 0.8 --snps snpfile.txt --out outfile --in /Path/To/TriTyperReference/ [--hwep 0.001] [--maf 0.05] [--cr 0.95]");
                        } else {
                            LDCalculator.proxyLookUpInReferenceDataset(in, snpfile, maf, hwep, cr, r2, out, dist);
                        }
                        break;
                    case MERGE:

                        if (in == null || region == null) {
                            System.out.println("USAGE: --merge --in dataset --in2 dataset2 --out outdir [--snps snpfile]");
                        } else {
                            GenotypeDataMerger m = new GenotypeDataMerger();
                            m.merge(in, in2, out, snpfile);
                        }
                        break;
                    case GETMAF:

                        if (in == null || region == null) {
                            System.out.println("USAGE: --getmaf snplistfile --in dataset");
                        } else {
                            GenotypeDataQuery dq = new GenotypeDataQuery();
                            dq.getSNPMAF(in, region);
                        }
                        break;
                    case GETSNPSTATS:

                        if (in == null || region == null) {
                            System.out.println("USAGE: --in dataset");
                        } else {
                            GenotypeDataQuery dq = new GenotypeDataQuery();
                            if (in2 != null) {
                                dq.getSNPStatsForAllSNPs(in, in2);
                            } else {
                                dq.getSNPStatsForAllSNPs(in);
                            }
                        }
                        break;
                    case SORTFILE:
                        if (in == null) {
                            System.out.println("USAGE: --in eQTLFile --out eQTLFile");
                        } else {
                            EQTLFileSorter f = new EQTLFileSorter();
                            f.run(in, out);
                        }
                        break;
                    case GETSNPSFROMREGION:
                        if (in == null || region == null) {
                            System.out.println("To use --getsnpsfromregion, please use --in to point to the genotype data and supply a region to query.");
                            printUsage();
                        } else {
                            int chr = -1;
                            int chrposA = -1;
                            int chrposB = -1;
                            GenotypeDataQuery q = new GenotypeDataQuery();
                            try {
                                String[] elems = region.split(":");
                                chr = ChrAnnotation.parseChr(elems[0]);
                                elems = elems[1].split("-");
                                chrposA = Integer.parseInt(elems[0]);
                                chrposB = Integer.parseInt(elems[1]);
                            } catch (Exception e) {
                                System.err.println("Error: malformed query: " + region);;
                            }
                            q.getSNPsInRegion(in, chr, chrposA, chrposB);
                        }

                        break;

                    case GETSNPSINPROBEREGION:
                        if (snpfile == null || inexpannot == null || probefile == null) {
                            System.out.println("To use --getsnpsinproberegion, please use --snps, --probes, and --inexpannot");
                            printUsage();
                        } else {
                            ProbeSNPMapper psm = new ProbeSNPMapper();
                            psm.mapprobes(snpfile, inexpannot, probefile);
                        }

                        break;
                    case FDR:
                        if (in == null || threshold == null || nreqtls == null) {
                            System.out.println("To use --fdr, please use --in, --threshold, and --perm and --nreqtls");
                            printUsage();
                        } else {
                            FDR.calculateFDR(in, perm, nreqtls, threshold, createQQPlot, null, null);
                        }

                        break;
                    case META:
                        if (in == null || out == null) {
                            System.out.println("To use --meta, please use --in, and --out");
                            printUsage();
                        } else {
                            FixedEffectMetaAnalysis f = new FixedEffectMetaAnalysis();
                            f.run(in, out, minnrdatasets, minnrsamples);
                        }

                        break;
                    case DOTPLOT:
                        if (in == null) {
                            System.out.println("Usage: --dotplot --in /path/to/file.txt");
                        } else {
                            eQTLDotPlotter d = new eQTLDotPlotter();
                            d.plot(in);
                        }
                        break;

                }

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(-1);
            }
        }


    }

    private void printUsage() {
        System.out.print("\tUtil\n" + ConsoleGUIElems.LINE);
        System.out.println("Util contains small utilities.");

        System.out.println("");
        System.out.print("Available Utilities:\n" + ConsoleGUIElems.LINE);

        System.out.println("--getsnpsinregion\t\tGet SNPs in a certain region: chr positionA positionB: Y:12000-13000 would get all SNPs on chr Y between 12000 and 13000 bp\n"
                + "--getsnpsinproberegion\t\tGet SNPs in a certain set of probes (specify with --probes)\n"
                + "--getmaf\t\t\tgets maf for snp\n"
                + "--merge\t\t\t\tmerges two datasets\n"
                + "--snpstats\t\t\tGets HWE, MAF, and CR for all SNPs\n"
                + "--findproxy\t\t\tSearches for a proxy given a list of SNPs\n"
                + "--dotplot\t\t\tCreates dotplot from eQTL result file\n"
                + "--regress\t\t\tRemoves eQTL effects from gene expression data.\n"
                + "--convertbinarymatrix\t\t\tConverts binary matrix to text\n");
        System.out.println("");

        System.out.print("Command line options:\n" + ConsoleGUIElems.LINE);
        System.out.println("--in\t\t\tdir\t\tLocation of the genotype data\n"
                + "--out\t\t\tdir\t\tLocation where the output should be stored\n"
                + "--inexp\t\t\tstring\t\tLocation of expression data\n"
                + "--inexpplatform\t\tstring\t\tGene expression platform\n"
                + "--inexpannot\t\tstring\t\tLocation of annotation file for gene expression data\n"
                + "--gte\t\t\tstring\t\tLocation of genotype to expression coupling file\n"
                + "--snps\t\t\tstring\t\tLocation of snp file\n"
                + "--probes\t\tstring\t\tLocation of probe file\n");

        System.out.println("");
    }
}
