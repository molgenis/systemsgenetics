package betaqtl;

import org.apache.commons.cli.*;

import java.io.IOException;

public class Main {

    // todo: clean up options
    public static void main(String[] args) {

        Options options = new Options();
        options.addOption(OptionBuilder.withLongOpt("mode")
                .withDescription("Mode: [metaqtl|betaqtl|betaqtlsingleds|betaqtlplot|regressqtl|sortfile|determineld|determineldgwas|concatconditional]")
                .isRequired()
                .hasArg()
                .withArgName("STRING")
                .create("m"));

        options.addOption(OptionBuilder.withLongOpt("vcf")
                .withDescription("Tabix indexed VCF")
                .hasArg()
                .withArgName("PATH")
                .create("v"));
        options.addOption(OptionBuilder.withLongOpt("exp")
                .withDescription("Expression matrix (can be gzipped)")
                .hasArg()
                .withArgName("PATH")
                .create("e"));
        options.addOption(OptionBuilder.withLongOpt("chr")
                .withDescription("Chromosome number")
                .hasArg()
                .withArgName("INT")
                .create());
        options.addOption(OptionBuilder.withLongOpt("gte")
                .withDescription("Genotype to expression to dataset linkfile (tab separated)")
                .hasArg()
                .withArgName("PATH")
                .create("g"));

        options.addOption(OptionBuilder.withLongOpt("snpgenelimit")
                .withDescription("SNP-gene limit file (one line per snp-gene combination, tab separated)")
                .hasArg()
                .withArgName("PATH")
                .create("sgl"));
        options.addOption(OptionBuilder.withLongOpt("snplimit")
                .withDescription("SNP limit file (one line per gene ID)")
                .hasArg()
                .withArgName("PATH")
                .create("sl"));

        options.addOption(OptionBuilder.withLongOpt("genelimit")
                .withDescription("Gene limit file (one line per gene ID)")
                .hasArg()
                .withArgName("PATH")
                .create("gl"));
        options.addOption(OptionBuilder.withLongOpt("annotation")
                .withDescription("Gene annotation file")
                .hasArg()
                .withArgName("PATH")
                .create("a"));
        options.addOption(OptionBuilder.withLongOpt("out")
                .withDescription("Output prefix")
                .hasArg()
                .withArgName("PATH")
                .create("o"));

        options.addOption(OptionBuilder.withLongOpt("seed")
                .withDescription("Random seed [default: 123456789]")
                .hasArg()
                .withArgName("LONG")
                .create());

        options.addOption(OptionBuilder.withLongOpt("perm")
                .withDescription("Number of permutations [default: 1000]")
                .hasArg()
                .withArgName("INT")
                .create());

        options.addOption(OptionBuilder.withLongOpt("ciswindow")
                .withDescription("Cis window size [default: 1mb]")
                .withArgName("INT")
                .hasArg()
                .create());

        options.addOption(OptionBuilder.withLongOpt("maf")
                .withDescription("Minor allele frequency threshold [default: 0.01]")
                .hasArg()
                .withArgName("FLOAT")
                .create());
        options.addOption(OptionBuilder.withLongOpt("cr")
                .withDescription("Call-rate threshold [default: 0.95]")
                .hasArg()
                .withArgName("FLOAT")
                .create());
        options.addOption(OptionBuilder.withLongOpt("hwep")
                .withDescription("Hardy-Weinberg p-value threshold [default: 0.0001]")
                .hasArg()
                .withArgName("FLOAT")
                .create());
        options.addOption(OptionBuilder.withLongOpt("snpannotation")
                .withDescription("SNP annotation file, tab separated: SNPID chr pos")
                .hasArg()
                .withArgName("FILE")
                .create());
        options.addOption(OptionBuilder.withLongOpt("minobservations")
                .withDescription("Require at least this many observations per dataset (i.e. non-NaN genotypes/phenotypes) [default: 10]")
                .hasArg()
                .withArgName("FLOAT")
                .create());
        options.addOption(OptionBuilder.withLongOpt("norank")
                .withDescription("Do not rank expression data")
                .create());
        options.addOption(OptionBuilder.withLongOpt("outputall")
                .withDescription("Output all associations, not just top association per gene")
                .create());
        options.addOption(OptionBuilder.withLongOpt("outputallpermutations")
                .withDescription("Output all permuted associations, not just top association per gene")
                .create());
        options.addOption(OptionBuilder.withLongOpt("snplog")
                .withDescription("Output SNP summary stats per snp/gene pair.")
                .create());
        options.addOption(OptionBuilder.withLongOpt("nrdatasets")
                .withDescription("Minimum number of datasets required in meta-analysis [default: 2]")
                .hasArg()
                .withArgName("INT")
                .create());
        options.addOption(OptionBuilder.withLongOpt("replacemissinggenotypes")
                .withDescription("Replace missing genotypes with average genotype: use this when both genotypes and expression data have missing values and perm > 0")
                .create());

        options.addOption(OptionBuilder.withLongOpt("input")
                .withDescription("Input file")
                .hasArg()
                .create());
        options.addOption(OptionBuilder.withLongOpt("input2")
                .withDescription("Input file 2")
                .hasArg()
                .create());
        options.addOption(OptionBuilder.withLongOpt("nriters")
                .withDescription("Nr iters to concatenate using concatconditional")
                .hasArg()
                .create());
        options.addOption(OptionBuilder.withLongOpt("sortbyz")
                .withDescription("Sort by Z-score")
                .create());


        options.addOption(OptionBuilder.withLongOpt("eqtlset")
                .withDescription("[determineldgwas] - List of eQTL snps to test - txt.gz file")
                .hasArg()
                .withArgName("FILE")
                .create());

        options.addOption(OptionBuilder.withLongOpt("eqtlfile")
                .withDescription("[determineldgwas] - eQTL file - txt.gz file")
                .hasArg()
                .withArgName("FILE")
                .create());

        options.addOption(OptionBuilder.withLongOpt("gwasset")
                .withDescription("[determineldgwas] - List of GWAS snps to test - txt.gz file")
                .hasArg()
                .withArgName("FILE")
                .create());

        options.addOption(OptionBuilder.withLongOpt("gwasassocfile")
                .withDescription("[determineldgwas] - GWAS association file - txt.gz file")
                .hasArg()
                .withArgName("FILE")
                .create());

        options.addOption(OptionBuilder.withLongOpt("gwaslistfile")
                .withDescription("[determineldgwas] - GWAS study file - txt.gz file")
                .hasArg()
                .withArgName("FILE")
                .create());

        options.addOption(OptionBuilder.withLongOpt("prunedistance")
                .withDescription("[determineldgwas] - Pruning distance [default: 1mb]")
                .hasArg()
                .withArgName("INT")
                .create());

        options.addOption(OptionBuilder.withLongOpt("prunethreshold")
                .withDescription("[determineldgwas] - Pruning LD threshold [default: 0.2]")
                .hasArg()
                .withArgName("FLOAT")
                .create());

        options.addOption(OptionBuilder.withLongOpt("ldthreshold")
                .withDescription("[determineldgwas] - LD threshold [default: 0.8]")
                .hasArg()
                .withArgName("FLOAT")
                .create());

        options.addOption(OptionBuilder.withLongOpt("skipchr6")
                .withDescription("[determineldgwas] - Skip chr6 when calculating LD")
                .create());

        options.addOption(OptionBuilder.withLongOpt("matchbyrsid")
                .withDescription("[determineldgwas] - Match by RsId in stead of full variant id")
                .create());


        try {
            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, args);

            String mode = cmd.getOptionValue("mode");

            String vcf = null;
            if (cmd.hasOption("vcf")) {
                vcf = cmd.getOptionValue("vcf");
            }

            int chrom = -1;
            if (cmd.hasOption("chr")) {
                chrom = Integer.parseInt(cmd.getOptionValue("chr"));
            }

            String linkfile = null;
            if (cmd.hasOption("gte")) {
                linkfile = cmd.getOptionValue("gte");
            }

            String genelimit = null;
            if (cmd.hasOption("genelimit")) {
                genelimit = cmd.getOptionValue("genelimit");
            }

            String snplimit = null;
            if (cmd.hasOption("snplimit")) {
                snplimit = cmd.getOptionValue("snplimit");
            }
            String snpgenelimit = null;
            if (cmd.hasOption("snpgenelimit")) {
                snpgenelimit = cmd.getOptionValue("snpgenelimit");
            }

            String genexpression = null;
            if (cmd.hasOption("exp")) {
                genexpression = cmd.getOptionValue("exp");
            }
            String geneannotation = null;
            if (cmd.hasOption("annotation")) {
                geneannotation = cmd.getOptionValue("annotation");
            }
            String output = null;
            if (cmd.hasOption("out")) {
                output = cmd.getOptionValue("out");
            }

            String input = null;
            if (cmd.hasOption("input")) {
                input = cmd.getOptionValue("input");
            }

            System.out.println(mode);
            switch (mode) {
                case "concatconditional":
                    Integer nriters = null;
                    if (cmd.hasOption("nriters")) {
                        nriters = Integer.parseInt(cmd.getOptionValue("nriters"));
                    }
                    if (input == null || vcf == null || linkfile == null || nriters == null || output == null) {
                        System.out.println("Usage: --input qtlfilestr --vcf vcffile.vcf.gz --out output --gte samplelinkfile.txt --nriters nriters ");
                        System.out.println("Input: " + input);
                        System.out.println("VCF: " + vcf);
                        System.out.println("Output: " + output);
                        System.out.println("SampleLink: " + linkfile);
                        System.out.println("Nr Iters: " + nriters);
                    } else {
                        ConcatenateConditionalResults c = new ConcatenateConditionalResults();
                        c.run(input, vcf, linkfile, nriters, output);
                    }
                    break;
                case "determineld":
                    System.out.println("Determine LD between SNPs in two QTL files");
                    String input2 = null;
                    if (cmd.hasOption("input2")) {
                        input2 = cmd.getOptionValue("input2");
                    }
                    if (input == null || input2 == null || output == null || vcf == null) {
                        System.out.println("Usage: --input qtlfile1 --input2 qtlfile2 --vcf vcffile.vcf.gz --out output [--gte samplelinkfile.txt]");
                        System.out.println("Input: " + input);
                        System.out.println("Input2 " + input2);
                        System.out.println("Output: " + output);
                        System.out.println("VCF: " + vcf);
                        System.out.println("GTE: " + linkfile);
                    } else {
                        LDCalculator calc = new LDCalculator();
                        calc.run(input, input2, vcf, linkfile, output);
                    }
                    break;
                case "determineldgwas":
                    System.out.println("Determine LD between QTL SNPs and GWAS SNPs");
                    String eqtlset = null;
                    String eqtlfile = null;
                    String gwasset = null;
                    String gwasassocfile = null;
                    String gwaslistfile = null;
                    int prunedistance = 1000000;
                    double prunethreshold = 0.2;
                    double ldthreshold = 0.8;
                    boolean skipchr6 = false;
                    boolean matchbyrsid = false;

                    if (cmd.hasOption("eqtlset")) {
                        eqtlset = cmd.getOptionValue("eqtlset");
                    }
                    if (cmd.hasOption("eqtlfile")) {
                        eqtlfile = cmd.getOptionValue("eqtlfile");
                    }
                    if (cmd.hasOption("gwasset")) {
                        gwasset = cmd.getOptionValue("gwasset");
                    }
                    if (cmd.hasOption("gwasassocfile")) {
                        gwasassocfile = cmd.getOptionValue("gwasassocfile");
                    }
                    if (cmd.hasOption("gwaslistfile")) {
                        gwaslistfile = cmd.getOptionValue("gwaslistfile");
                    }
                    if (cmd.hasOption("prunedistance")) {
                        prunedistance = Integer.parseInt(cmd.getOptionValue("prunedistance"));
                    }
                    if (cmd.hasOption("prunethreshold")) {
                        prunethreshold = Double.parseDouble(cmd.getOptionValue("prunethreshold"));
                    }
                    if (cmd.hasOption("ldthreshold")) {
                        ldthreshold = Double.parseDouble("ldthreshold");
                    }
                    if (cmd.hasOption("skipchr6")) {
                        skipchr6 = true;
                    }
                    if (cmd.hasOption("matchbyrsid")) {
                        matchbyrsid = true;
                    }

                    if (eqtlset == null || eqtlfile == null
                            || gwasset == null || gwasassocfile == null || gwaslistfile == null
                            || vcf == null || gwasassocfile == null || gwasassocfile == null || gwasassocfile == null) {
                        System.out.println("Usage: --eqtlset eqtlset.txt.gz " +
                                "--eqtlfile eqtlfile.txt.gz " +
                                "--gwasset gwasset.txt.gz " +
                                "--gwasassocfile gwasassocfile.txt.gz " +
                                "--gwaslistfile gwaslist.txt.gz " +
                                "--vcf vcfCHR.vcf.gz" +
                                "[--gte samplelinkfile.txt] " +
                                "[--prunedistance 1E6] " +
                                "[--prunethreshold 0.2] " +
                                "[--ldthreshold 0.8] " +
                                "--out output.txt.gz " +
                                "[--skipchr6] " +
                                "[--matchbyrsid]");
                        System.out.println("eqtlset: " + eqtlset);
                        System.out.println("eqtlfile: " + eqtlfile);
                        System.out.println("gwasset: " + gwasset);
                        System.out.println("gwasassocfile: " + gwasassocfile);
                        System.out.println("gwaslistfile: " + gwaslistfile);
                        System.out.println("vcf: " + vcf);
                        System.out.println("gte: " + linkfile);
                        System.out.println("prunedistance: " + prunedistance);
                        System.out.println("prunethreshold: " + prunethreshold);
                        System.out.println("ldthreshold: " + ldthreshold);
                        System.out.println("out: " + output);
                        System.out.println("skipchr6: " + skipchr6);
                        System.out.println("matchbyrsid: " + matchbyrsid);

                    } else {
                        GWASLDCalculator calc = new GWASLDCalculator();

                        calc.runPerGWAS3(eqtlset, eqtlfile, gwasset, gwasassocfile, gwaslistfile, vcf, linkfile, prunedistance, prunethreshold, ldthreshold, output, skipchr6, matchbyrsid);

                    }
                    break;
                case "sortfile":
                    System.out.println("QTL file sorter sorts by position by default; use --sortz to sort by Z-score");
                    QTLFileSorter sorter = new QTLFileSorter();
                    if (input == null || output == null) {
                        System.out.println("Usage: --input infile.txt[.gz] --out outfile.txt[.gz] [--sortbyz|--sortbysnppos|--sortbygenepos]");
                        System.out.println("Input: " + input);
                        System.out.println("Output: " + output);
                        System.out.println("Sort by Z: " + cmd.hasOption("sortbyz"));
                    } else {
                        if (cmd.hasOption("sortbyz")) {
                            sorter.run(input, output, QTLFileSorter.SORTBY.Z);
                        } else if (cmd.hasOption("sortbysnppos")) {
                            sorter.run(input, output, QTLFileSorter.SORTBY.SNPPOS);
                        } else {
                            sorter.run(input, output, QTLFileSorter.SORTBY.GENEPOS);
                        }
                    }
                    break;
                case "regressqtl":
                    if (vcf == null || linkfile == null || geneannotation == null || genexpression == null || output == null || snpgenelimit == null) {
                        System.err.println("Usage: -m regressqtl --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ --snpgenelimit snpgenecombos.txt.gz");
                        System.err.println("Optional: --snpannotation snpannotation.txt.gz --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6 --nrdatasets 2 --genelimit");
                        System.out.println("VCF: " + vcf);
                        System.out.println("Linkfile: " + linkfile);
                        System.out.println("Annotation: " + geneannotation);
                        System.out.println("Exp: " + genexpression);
                        System.out.println("Out: " + output);
                        System.out.println("SNP/Gene combos: " + snpgenelimit);
                        System.out.println("Genelimit: " + genelimit);
                    } else {
                        QTLRegression qtlr = new QTLRegression(vcf, chrom, linkfile, snplimit, genelimit, snpgenelimit, genexpression, geneannotation, output);
                        if (cmd.hasOption("replacemissinggenotypes")) {
                            qtlr.setReplaceMissingGenotypes(true);
                        }

                        if (cmd.hasOption("norank")) {
                            qtlr.setRankData(false);
                        }

                        if (cmd.hasOption("minobservations")) {
                            int t = Integer.parseInt(cmd.getOptionValue("minobservations"));
                            qtlr.setMinObservations(t);
                        }

                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            qtlr.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            qtlr.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            qtlr.setHwepthreshold(t);
                        }

                        if (cmd.hasOption("ciswindow")) {
                            int t = Integer.parseInt(cmd.getOptionValue("ciswindow"));
                            qtlr.setCisWindow(t);
                        }

                        if (cmd.hasOption("nrdatasets")) {
                            int t = Integer.parseInt(cmd.getOptionValue("nrdatasets"));
                            qtlr.setMinNumberOfDatasets(t);
                        }
                        if (cmd.hasOption("snpannotation")) {
                            qtlr.loadSNPAnnotation(cmd.getOptionValue("snpannotation"));
                        }
                        qtlr.run();
                    }
                    break;
                case "betaqtlplot":
                    if (vcf == null || chrom == -1 || linkfile == null || geneannotation == null || genexpression == null || output == null) {
                        System.err.println("Required: --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ ");
                        System.err.println("Optional: --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6 --nrdatasets 2");
                        System.out.println("VCF: " + vcf);
                        System.out.println("Chrom: " + chrom);
                        System.out.println("GTE:" + linkfile);
                        System.out.println("Gene annotation: " + geneannotation);
                        System.out.println("Gene expression: " + genexpression);
                        System.out.println("Output: " + output);
                    } else {
                        BetaQTLPlot bpp = new BetaQTLPlot(vcf, chrom, linkfile, snplimit, genelimit, snpgenelimit, genexpression, geneannotation, output);
                        if (cmd.hasOption("replacemissinggenotypes")) {
                            bpp.setReplaceMissingGenotypes(true);
                        }
                        if (cmd.hasOption("norank")) {
                            bpp.setRankData(false);
                        }
                        // bQTL.setMinObservations();
                        if (cmd.hasOption("outputall")) {
                            bpp.setOutputAll(true);
                        }
                        if (cmd.hasOption("minobservations")) {
                            int t = Integer.parseInt(cmd.getOptionValue("minobservations"));
                            bpp.setMinObservations(t);
                        }

                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            bpp.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            bpp.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            bpp.setHwepthreshold(t);
                        }

                        if (cmd.hasOption("nrdatasets")) {
                            int t = Integer.parseInt(cmd.getOptionValue("nrdatasets"));
                            bpp.setMinNumberOfDatasets(t);
                        }
                        bpp.plot();
                    }
                    break;
                case "betaqtl":
                    if (vcf == null || chrom == -1 || linkfile == null || geneannotation == null || genexpression == null || output == null) {
                        System.err.println("Usage: --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ ");
                        System.err.println("Optional: --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6 --nrdatasets 2");
                        System.err.println("Optional: --seed 123456789 --outputall --perm 1000 --snplog --outputallpermutations");

                        System.out.println("You've set the following:");
                        System.out.println("VCF: " + vcf);
                        System.out.println("Chrom: " + chrom);
                        System.out.println("GTE:" + linkfile);
                        System.out.println("Gene annotation: " + geneannotation);
                        System.out.println("Gene expression: " + genexpression);
                        System.out.println("Output: " + output);
                    } else {
                        BetaQTL2ParallelCis bQTL = new BetaQTL2ParallelCis(vcf, chrom, linkfile, snplimit, genelimit, snpgenelimit, genexpression, geneannotation, output);
                        if (cmd.hasOption("snpannotation")) {
                            bQTL.loadSNPAnnotation(cmd.getOptionValue("snpannotation"));
                        }
                        if (cmd.hasOption("replacemissinggenotypes")) {
                            bQTL.setReplaceMissingGenotypes(true);
                        }
                        if (cmd.hasOption("norank")) {
                            bQTL.setRankData(false);
                        }
                        // bQTL.setMinObservations();
                        if (cmd.hasOption("outputall")) {
                            bQTL.setOutputAll(true);
                        }
                        if (cmd.hasOption("outputallpermutations")) {
                            bQTL.setOutputAllPermutations(true);
                        }
                        if (cmd.hasOption("snplog")) {
                            bQTL.setOutputSNPLog(true);
                        }
                        if (cmd.hasOption("minobservations")) {
                            int t = Integer.parseInt(cmd.getOptionValue("minobservations"));
                            bQTL.setMinObservations(t);
                        }

                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            bQTL.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            bQTL.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            bQTL.setHwepthreshold(t);
                        }

                        if (cmd.hasOption("ciswindow")) {
                            int t = Integer.parseInt(cmd.getOptionValue("ciswindow"));
                            bQTL.setCisWindow(t);
                        }
                        if (cmd.hasOption("perm")) {
                            int t = Integer.parseInt(cmd.getOptionValue("perm"));
                            bQTL.setNrPermutations(t);
                        }
                        if (cmd.hasOption("seed")) {
                            long seed = Long.parseLong(cmd.getOptionValue("seed"));
                            bQTL.setRandomSeed(seed);
                        }
                        if (cmd.hasOption("nrdatasets")) {
                            int t = Integer.parseInt(cmd.getOptionValue("nrdatasets"));
                            bQTL.setMinNumberOfDatasets(t);
                        }
                        bQTL.run();
                    }
                    break;
                case "betaqtlsingleds":
                    if (vcf == null || chrom == -1 || linkfile == null || geneannotation == null || genexpression == null || output == null) {
                        System.err.println("Usage: --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ ");
                        System.err.println("Optional: --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6");
                        System.err.println("Optional: --seed 123456789 --outputall --perm 1000 ");
                    } else {

                        BetaQTLSingleDataset ds = new BetaQTLSingleDataset(vcf, chrom, linkfile, genelimit, genexpression, geneannotation, output);

                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            ds.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            ds.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            ds.setHwepthreshold(t);
                        }
                        if (cmd.hasOption("ciswindow")) {
                            int t = Integer.parseInt(cmd.getOptionValue("ciswindow"));
                            ds.setCisWindow(t);
                        }
                        if (cmd.hasOption("perm")) {
                            int t = Integer.parseInt(cmd.getOptionValue("perm"));
                            ds.setNrPermutations(t);
                        }
                        ds.fastqtlclone2();
                    }
                    break;
                case "metaqtl":
                    System.out.println("Not implemented yet.");
                    break;
                default:

            }


        } catch (ParseException e) {
            System.out.println("Command line parse exception: " + e.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("betaqtl.jar", options);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
