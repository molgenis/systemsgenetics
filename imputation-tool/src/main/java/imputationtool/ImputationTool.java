/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputationtool;

import imputationtool.postprocessing.BeagleImputationQuality;
import imputationtool.postprocessing.TriTyperDatasetCorrelator;
import imputationtool.postprocessing.TriTypertoPedMapDatExcludedSNPAnalyzer;
import imputationtool.preprocessing.BatchGenerator;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.regex.Pattern;
import umcg.genetica.io.trityper.converters.BeagleImputedToTriTyper;
import umcg.genetica.io.trityper.converters.FinalReportToTriTyper;
import umcg.genetica.io.trityper.converters.ImputeImputedToTriTyperV2;
import umcg.genetica.io.trityper.converters.MachImputedToTriTyper;
import umcg.genetica.io.trityper.converters.MinimacImputedToTriTyper;
import umcg.genetica.io.trityper.converters.PedAndMapToTriTyper;
import umcg.genetica.io.trityper.converters.TriTyperReferenceConcordantPedAndMapExporter;
import umcg.genetica.io.trityper.converters.TriTyperToPedAndMapConverter;
import umcg.genetica.io.trityper.converters.TriTyperToPlinkDosage;
import umcg.genetica.io.trityper.converters.TriTyperToVCF;
import umcg.genetica.io.trityper.converters.VCFToTriTyper;
import umcg.genetica.io.trityper.util.TriTyperConcatDatasets;
import umcg.genetica.io.trityper.util.TriTyperGenotypeDataMerger;

/**
 *
 * @author harm-jan
 */
public class ImputationTool {

    public static final Pattern SEMI_COLON_PATTERN = Pattern.compile(";");
    public static final String VERSION = ImputationTool.class.getPackage().getImplementationVersion();

    public static void main(String[] args) {

        System.out.println("ImputationTool " + VERSION + "\n\n");

        try {
            ImputationTool t = new ImputationTool();

            if (args.length == 0) {
                t.printUsage();
            } else {
                t.processArgs(args);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.exit(0);
    }

    private void processArgs(String[] args) throws IOException, Exception {

        String mode = null;
        String out = null;
        String in = null;
        String inName = null;
        String in2 = null;
        String in2Name = null;
        String fam = null;
        String hap = null;
        String tmp = null;	    // template
        String ext = null;	    // extension
        Integer chr = null;	    // chromosome
        Integer size = null;	    // sample size
        String batch = null;	    // batch file
        String beagle = null;
        String batchdesc = null;    // batch description
        String snps = null;	    // snp list
        Double threshold = null;    // correlation threshold
        String exclude = null;
        Integer chrStart = 1;
        Integer chrEnd = 22;
        Integer nrSamples = null;
        String sampleFile = null;
        String sampleFileToInclude = null;
        String pattern = null;
        String fileMatchRegex = null;

        boolean splitbychromosome = false;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--mode")) {
                mode = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--name")) {
                inName = val;
            } else if (arg.equals("--pattern")) {
                pattern = val;
            } else if (arg.equals("--in2")) {
                in2 = val;
            } else if (arg.equals("--name2")) {
                in2Name = val;
            } else if (arg.equals("--hap")) {
                hap = val;
            } else if (arg.equals("--fam")) {
                fam = val;
            } else if (arg.equals("--tpl")) {
                tmp = val;
            } else if (arg.equals("--chr")) {
                try {
                    chr = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value supplied for --chr is not an integer");
                    System.exit(-1);
                }
                chrStart = chr;
                chrEnd = chr;
            } else if (arg.equals("--batch-file")) {
                batch = val;
            } else if (arg.equals("--size")) {
                try {
                    size = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value supplied for --size is not an integer");
                    System.exit(-1);
                }

            } else if (arg.equals("--exclude")) {
                exclude = val;
            } else if (arg.equals("--beagle")) {
                beagle = val;
            } else if (arg.equals("--batchdesc")) {
                batchdesc = val;
            } else if (arg.equals("--snps")) {
                snps = val;
            } else if (arg.equals("--split")) {
                splitbychromosome = true;
            } else if (arg.equals("--chrstart")) {

                try {
                    chrStart = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value supplied for --chrStart is not an integer");
                    System.exit(-1);
                }
            } else if (arg.equals("--chrend")) {
                try {
                    chrEnd = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value supplied for --chrEnd is not an integer");
                    System.exit(-1);
                }
            } else if (arg.equals("--threshold")) {
                try {
                    threshold = Double.parseDouble(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value supplied for --threshold is not an double");
                    System.exit(-1);
                }
            } else if (arg.equals("--nrSamples")) {
                try {
                    nrSamples = Integer.valueOf(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value supplied for --nrSamples is not an integer");
                    System.exit(-1);
                }
            } else if (arg.equals("--samples")) {
                sampleFile = val;
            } else if (arg.equals("--samplestoinclude")) {
                sampleFileToInclude = val;
            } else if (arg.equals("--fileMatchRegex")) {
                fileMatchRegex = val;
            }
        }

        if (mode == null) {
            System.out.println("Mode is null");

            printUsage();

        } else if (mode.equals("btt")) {

            convertImputedBeagleToTriTyper(in, tmp, ext, out, fam);
            // --mode btt --in indir --tpl template --ext ext --out outdir --fam famfile
        } else if (mode.equals("bttb")) {
            System.out.println("in:\t" + in);
            System.out.println("tpl:\t" + tmp);
            System.out.println("size:\t" + size);
            System.out.println("out:\t" + out);
            convertImputedBeagleToTriTyperFromBatches(in, tmp, out, size, chrStart, chrEnd);
            // --mode bttb --in indir --tpl template --out outdir --size samplesize
        } else if (mode.equals("ttpd")) {
            convertTriTyperToPlinkDosage(in, out, fam, splitbychromosome);
            // --mode ttpd --in indir --beagle beagledir --tpl template --batchdescriptor batchdescriptor --out outdir --fam famfile
        } else if (mode.equals("ttpm")) {
            convertTriTyperToPedAndMap(in, out, fam, snps, splitbychromosome);

            // --mode ttpm --in indir --out outdir --fam famfile
        } else if (mode.equals("ftt")) {
            convertFinalReportToTriTyper(in, out);
            // --mode ttpm --in indir --out outdir --fam famfile
        } else if (mode.equals("ttpmh")) {
            convertTriTyperToPedAndMapHapMap(in, hap, out, fam, batch, chr, exclude, false);
            // --mode ttpmh --in indir --hap hapmapdir --out outdir --fam famfile --batch-file batchfile --chr chromosome
        } else if (mode.equals("corr")) {
            System.out.println("in:\t" + in);
            System.out.println("inname:\t" + inName);
            System.out.println("in2:\t" + in2);
            System.out.println("in2name:\t" + in2Name);
            System.out.println("out:\t" + out);

            correlateDatasets(in, inName, in2, in2Name, out, snps);
            // --mode corr --in indir --inname name --in2 indir2 --in2name inname2 --out outdir --snps snplist
        } else if (mode.equals("pmtt")) {
            convertPedAndMapToTriTyper(in, out);
            // --mode pmtt --in indir --out outdir
        } else if (mode.equals("ecra")) {
            getExcludedSNPsWithCallRate(in, threshold);
            // --mode ecra --in indir --threshold threshold
        } else if (mode.equals("r2dist")) {
            createR2Distributions(in, tmp, out, size);
            // --mode r2dist --in indir --tpl template --out outdir --size numbatches
        } else if (mode.equals("batch")) {
            createBatches(in, out, size);
            // --mode batch --in indir --tpl template --out outdir --size batchsize
        } else if (mode.equals("corrb")) {
            System.out.println("in:\t" + in);
            System.out.println("inname:\t" + inName);
            System.out.println("in2:\t" + in2);
            System.out.println("in2name:\t" + in2Name);
            System.out.println("out:\t" + out);

            correlateDatasets(in, inName, in2, in2Name, out, snps, beagle, tmp, size);
        } else if (mode.equals("merge")) {
            System.out.println("in");
            System.out.println("in2");
            System.out.println("out");
            concatDatasets(in, in2, out, snps);
        } else if (mode.equals("concat")) {
            System.out.println("in:\t" + in);
            System.out.println("out:\t" + out);
            concatDatasets(in, out);
        } else if (mode.equals("itt")) {
            System.out.println("in:\t" + in);
            System.out.println("out:\t" + out);
            System.out.println("nrSamples:\t" + nrSamples);

            convertImputeToTriTyper(in, out, nrSamples, sampleFile, sampleFileToInclude, fileMatchRegex, snps);
        } else if (mode.equals("ttvcf")) {
            System.out.println("in:\t" + in);
            System.out.println("out:\t" + out);
            convertTriTyperToVcf(in, out);

        } else if (mode.equals("vcftt")) {
            System.out.println("in:\t" + in);
            System.out.println("out:\t" + out);
            System.out.println("pattern:\t" + pattern);
            convertVCFToTriTyper(in, out, pattern);

        } else if (mode.equals("mmtt")) {
            System.out.println("in:\t" + in);
            System.out.println("out:\t" + out);
            convertMinimacToTriTryper(in, out);

        } else if (mode.equals("mtt")) {
            System.out.println("in:\t" + in);
            System.out.println("out:\t" + out);
            convertMachToTriTyper(in, out);

        } else if (mode.equals("mmtt2")) {
            System.out.println("in:\t" + in);
            System.out.println("out:\t" + out);
            convertMinimacToTriTryper2(in, out);

        } else {
            printUsage();
        }
//        else if (mode.equals("pmbg")) {
//            System.out.println("in:\t" + in);
//            System.out.println("batch:\t" + batch);
//            convertPedAndMapToBeagle(in, batch);
//        } 
    }

    public void printUsage() {

        System.out.println("------------------------\nPreProcessing\n------------------------\n");
        System.out.println("# Create random batches of cases and controls from a TriTyper dataset. Creates a file called batches.txt in outdir.\n"
                + "--mode batch --in TriTyperdir --out outdir --size batchsize\n");

        System.out.println("------------------------\nIllumina FinalReport files\n------------------------\n");
        System.out.println("# Convert a dir with MACH imputed data into TriTyper\n"
                + "--mode ftt --in finalreportfile --out TriTyperDir");

        System.out.println("------------------------\nImputation\n------------------------\n");
        System.out.println("# Convert Impute Imputed data into TriTyper\n"
                + "--mode itt --in ImputeDir --out TriTyperDir --nrSamples numberOfSamplesInImputedData [--samples samplelistfile.txt] [--samplestoinclude samplelistfiletoinclude.txt] [--fileMatchRegex pattern] [--snps snpfile.txt]");

        System.out.println("# Convert a dir with MACH imputed data into TriTyper\n"
                + "--mode mtt --in Imputation restult dir --out TriTyperDir");

        System.out.println("# Convert a dir with Minimac imputed data into TriTyper\n"
                + "--mode mmtt --in Imputation restult dir --out TriTyperDir");

        System.out.println("# Convert a single Minimac imputed file into TriTyper\n"
                + "--mode mmtt2 --in Imputation restult file prefix (.info and .dose should be present) --out TriTyperDir");

        System.out.println();

        System.out.println("------------------------\nBeagle\n------------------------\n");
        System.out.println("# Convert beagle files (one file/chromosome) to TriTyper. Filetemplate is a template for the batch filenames, The text CHROMOSOME will be replaced by the chromosome number.\n"
                + "--mode btt --in BeagleDir --tpl template --ext ext --out TriTyperDir [--fam famfile]\n");

        System.out.println("# Convert batches of beagle files (multiple files / chromosome) to trityper files. Filetemplate is a template for the batch filenames, The text CHROMOSOME will be replaced by the chromosome number, BATCH by the batchname.\n"
                + "--mode bttb --in BeagleDirdir --tpl template --out TriTyperDir --size numbatches [--chr chromosome] [--chrstart startchr] [--chrend endchr]\n");

        System.out.println("------------------------\nPed+Map (Plink files)\n------------------------\n");

        System.out.println("# Converts Ped and Map files created by ttpmh to Beagle format\n"
                + "--mode pmbg --in indir --batch-file batches.txt\n");

        System.out.println("# Converts TriTyper file to Plink Dosage format. Filetemplate is a template for the batch filenames, The text CHROMOSOME will be replaced by the chromosome number, BATCH by the batchname.\n"
                + "--mode ttpd --in indir --out outdir --fam famfile [--split]\n");

        System.out.println("# Converts PED and MAP files to TriTyper.\n"
                + "--mode pmtt --in Ped+MapDir --out TriTyperDir\n");

        System.out.println("# Converts TriTyper file to PED and MAP files. The FAM file is optional. --split splits the ped and map files per chromosome. Providing a SNP file will override --split\n"
                + "--mode ttpm --in indir --out outdir [--fam famfile] [--split] [--snps snpfile]\n");

        System.out.println("# Converts TriTyper dataset to Ped+Map concordant to reference (hap) dataset. Supply a batchfile if you want to export in batches. Supply a chromosome if you want to export a certain chromosome.\n"
                + "--mode ttpmh --in TriTyperDir --hap TriTyperReferenceDir --out outdir [--fam famfile] [--batch-file batchfile] [--chr chromosome] [--exclude fileName]\n");

        System.out.println("---------------------\nPostProcessing\n---------------------\n");
        System.out.println("# Correlates genotypes of imputed vs non-imputed datasets. Saves a file called correlationOutput.txt in outdir, containing correlation per chromosome as well as correlation distribution.\n"
                + "--mode corr --in TriTyperDir --name datasetname --in2 TriTyperDir2 --name2 datasetname2 --out outdir [--snps snplist]\n");

        System.out.println("# Correlates genotypes of imputed vs non-imputed datasets. Also take Beagle imputation score (R2) into account. Saves a file called correlationOutput.txt in outdir, containing correlation per chromosome as well as correlation distribution.\n"
                + "--mode corrb --in TriTyperDir --name datasetname --in2 TriTyperDir2 --name2 datasetname2 --out outdir --beagle beagleDir --tpl template --size numBatches \n");

        System.out.println("# Gets all the excluded snps from chrx.excludedsnps.txt with a certain call-rate threshold (0 < threshold < 1.0)\n"
                + "--mode ecra --in TriTyperDir --threshold threshold\n");

        System.out.println("# Generates R2 distribution (beagle quality score) for each batch and chromosome, and tests each batch against chromosome R2 distribution, using WilcoxonMannWhitney test\n"
                + "--mode r2dist --in BeagleDir --template template --out outdir --size numbatches\n");

        System.out.println("# Merge two TriTyper datasets\n"
                + "--mode merge --in TriTyper1Dir --in2 TryTyper2Dir --out outdir\n");

        System.out.println("# Concatinate TriTyper datasets. All should contain same number of individuals in the same order.\n"
                + "--mode concat --in \"TriTyperDir1;TriTyperDir2;TriTyperDirN\" --out outdir\n");

        System.out.println("---------------------\nVCF\n---------------------\n");
        System.out.println("# Convert TriTyper to VCF and vice versa\n"
                + "--mode ttvcf --in indir --out outdir \n"
                + "--mode vcftt --in indir --out outdir [--pattern pattern]\n");

        System.exit(0);
    }

    private void convertTriTyperToPedAndMapHapMap(String baseDir, String hapmapDir, String outputDir, String famfile, String batchFile, Integer chromosome, String excludesnps, boolean singlefile) throws IOException {
        if (baseDir == null || hapmapDir == null || outputDir == null) {
            System.out.println("Please supply values for TriTyper inputdir (--in), HapMap dir (--hap) and outputdir (--out) when to running --mode ttpmh\n\n");
            // printUsage();
            System.exit(0);
        }

        TriTyperReferenceConcordantPedAndMapExporter ttb = new TriTyperReferenceConcordantPedAndMapExporter();
        ttb.export(hapmapDir, baseDir, batchFile, excludesnps, outputDir);
    }

    private void convertTriTyperToPedAndMap(String baseDir, String outputDir, String famfile, String snpList, boolean splitbychromosome) throws IOException, Exception {
        if (baseDir == null || outputDir == null) {
            System.out.println("Please supply values for TriTyper inputdir (--in) and outputdir (--out) when to running --mode ttpm\n\n");
            // printUsage();
            System.exit(0);
        }
        TriTyperToPedAndMapConverter ttpm = new TriTyperToPedAndMapConverter();
        if (snpList != null) {
            ttpm.exportSubsetOfSNPs(baseDir, outputDir, snpList, null);
        } else {
            ttpm.exportAllSNPs(baseDir, outputDir, splitbychromosome);
        }
    }

    private void convertImputedBeagleToTriTyper(String inputLocation, String descriptor, String extension, String trityperoutputdir, String famfile) throws IOException {
        if (inputLocation == null || descriptor == null || extension == null || trityperoutputdir == null) {
            System.out.println("Please supply values for Beagle inputdir (--beagle), file template (--tpl), batch descriptor (--batchdesc) and TriTyper outputdir (--out) when to running --mode btt\n\n");
            // printUsage();
            System.exit(0);
        }

        BeagleImputedToTriTyper btt = new BeagleImputedToTriTyper();
        if (famfile != null) {
            btt.loadFamFile(famfile);
        }

        btt.importImputedDataWithDosageInformationBeagle(inputLocation, descriptor, extension, trityperoutputdir);
    }

    private void convertImputedBeagleToTriTyperFromBatches(String inputLocation, String descriptor, String outputLocation, Integer numBatches, Integer chrStart, Integer chrEnd) throws IOException {
        if (inputLocation == null || descriptor == null || outputLocation == null || numBatches == null) {
            System.out.println("Please supply values for Beagle inputdir (--beagle), file template (--tpl), number of batches (--size) and TriTyper outputdir (--out) when to running --mode bttb\n\n");
            System.exit(0);
            // printUsage();
        }

        BeagleImputedToTriTyper btt = new BeagleImputedToTriTyper();
        btt.importImputedDataWithDosageInformationBeagleBatches(inputLocation, descriptor, numBatches, outputLocation, chrStart, chrEnd);
    }

    private void correlateDatasets(String dataset1, String dataset1Name, String dataset2, String dataset2Name, String output, String snpList) throws IOException {
        if (dataset1 == null || dataset1Name == null || dataset2 == null || dataset2Name == null || output == null) {
            System.out.println("Please supply values for dataset location (--in), datasetname (--name), dataset 2 location (--in2), dataset 2 name (--name2), and output location (--out)\n\n");
            System.exit(0);
        }

        TriTyperDatasetCorrelator c = new TriTyperDatasetCorrelator(dataset1, dataset1Name, dataset2, dataset2Name);
        if (snpList != null) {
            c.confineToSNPs(snpList);
        }
        c.run(output);
    }

    private void getExcludedSNPsWithCallRate(String location, Double threshold) throws IOException {
        if (location == null || threshold == null) {
            System.out.println("Please supply values for input location (--in) and threshold (--threshold)\n\n");
            System.exit(0);
        }

        TriTypertoPedMapDatExcludedSNPAnalyzer q = new TriTypertoPedMapDatExcludedSNPAnalyzer();
        q.printSNPsWithCallRateHigherThan(threshold, location);
    }

    private void convertPedAndMapToTriTyper(String dataLocation, String outputLocation) {
        if (dataLocation == null || outputLocation == null) {
            System.out.println("Please supply values for input location (--in) and output (--out)");
            System.exit(0);
        }
        PedAndMapToTriTyper p = new PedAndMapToTriTyper();
        try {
            p.importPEDFile(dataLocation, outputLocation);
        } catch (IOException e) {
            e.printStackTrace();

        }
    }

//    private void convertPedAndMapToTriTyper(String dataLocation, String outputLocation, String mapdelimiter, int chrcol, int chrpos, int snpcol, String peddelimiter) {
//	PedAndMapToTriTyper p = new PedAndMapToTriTyper();
//	String[] casesToInclude = {"2"};
//	p.importPEDFile(dataLocation, mapdelimiter, chrcol, chrpos, snpcol, peddelimiter, outputLocation, casesToInclude);
//    }
    // ttpd trityperdir beagledir datasetdescriptor batchdescriptor outputdir
    private void convertTriTyperToPlinkDosage(String trityperdir, String outputdir, String famFile, boolean splitperchr) throws IOException {
        // if(trityperdir == null || beagledir == null || datasetdescriptor == null || batchdescriptor == null || outputdir == null){
        if (trityperdir == null || outputdir == null) {
            System.out.println("Please supply values for TriTyper input (--in), beagle dir (--beagle), dataset descriptor (--tpl), batch descriptor (--batchdesc), and output dir (--out)\n\n");
            System.exit(0);
        }
        TriTyperToPlinkDosage p = new TriTyperToPlinkDosage();

        p.outputDosageInformation(trityperdir, outputdir, famFile, splitperchr);
        // String inputFile, String datasetDir, String datasetDescription, String batchname,  String outputDir
        // p.splitDosageInformationPerChromosome(outputdir, beagledir, datasetdescriptor, batchdescriptor, outputdir);
    }

    private void createR2Distributions(String inputLocation, String descriptor, String outputLocation, Integer numBatches) throws IOException {
        if (inputLocation == null || descriptor == null || outputLocation == null || numBatches == null) {
            System.out.println("Please supply values for Beagle input dir (--in), dataset template (--tpl), outputlocation (--out) and num batches (--size)\n\n");
            System.exit(0);
        }
        BeagleImputationQuality b = new BeagleImputationQuality();
        b.determineImputationQualityDistribution(inputLocation, descriptor, numBatches, outputLocation);
    }

    private void createBatches(String in, String out, Integer size) throws IOException {
        if (in == null || out == null || size == null) {
            System.out.println("Please supply values for TriTyper input directory (--in), output directory (--out) and requested batch size (--size)\n\n");
            System.exit(0);
        }
        BatchGenerator b = new BatchGenerator();
        b.generateBatchesFromTriTyper(in, out, size);

    }

    private void correlateDatasets(String dataset1, String dataset1Name, String dataset2, String dataset2Name, String output, String snpList, String beagle, String tmp, Integer size) throws IOException {
        if (dataset1 == null || dataset1Name == null || dataset2 == null || dataset2Name == null || output == null) {
            System.out.println("Please supply values for dataset location (--in), datasetname (--name), dataset 2 location (--in2), dataset 2 name (--name2), and output location (--out)\n\n");
            System.exit(0);
        }

        TriTyperDatasetCorrelator c = new TriTyperDatasetCorrelator(dataset1, dataset1Name, dataset2, dataset2Name, beagle, tmp, size);
        if (snpList != null) {
            c.confineToSNPs(snpList);
        }
        c.run(output);
    }

    private void concatDatasets(String in, String in2, String out, String snps) throws IOException {

        if (in != null && in2 != null && out != null) {
            TriTyperGenotypeDataMerger merger = new TriTyperGenotypeDataMerger();
            merger.combinePrioritizerDatasetsMergeCommonSNPs(in, in2, out, snps);
        } else {
            System.out.println("Please provide: --in, --in2 and --out for --mode merge");
            System.exit(0);
        }
    }

    private void convertImputeToTriTyper(String in, String out, Integer nrSamples, String samplesFile, String samplesToIncludeFile, String fileMatchRegex, String snps) throws IOException, Exception {
        if (in == null || out == null || nrSamples == null) {
            System.out.println("Please provide: --in, --nrSamples and --out for --mode itt");
            System.exit(0);
        } else {
            ImputeImputedToTriTyperV2 t = new ImputeImputedToTriTyperV2();
            t.importImputedDataWithProbabilityInformationImpute(in, out, nrSamples, samplesFile, samplesToIncludeFile, fileMatchRegex, snps);
        }

    }

    private void convertTriTyperToVcf(String in, String out) throws IOException, Exception {
        if (in == null || out == null) {
            System.out.println("Please provide: --in and --out for --mode ttvcf");
            System.exit(0);
        } else {
            System.out.println("Warning! All genotypes are exported as if phased");
            TriTyperToVCF ttvcfConvertor = new TriTyperToVCF();
            ttvcfConvertor.convert(in, out, null);
        }

    }

    private void convertVCFToTriTyper(String in, String out, String pattern) throws IOException, Exception {
        if (in == null || out == null) {
            System.out.println("Please provide: --in and --out for --mode vcftt");
            System.exit(0);
        } else {

            VCFToTriTyper tt = new VCFToTriTyper();
            tt.parse(in, out, pattern);
        }

    }

    private void convertFinalReportToTriTyper(String in, String out) throws IOException {
        boolean isIlluminaFinalReportFile = true;
        String delimiter = "\t";
        String decimalSeparator = ".";
//        for (int a=0; a<args.length; a++) {
//            String argument = args[a].trim();
//            if (argument.startsWith("-commadelimited")) delimiter = ",";
//            if (argument.startsWith("-tabdelimited")) delimiter = "\t";
//            if (argument.startsWith("-spacedelimited")) delimiter = " ";
//            if (argument.startsWith("-semicolondelimited")) delimiter = ";";
//            if (argument.startsWith("-finalreportformat")) isIlluminaFinalReportFile = true;
//            if (argument.startsWith("-decimalseparatoriscomma")) decimalSeparator = ",";
//        }
        String inputFile = in.trim();
        String outputDir = out.trim();
        if (!outputDir.endsWith("/")) {
            outputDir += "/";
        }

        new FinalReportToTriTyper(inputFile, outputDir, isIlluminaFinalReportFile, delimiter, decimalSeparator);
    }
//    private void convertPedAndMapToBeagle(String in, String batch) {
//        PedAndMapToBeagle pmbg = new PedAndMapToBeagle(in, batch);
//    }

    private void convertMinimacToTriTryper(String in, String out) throws IOException {
        if (in == null || out == null) {
            System.out.println("Please provide: --in and --out for --mode mmtt");
            System.exit(-1);
        } else {
            MinimacImputedToTriTyper.convertMinimacImputedToTriTyper(in, out);
        }
    }

    private void convertMinimacToTriTryper2(String in, String out) throws IOException {
        if (in == null || out == null) {
            System.out.println("Please provide: --in and --out for --mode mmtt2");
            System.exit(-1);
        } else {
            MinimacImputedToTriTyper.convertSingleFileMinimacImputedToTriTyper(in, out);
        }
    }

    private void concatDatasets(String in, String out) throws IOException, FileNotFoundException, Exception {
        if (in == null || out == null) {
            System.out.println("Please provide: --in and --out for --mode concat");
            System.exit(-1);
        } else {

            TriTyperConcatDatasets triTyperConcatDatasets = new TriTyperConcatDatasets(out);

            String[] inputDatasets = SEMI_COLON_PATTERN.split(in);

            for (String inputDataset : inputDatasets) {
                if (inputDataset.length() == 0) {
                    //allow ; at begin of end of the string
                    continue;
                }
                triTyperConcatDatasets.addDatasetToConcat(inputDataset);
            }

            triTyperConcatDatasets.writeConcatedDataset();

        }
    }

    private void convertMachToTriTyper(String in, String out) throws IOException {
        if (in == null || out == null) {
            System.out.println("Please provide: --in and --out for --mode mmtt");
            System.exit(-1);
        } else {
            MachImputedToTriTyper m = new MachImputedToTriTyper(in, out);
        }
    }
}
