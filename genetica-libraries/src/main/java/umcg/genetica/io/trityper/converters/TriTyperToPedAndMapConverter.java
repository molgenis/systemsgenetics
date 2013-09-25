/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.converters;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class TriTyperToPedAndMapConverter {

    public void exportSubsetOfSNPs(String datadir, String outputDir, String snpSubsetFile, String indsToExport) throws IOException, Exception {
        TriTyperGenotypeData dataGenotypeDataset = new TriTyperGenotypeData();
        dataGenotypeDataset.load(datadir);
        String[] snps = dataGenotypeDataset.getSNPs();
        int numsnps = snps.length;

        HashSet<String> hashInds = null;
        if (indsToExport != null) {
            TextFile infile = new TextFile(indsToExport, TextFile.R);
            hashInds = (HashSet<String>) infile.readAsSet(0, TextFile.tab);
            infile.close();
        }


        if (outputDir.endsWith("/")) {
            outputDir += "/";
        }
        if (!Gpio.isDir(outputDir)) {
            Gpio.createDir(outputDir);
        }
        TextFile log = new TextFile(outputDir + "log.txt", TextFile.W);
        TextFile tf = new TextFile(snpSubsetFile, TextFile.R);
        ArrayList<String> requestedSNPs = tf.readAsArrayList();
        tf.close();

        ArrayList<Integer> snpids = new ArrayList<Integer>();

        HashSet<String> availableSNPs = new HashSet<String>(snps.length);
        availableSNPs.addAll(Arrays.asList(snps));

        for (String s : requestedSNPs) {
            if (!availableSNPs.contains(s)) {
                System.out.println(s + "\tnot detected in dataset!");
                log.writeln(s + "\tnot detected in dataset!");
            }
        }

        HashSet<String> requestedSNPHash = new HashSet<String>();
        requestedSNPHash.addAll(requestedSNPs);
        System.out.println(requestedSNPHash.size() + " unique snps loaded from file");

        for (int snpid = 0; snpid < snps.length; snpid++) {
            if (requestedSNPHash.contains(snps[snpid])) {
                snpids.add(snpid);
            }
        }
        System.out.println(snpids.size() + " snps detected from " + snpSubsetFile);

        boolean[] snppassesQC = new boolean[snpids.size()];
        for (int i = 0; i < snppassesQC.length; i++) {
            snppassesQC[i] = true;
        }

        TextFile pedFile = new TextFile(outputDir + "/output.ped", true);
        exportSetOfSNPsToPedFile(snppassesQC, snpids, dataGenotypeDataset, pedFile, log, hashInds);
        pedFile.close();

        TextFile mapfile = new TextFile(outputDir + "/output.map", TextFile.W);
        exportSetOfSNPsToMapFile(snppassesQC, snpids, dataGenotypeDataset, mapfile);
        snppassesQC = null;
        mapfile.close();
        exportFamFile(dataGenotypeDataset, outputDir);
        log.close();
    }

    public void exportAllSNPs(String datadir, String outputDir, boolean splitByChromosome) throws IOException {
        TriTyperGenotypeData dataGenotypeDataset = new TriTyperGenotypeData();
        dataGenotypeDataset.load(datadir);
        String[] snps = dataGenotypeDataset.getSNPs();
        int numsnps = snps.length;

        if (outputDir.endsWith("/")) {
            outputDir += "/";
        }
        if (!Gpio.isDir(outputDir)) {
            Gpio.createDir(outputDir);
        }

        TextFile log = new TextFile(outputDir + "log.txt", TextFile.W);

        if (splitByChromosome) {
            for (int chr = 1; chr < 23; chr++) {
                ArrayList<Integer> snpids = new ArrayList<Integer>();
                for (int i = 0; i < snps.length; i++) {
                    if (dataGenotypeDataset.getChr(i) == chr) {
                        snpids.add(i);
                    }
                }

                boolean[] snppassesQC = new boolean[snpids.size()];
                for (int i = 0; i < snppassesQC.length; i++) {
                    snppassesQC[i] = true;
                }

                TextFile pedFile = new TextFile(outputDir + "/output." + chr + ".ped", true);
                exportSetOfSNPsToPedFile(snppassesQC, snpids, dataGenotypeDataset, pedFile, log, null);
                pedFile.close();


                TextFile datFile = new TextFile(outputDir + "/output." + chr + ".dat", true);
                exportSetOfSNPsToDatFile(snppassesQC, snpids, dataGenotypeDataset, datFile);
                datFile.close();

                TextFile mapfile = new TextFile(outputDir + "/output." + chr + ".map", TextFile.W);
                exportSetOfSNPsToMapFile(snppassesQC, snpids, dataGenotypeDataset, mapfile);
                snppassesQC = null;
                mapfile.close();
            }
        } else {
            ArrayList<Integer> snpids = new ArrayList<Integer>();
            for (int i = 0; i < snps.length; i++) {
                snpids.add(i);
            }
            boolean[] snppassesQC = new boolean[snpids.size()];
            for (int i = 0; i < snppassesQC.length; i++) {
                snppassesQC[i] = true;
            }

            TextFile pedFile = new TextFile(outputDir + "/output.ped", true);
            exportSetOfSNPsToPedFile(snppassesQC, snpids, dataGenotypeDataset, pedFile, log, null);
            pedFile.close();


            TextFile datFile = new TextFile(outputDir + "/output.dat", true);
            exportSetOfSNPsToDatFile(snppassesQC, snpids, dataGenotypeDataset, datFile);
            datFile.close();

            TextFile mapfile = new TextFile(outputDir + "/output.map", TextFile.W);
            exportSetOfSNPsToMapFile(snppassesQC, snpids, dataGenotypeDataset, mapfile);
            snppassesQC = null;
            mapfile.close();
        }


        exportFamFile(dataGenotypeDataset, outputDir);
        log.close();
        // write snp data
    }

    private void exportSetOfSNPsToDatFile(boolean[] snppassesQC, ArrayList<Integer> snpids, TriTyperGenotypeData dataGenotypeDataset, TextFile datfile) throws IOException {
        int numsnps = snpids.size();
        datfile.writeln("A Status");
        for (int s = 0; s < numsnps; s++) {
            if (snppassesQC[s]) {
                int snpid = snpids.get(s);
                SNP snpObj = dataGenotypeDataset.getSNPObject(snpid);

                datfile.writeln("M " + snpObj.getName());
            }
        }
    }

    private void exportSetOfSNPsToMapFile(boolean[] snppassesQC, ArrayList<Integer> snpids, TriTyperGenotypeData dataGenotypeDataset, TextFile mapfile) throws IOException {
        int numsnps = snpids.size();
        for (int s = 0; s < numsnps; s++) {
            if (snppassesQC[s]) {
                int snpid = snpids.get(s);
                SNP snpObj = dataGenotypeDataset.getSNPObject(snpid);

                mapfile.writeln(ChrAnnotation.parseByte(snpObj.getChr()) + " " + snpObj.getName() + " " + 0 + " " + snpObj.getChrPos());
            }
        }
    }

    public void exportFamFile(TriTyperGenotypeData ds, String outdir) throws IOException {
        TextFile famout = new TextFile(outdir + "output.fam", TextFile.W);
        String[] individuals = ds.getIndividuals();

        Boolean[] indIsIncluded = ds.getIsIncluded();
        Boolean[] indIsFemale = ds.getIsFemale();
        Boolean[] indIsCase = ds.getIsCase();

        for (int i = 0; i < individuals.length; i++) {
            if (indIsIncluded[i] != null && indIsIncluded[i]) {
                System.out.println(individuals[i]);
                String sex = "-9";
                String condition = "-9";

                if (indIsFemale[i] == null) {
                } else if (indIsFemale[i]) {
                    sex = "2";
                } else {
                    sex = "1";
                }
                if (indIsCase[i] == null) {
                } else if (indIsCase[i]) {
                    condition = "2";
                } else {
                    condition = "1";
                }

                StringBuilder sb = new StringBuilder();
                sb.append("1");
                sb.append(" ");
                sb.append(individuals[i]);

                sb.append(" ");
                sb.append("0");
                sb.append(" ");
                sb.append("0");
                sb.append(" ");
                sb.append(sex);
                sb.append(" ");
                sb.append(condition);

                sb.append("\n");
                famout.write(sb.toString());
            }
        }

        famout.close();
    }

    private void exportSetOfSNPsToPedFile(boolean[] snppassesQC, ArrayList<Integer> snpids, TriTyperGenotypeData dataGenotypeDataset, TextFile pedFile, TextFile log, HashSet<String> indsToExport) throws IOException {
        String[] individuals = dataGenotypeDataset.getIndividuals();
        int numDatasetIndividuals = individuals.length;
        SNPLoader loader = dataGenotypeDataset.createSNPLoader();
        int batchsize = 1000;

        double requiredmemorypersample = (double) (2 * snpids.size()) / 1048576; // in bytes, per sample

        double freeMem = (double) Runtime.getRuntime().totalMemory() / 1048576;

        System.out.println("Free memory: " + freeMem + " Mb, taking 50%: " + (freeMem * 0.50) + " Mb");
        freeMem *= 0.50;
        System.out.println("Req: " + requiredmemorypersample + " mb/sample");

        batchsize = (int) Math.floor(freeMem / requiredmemorypersample);


        if (batchsize > numDatasetIndividuals) {
            batchsize = numDatasetIndividuals;
        }

        System.out.println("Batch size: " + batchsize);

        int included = 0;

        int numsnps = snpids.size();
        int indsprocessed = 0;
        while (indsprocessed < numDatasetIndividuals) {
            if (indsprocessed + batchsize > numDatasetIndividuals) {
                batchsize = numDatasetIndividuals - indsprocessed;
            }
            byte[][] indSNPData = new byte[batchsize][snpids.size() * 2];

            System.out.println("");
            ProgressBar pb = new ProgressBar(snpids.size());
            int pos = 0;

            for (int s = 0; s < numsnps; s++) {
                if (snppassesQC[s]) {
                    int snpid = snpids.get(s);
                    SNP snpObj = dataGenotypeDataset.getSNPObject(snpid);
                    loader.loadGenotypes(snpObj);
                    if (snpObj.getMAF() > 0) {
                        int batchind = 0;

                        byte[] snpAlleles1 = snpObj.getAllele1();
                        byte[] snpAlleles2 = snpObj.getAllele2();

                        for (int i = indsprocessed; i < indsprocessed + batchsize; i++) {
                            byte allele1 = snpAlleles1[i];
                            byte allele2 = snpAlleles2[i];
                            indSNPData[batchind][s] = allele1; // gt1
                            indSNPData[batchind][numsnps + s] = allele2; // gt2
                            batchind++;
                        }
                    } else {
                        snppassesQC[s] = false;
                        log.writeln(snpObj.getName() + "\tChr: " + snpObj.getChr() + "\tChrPos: " + snpObj.getChrPos() + "\tdoes not pass MAF QC: " + snpObj.getMAF() + ", HWE: " + snpObj.getHWEP());
                    }
                    snpObj.clearGenotypes();
                } else {
                    // do nothing.
                }

                pb.iterate();

            }

            pb.close();
            System.out.println("");
            System.out.println("Now writing batch");

            int batchind = 0;
            ProgressBar pb2 = new ProgressBar(batchsize);
            for (int i = indsprocessed; i < indsprocessed + batchsize; i++) {

                if (dataGenotypeDataset.getIsIncluded()[i] == null) {
                    System.err.println("ERROR: " + dataGenotypeDataset.getIndividuals()[i] + " has no inclusion information?");
                } else if (dataGenotypeDataset.getIsIncluded()[i]) {

                    String sex = "-9";
                    String condition = "-9";

                    StringBuilder sb = new StringBuilder();
                    if (dataGenotypeDataset.getIsFemale()[i]) {
                        sex = "2";
                    } else {
                        sex = "1";
                    }

                    if (dataGenotypeDataset.getIsCase()[i]) {
                        condition = "2";
                    } else {
                        condition = "1";
                    }


                    sb.append("1");
                    sb.append(" ");
                    sb.append(individuals[i]);

                    sb.append(" ");
                    sb.append("0");
                    sb.append(" ");
                    sb.append("0");
                    sb.append(" ");
                    sb.append(sex);
                    sb.append(" ");
                    sb.append(condition);


                    byte[] genotypes = indSNPData[batchind];
                    for (int s = 0; s < numsnps; s++) {

                        if (snppassesQC[s]) {
                            byte value1 = genotypes[s];
                            byte value2 = genotypes[s + numsnps];
                            sb.append(" ");
                            String variant1 = BaseAnnot.toString(value1) == null ? "0" : BaseAnnot.toString(value1);
                            sb.append(variant1);
                            sb.append(" ");
                            String variant2 = BaseAnnot.toString(value2) == null ? "0" : BaseAnnot.toString(value2);
                            sb.append(variant2);
                        }
                    }

                    if (indsToExport == null || indsToExport.contains(individuals[i])) {
                        pedFile.writeln(sb.toString());
                    }



                    sb = null;
                    included++;
                }

                batchind++;
                pb2.iterate();
            }


            for (int q = 0; q < indSNPData.length; q++) {
                indSNPData[q] = null;
            }
            indSNPData = null;

            pb2.close();
            System.out.println(indsprocessed + " individuals processed, " + included + " included");
            indsprocessed += batchsize;

        }
        loader.close();
    }
}
