package umcg.genetica.io.trityper;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.util.*;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import umcg.genetica.io.Gpio;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.util.RankArray;

/**
 *
 * @author Lude, Marc Jan
 */
public class ConvertDoubleMatrixDataToTriTyper {

    private static Pattern SPLIT_ON_TAB = Pattern.compile("\t");

    public static void main(String[] args) throws IOException {

        String mappingFile = args[0];
        String dataMatrix = args[1];
        String outputFolder = args[2];
        String option1 = args[3];

        if (!(new File(outputFolder).exists())) {
            Gpio.createDir(outputFolder);
        }

        HashSet<String> hashCpGSites = new HashSet<String>();
        try {
            System.out.println("Writing SNPMappings.txt to file:");
            int nrSites = 0;
            java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(mappingFile)));
            java.io.BufferedWriter outSNPMappings = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputFolder + "/SNPMappings.txt")));
            String str;
            in.readLine();
            while ((str = in.readLine()) != null) {
                String[] data = SPLIT_ON_TAB.split(str);
                hashCpGSites.add(data[1]);
                outSNPMappings.write(data[3] + '\t' + data[4] + '\t' + data[1] + '\n');
                nrSites++;
            }
            System.out.println("Number of sites in mapping file:\t" + nrSites);
            outSNPMappings.close();
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }

        DoubleMatrixDataset<String, String> dataset = null;
        try {
            dataset = DoubleMatrixDataset.loadSubsetOfTextDoubleData(dataMatrix, "\t", hashCpGSites, null);
        } catch (IOException ex) {
            Logger.getLogger(ConvertDoubleMatrixDataToTriTyper.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(0);
        }
        if (dataset != null && !dataset.getHashCols().isEmpty() && !dataset.getHashRows().isEmpty()) {
//            if(args.length > 3 && args[3].equals("scale")){
//                ConvertBetaAndMvalues.transformMToBetavalue(dataset.getMatrix());
//            }
            if (args[3].equals("rank")) {
                rankRows(dataset.getMatrix());
            }
            rescaleValue(dataset.getMatrix(), null);

            try {
                System.out.println("\nWriting SNPs.txt to file:");
                BufferedWriter out = new BufferedWriter(new FileWriter(outputFolder + "/SNPs.txt"));
                for (String snp : dataset.getRowObjects()) {
                    out.write(snp + "\n");
                }
                out.close();
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }

            try {
                System.out.println("\nWriting Individuals.txt and Phenotype.txt to file:");
                BufferedWriter outIndNew = new BufferedWriter(new FileWriter(outputFolder + "/Individuals.txt"));
                BufferedWriter outPhenoNew = new BufferedWriter(new FileWriter(outputFolder + "/PhenotypeInformation.txt"));
                for (String ind : dataset.getColObjects()) {
                    outIndNew.write(ind + '\n');
                    outPhenoNew.write(ind + "\tcontrol\tinclude\tmale\n");
                }
                outIndNew.close();
                outPhenoNew.close();
            } catch (Exception e) {
                System.out.println("Error:\t" + e.getMessage());
                e.printStackTrace();
            }

            int nrSNPs = dataset.rows();
            int nrSamples = dataset.columns();

            WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, new File(outputFolder + "/GenotypeMatrix.dat"), false);
            WGAFileMatrixImputedDosage fileMatrixDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamples, new File(outputFolder + "/ImputedDosageMatrix.dat"), false);
            byte[] alleles = new byte[2];
            alleles[0] = 84;
            alleles[1] = 67;

            for (int snp = 0; snp < nrSNPs; snp++) {
                DoubleMatrix1D snpRow = dataset.getMatrix().viewRow(snp);
                byte[] allele1 = new byte[nrSamples];
                byte[] allele2 = new byte[nrSamples];
                byte[] dosageValues = new byte[nrSamples];
                for (int ind = 0; ind < nrSamples; ind++) {
                    if (snpRow.get(ind) > 100) {
                        allele1[ind] = alleles[1];
                        allele2[ind] = alleles[1];
                    } else {
                        allele1[ind] = alleles[0];
                        allele2[ind] = alleles[0];
                    }

                    int dosageInt = (int) Math.round(snpRow.get(ind));
                    byte value = (byte) (Byte.MIN_VALUE + dosageInt);
                    dosageValues[ind] = value;
                }
                fileMatrixGenotype.setAllele1(snp, 0, allele1);
                fileMatrixGenotype.setAllele2(snp, 0, allele2);
                fileMatrixDosage.setDosage(snp, 0, dosageValues);
            }
            fileMatrixGenotype.close();
            fileMatrixDosage.close();
        }
        System.out.println("Finished.");
    }

    public static void rescaleValue(DoubleMatrix2D matrix, Double multiplier) {
        if (multiplier != null) {
            for (int p = 0; p < matrix.rows(); p++) {
                double min = matrix.viewRow(p).getMinLocation()[0];
                double denominator = (matrix.viewRow(p).getMaxLocation()[0] - min) * (1 / multiplier);
                for (int s = 0; s < matrix.columns(); s++) {
                    matrix.setQuick(p, s, ((matrix.getQuick(p, s) - min) / denominator));
                }
            }
        } else {
            for (int p = 0; p < matrix.rows(); p++) {
                double min = matrix.viewRow(p).getMinLocation()[0];
                double denominator = matrix.viewRow(p).getMaxLocation()[0] - min;
                for (int s = 0; s < matrix.columns(); s++) {
                    matrix.setQuick(p, s, ((matrix.getQuick(p, s) - min) / denominator));
                }
            }
        }

    }

    private static void rankRows(DoubleMatrix2D matrix) {
        RankArray rda = new RankArray();
        for (int p = 0; p < matrix.rows(); p++) {
            double[] rankedValues = rda.rank(matrix.viewRow(p).toArray(), false);
            for (int s = 0; s < matrix.columns(); s++) {
                matrix.setQuick(p, s, rankedValues[s]);
            }
        }
    }

}
