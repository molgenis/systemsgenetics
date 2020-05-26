package umcg.genetica.io.trityper;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.util.*;
import java.io.*;
import java.util.regex.Pattern;
import umcg.genetica.io.Gpio;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.util.RankArray;

/**
 *
 * @author Lude, Marc Jan
 */
public class ConvertHLAHaplotypesToTriTyper {

    private static final Pattern SPLIT_ON_TAB = Pattern.compile("\t");
    private static final Pattern SPLIT_ON_SLASH = Pattern.compile("/");

    public static void main(String[] args) throws IOException {

        String inputHlaHaplotypes = "D:\\UMCG\\Projects\\LLD_Phenotypes\\HLA_Phenotypes.txt";
        String outputFolder = "D:\\UMCG\\Projects\\LL-DeepMGS\\RiskScores\\HaplotypesLLD\\";
        String option1 = "";

        if (!(new File(outputFolder).exists())) {
            Gpio.createDir(outputFolder);
        }
        
        
        DoubleMatrixDataset<String, String> dataset = readHaplotypeFile(inputHlaHaplotypes);

        System.out.println("Faking SNPMappings.txt and SNPs.txt to file:");

        HashSet<String> hashCpGSites = new HashSet<String>();
        try {
            System.out.println("Writing SNPMappings.txt to file:");
            
            java.io.BufferedWriter outSNPMappings = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputFolder + "/SNPMappings.txt")));
            java.io.BufferedWriter SNPs = new java.io.BufferedWriter(new java.io.FileWriter(new File(outputFolder + "/SNPs.txt")));
            
            int pos = 1000;
            for(String haplotype : dataset.getRowObjects()) {
                hashCpGSites.add(haplotype);
                outSNPMappings.write("6\t" + pos+ '\t' + haplotype + '\n');
                SNPs.write(haplotype + '\n');
                pos+=1000;
            }
            
           
            outSNPMappings.close();
            SNPs.close();
        } catch (Exception e) {
            System.out.println("Error:\t" + e.getMessage());
            e.printStackTrace();
        }

        
        if (dataset != null && !dataset.getHashCols().isEmpty() && !dataset.getHashRows().isEmpty()) {
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
                    if (snpRow.get(ind) < 100) {
                        allele1[ind] = alleles[1];
                        allele2[ind] = alleles[1];
                    } else if(snpRow.get(ind) == 100){
                        allele1[ind] = alleles[1];
                        allele2[ind] = alleles[0];
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

    public static DoubleMatrix2D rescaleValue(DoubleMatrix2D matrix, Double multiplier) {
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
        return matrix;
    }

    private static DoubleMatrix2D rankRows(DoubleMatrix2D matrix) {
        RankArray rda = new RankArray();
        for (int p = 0; p < matrix.rows(); p++) {
            double[] rankedValues = rda.rank(matrix.viewRow(p).toArray(), false);
            for (int s = 0; s < matrix.columns(); s++) {
                matrix.setQuick(p, s, rankedValues[s]);
            }
        }
        return matrix;
    }

    private static DoubleMatrixDataset<String, String> readHaplotypeFile(String inputHlaHaplotypes) throws IOException {
        java.io.BufferedReader in = new java.io.BufferedReader(new java.io.FileReader(new File(inputHlaHaplotypes)));
        String str;
        
        HashSet<String> uniqueHaploTypeCombinations = new HashSet<>();
        
        int sampleCounter = 0;
        while ((str = in.readLine()) != null) {
            String[] parts = SPLIT_ON_TAB.split(str);
            uniqueHaploTypeCombinations.add(parts[1]);
            sampleCounter++;
        }
        in.close();
        
        System.out.println("Number of unique Haplotype combinations: "+uniqueHaploTypeCombinations.size());
        System.out.println("Number of samples: "+sampleCounter);
        
        DoubleMatrixDataset<String, String> haplotypesDataset = new DoubleMatrixDataset<String, String>(sampleCounter, uniqueHaploTypeCombinations.size());
        
        //Zero-types        haplotype 0;
        //One-types         haplotype 1;
        //Two-types         haplotype 2;
        
        int columns = 0;
        
        for(String col : uniqueHaploTypeCombinations){
            haplotypesDataset.getHashCols().put(col, columns);
            columns++;
        }
        
        HashMap<String, HashSet<Integer>> uniqueHaploTypes = new HashMap<>();
        
        for(String col : uniqueHaploTypeCombinations){
            int colNumer = haplotypesDataset.getHashCols().get(col);
        
            String[] haplotypes = SPLIT_ON_SLASH.split(col);
            
            for(String h : haplotypes){
                if(!h.equals("-")){
                    if(!uniqueHaploTypes.containsKey(h)){
                        uniqueHaploTypes.put(h, new HashSet<Integer>());
                    } 
                    uniqueHaploTypes.get(h).add(colNumer);
                }
            }
        }
        
        System.out.println("Number of unique Haplotypes: "+uniqueHaploTypes.size());
        
        
//        System.out.println(haplotypesDataset.getHashCols().size());
        
        in = new java.io.BufferedReader(new java.io.FileReader(new File(inputHlaHaplotypes)));
        
        sampleCounter = 0;
        while ((str = in.readLine()) != null) {
            String[] parts = SPLIT_ON_TAB.split(str);
            haplotypesDataset.getHashRows().put(parts[0], sampleCounter);
            
            String[] haploTypesSample = SPLIT_ON_SLASH.split(parts[1]);
            int currentNumber = haplotypesDataset.getHashCols().get(parts[1]);
//            System.out.println(currentNumber);
//            System.out.println(sampleCounter);
            
            haplotypesDataset.getMatrix().setQuick(sampleCounter, currentNumber, 200);
            
//            System.out.println(haploTypesSample[0]);
//            System.out.println("\t"+haploTypesSample[1]);
            if(uniqueHaploTypes.containsKey(haploTypesSample[0])){
//                System.out.println(uniqueHaploTypes.get(haploTypesSample[0]).toString());
                for(int otherPos : uniqueHaploTypes.get(haploTypesSample[0]) ){
                    if(otherPos!=currentNumber){
                        haplotypesDataset.getMatrix().setQuick(sampleCounter, otherPos, 100);
                    }
                }
            } else {
//                System.out.println("-");
            } 
            
            if(uniqueHaploTypes.containsKey(haploTypesSample[1])){
//                System.out.println("\t"+uniqueHaploTypes.get(haploTypesSample[1]).toString());
                if(haploTypesSample[1].equals(haploTypesSample[0])){
                    for(int otherPos : uniqueHaploTypes.get(haploTypesSample[1]) ){
                        if(otherPos!=currentNumber){
                            haplotypesDataset.getMatrix().setQuick(sampleCounter, otherPos, 100);
                        }
                    }
                }
            } else {
//                System.out.println("-");
            }
            sampleCounter++;
        }
        in.close();
        return haplotypesDataset.viewDice();
    }

}
