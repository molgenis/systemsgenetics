/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputationtool.postprocessing;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import org.apache.commons.collections.primitives.ArrayDoubleList;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class CompareImputedSNPsToRealSNPs {

    private final TriTyperGenotypeData ggDataset1;
    private final TriTyperGenotypeData ggDataset2;  
    
    public static void main(String[] args) throws IOException {
        //LiverOmni  vs LiverCyto
        //0.9997752239290337
        //CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverOmni\\", "Unimputed", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverCyto\\", "Imputed");
        
        //unimputed merged vs LiverCyto
        //0.6492976314532056
        //CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverMerged\\", "Unimputed", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverCyto\\", "Imputed");
        
        //unimputed merged vs LiverOmni 
        //0.3528032918756224
        //CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverMerged\\", "Unimputed", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverOmni\\", "Imputed");
        
        //LiverCyto vs original imputed
        //0.9961702580911673
        //CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverCyto\\", "Unimputed", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\", "Imputed");
        
        //LiverOmni vs original imputed
        //0.9982103900136859
        //CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverOmni\\", "Unimputed", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\", "Imputed");
        
        //LiverCyto  vs omni imputed 100g
        // 0.9956201693043798
        //CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverCyto\\", "Unimputed", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverImputed100G\\", "Imputed");
        
        //LiverOmni vs omni imputed 100g
        // 0.9959821144437336
        //CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverOmni\\", "Unimputed", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverImputed100G\\", "Imputed");

        //original imputed vs omni imputed 100g
        // 0.9893960960146694
        CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\", "Imputed_1", "D:\\UMCG\\SAT-VAT-Liver-Muscle-ImputeTriTyper\\Liver\\LiverImputed100G\\", "Imputed");
    }
    
    public void run(String[] args) throws IOException {
        CompareImputedSNPsToRealSNPs c = new CompareImputedSNPsToRealSNPs("/Data/GeneticalGenomicsDatasets/RotterdamStudy/TriTyper773Samples-GenotypesQCedRotterdam/", "Unimputed", "/Data/GeneticalGenomicsDatasets/RotterdamStudy/MachImputed/TriTyperFixed/", "Imputed");
    }

    public CompareImputedSNPsToRealSNPs(String dataset1, String dataset1Name, String dataset2, String dataset2Name) throws IOException {

        ggDataset1 = new TriTyperGenotypeData();
        ggDataset1.load(dataset1);

        ggDataset2 = new TriTyperGenotypeData();
        ggDataset2.load(dataset2);

        //Added because of potentialy missing snps in both.
        ArrayList<String> snpsTmp = new ArrayList<String>();
        if (ggDataset1.getSNPs().length < ggDataset2.getSNPs().length) {
            
            List<String> list = Arrays.asList(ggDataset2.getSNPs());
            HashSet<String> testSnps = new HashSet<String>(list);
            list = null;
            String[] snps = ggDataset1.getSNPs();
            
            for(int i=0; i<snps.length; ++i){
                if(testSnps.contains(snps[i])){
                    snpsTmp.add(snps[i]);
                }
            }
            testSnps = null;
            
        } else {
            
            List<String> list = Arrays.asList(ggDataset1.getSNPs());
            HashSet<String> testSnps = new HashSet<String>(list);
            list = null;
            String[] snps = ggDataset2.getSNPs();
            
            for(int i=0; i<snps.length; ++i){
                if(testSnps.contains(snps[i])){
                    snpsTmp.add(snps[i]);
                }
            }
            testSnps = null;
            
        }
        System.out.println("Numer of snps shared: "+ snpsTmp.size());
        String[] snps = snpsTmp.toArray(new String[snpsTmp.size()]);
        snpsTmp = null;
        
        
        String[] ind1 = ggDataset1.getIndividuals();
        //String[] ind2 = ggDataset2.getIndividuals();
        int[] ids1 = new int[ind1.length];
//        int[] ids2 = new int[ind1.length];
        for (int i = 0; i < ids1.length; i++) {
            String name1 = ind1[i];
            Integer id2 = ggDataset2.getIndividualId(name1);

            if (id2 != -9) {
                ids1[i] = id2;
            } else {
                ids1[i] = -1;
            }
        }

        long numTotalConcordant = 0;
        long numTotalCalled = 0;
        ArrayDoubleList maf1 = new ArrayDoubleList();
        ArrayDoubleList maf2 = new ArrayDoubleList();
        ArrayDoubleList hwep1 = new ArrayDoubleList();
        ArrayDoubleList hwep2 = new ArrayDoubleList();
        
        int q = 0;
        int concordant = 0;
        int called = 0;
        int nrIncompatible = 0;
        int nrAT = 0;


        SNPLoader loader1 = ggDataset1.createSNPLoader();
        SNPLoader loader2 = ggDataset2.createSNPLoader();

        for (String snp : snps) {

            Integer snp1Id = ggDataset1.getSnpToSNPId().get(snp);
            Integer snp2Id = ggDataset2.getSnpToSNPId().get(snp);

            if (snp1Id != -9 && snp2Id != -9) {
                SNP snp1 = ggDataset1.getSNPObject(snp1Id);
                SNP snp2 = ggDataset2.getSNPObject(snp2Id);

                loader1.loadGenotypes(snp1);
                loader2.loadGenotypes(snp2);
                
                int[] genotype1 = new int[ids1.length];
                int[] genotype2 = new int[ids1.length];

                byte[] alleles1 = snp1.getAlleles();
                byte[] alleles2 = snp2.getAlleles();
                
                if(alleles1 != null && alleles2 != null){
//                    System.out.println("soutje");


                    String str1allele1 = BaseAnnot.toString(alleles1[0]);
                    String str1allele2 = BaseAnnot.toString(alleles1[1]);


                    String str2allele1 = BaseAnnot.toString(alleles2[0]);
                    String str2allele2 = BaseAnnot.toString(alleles2[1]);

                    byte minor1 = snp1.getMinorAllele();
                    byte minor2 = snp2.getMinorAllele();

                    String strminor1 = BaseAnnot.toString(minor1);
                    String strminor2 = BaseAnnot.toString(minor2);

                    boolean atOrCG = false;
                    boolean inverted = false;


                    boolean incompatible = false;

                    if (BaseAnnot.getComplement(alleles1[0]) == alleles1[1]) {
                        // AT or CG snp
                        atOrCG = true;
                        nrAT++;
                    } else {
                        if (alleles1[0] == alleles2[0] && alleles1[1] == alleles2[1]) {
                            inverted = false;
                        } else if (alleles1[0] == alleles2[1] && alleles1[1] == alleles2[0]) {
                            inverted = true;
                        } else if (BaseAnnot.getComplement(alleles1[0]) == alleles2[0] && BaseAnnot.getComplement(alleles1[1]) == alleles2[1]) {
                            inverted = false;
                        } else if (BaseAnnot.getComplement(alleles1[0]) == alleles2[1] && BaseAnnot.getComplement(alleles1[1]) == alleles2[0]) {
                            inverted = true;
                        } else {
                            incompatible = true;
                            nrIncompatible++;
                        }
                    }

                    if (!atOrCG && !incompatible) {

                        for (int i = 0; i < ids1.length; i++) {
                            if (ids1[i] != -1) {
                                genotype1[i] = snp1.getGenotypes()[i];
                                genotype2[i] = snp2.getGenotypes()[ids1[i]];

                                int gt1 = genotype1[i];
                                int gt2 = genotype2[i];

                                if (inverted) {
                                    if (gt2 == 0) {
                                        gt2 = 2;
                                    } else if (gt2 == 2) {
                                        gt2 = 0;
                                    }
                                }

                                if (gt1 != -1 && gt2 != -1) {
                                    if (gt2 == gt1) {
                                        numTotalConcordant++;
                                    }
                                    numTotalCalled++;
                                }
                            }
                        }

//                    if(called > 0 && concordant == 0){
//                        System.out.println("0 Concordant!");
//                    }



                    }

                    
                    if(!Double.isNaN(snp1.getMAF())&& !Double.isNaN(snp2.getMAF())){
                        maf1.add(snp1.getMAF());
                        maf2.add(snp2.getMAF());
                    }
                    
                    
                    hwep1.add(snp1.getHWEP());
                    hwep2.add(snp2.getHWEP());
                    
//                  System.out.println("Num shared: "+ concordant +" / "+ called);

                    snp1.clearGenotypes();
                    snp2.clearGenotypes();

                }
            }

            q++;

            if (q % 10000 == 0) {
                System.out.println(q + "snps parsed");
            }
        }

        loader1.close();
        loader2.close();
        System.out.println("Total called: " + numTotalCalled + "\t" + numTotalConcordant + "\t" + ((double) numTotalConcordant / numTotalCalled) + "\t" + nrIncompatible + "\t" + nrAT);
        System.out.println("Correlations between minor allel frequencies: "+JSci.maths.ArrayMath.correlation(maf1.toArray(), maf2.toArray()));
        System.out.println("Correlations between hwep "+JSci.maths.ArrayMath.correlation(hwep1.toArray(), hwep2.toArray()));
    }
}
