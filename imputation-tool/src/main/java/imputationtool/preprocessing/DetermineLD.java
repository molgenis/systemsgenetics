///*
// * DetermineLD.java
// *
// * Created on August 8, 2006, 4:32 PM
// *
// * To change this template, choose Tools | Template Manager
// * and open the template in the editor.
// */
//
//package imputationtool.preprocessing;
//
//import java.lang.Math;
//import metaqtl2.Containers.GenDataset;
//import metaqtl2.Containers.SNP;
//
///**
// *
// * @author lude
// */
//public class DetermineLD {
//
//    public final static int INCLUDE_CASES_AND_CONTROLS = 1;
//    public final static int INCLUDE_CASES = 2;
//    public final static int INCLUDE_CONTROLS = 3;
//    
//    public final static int RETURN_R_SQUARED = 4;
//    public final static int RETURN_D_PRIME = 5;
//    
//    public double h11 = 0;
//    public double h12 = 0;
//    public double h21 = 0;
//    public double h22 = 0;
//    
//    //Triallelic LD calculations:
//    public double h31 = 0;
//    public double h32 = 0;
//    
//    public int nrCalledGenotypes = 0;
//    
//    /** Creates a new instance of DetermineLD */
//    public DetermineLD(){
//    }
//    
//    public synchronized double getRSquared(SNP snpX, SNP snpY, GenDataset genotypeData, int returnType, int individualsToInclude, boolean print) {
//
//        if (snpX==null||snpY==null) return 0;
//
//	int[] genotypesX = snpX.getGenotype();
//	int[] genotypesY = snpY.getGenotype();
//	
//	if(genotypesX.length!=genotypesY.length)return 0;
//        //Get genotypes:
//        int[][] genotypes = new int[3][3];
//        nrCalledGenotypes = 0;
//	
//
//        for (int ind=0; ind<genotypesX.length; ind++) {
//            if (individualsToInclude==INCLUDE_CASES_AND_CONTROLS||(genotypeData.isCase(ind)&&individualsToInclude==INCLUDE_CASES)||(!genotypeData.isCase(ind)&&individualsToInclude==INCLUDE_CONTROLS)) {
//		
//                int genotypeX = genotypesX[ind];
//                int genotypeY = genotypesY[ind];
//                if (genotypeX!=-1&&genotypeY!=-1) {
//                    genotypes[genotypeX][genotypeY]++;
//                    nrCalledGenotypes++;
//                }
//            }
//        }
//        //System.out.println("NrCalledGenotypes:\t" + nrCalledGenotypes);
//
//        //Determine genotype frequencies:
//        double[][] genotypesFreq = new double[3][3];
//        for (int x=0; x<3; x++) {
//            for (int y=0; y<3; y++) {
//                genotypesFreq[x][y] = (double) genotypes[x][y] / (double) nrCalledGenotypes;
//                if (print) System.out.print(genotypes[x][y] + "\t");
//            }
//            if (print) System.out.println("");
//        }
//
//        //Determine alle frequencies:
//        double[][] alleleFreq = new double[2][2];
//        //SNP X:
//        alleleFreq[0][0] = (genotypesFreq[0][0] + genotypesFreq[0][1] + genotypesFreq[0][2]) + (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
//        alleleFreq[0][1] = (genotypesFreq[2][0] + genotypesFreq[2][1] + genotypesFreq[2][2]) + (genotypesFreq[1][0] + genotypesFreq[1][1] + genotypesFreq[1][2]) / 2d;
//        //SNP Y:
//        alleleFreq[1][0] = (genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[2][0]) + (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;
//        alleleFreq[1][1] = (genotypesFreq[0][2] + genotypesFreq[1][2] + genotypesFreq[2][2]) + (genotypesFreq[0][1] + genotypesFreq[1][1] + genotypesFreq[2][1]) / 2d;
//
//        if (print) {
//            System.out.println("Allele freq:");
//            System.out.println(alleleFreq[0][0]);
//            System.out.println(alleleFreq[0][1]);
//            System.out.println(alleleFreq[1][0]);
//            System.out.println(alleleFreq[1][1]);
//        }
//        
//        //Precalculate triangles of non-double heterozygote:
//        double[][] genotypesTriangleFreq = new double[3][3];
//        genotypesTriangleFreq[0][0] = 2d * genotypesFreq[0][0] + genotypesFreq[1][0] + genotypesFreq[0][1];
//        genotypesTriangleFreq[2][0] = 2d * genotypesFreq[2][0] + genotypesFreq[1][0] + genotypesFreq[2][1];
//        genotypesTriangleFreq[0][2] = 2d * genotypesFreq[0][2] + genotypesFreq[0][1] + genotypesFreq[1][2];
//        genotypesTriangleFreq[2][2] = 2d * genotypesFreq[2][2] + genotypesFreq[1][2] + genotypesFreq[2][1];
//        if (print) {
//            System.out.println("Triangle freq:");
//            System.out.println(genotypesTriangleFreq[0][0]);
//            System.out.println(genotypesTriangleFreq[0][2]);
//            System.out.println(genotypesTriangleFreq[2][0]);
//            System.out.println(genotypesTriangleFreq[2][2]);
//        }
//
//        //Calculate expected genotypes, assuming equilibrium, take this as start:
//        h11 = alleleFreq[0][0] * alleleFreq[1][0]; 
//        h12 = alleleFreq[0][0] * alleleFreq[1][1]; 
//        h21 = alleleFreq[0][1] * alleleFreq[1][0]; 
//        h22 = alleleFreq[0][1] * alleleFreq[1][1];
//
//        //Calculate the frequency of the two double heterozygotes:
//        double x12y12 = h11 * h22 / (h11*h11 + h12*h21) * genotypesFreq[1][1];
//        double x12y21 = h12 * h21 / (h11*h11 + h12*h21) * genotypesFreq[1][1];
//        
//        if (print) {
//            System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22 + "\t\t" + x12y12 + "\t" + x12y21);
//        }
//
//        //Perform iterations using EM algorithm:
//        for (int itr=0; itr<25; itr++) {
//
//            h11 = (x12y12 + genotypesTriangleFreq[0][0]) / 2;
//            h12 = (x12y21 + genotypesTriangleFreq[2][0]) / 2;
//            h21 = (x12y21 + genotypesTriangleFreq[0][2]) / 2;
//            h22 = (x12y12 + genotypesTriangleFreq[2][2]) / 2;
//
//            x12y12 = h11 * h22 / (h11*h22 + h12*h21) * genotypesFreq[1][1];
//            x12y21 = h12 * h21 / (h11*h22 + h12*h21) * genotypesFreq[1][1];
//
//            if (print) {
//                System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22 + "\t\t" + x12y12 + "\t" + x12y21);
//            }
//
//        }
//
//        double d = h11  - (alleleFreq[0][0] * alleleFreq[1][0]);
//        
//        if (returnType==RETURN_R_SQUARED) {
//            double rSquared = d * d / (alleleFreq[0][0]*alleleFreq[0][1]*alleleFreq[1][0]*alleleFreq[1][1]);
//            return rSquared;
//        } else {
//            double dMax = 0;
//            if (d<0) {
//                double a = alleleFreq[0][1]*alleleFreq[1][1];
//                if (alleleFreq[0][0]>alleleFreq[1][0]) a = alleleFreq[0][0]*alleleFreq[1][0];
//                double b = alleleFreq[0][0]*alleleFreq[1][0];
//                if (alleleFreq[0][0]>alleleFreq[1][0]) b = alleleFreq[0][1]*alleleFreq[1][1];
//                dMax = Math.min(a,b);
//            } else {
//                double a = alleleFreq[0][1]*alleleFreq[1][0];
//                if (alleleFreq[0][0]>alleleFreq[1][0]) a = alleleFreq[0][0]*alleleFreq[1][1];
//                double b = alleleFreq[0][0]*alleleFreq[1][1];
//                if (alleleFreq[0][0]>alleleFreq[1][0]) b = alleleFreq[0][1]*alleleFreq[1][0];
//                dMax = Math.min(a,b);
//            }
//            double dPrime = Math.abs(d / dMax);
//            /*
//            if (dPrime>1.01) {
//                System.out.println("");
//                System.out.println(genotypes[0][0] + "\t" + genotypes[0][1] + "\t" + genotypes[0][2]);
//                System.out.println(genotypes[1][0] + "\t" + genotypes[1][1] + "\t" + genotypes[1][2]);
//                System.out.println(genotypes[2][0] + "\t" + genotypes[2][1] + "\t" + genotypes[2][2]);
//                System.out.println(alleleFreq[0][0] + "\t" + alleleFreq[0][1] + "\t" + alleleFreq[1][0] + "\t" + alleleFreq[1][1]);
//                System.out.println(h11 + "\t" + h12 + "\t" + h21 + "\t" + h22);
//                System.out.println(d + "\t" + dMax + "\t" + dPrime);
//            }
//             */
//            return Math.min(1,dPrime);
//        }
//    }
//    
//
//    
//    
//}
