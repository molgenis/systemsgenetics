package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author adriaan
 * Used as a class to store all variables in that are being accessed by many functions. 
 * 
 * 
 */
public class GlobalVariables {
    
    //How much will we output in stdout.
    //Want to implement the following settings:
    //   0: fully quiet and errors
    //   1: Shows start, end and warnings
    //  10: Shows human readable statistics about what's been done.
    // 100: Shows very big debug statistics
    
    public static int verbosity = 10;
    
    
    public static int maximumIterations = 2000;
    
    public static double simplexThreshold  = 1e-8;
    
    public static int numberOfTestPerformed = 0;
    
    
    
    //BETA BINOMIAL SPECIFIC PARAMETERS
    
    
    //TESTING SPECIFIC PARAMETERS
    public static int minReads = 1;
    public static int minHets  = 1;

    public static double minHetReads = 0.0; 
            
    
    
    //I have some fully phased data and this is what wasp comes up with.
    public static double hetProb  = 1.0;
    public static double seqError = 0.0;
    
    //vcf minimum probability:
    public static double variantProb = 0.99;
    
    //Make scatter plots
    public static String plotDir = "";
    
    
}
