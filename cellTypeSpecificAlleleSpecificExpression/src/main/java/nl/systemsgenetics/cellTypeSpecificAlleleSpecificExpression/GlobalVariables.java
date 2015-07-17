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
    
    public static int verbosity = 0;
    
    public static int minReads = 1;
    public static int minHets  = 1;
    
    public static int maximumIterations = 20000;
    
    public static double simplexThreshold  = 1e-8;
    
    public static int numberOfTestPerformed = 0;
    
    
    
    //BETA BINOMIAL SPECIFIC PARAMETERS
    
    //I have some fully phased data and this is what wasp comes up with.
    public static double hetProb  = 0.980198;
    public static double seqError = 0.005;
}
