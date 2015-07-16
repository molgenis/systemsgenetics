package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author adriaan
 * Used as a class to store all variables in. 
 * Needs to be implemented soon, but the code needs quite some refactoring  
 * 
 */
public class GlobalVariables {
    
    public static int verbosity = 0;
    
    public static int minReads = 1;
    public static int minHets  = 1;
    
    public static int maximumIterations = 20000;
    public static double simplexPrecision  = 1e7;
    
    public static int numberOfTestPerformed = 0;
    
}
