package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

/*
 * Copyright (C) 2015 Adriaan van der Graaf
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
    
    
    public static int maximumIterations = 20000;
    
    public static double simplexThreshold  = 1e-8;
    
    public static int numberOfTestPerformed = 0;
    
    
    
    //BETA BINOMIAL SPECIFIC PARAMETERS
    
    
    //TESTING SPECIFIC PARAMETERS
    public static int minReads = 1;
    public static int minHets  = 1;
    
    
    
    //I have some fully phased data and this is what wasp comes up with.
    public static double hetProb  = 1.0;
    public static double seqError = 0.0;
    
    //vcf minimum probability:
    public static double variantProb = 0.99;
}
