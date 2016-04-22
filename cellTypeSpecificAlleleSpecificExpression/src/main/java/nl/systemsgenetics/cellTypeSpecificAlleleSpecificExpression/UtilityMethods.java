/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.cellTypeSpecificAlleleSpecificExpression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author adriaan
 * This class is used as a combination of static methods, so we can use a 
 * Utility list.
 */
public class UtilityMethods {
    

    public static ArrayList<String> readFileIntoStringArrayList(String filenames_file) throws IOException {
        /**
         * As the name implies, this will read files into an ArrayListString
         * Probably there is a more efficient way to do this, but currently it works.
         */
        
        
        ArrayList<String> lines = new ArrayList<String>();

        File file = new File(filenames_file);
        FileReader fr = new FileReader(file);
        BufferedReader br = new BufferedReader(fr);
        String line;

        while ((line = br.readLine()) != null) {
            
            lines.add(line);
            
        }
        br.close();
        fr.close();
        
        
        return lines;
    }
    
    static ArrayList<IndividualSnpData> parse_AS_lines(ArrayList<String> all_line_data, ArrayList<File> all_files) {

        ArrayList<IndividualSnpData> all_individuals;
        all_individuals = new ArrayList<IndividualSnpData>();
        int i = 0;

        for (String iLine : all_line_data) {

            all_individuals.add(new IndividualSnpData(all_files.get(i).getAbsolutePath(), iLine));

            i++;
        }

        return all_individuals;
    }
    
    //This 'function' will look for heterozygotes that also have enough reads to start testing
    static ArrayList<IndividualSnpData> isolateValidHeterozygotesFromIndividualSnpData(ArrayList<IndividualSnpData> all_individuals) {
        
        ArrayList<IndividualSnpData> hets;
        hets = new ArrayList<IndividualSnpData>();
        
        for(IndividualSnpData sample : all_individuals){
            
            String genotype = sample.getGenotype();
            
            boolean validReads = valid_heterozygote(sample);
            
            //assuming the genotype is formatted as: "[C, A]"
            if(isGenoHeterozygote(genotype) && validReads){
                hets.add(sample);
            }       
        }
        
        return hets;
    }
    
    static boolean valid_heterozygote(IndividualSnpData sample){            
        double thisHetReadProp = (double)sample.getRefNum()  / ((double)sample.getRefNum() + (double)sample.getAltNum());
        return ((thisHetReadProp <= (1 - GlobalVariables.minHetReads)) && (thisHetReadProp >= (GlobalVariables.minHetReads)));
    }
    
    
    //This 'function' will only look for heterozygote snps for phased entry
    static ArrayList<IndividualSnpData> isolateOnlyHeterozygotesFromIndividualSnpData(ArrayList<IndividualSnpData> all_individuals) {
        
        ArrayList<IndividualSnpData> hets;
        hets = new ArrayList<IndividualSnpData>();
        
        for(IndividualSnpData sample : all_individuals){
            
            String genotype = sample.getGenotype();
            
            //assuming the genotype is formatted as: "[C, A]"
            
            if(isGenoHeterozygote(genotype)){
                hets.add(sample);
            }       
        }
        
        return hets;
    }
    
    static boolean isGenoHeterozygote(String genotype){
        return (genotype.charAt(1) != genotype.charAt(4));
    }
    
    static void printFinalTestStats(){
        if(GlobalVariables.verbosity >= 1){
            System.out.println("\n-------------------------------------------");
            System.out.println("Finished " + GlobalVariables.numberOfTestPerformed + " tests, now closing.");
            System.out.println("-------------------------------------------");
        }
    }
    
}
