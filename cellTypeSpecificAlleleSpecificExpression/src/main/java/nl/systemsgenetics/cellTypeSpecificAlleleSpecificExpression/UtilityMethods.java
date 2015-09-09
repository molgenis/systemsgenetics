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
    
    static ArrayList<IndividualSnpData> isolateHeterozygotesFromIndividualSnpData(ArrayList<IndividualSnpData> all_individuals) {
        
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
