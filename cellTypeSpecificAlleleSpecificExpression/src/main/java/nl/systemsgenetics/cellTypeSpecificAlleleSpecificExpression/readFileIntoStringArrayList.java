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
 */
public class readFileIntoStringArrayList {
    
    private ArrayList<String> lines;
    
    readFileIntoStringArrayList(String filenames_file) throws IOException {
        this.lines = new ArrayList<String>();

        File file = new File(filenames_file);
        FileReader fr = new FileReader(file);
        BufferedReader br = new BufferedReader(fr);
        String line;

        while ((line = br.readLine()) != null) {
            
            lines.add(line);
            
        }
        br.close();
        fr.close();
        
    }

    /**
     * @return the lines
     */
    public ArrayList<String> getLines() {
        return lines;
    }
    
    
}
