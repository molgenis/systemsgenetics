/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.chrContacts;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;
import static umcg.genetica.io.chrContacts.SortInterChrContacts.writeRawInterContactInformation;

/**
 *
 * @author MaKKie_Admin
 */
public class SortIntraChrContacts {
    
    public static void readNonSortedWriteSorted(String fileToReads, String fileToWrite){
        ArrayList<ChrContact> contacts = null;
        try {
            contacts = readRawIntraContactInformation(fileToReads);
        } catch (IOException ex) {
            Logger.getLogger(SortIntraChrContacts.class.getName()).log(Level.SEVERE, null, ex);
        }
        Collections.sort(contacts);
        
        try {
            writeRawInterContactInformation(contacts, fileToWrite);
        } catch (IOException ex) {
            Logger.getLogger(SortIntraChrContacts.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    private static ArrayList<ChrContact> readRawIntraContactInformation(String fileToReads) throws IOException {
        ArrayList<ChrContact> chrContactInfo = new ArrayList<ChrContact>();
        HashSet<String> obsChrContactInfo = new HashSet<String>();

        BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(fileToReads), "UTF-8"));

        String row;

        while ((row = input.readLine()) != null) {
            String[] parts = StringUtils.split(row, '\t');

            int posChrX = Integer.parseInt(parts[0]);
            int posChrY = Integer.parseInt(parts[1]);
            
            
            int posChr1;
            int posChr2;
            
            if(posChrX<posChrY){
                posChr1 = posChrX;
                posChr2 = posChrY;
            } else {
                posChr1 = posChrY;
                posChr2 = posChrX;
            }
            
            if(!obsChrContactInfo.contains(posChr1+"-"+posChr2)){
                double contact = Double.parseDouble(parts[2]);
                chrContactInfo.add(new ChrContact(posChr1, posChr2, contact));
                obsChrContactInfo.add(posChr1+"-"+posChr2);
            }
            
        }
        input.close();
        return chrContactInfo;

    }

}
