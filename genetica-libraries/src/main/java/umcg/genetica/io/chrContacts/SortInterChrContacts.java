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
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author MaKKie_Admin
 */
public class SortInterChrContacts {
    
    public static void readNonSortedWriteSorted(String fileToReads, String fileToWrite){
        ArrayList<ChrContact> contacts = null;
        try {
            contacts = readRawInterContactInformation(fileToReads);
        } catch (IOException ex) {
            Logger.getLogger(SortInterChrContacts.class.getName()).log(Level.SEVERE, null, ex);
        }
        Collections.sort(contacts);
        
        try {
            writeRawInterContactInformation(contacts, fileToWrite);
        } catch (IOException ex) {
            Logger.getLogger(SortInterChrContacts.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    private static ArrayList<ChrContact> readRawInterContactInformation(String fileToReads) throws IOException {
        ArrayList<ChrContact> chrContactInfo = new ArrayList<ChrContact>();

        BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(fileToReads), "UTF-8"));

        String row;

        while ((row = input.readLine()) != null) {
            String[] parts = StringUtils.split(row, '\t');

            int posChr1 = Integer.parseInt(parts[0]);
            int posChr2 = Integer.parseInt(parts[1]);
            double contact = Double.parseDouble(parts[2]);
            chrContactInfo.add(new ChrContact(posChr1, posChr2, contact));
        }
        input.close();
        return chrContactInfo;

    }
    
    public static void writeRawInterContactInformation(ArrayList<ChrContact> contacts, String fileToWrite) throws IOException {

        TextFile outWriter = new TextFile(fileToWrite, TextFile.W);

        for(ChrContact contact : contacts){
            outWriter.writeln(contact.getChrLocationSmaller()+"\t"+contact.getChrLocationLarger()+"\t"+contact.getContactValue());
        }
        outWriter.close();

    }
}
