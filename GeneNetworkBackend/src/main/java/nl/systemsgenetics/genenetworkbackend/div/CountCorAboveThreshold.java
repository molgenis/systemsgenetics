/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.div;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author patri
 */
public class CountCorAboveThreshold {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		File corFile = new File(args[0]);
		double minCor = Double.parseDouble(args[1]);
		
		System.err.println("File: " + corFile.getAbsolutePath());
		System.err.println("Min cor: " + minCor);
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(corFile))));
		
		reader.readLine();//skip header
		
		String line; 
		int countAboveCor;
		while( (line = reader.readLine()) != null){
			
			countAboveCor = 0;
			
			String[] data = StringUtils.splitPreserveAllTokens(line, '\t');
			
			for(int i = 1; i < data.length ; ++i){
				
				if(Double.parseDouble(data[i]) >= minCor){
					countAboveCor++;
				}
				
			}
			
			System.out.println(data[0] + "\t" + countAboveCor);
			
		}
		
		reader.close();
		
	}
	
}
