/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.ucsc.UCSCDataObject;
import umcg.genetica.io.ucsc.WigFile;

/**
 *
 * @author Matthieu
 */
public class WriteWigToText {
	public WriteWigToText()throws IOException{
		
	}
	
	/**
	 * Writes the contents of a specified .wig file to a specified .txt file as: chr, start, stop, score. Each field is separated by a tab.
	 * @param wigFile
	 * @param textFileOutLocation
	 * @throws IOException 
	 */
	public void writeWigToText(String wigFile, String textFileOutLocation)throws IOException{
		WigFile wf = new WigFile("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\Gm12864Ctcf.wig", false);
		
		//TextFile outputFile = new TextFile(textFileOut, true);
		TextFile outputFile = new TextFile(textFileOutLocation, true);
		long totalN = wf.countLines();
		long n = 0;
		while(n < totalN){
			UCSCDataObject ucscdo = wf.parseLn();
			if(ucscdo != null){
				outputFile.write(ucscdo.getChr() + "\t" + ucscdo.getPositionStart() + "\t" + ucscdo.getPositionEnd() + "\t" + ucscdo.getValue());
			}
			n++;
		}
		outputFile.close();
		wf.close();
	}
	
	/**
	 * Copies the exact content of a specified .wig file to a specified .txt file.
	 * @param wigFile
	 * @param copyFileLocation
	 * @throws IOException 
	 */
	public void copyWigToText(String wigFile, String copyFileLocation)throws IOException{
		WigFile wf = new WigFile(wigFile, false);
		
		String fileLine;
		TextFile copyFile = new TextFile(copyFileLocation, true);
		while( (fileLine=wf.readLine())!=null ){
			copyFile.write(fileLine);
		}
		
		copyFile.close();
		wf.close();
	}
}
