/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author Matthieu
 */
public class SplitWigTextDataOnChromosome {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	private static final Pattern PATH_PATTERN = Pattern.compile("\\\\");
	private static final Pattern FILENAME_PATTERN = Pattern.compile("\\.");
	
	public static void main(String[] args)throws IOException{
		SplitWigTextDataOnChromosome swtdoc = new SplitWigTextDataOnChromosome("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\Gm12864Ctcf.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\Gm12864Ctcf\\");
	}
	
	
	public SplitWigTextDataOnChromosome(String wigTextFile, String outputFolder)throws IOException{
		String fileLine;
		String[] fileLineData;
		String previousChr = "#";
		String fileName = getFileName(wigTextFile);
		
		ArrayList<String> dataToWriteToFile = new ArrayList<String>();
		TextFile tf = new TextFile(wigTextFile, false);
		while((fileLine=tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			String chr = new String(fileLineData[0]);
			String start = new String(fileLineData[1]);
			String stop = new String(fileLineData[2]);
			String score = new String(fileLineData[3]);
			
			
			if(chr.equals(previousChr)){
				dataToWriteToFile.add(fileLine);
			}
			
			else{
				//Write data to file.
				if(!(previousChr.equals("#"))){
					writeToTextFile(dataToWriteToFile, outputFolder, fileName, previousChr);
				
					//Change the saved data to hold the next batch of data.
					dataToWriteToFile = new ArrayList<String>();
					dataToWriteToFile.add(fileLine);
					//previousChr = chr;
				}
				previousChr = chr;
			}
			
		}
		tf.close();
	}
	
	
	private void writeToTextFile(ArrayList<String> dataToWrite, String outputFolder, String mainFileName, String chromosome)throws IOException{
		String outputFile = outputFolder+mainFileName+"_chr"+chromosome+".txt";
		PrintWriter pw = new PrintWriter( new FileWriter(outputFile) );
		for(String entry : dataToWrite){
			pw.println(entry);
		}
		pw.close();
	}
	
	
	private String getFileName(String filePath){
		String[] filePathData = PATH_PATTERN.split(filePath);
		int size = filePathData.length;
		String[] fileNameData = FILENAME_PATTERN.split(filePathData[size-1]);
		return new String(fileNameData[0]);
	}
}
