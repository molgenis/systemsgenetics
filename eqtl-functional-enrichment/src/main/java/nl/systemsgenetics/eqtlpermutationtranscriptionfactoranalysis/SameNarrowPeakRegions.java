/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author Matthieu
 */
public class SameNarrowPeakRegions {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	
	public static void main(String[] args)throws IOException{
		SameNarrowPeakRegions snpr = new SameNarrowPeakRegions("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\narrowPeak\\wgEncodeUwHistoneGm06990H3k27me3StdPkRep1.B36.txt",
				"C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\narrowPeak\\wgEncodeUwHistoneGm06990H3k27me3StdPkRep2.B36.txt");
	}
	
	public SameNarrowPeakRegions(String firstFileLoc, String secondFileLoc)throws IOException{
		HashMap<String, TreeMap<Integer, Integer>> readDataset = readDataset(firstFileLoc);
		searchWithSetTwoThroughSetOne(secondFileLoc, readDataset);
	}
	
	public HashMap<String, TreeMap<Integer, Integer>> readDataset(String datasetFileLocation)throws IOException{
		
		HashMap<String, TreeMap<Integer, Integer>> data = new HashMap<String, TreeMap<Integer, Integer>>();
		TreeMap<Integer, Integer> tmp;
		
		String fileLine;
		String[] fileLineData;
		TextFile tf = new TextFile(datasetFileLocation, false);
		while((fileLine=tf.readLine())!=null){
			
			fileLineData = TAB_PATTERN.split(fileLine);
			String chr = new String(fileLineData[0]);
			int startPos = Integer.parseInt(fileLineData[1]);
			int stopPos = Integer.parseInt(fileLineData[2]);
			
			if(data.containsKey(chr)){
				tmp = data.get(chr);
				tmp.put(startPos, stopPos);
			}
			
			else{
				tmp = new TreeMap<Integer, Integer>();
				tmp.put(startPos, stopPos);
				data.put(chr, tmp);
			}
		}
		tf.close();
		return data;
	}
	
	
	public void searchWithSetTwoThroughSetOne(String otherDataFileLocation, HashMap<String, TreeMap<Integer, Integer>> otherDataset)throws IOException{
		int n = 0;
		String fileLine;
		String[] fileLineData;
		TextFile tf = new TextFile(otherDataFileLocation, false);
		while((fileLine=tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			String chr = new String(fileLineData[0]);
			int startPos = Integer.parseInt(fileLineData[1]);
			int stopPos = Integer.parseInt(fileLineData[2]);
			
			if(otherDataset.containsKey(chr)){
				
				TreeMap<Integer, Integer> chromosomeEntries = otherDataset.get(chr);
				if(chromosomeEntries.containsKey(startPos)){
					int otherStop = chromosomeEntries.get(startPos);
					n++;
					System.out.println(chr + "\t" + startPos + "\t" + stopPos + "\t" + startPos + "\t" + otherStop);
				}
				
			}
		}
		System.out.println(n);
		tf.close();
	}
}
