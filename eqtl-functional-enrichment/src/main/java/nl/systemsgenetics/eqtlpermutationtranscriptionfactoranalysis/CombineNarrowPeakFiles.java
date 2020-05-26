/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author Matthieu
 */
public class CombineNarrowPeakFiles {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public static void main(String[] args)throws IOException{
		CombineNarrowPeakFiles cnpf = new CombineNarrowPeakFiles();
	}
	
	public CombineNarrowPeakFiles()throws IOException{
		HashMap<String, TreeMap<Integer, NarrowPeakElement>> data = readBaseFile("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\narrowPeak\\wgEncodeUwHistoneK562H3k4me3StdHotspotsRep1.broadPeak.B36.txt");
		addNewRegionsFromOtherFile("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\narrowPeak\\wgEncodeUwHistoneK562H3k4me3StdHotspotsRep2.broadPeak.B36.txt", data);
		writeCombinedEntriesTopFile(data, "C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\3.histones\\narrowPeak\\wgEncodeUwHistoneK562H3k4me3StdHotspots.broadPeak.B36.txt");
	}
	
	public HashMap<String, TreeMap<Integer, NarrowPeakElement>> readBaseFile(String baseFileLocation)throws IOException{
		HashMap<String, TreeMap<Integer, NarrowPeakElement>> baseData = new HashMap<String, TreeMap<Integer, NarrowPeakElement>>();
		TreeMap<Integer, NarrowPeakElement> baseDataSubMap;
		String fileLine;
		String[] fileLineData;
		
		TextFile tf = new TextFile(baseFileLocation, false);
		while((fileLine=tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			String chr = new String(fileLineData[0]);
			int startPos = Integer.parseInt(new String(fileLineData[1]));
			int stopPos = Integer.parseInt(new String(fileLineData[2]));
			double peakValue = Double.valueOf(new String(fileLineData[6]));
			double log10pvalue = Double.valueOf(new String(fileLineData[7]));
			
			NarrowPeakElement npe = new NarrowPeakElement(chr, startPos, stopPos, peakValue, log10pvalue, fileLine);
			
			if(baseData.containsKey(chr)){
				baseDataSubMap = baseData.get(chr);
				baseDataSubMap.put(startPos, npe);
			}
			else{
				baseDataSubMap = new TreeMap<Integer, NarrowPeakElement>();
				baseDataSubMap.put(startPos, npe);
				baseData.put(chr, baseDataSubMap);
			}
		}
		tf.close();
		return baseData;
	}
	
	
	public void addNewRegionsFromOtherFile(String otherFileLocation, HashMap<String, TreeMap<Integer, NarrowPeakElement>> baseData)throws IOException{
		TreeMap<Integer, NarrowPeakElement> tmp;
		String fileLine;
		String[] fileLineData;
		
		TextFile tf = new TextFile(otherFileLocation, false);
		while((fileLine=tf.readLine())!=null){
			fileLineData = TAB_PATTERN.split(fileLine);
			String chr = new String(fileLineData[0]);
			int startPos = Integer.parseInt(new String(fileLineData[1]));
			int stopPos = Integer.parseInt(new String(fileLineData[2]));
			double peakValue = Double.valueOf(new String(fileLineData[6]));
			double log10pvalue = Double.valueOf(new String(fileLineData[7]));
			NarrowPeakElement npe = new NarrowPeakElement(chr, startPos, stopPos, peakValue, log10pvalue, fileLine);
			
			if(baseData.containsKey(chr)){
				tmp = baseData.get(chr);
				
				if(!(tmp.containsKey(startPos))){
					tmp.put(startPos, npe);
				}
			}
		}
		tf.close();
	}
	
	
	public void writeCombinedEntriesTopFile(HashMap<String, TreeMap<Integer, NarrowPeakElement>> combinedData, String outputFileLocation)throws IOException{
		int n = 0;
		TextFile tf = new TextFile(outputFileLocation, true);
		Iterator<Entry<String, TreeMap<Integer, NarrowPeakElement>>> dataIterator = combinedData.entrySet().iterator();
		while(dataIterator.hasNext()){
			Entry<String, TreeMap<Integer, NarrowPeakElement>> chromosomeEntries = dataIterator.next();
			TreeMap<Integer, NarrowPeakElement> chrPosEntries = chromosomeEntries.getValue();
			
			Iterator<Entry<Integer, NarrowPeakElement>> chrPosEntriesIterator = chrPosEntries.entrySet().iterator();
			while(chrPosEntriesIterator.hasNext()){
				Entry<Integer, NarrowPeakElement> chrPosElement = chrPosEntriesIterator.next();
				tf.write(chrPosElement.getValue().getInfoLine() + "\n");
				n++;
			}
		}
		System.out.println(n);
		tf.close();
	}
}
