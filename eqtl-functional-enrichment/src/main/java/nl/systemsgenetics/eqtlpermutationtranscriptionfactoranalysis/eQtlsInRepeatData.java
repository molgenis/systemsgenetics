/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;
import umcg.genetica.genomicboundaries.GenomicBoundaries;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author Matthieu
 */
public class eQtlsInRepeatData {
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	
	public static void main(String[] args) throws IOException{
		eQtlsInRepeatData aap = new eQtlsInRepeatData();
	}
	
	
	public eQtlsInRepeatData()throws IOException{
		System.out.println("[O]:> Read the repeat data.");
		readRepeatDataV3("C:\\Users\\Matthieu\\Documents\\Afstudeerstage\\Pilot\\4.eQtlsInRepeats\\repeats.txt");
		System.out.println("[O]:> Done reading the repeat data.");
	}
	
	
	/*
	 * =========================================================================
	 * = START: EQTL PROCESSING METHODS.
	 * =========================================================================
	 */
	public EQTL[] readEQtlData(String eqtlFileLocation) throws IOException{
		eQTLTextFile eqtlData = new eQTLTextFile(eqtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	public HashMap<String, EQTL> getTopEQtlEffects(EQTL[] eqtls) throws IOException{
		HashMap<String, EQTL> eqtlMap = new HashMap<String, EQTL>();
		for(EQTL eqtl : eqtls){
			
			if( eqtlMap.containsKey(eqtl.getProbe()) ){
				EQTL tmpEqtl = eqtlMap.get(eqtl.getProbe());
				
				if( Math.abs(eqtl.getZscore()) > Math.abs(tmpEqtl.getZscore())){
					eqtlMap.put(eqtl.getProbe(), eqtl);
				}
			}
			
			else{
				eqtlMap.put(eqtl.getProbe(), eqtl);
			}
		}
		return eqtlMap;
	}
	
	
	public Set<String> makeRsIdList(EQTL[] eqtls){
		Set<String> rsIdList = new HashSet<String>();
		for(EQTL eqtl : eqtls){
			rsIdList.add(eqtl.getRsName());
		}
		return rsIdList;
	}
	
	
	public HashMap<String, HashSet<Integer>> getNonTopEqtlEffects(EQTL[] eqtls, HashMap<String, EQTL> topEffects){
		HashMap<String, HashSet<Integer>> otherEffects = new HashMap<String, HashSet<Integer>>();
		HashSet<Integer> tmp;
		int n = 0;
		
		for(EQTL eqtl : eqtls){
			String probe = eqtl.getProbe();
			int pos = eqtl.getRsChrPos();
			
			if(topEffects.containsKey(probe)){
				EQTL topEqtl = topEffects.get(probe);
				int topPos = topEqtl.getRsChrPos();
				
				if(pos != topPos){
					if(otherEffects.containsKey(probe)){
						tmp = otherEffects.get(probe);
						tmp.add(pos);
					}
					else{
						tmp = new HashSet<Integer>();
						tmp.add(pos);
						otherEffects.put(probe, tmp);
					}
					n++;
				}
			}
		}
		return otherEffects;
	}
	/*
	 * =========================================================================
	 * = END: EQTL PROCESSING METHODS.
	 * =========================================================================
	 */
	
	
	
	/* 
	 * =========================================================================
	 * = START: READ REPEAT DATA METHODS.
	 * =========================================================================
	 */
	public void readRepeatData(String repeatFile)throws IOException{
		HashMap<String, TreeMap<Integer, RepeatElement>> repeatData = new HashMap<String, TreeMap<Integer, RepeatElement>>();
		TreeMap<Integer, RepeatElement> chromosomeEntries;
		
		String fileLine;
		String[] fileLineData;
		
		BufferedReader repeatReader = new BufferedReader(new FileReader(repeatFile));
		while( (fileLine = repeatReader.readLine())!=null ){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			String chr = new String(fileLineData[0]);
			int start = Integer.parseInt(fileLineData[1]);
			int stop = Integer.parseInt(fileLineData[2]);
			String repeatClass = new String(fileLineData[3]);
			String repeatFamily = new String(fileLineData[4]);
			String repeatSubFamily = new String(fileLineData[5]);
			
			if(repeatData.containsKey(chr)){
				chromosomeEntries = repeatData.get(chr);
				chromosomeEntries.put(start, new RepeatElement(chr, start, stop, repeatClass, repeatFamily, repeatSubFamily));
			}
			else{
				chromosomeEntries = new TreeMap<Integer, RepeatElement>();
				chromosomeEntries.put(start, new RepeatElement(chr, start, stop, repeatClass, repeatFamily, repeatSubFamily));
				repeatData.put(chr, chromosomeEntries);
			}
		}
		repeatReader.close();
	}
	
	
	public void readRepeatDataV2(String repeatFile)throws IOException{
		HashMap<String, TreeMap<Integer, TreeMap<Integer, RepeatElement>>> repeatData  = new HashMap<String, TreeMap<Integer, TreeMap<Integer, RepeatElement>>>();
		TreeMap<Integer, TreeMap<Integer, RepeatElement>> startPosEntries;
		TreeMap<Integer, RepeatElement> stopPosEntries;
		
		String fileLine;
		String[] fileLineData;
		
		BufferedReader repeatReader = new BufferedReader(new FileReader(repeatFile));
		while( (fileLine = repeatReader.readLine())!=null ){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			String chr = new String(fileLineData[0]);
			int start = Integer.parseInt(fileLineData[1]);
			int stop = Integer.parseInt(fileLineData[2]);
			String repeatClass = new String(fileLineData[3]);
			String repeatFamily = new String(fileLineData[4]);
			String repeatSubFamily = new String(fileLineData[5]);
			
			
			if(repeatData.containsKey(chr)){
				startPosEntries = repeatData.get(chr);
				
			}
		}
		repeatReader.close();
	}
	
	
	
	public void readRepeatDataV3(String repeatFile)throws IOException{
		GenomicBoundaries<Object> repeatRegions = new GenomicBoundaries();
		
		String fileLine;
		String[] fileLineData;
		
		BufferedReader repeatReader = new BufferedReader(new FileReader(repeatFile));
		while( (fileLine = repeatReader.readLine())!=null ){
			fileLineData = TAB_PATTERN.split(fileLine);
			
			String chr = new String(fileLineData[0]);
			int start = Integer.parseInt(fileLineData[1]);
			int stop = Integer.parseInt(fileLineData[2]);
			String repeatClass = new String(fileLineData[3]);
			String repeatFamily = new String(fileLineData[4]);
			String repeatSubFamily = new String(fileLineData[5]);
			String supportData = repeatClass+"|"+repeatFamily+"|"+repeatSubFamily;
			
			repeatRegions.addBoundary(chr, start, stop, supportData);
		}
	}
	/*
	 * =========================================================================
	 * = END: READ REPEAT DATA METHODS.
	 * =========================================================================
	 */
}