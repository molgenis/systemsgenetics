/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlpermutationtranscriptionfactoranalysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author Matthieu
 */
public class eQtlDataParser {
	
	
	public eQtlDataParser(){
		
	}
	
	
	public EQTL[] readEQtlData(String eqtlFileLocation) throws IOException{
		QTLTextFile eqtlData = new QTLTextFile(eqtlFileLocation, false);
		return eqtlData.read();
	}
	
	
	
	public HashMap<String, EQTL> getTopEffects(EQTL[] eqtls, HashSet<String> sharedProbes) throws IOException{
		HashMap<String, EQTL> eqtlMap = new HashMap<String, EQTL>();
		for(EQTL eqtl : eqtls){
			
			String eqtlProbe = eqtl.getProbe();
			if(sharedProbes.contains(eqtlProbe)){
				
				//Check if the eQTL for the shared probe has the highest absolute Z-Score.
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
		}
		return eqtlMap;
	}
	
	
	
	public HashMap<String, HashSet<Integer>> getNonTopEffects(EQTL[] eqtls, HashMap<String, EQTL> topEffects){
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
	
	
	
	public Set<String> makeRsIdList(EQTL[] eqtls){
		Set<String> rsIdList = new HashSet<String>();
		for(EQTL eqtl : eqtls){
			rsIdList.add(eqtl.getRsName());
		}
		return rsIdList;
	}
	
}
