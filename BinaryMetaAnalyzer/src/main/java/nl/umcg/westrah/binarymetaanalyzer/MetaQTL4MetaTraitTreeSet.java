/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import java.util.NavigableSet;
import java.util.TreeSet;

/**
 *
 * @author Patrick Deelen
 */
public class MetaQTL4MetaTraitTreeSet extends TreeSet<MetaQTL4MetaTrait>{
	
	public NavigableSet<MetaQTL4MetaTrait> getTraitByRange(String startChr, int startPos, String stopChr, int stopPos){
		return subSet(
				new MetaQTL4MetaTrait(Integer.MIN_VALUE, null, startChr, startPos, startPos, null, null), true,
				new MetaQTL4MetaTrait(Integer.MAX_VALUE, null, stopChr, stopPos, stopPos, null, null), true);
	}
	
	public NavigableSet<MetaQTL4MetaTrait> getTraitByRange(String chr, int startPos, int stopPos){
		return getTraitByRange(chr, startPos, chr, stopPos);
	}
	
	public NavigableSet<MetaQTL4MetaTrait> getTraitInWindow(String chr, int windowCenter, int windowHalfSize){
		return getTraitByRange(chr, windowCenter - windowHalfSize, chr, windowCenter + windowHalfSize);
	}
	
}
