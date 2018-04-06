package umcg.genetica.collections;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import edu.ufl.cise.colamd.tdouble.Colamd_Col;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.NavigableMap;
import java.util.TreeMap;
import org.molgenis.genotype.util.ChrPos;

/**
 *
 * @author Patrick Deelen
 */
public class ChrPosTreeMap<E> implements Iterable<E>{

	private final HashMap<String, TreeMap<Integer, E>> data;

	public ChrPosTreeMap() {
		data = new HashMap<String, TreeMap<Integer, E>>();
	}
	
	public void put(String chr, Integer pos, E element){
		TreeMap<Integer, E> chrElements = data.get(chr);
		if(chrElements == null){
			chrElements = new TreeMap<Integer, E>();
			data.put(chr, chrElements);
		}
		chrElements.put(pos, element);
	}
	
	/**
	 * 
	 * @param chr
	 * @param pos
	 * @return null if not found
	 */
	public E get(String chr, Integer pos){
		TreeMap<Integer, E> chrElements = data.get(chr);
		if(chrElements == null){
			return null;
		} else {
			return chrElements.get(pos);
		}
	}
	
	/**
	 * 
	 * @param chr
	 * @param pos
	 * @return the removed element or null
	 */
	public E remove(String chr, Integer pos){
		
		TreeMap<Integer, E> chrElements = data.get(chr);
		if(chrElements == null){
			return null;
		} else {
			return chrElements.remove(pos);
		}
		
	}
	
	public Iterator<E> chrIterator(String chr) {
		
		TreeMap<Integer, E> chrElements = data.get(chr);
		if(chrElements == null){
			return null;
		} else {
			return chrElements.values().iterator();
		}
		
	}
	
	public Iterable<String> getChrs(){
		return data.keySet();
	}
	
	public Iterable<Integer> getChrPositions(String chr){
		TreeMap<Integer, E> chrElements = data.get(chr);
		if(chrElements == null) {
			return null;
		} else {
			return chrElements.keySet();
		}
	}
	
	
	public NavigableMap<Integer,E> getChrRange(String chr, Integer fromKey, boolean fromInclusive, Integer toKey, boolean toInclusive){
		
		TreeMap<Integer, E> chrElements = data.get(chr);
		if(chrElements == null){
			return Collections.emptyNavigableMap();
		} else {
			return chrElements.subMap(fromKey, fromInclusive, toKey, toInclusive);
		}
		
	}
	
	/**
	 * not tested
	 * 
	 * @return 
	 */
	public Iterator<ChrPos> getChrPosIterator() {
		
		Iterator[] chrIterators = new Iterator[data.size()];
		
		int i = 0;
		for(String chr : data.keySet()){
			
			chrIterators[i] = Iterators.transform(data.get(chr).keySet().iterator(), new CreateChrPos(chr));
		
			++i;
		}
		
		return Iterators.concat(chrIterators);
		
	}

	@Override
	public Iterator<E> iterator() {
		
		Iterator[] chrIterators = new Iterator[data.size()];
		
		int i = 0;
		for(TreeMap<Integer, E> chrData : data.values()){
			chrIterators[i] = chrData.values().iterator();
			++i;
		}
		return Iterators.concat(chrIterators);
		
	}
	
	public int size(){
		int count = 0;

		for(TreeMap<Integer, E> chrResults : data.values()){
			count += chrResults.size();
		}
		return count;
	}
	
	private class CreateChrPos implements Function<Integer, ChrPos>{

		private final String chr;

		public CreateChrPos(String chr) {
			this.chr = chr;
		}

		@Override
		public ChrPos apply(Integer input) {
			return new ChrPos(chr, input);
		}
		
	}
	
}
