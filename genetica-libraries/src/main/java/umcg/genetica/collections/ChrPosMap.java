package umcg.genetica.collections;

import com.google.common.collect.Iterators;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashMap;
import java.util.Iterator;

/**
 *
 * @author Patrick Deelen
 */
public class ChrPosMap<E> implements Iterable<E>{

	private final HashMap<String, TIntObjectMap<E>> data;

	public ChrPosMap() {
		data = new HashMap<String, TIntObjectMap<E>>();
	}
	
	public void put(String chr, int pos, E element){
		TIntObjectMap<E> chrElements = data.get(chr);
		if(chrElements == null){
			chrElements = new TIntObjectHashMap<E>();
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
	public E get(String chr, int pos){
		TIntObjectMap<E> chrElements = data.get(chr);
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
	public E remove(String chr, int pos){
		
		TIntObjectMap<E> chrElements = data.get(chr);
		if(chrElements == null){
			return null;
		} else {
			return chrElements.remove(pos);
		}
		
	}
	
	public Iterator<E> chrIterator(String chr) {
		
		TIntObjectMap<E> chrElements = data.get(chr);
		if(chrElements == null){
			return null;
		} else {
			return chrElements.valueCollection().iterator();
		}
		
	}

	@Override
	public Iterator<E> iterator() {
		
		Iterator[] chrIterators = new Iterator[data.size()];
		
		int i = 0;
		for(TIntObjectMap<E> chrData : data.values()){
			chrIterators[i] = chrData.valueCollection().iterator();
			++i;
		}
		return Iterators.concat(chrIterators);
		
	}
	
	public int size(){
		int count = 0;

		for(TIntObjectMap<E> chrResults : data.values()){
			count += chrResults.size();
		}
		return count;
	}
	
}
