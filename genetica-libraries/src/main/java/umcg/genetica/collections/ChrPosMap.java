package umcg.genetica.collections;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashMap;

/**
 *
 * @author Patrick Deelen
 */
public class ChrPosMap<E> {

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
	
}
