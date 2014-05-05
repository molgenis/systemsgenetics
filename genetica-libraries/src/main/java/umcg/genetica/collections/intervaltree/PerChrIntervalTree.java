package umcg.genetica.collections.intervaltree;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author Patrick Deelen
 */
public class PerChrIntervalTree <E extends Range> {

	private final Class<E> classE;
	private final HashMap<String, IntervalTree<E>> perChrIntervalTree;

	public PerChrIntervalTree(Class<E> classE){
		this.classE = classE;
		this.perChrIntervalTree = new HashMap<String, IntervalTree<E>>();
	}
	
	public void addChrElements(String chr, List<E> elements) throws Exception{
		if(perChrIntervalTree.containsKey(chr)){
			throw new Exception("Cannot create " + chr + " interval tree. Already exisits. Appending not possible");
		}
		perChrIntervalTree.put(chr, new IntervalTree<E>(elements, classE));
	}
	
	public List<E> searchPosition(String chr, int pos){
		IntervalTree<E> chrIntervalTree = perChrIntervalTree.get(chr);
		if(chrIntervalTree == null){
			return Collections.emptyList();
		}
		return chrIntervalTree.getElementsOverlappingQuery(pos);
	}
	
	public int size(){
		int size = 0;
		for(IntervalTree<E> chrIntervalTree : perChrIntervalTree.values()){
			size += chrIntervalTree.size();
		}
		return size;
	}
	
}
