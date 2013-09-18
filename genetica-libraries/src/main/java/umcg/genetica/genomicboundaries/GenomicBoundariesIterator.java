package umcg.genetica.genomicboundaries;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;

/**
 *
 * @author PatrickDeelen
 */
public class GenomicBoundariesIterator<V> implements Iterator<GenomicBoundary<V>>{

	private Iterator<TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> chromosomeIterator;
	private Iterator<ArrayList<GenomicBoundary<V>>> currentChromosomePositionBoundariesIterator;
	private Iterator<GenomicBoundary<V>> currentChromosomeBoundariesIterator;
	private ArrayList<GenomicBoundary<V>> currentChromosomePosBoundaries;
	
	private boolean end = false;
	
	
	public GenomicBoundariesIterator(HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> genomicsBoundaries) {
		chromosomeIterator = genomicsBoundaries.values().iterator();
		
		if(chromosomeIterator.hasNext()){
			currentChromosomePositionBoundariesIterator = chromosomeIterator.next().values().iterator();
			
			if(currentChromosomePositionBoundariesIterator.hasNext()){
				currentChromosomePosBoundaries = currentChromosomePositionBoundariesIterator.next();
				currentChromosomeBoundariesIterator = currentChromosomePosBoundaries.iterator();
			}
		} else {
			end = true;
		}
	}
	
	public boolean hasNext() {

		if(end){
			return false;
		}

		if(currentChromosomeBoundariesIterator.hasNext()){
			return true;
		} else {
			
			if(currentChromosomePositionBoundariesIterator.hasNext()){
				currentChromosomePosBoundaries = currentChromosomePositionBoundariesIterator.next();
				currentChromosomeBoundariesIterator = currentChromosomePosBoundaries.iterator();
				if(currentChromosomeBoundariesIterator.hasNext()){
					return true;
				}
			}
			
			else{
				while(chromosomeIterator.hasNext()){
					currentChromosomePositionBoundariesIterator = chromosomeIterator.next().values().iterator();

					if(currentChromosomePositionBoundariesIterator.hasNext()){
						currentChromosomePosBoundaries = currentChromosomePositionBoundariesIterator.next();
						currentChromosomeBoundariesIterator = currentChromosomePosBoundaries.iterator();
						if(currentChromosomeBoundariesIterator.hasNext()){
							return true;
						}
					}


				}
			}
		}
		end = true;
		return false;
	}
	
	public GenomicBoundary<V> next() {
		if(currentChromosomeBoundariesIterator.hasNext()){
			return currentChromosomeBoundariesIterator.next();
		} else {
			
			if(currentChromosomePositionBoundariesIterator.hasNext()){
				currentChromosomePosBoundaries = currentChromosomePositionBoundariesIterator.next();
				currentChromosomeBoundariesIterator = currentChromosomePosBoundaries.iterator();
				if(currentChromosomeBoundariesIterator.hasNext()){
					return currentChromosomeBoundariesIterator.next();
				}
			}
			else{
				while(chromosomeIterator.hasNext()){

					currentChromosomePositionBoundariesIterator = chromosomeIterator.next().values().iterator();

					if(currentChromosomePositionBoundariesIterator.hasNext()){
						currentChromosomePosBoundaries = currentChromosomePositionBoundariesIterator.next();
						currentChromosomeBoundariesIterator = currentChromosomePosBoundaries.iterator();
						if(currentChromosomeBoundariesIterator.hasNext()){
							return currentChromosomeBoundariesIterator.next();
						}
					}
				}
			}
		}
		throw new NoSuchElementException();
	}

	public void remove() {
		
		currentChromosomeBoundariesIterator.remove();
		if(currentChromosomePosBoundaries.isEmpty()){
			currentChromosomePositionBoundariesIterator.remove();
		}
		
	}

}
