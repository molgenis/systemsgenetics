/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gtf;

import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 *
 * @author Patrick Deelen
 */
public class GffElementIterator implements Iterator<GffElement> {

	private final GffReaderInterface gffReader;
	
	private boolean atNext = false;
	private boolean atEnd = false;
	private GffElement next;
	
	public GffElementIterator(GffReaderInterface gffReader) {
		this.gffReader = gffReader;
	}

	@Override
	public boolean hasNext(){
		
		if(atEnd){
			return false;
		}
		
		if(atNext){
			return true;
		}
		
		try {
			next = gffReader.nextElement();
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
		
		if(next == null){
			atEnd = true;
			try {
				gffReader.close();
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
			return false;
		} else {
			atNext = true;
			return true;
		}
		
	}

	@Override
	public GffElement next() throws NoSuchElementException {
		if(!hasNext()){
			throw new NoSuchElementException();
		}
		atNext = false;
		return next;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not supported ever.");
	}
	
}
