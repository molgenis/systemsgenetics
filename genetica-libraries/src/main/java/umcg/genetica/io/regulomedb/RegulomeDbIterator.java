/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.regulomedb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbIterator implements Iterator<RegulomeDbEntry> {
	
	private Iterator<RegulomeDbFile> regulomeDbFilesIterator;
	private BufferedReader regulomeDbFileReader = null;
	private boolean hasNext;
	private RegulomeDbEntry next;
	
	public RegulomeDbIterator(RegulomeDbFile regulomeDbFile) {
		
		ArrayList<RegulomeDbFile> regulomeDbFiles = new ArrayList<RegulomeDbFile>(1);
		regulomeDbFiles.add(regulomeDbFile);
		regulomeDbFilesIterator = regulomeDbFiles.iterator();
		
		try {
			
			if(regulomeDbFilesIterator.hasNext()){
				regulomeDbFileReader = new BufferedReader(new FileReader(regulomeDbFilesIterator.next().getRegulomeDbFile()));
				loadNext();
			} else {
				hasNext = false;
				next = null;
			}
			
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
		
	}

	public RegulomeDbIterator(Collection<RegulomeDbFile> regulomeDbFiles) {
		
		this.regulomeDbFilesIterator = regulomeDbFiles.iterator();
		
		try {
			
			if(regulomeDbFilesIterator.hasNext()){
				regulomeDbFileReader = new BufferedReader(new FileReader(regulomeDbFilesIterator.next().getRegulomeDbFile()));
				loadNext();
			} else {
				hasNext = false;
				next = null;
			}
			
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
		
	}
	
	

	@Override
	public boolean hasNext() {
		return hasNext;
	}

	@Override
	public RegulomeDbEntry next() {
		
		RegulomeDbEntry current = next;
		
		try {
			loadNext();
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
		
		return current;
		
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not supported ever.");
	}

	private void loadNext() throws IOException,Exception{
		
		String line = regulomeDbFileReader.readLine();
		
		if(line == null){
			while(regulomeDbFilesIterator.hasNext()){
				regulomeDbFileReader = new BufferedReader(new FileReader(regulomeDbFilesIterator.next().getRegulomeDbFile()));
				line = regulomeDbFileReader.readLine();
				if(line == null){
					continue;
				} else {
					hasNext = true;
					next = new RegulomeDbEntry(line);
					return;
				}
			}
			hasNext = false;
			next = null;
		} else {
			hasNext = true;
			next = new RegulomeDbEntry(line);
		}
		
	}
	
}
