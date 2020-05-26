/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.regulomedb;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbFiles implements Iterable<RegulomeDbEntry> {
	
	private ArrayList<RegulomeDbFile> regulomeDbFiles;

	public RegulomeDbFiles(ArrayList<RegulomeDbFile> regulomeDbFiles) {
		this.regulomeDbFiles = regulomeDbFiles;
	}

	public RegulomeDbFiles(Collection<File> regulomeDbFiles) throws FileNotFoundException, IOException {
		this.regulomeDbFiles = new ArrayList<RegulomeDbFile>(regulomeDbFiles.size());
		for(File regulomeDbFile : regulomeDbFiles){
			this.regulomeDbFiles.add(new RegulomeDbFile(regulomeDbFile));
		}
	}
	
	public void addRegulomeDbFile(RegulomeDbFile regulomeDbFile){
		regulomeDbFiles.add(regulomeDbFile);
	}

	@Override
	public Iterator<RegulomeDbEntry> iterator() {
		return new RegulomeDbIterator(regulomeDbFiles);
	}
	
	
	
}
