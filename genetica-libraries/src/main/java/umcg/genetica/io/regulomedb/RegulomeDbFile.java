/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.regulomedb;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;

/**
 *
 * @author Patrick Deelen
 */
public class RegulomeDbFile implements Iterable<RegulomeDbEntry> {
	
	private final File regulomeDbFile;

	public RegulomeDbFile(File regulomeDbFile) throws FileNotFoundException, IOException {
		if(!regulomeDbFile.exists()){
			throw new FileNotFoundException("RegulomeDB file not found: " + regulomeDbFile.getAbsolutePath());
		}
		if(!regulomeDbFile.canRead()){
			throw new IOException("Can not read RegulomeDB file: " + regulomeDbFile.getAbsolutePath());
		}
		this.regulomeDbFile = regulomeDbFile;
	}

	@Override
	public Iterator<RegulomeDbEntry> iterator() {
		return new RegulomeDbIterator(this);
	}

	public File getRegulomeDbFile() {
		return regulomeDbFile;
	}
	
	
	
}
