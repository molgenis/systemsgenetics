package org.molgenis.genotype.tabix;

import java.io.FileNotFoundException;

/**
 *
 * @author Patrick Deelen
 */
public class TabixFileNotFoundException extends FileNotFoundException {

	private String path;

	public TabixFileNotFoundException(String path, String s) {
		super(s);
		this.path = path;
	}

	public String getPath() {
		return path;
	}
	
	
}
