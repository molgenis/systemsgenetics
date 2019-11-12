/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.io.File;

/**
 *
 * @author patri
 */
public class PcaData {
	
	private final String name;
	private final String path;

	public PcaData(String name, String path) {
		this.name = name;
		this.path = path;
	}
	
	public File getEigenvectorPath(){
		return new File(path);
	}
	
}
