/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.pathway;

import java.io.File;
import java.util.Objects;

/**
 *
 * @author patri
 */
public class PathwayDatabase {
	
	private final String name;
	private final String location;

	public PathwayDatabase(String name, String location) {
		this.name = name;
		this.location = location;
	}
	
	public String getName() {
		return name;
	}

	public String getLocation() {
		return location;
	}

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 53 * hash + Objects.hashCode(this.name);
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final PathwayDatabase other = (PathwayDatabase) obj;
		if (!Objects.equals(this.name, other.name)) {
			return false;
		}
		if (!Objects.equals(this.location, other.location)) {
			return false;
		}
		return true;
	}
	
	public boolean exist(){
		return new File(location + ".dat").canRead();
	}
	
}
