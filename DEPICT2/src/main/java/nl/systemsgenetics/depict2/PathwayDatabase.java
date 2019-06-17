/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import java.util.Objects;

/**
 *
 * @author patri
 */
public class PathwayDatabase {
	
	private final String name;
	private final String location;
	private final boolean textBasedMatrix;//Only used for testing

	public PathwayDatabase(String name, String location) {
		this.name = name;
		this.location = location;
		this.textBasedMatrix = false;
	}

	protected PathwayDatabase(String name, String location, boolean textBasedMatrix) {
		this.name = name;
		this.location = location;
		this.textBasedMatrix = textBasedMatrix;
	}	

	public String getName() {
		return name;
	}

	public String getLocation() {
		return location;
	}

	protected boolean isTextBasedMatrix() {
		return textBasedMatrix;
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
	
	
	
}
