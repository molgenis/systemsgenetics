/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.collections.intervaltree;

import java.util.Objects;
import umcg.genetica.variantAnnotator.GenomicRange;

/**
 *
 * @author patri
 */
public class NamedGenomicRange implements GenomicRange{

	private final String name;
	private final String sequence;
	private final int start;
	private final int end;

	public NamedGenomicRange(String name, String sequence, int start, int end) {
		this.name = name;
		this.sequence = sequence;
		this.start = start;
		this.end = end;
	}
	
	@Override
	public String getSeqname() {
		return sequence;
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getEnd() {
		return end;
	}

	public String getName() {
		return name;
	}

	@Override
	public int hashCode() {
		int hash = 3;
		hash = 59 * hash + Objects.hashCode(this.name);
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
		final NamedGenomicRange other = (NamedGenomicRange) obj;
		if (this.start != other.start) {
			return false;
		}
		if (this.end != other.end) {
			return false;
		}
		if (!Objects.equals(this.name, other.name)) {
			return false;
		}
		if (!Objects.equals(this.sequence, other.sequence)) {
			return false;
		}
		return true;
	}

	
	
}
