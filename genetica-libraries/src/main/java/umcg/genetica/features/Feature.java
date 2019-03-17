/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.features;


import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;

import java.util.ArrayList;
import java.util.Objects;

/**
 * @author hwestra
 */
public class Feature {
	
	protected Feature parent;
	
	protected Chromosome chromosome;
	protected String name;
	protected Strand strand;
	protected int start = Integer.MAX_VALUE;
	protected int stop = -Integer.MAX_VALUE;
	protected int nrAT;
	protected int nrGC;
	protected int nrN;
	
	
	public Feature(Chromosome chromosome, int alignmentStart, int alignmentEnd) {
		this.chromosome = chromosome;
		this.start = alignmentStart;
		this.stop = alignmentEnd;
	}
	
	public Feature() {
	
	}
	
	public Feature(Feature feature) {
		this.chromosome = feature.chromosome;
		this.name = feature.getName();
		this.start = feature.start;
		this.stop = feature.stop;
		this.nrAT = feature.nrAT;
		this.nrGC = feature.nrGC;
		this.nrN = feature.nrN;
	}
	
	
	public static Feature parseFeature(String regionStr) {
		String[] elems = regionStr.split("_");
		if (elems.length == 2) {
			Chromosome chr = Chromosome.parseChr(elems[0]);
			String posStr = elems[1];
			String[] posStrElems = posStr.split("-");
			if (posStrElems.length == 2) {
				int start = Integer.parseInt(posStrElems[0]);
				int stop = Integer.parseInt(posStrElems[1]);
				return new Feature(chr, start, stop);
			} else {
				return null;
			}
		} else {
			return null;
		}
	}
	
	public Feature getParent() {
		return parent;
	}
	
	public void setParent(Feature p) {
		this.parent = p;
	}
	
	public Chromosome getChromosome() {
		return chromosome;
	}
	
	public void setChromosome(Chromosome chromosome) {
		this.chromosome = chromosome;
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public Strand getStrand() {
		return strand;
	}
	
	public void setStrand(Strand strand) {
		this.strand = strand;
	}
	
	public int getStart() {
		return start;
	}
	
	public void setStart(int start) {
		this.start = start;
	}
	
	public int getStop() {
		return stop;
	}
	
	public void setStop(int stop) {
		this.stop = stop;
	}
	
	@Override
	public int hashCode() {
		int hash = 5;
		hash = 67 * hash + Objects.hashCode(this.chromosome);
		if (useNameForEquals) {
			hash = 67 * hash + Objects.hashCode(this.name);
		}
		hash = 67 * hash + Objects.hashCode(this.strand);
		hash = 67 * hash + this.start;
		hash = 67 * hash + this.stop;
		return hash;
	}
	
	boolean useNameForEquals = true;
	
	public void useNameForComparison(boolean b) {
		useNameForEquals = b;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final Feature other = (Feature) obj;
		if (this.chromosome != other.chromosome) {
			return false;
		}
		if (useNameForEquals) {
			if (this.name != null && other.name != null) {
				if (!this.name.equals(other.name)) {
					return false;
				}
			}
		}
		
		if (this.strand != other.strand) {
			return false;
		}
		if (this.start != other.start) {
			return false;
		}
		if (this.stop != other.stop) {
			return false;
		}
		return true;
	}
	
	public boolean overlaps(Feature that) {
		if (this.equals(that)) {
			return true;
			
		}
		if (this.chromosome == null || that.chromosome == null) {
			return false;
		}
		if (this.chromosome.getNumber() != that.chromosome.getNumber()) {
//            System.out.println("chr diff");
			return false;
		}
		
		if (this.start >= that.start && this.stop <= that.stop) {
			return true;
			//   this     |-|
			//   that    |---|
			
		}
		if (that.start >= this.start && that.stop <= this.stop) {
			return true;
			//   this    |---|
			//   that     |-|
		}
		
		if (this.start >= that.start && this.start < that.stop) {
			return true;
			//   this      |---|
			//   that    |---|
		}
		if (this.start < that.start && this.stop > that.start) {
			return true;
			//   this  |---|
			//   that    |---|
		}
		
		
		return false;
	}
	
	public int determineBaseOverlap(Feature feat2) {
		
		if (!this.overlaps(feat2)) {
			return 0;
		}
		
		int feat2start = feat2.getStart();
		int feat2stop = feat2.getStop();
		int feat1start = this.getStart();
		int feat1stop = this.getStop();
		int overlap = 0;
		if (feat1start >= feat2start && feat1stop <= feat2stop) {
			// the overlap is the size of feat1
			overlap = feat1stop - feat1start;
		} else if (feat2start >= feat1start && feat2stop <= feat1stop) {
			// feat1   ----------
			// feat2     ------
			// the overlap is the size of feat2
			overlap = feat2stop - feat2start;
		} else if (feat1start >= feat2start) {
			// feat1       ----------
			// feat2     ------
			// overlap starts at feat1start, ends at feat2end
			overlap = feat2stop - feat1start;
		} else if (feat2start > feat1start) {
			// feat1    ----------
			// feat2            ------
			// overlap starts at feat2start, ends at feat1end
			
			overlap = feat1stop - feat2start;
		}
		
		return overlap;
	}
	
	public void setBaseProperties(int nrAT, int nrGC, int nrN) {
		this.nrAT = nrAT;
		this.nrGC = nrGC;
		this.nrN = nrN;
	}
	
	public int getNrAT() {
		return nrAT;
	}
	
	public void setNrAT(int nrAT) {
		this.nrAT = nrAT;
	}
	
	public int getNrGC() {
		return nrGC;
	}
	
	public void setNrGC(int nrGC) {
		this.nrGC = nrGC;
	}
	
	public int getNrN() {
		return nrN;
	}
	
	public void setNrN(int nrN) {
		this.nrN = nrN;
	}
	
	public double getGCContent() {
		return (double) nrGC / (nrGC + nrAT + nrN);
	}
	
	public int getBaseSum() {
		return (nrGC + nrAT + nrN);
	}
	
	public String toBedString() {
		return getChromosome().getName() + "\t" + getStart() + "\t" + getStop();
	}
	
	@Override
	public String toString() {
		return getChromosome().getName() + "_" + getStart() + "-" + getStop();
	}
	
	
	public int getSize() {
		return stop - start;
	}
	
	public boolean overlaps(ArrayList<Feature> regions) {
		if (regions == null) {
			System.out.println("QQC?");
			
			return false;
		}
		for (Feature f : regions) {
			if (this.overlaps(f)) {
				return true;
			}
		}
		return false;
	}
	
	public Feature newFeatureFromCoordinates() {
		return new Feature(this.chromosome, this.start, this.stop);
	}
}
