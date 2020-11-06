/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.gene;

import nl.systemsgenetics.downstreamer.summarystatistic.Locus;
import nl.systemsgenetics.downstreamer.summarystatistic.LocusUtils;
import nl.systemsgenetics.downstreamer.summarystatistic.OverlappableGenomicRange;
import umcg.genetica.collections.intervaltree.Range;

import java.util.Objects;

/**
 * @author patri
 */
public class Gene implements OverlappableGenomicRange, Range {

    private final String gene;
    private final String geneSymbol;
    private final String chr;
    private final int start;
    private final int end;
    private final String band;


    public Gene(String gene, String chr, int start, int stop, String band) {
        this.gene = gene;
        this.chr = chr.intern();
        this.start = start;
        this.end = stop;
        this.band = band;
        this.geneSymbol = null;
    }

    public Gene(String gene, String chr, int start, int stop, String band, String geneSymbol) {
        this.gene = gene;
        this.chr = chr.intern();
        this.start = start;
        this.end = stop;
        this.band = band;
        this.geneSymbol = geneSymbol;
    }

    public String getGene() {
        return gene;
    }

    public String getGeneSymbol() {
        return geneSymbol;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    @Override
    public String getSequenceName() {
        return getChr();
    }

    @Override
    public boolean isOverlapping(OverlappableGenomicRange other) {
        return LocusUtils.partialGenomicRangeOverlap(this, other);
    }

    @Override
    public boolean isOverlapping(OverlappableGenomicRange other, int window) {
        return LocusUtils.partialGenomicRangeOverlapWindow(this, other, window);
    }

    public String getBand() {
        return band;
    }

    public String getChrAndArm() {

        StringBuilder sb = new StringBuilder(chr);
        if (!this.band.equals("")) {
            sb.append('_');
            sb.append(band.charAt(0));
        }
        return sb.toString().intern();

    }

    public int getLength() {
        return getEnd() - getStart();
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 67 * hash + Objects.hashCode(this.gene);
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
        final Gene other = (Gene) obj;
        if (this.start != other.start) {
            return false;
        }
        if (this.end != other.end) {
            return false;
        }
        if (!Objects.equals(this.gene, other.gene)) {
            return false;
        }
        if (!Objects.equals(this.chr, other.chr)) {
            return false;
        }
        if (!Objects.equals(this.band, other.band)) {
            return false;
        }
        return true;
    }

}
