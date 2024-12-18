/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend;

import htsjdk.samtools.util.Interval;

import java.util.Objects;

/**
 * @author patri
 */
public class Gene extends Interval {

    private final String gene;
    private final String geneSymbol;
    private final String band;

    public Gene(String gene, String chr, int start, int stop, String band) {
        super(chr.intern(), start, stop);
        this.gene = gene;
        this.band = band;
        this.geneSymbol = null;
    }

    public Gene(String gene, String chr, int start, int stop, String band, String geneSymbol) {
        super(chr.intern(), start, stop);
        this.gene = gene;
        this.band = band;
        this.geneSymbol = geneSymbol;
    }

    public String getGene() {
        return gene;
    }

    public String getGeneSymbol() {
        return geneSymbol;
    }

    public String getBand() {
        return band;
    }

    public String getChrAndArm() {
        StringBuilder sb = new StringBuilder(getContig());
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
        if (getStart() != other.getStart()) {
            return false;
        }
        if (getStart() != other.getEnd()) {
            return false;
        }
        if (!Objects.equals(this.gene, other.gene)) {
            return false;
        }
        if (!Objects.equals(getContig(), other.getContig())) {
            return false;
        }
        if (!Objects.equals(this.band, other.band)) {
            return false;
        }
        return true;
    }

}
