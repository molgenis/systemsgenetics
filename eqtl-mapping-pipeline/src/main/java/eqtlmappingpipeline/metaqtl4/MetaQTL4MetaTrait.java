/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl4;

import org.molgenis.genotype.util.ChromosomeComparator;

/**
 *
 * @author harmjan
 */
public class MetaQTL4MetaTrait implements Comparable<Object> {

    static boolean mapsToChromosome() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private final int metaTraitId;
    private final String metaTraitName;
    private final String chr;
    private final int chrStart;
    private final int chrEnd;
    private final int chrMidpoint;
    private final String annotation;
    private final String[] platformIds;
    private static final ChromosomeComparator chrComparator = new ChromosomeComparator();
    private int currentMetaId;

    public MetaQTL4MetaTrait(int metaTraitId, String metaTraitName, String chr, int chrStart, int chrEnd, String annotation, String[] platformIds) {
        this.metaTraitId = metaTraitId;
        this.metaTraitName = metaTraitName;
        this.chr = chr;
        this.chrStart = chrStart;
        this.chrEnd = chrEnd;
        if (chrStart >= 0 && chrEnd > 1 && chrStart != chrEnd) {
            chrMidpoint = Math.abs((chrStart + chrEnd) / 2);
        } else if (chrStart == chrEnd) {
            chrMidpoint = chrStart;
        } else {
            chrMidpoint = -1;
        }

        this.annotation = annotation;
        this.platformIds = platformIds;
    }

    @Override
    public int compareTo(Object t) {
        if (this.equals(t)) {
            return 0;
        } else {
            MetaQTL4MetaTrait that = (MetaQTL4MetaTrait) t;
            if (!this.getChr().equals(that.getChr())) {
                return chrComparator.compare(chr, that.getChr());
            } else {
                if (this.getChrMidpoint() == that.getChrMidpoint()) {
                    //Assume id is unique
                    return this.metaTraitId - that.metaTraitId;
                } else {
                    return this.getChrMidpoint() - that.getChrMidpoint();
                }
            }
        }
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 59 * hash + this.metaTraitId;
        hash = 59 * hash + (this.chr != null ? this.chr.hashCode() : 0);
        hash = 59 * hash + this.chrMidpoint;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final MetaQTL4MetaTrait other = (MetaQTL4MetaTrait) obj;
        if (this.metaTraitId != other.metaTraitId) {
            return false;
        }
        if ((this.chr == null) ? (other.chr != null) : !this.chr.equals(other.chr)) {
            return false;
        }
        if (this.chrMidpoint != other.chrMidpoint) {
            return false;
        }
        return true;
    }

    public int getMetaTraitId() {
        return metaTraitId;
    }

    public String getMetaTraitName() {
        return metaTraitName;
    }

    public String getChr() {
        return chr;
    }

    public int getChrStart() {
        return chrStart;
    }

    public int getChrEnd() {
        return chrEnd;
    }

    public int getChrMidpoint() {
        return chrMidpoint;
    }

    public String getAnnotation() {
        return annotation;
    }

    public String[] getPlatformIds() {
        return platformIds;
    }

    void setMetaTraitId(int id) {
        this.currentMetaId = id;
    }
    public int getCurrentMetaId(){
        return currentMetaId;
    }
}
