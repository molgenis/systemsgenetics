/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.ucsc;

/**
 *
 * @author harmjan
 */
public class UCSCDataObject implements Comparable<UCSCDataObject> {

    private byte chr;
    private int positionStart;
    private int positionEnd;
    private double value;
    private SORTBY s = SORTBY.CHRPOS;

    

    public static enum SORTBY {

        CHRPOS, VALUE
    }

    public UCSCDataObject(byte chr, int posStart, int posEnd, double v, SORTBY sort) {
        this.chr = chr;
        this.positionStart = posStart;
        this.positionEnd = posEnd;
        this.value = v;
        this.s = sort;
    }

    @Override
    public int compareTo(UCSCDataObject o) {
        if (o == null) {
            throw new IllegalArgumentException("Error comparing UCSC data objects: including NULL objects in sort!");
        }
        if (s == SORTBY.CHRPOS) {
            if (this.equals(o)) {
                return 0;
            } else {
                if (this.chr > o.chr) {
                    return 1;
                } else if (this.chr < o.chr) {
                    return -1;
                } else if (this.positionStart > o.positionStart) {
                    return 1;
                } else {
                    return -1;
                }
            }
        } else {
            if (this.equals(o)) {
                return 0;
            } else {
                if (this.value > o.value) {
                    return 1;
                } else {
                    return -1;
                }
            }
        }
    }

    public boolean equals(UCSCDataObject o) {
        if (o == null) {
            return false;
        }

        if (o.s != this.s) {
            throw new IllegalArgumentException("Error: sort types are not identical!");
        }


        if (s == SORTBY.CHRPOS) {
            if (o.chr == this.chr && o.positionStart == positionStart && o.positionEnd == positionEnd) {
                return true;
            } else {
                return false;
            }
        } else {
            if (o.value == this.value) {
                return true;
            } else {
                return false;
            }
        }


    }

    @Override
    public String toString() {
        return "chr: " + chr + "\tstart: " + positionStart + "\tend: " + positionEnd + "\tvalue: " + value;
    }
    
    
    public byte getChr() {
        return chr;
    }

    public int getPositionEnd() {
        return positionEnd;
    }

    public int getPositionStart() {
        return positionStart;
    }
    
    public double getValue(){
        return value;
    }
}
