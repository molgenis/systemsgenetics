/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.containers;

/**
 *
 * @author harmjan
 */
public class SortableSNP implements Comparable<SortableSNP> {
    

    public static enum SORTBY {

        NAME, ID, CHR, CHRPOS, CHRANDCHRPOS, EFFECT
    }
    public SORTBY s;
    public byte chr = -1;
    public int chrpos = -1;
    public int id = -1;
    public String name = null;

    public final double effect;
    
    public SortableSNP(String name, int id, byte chr, int chrpos, SORTBY s) {
        this.s = s;
        this.name = name;
        this.chr = chr;
        this.chrpos = chrpos;
        this.id = id;
        this.effect = 0;
    }
    
    public SortableSNP(String name, int id, byte chr, int chrpos, double effect, SORTBY s) {
        this.s = s;
        this.name = name;
        this.chr = chr;
        this.chrpos = chrpos;
        this.id = id;
        this.effect = effect;
    }
    

    @Override
    public int compareTo(SortableSNP t) {
        int out = 0;
        switch (s) {
            case ID:
                if (this.id >= t.id) {
                    out = 1;
                } else {
                    out = -1;
                }
                break;
            case CHR:
                if (this.chr >= t.chr) {
                    out = 1;
                } else {
                    out = -1;
                }
                break;
            case CHRPOS:
                if (this.chrpos >= t.chrpos) {
                    out = 1;
                } else {
                    out = -1;
                }
                break;
            case CHRANDCHRPOS:
                if (this.chr > t.chr) {
                    out = 1;
                } else if (this.chr == t.chr) {
                    int output = 0;
                    if (this.chrpos > t.chrpos) {
                        output = 1;
                    } else {
                        output = -1;
                    }
                    if (this.chrpos > t.chrpos) {
                        return 1;
                    } else {
                        return -1;
                    }
                } else {
                    out = -1;
                }
                break;
            case EFFECT:
                if(this.effect > t.effect){
                    return 1;
                } else {
                    return -1;
                }
        }
        return out;
    }
}
