/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.qcpca;

/**
 *
 * @author harmjan
 */
public class SortableSNP implements Comparable<SortableSNP> {
    public static enum SORTBY {NAME, ID, CHR, CHRPOS}

    SORTBY s;
    
    public byte chr      = -1;
    public int chrpos   = -1;
    public int id       = -1;
    public String name = null;
    
    
    public SortableSNP(String name, int id, byte chr, int chrpos, SORTBY s) {
        this.s = s;
	this.name = name;
        this.chr = chr;
        this.chrpos = chrpos;
        this.id = id;
    }
    
    public int compareTo(SortableSNP t) {
        int out = 0;
        switch(s){
            case ID:
                if(this.id >= t.id){
                    out = 1;
                } else {
                    out = -1;
                }
            break;
            case CHR:
                if(this.chr >= t.chr){
                    out = 1;
                } else {
                    out = -1;
                }
            break;
            case CHRPOS:
                if(this.chrpos >= t.chrpos){
                    out = 1;
                } else {
                    out = -1;
                }
            break;
        }
        return out;
    }
}
