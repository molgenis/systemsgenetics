/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlannotation;

/**
 *
 * This format is used to provide called peaks of signal enrichment based on
 * pooled, normalized (interpreted) data. It is a BED6+4 format.
 *
 * chrom - Name of the chromosome (or contig, scaffold, etc.). chromStart - The
 * starting position of the feature in the chromosome or scaffold. The first
 * base in a chromosome is numbered 0. chromEnd - The ending position of the
 * feature in the chromosome or scaffold. The chromEnd base is not included in
 * the display of the feature. For example, the first 100 bases of a chromosome
 * are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
 * name - Name given to a region (preferably unique). Use '.' if no name is
 * assigned. score - Indicates how dark the peak will be displayed in the
 * browser (0-1000). If all scores were '0' when the data were submitted to the
 * DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the
 * average signalValue per base spread is between 100-1000. strand - +/- to
 * denote strand or orientation (whenever applicable). Use '.' if no orientation
 * is assigned. signalValue - Measurement of overall (usually, average)
 * enrichment for the region. pValue - Measurement of statistical significance
 * (-log10). Use -1 if no pValue is assigned. qValue - Measurement of
 * statistical significance using false discovery rate (-log10). Use -1 if no
 * qValue is assigned. peak - Point-source called for this peak; 0-based offset
 * from chromStart. Use -1 if no point-source called.
 *
 * @author MarcJan
 */
public class Bed implements Comparable<Bed> {

    final private String chrom;
    final private int chromStart;
    final private int chromEnd;


    Bed(String chr, int start, int stop) {
        chrom = chr;
        chromStart = start;
        chromEnd = stop;

    }

    Bed(String[] row) {
        chrom = row[0];
        chromStart = Integer.parseInt(row[1]);
        chromEnd = Integer.parseInt(row[2]);
    }


    public int getNumChrom() {
        return getNumericChr(chrom);
    }

    public String getChrom() {
        return chrom;
    }

    public int getChromStart() {;
        return chromStart;
    }

    public int getChromEnd() {
        return chromEnd;
    }

    @Override
    public int compareTo(Bed other) {
        if (other.getNumChrom() > this.getNumChrom()) {
            return -1;
        } else if (other.getNumChrom() < this.getNumChrom()) {
            return 1;
        } else {
            if (other.getChromStart() > this.getChromStart()) {
                return -1;
            } else if (other.getChromStart() < this.getChromStart()) {
                return 1;
            } else {
                if (other.getChromEnd() > this.getChromEnd()) {
                    return -1;
                } else if (other.getChromEnd() < this.getChromEnd()) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }
    }

    private int getNumericChr(String chrom) {
        String relevant = chrom.replace("chr", "");
        
        if(relevant.equals("X")){
            return(23);
        } else if(relevant.equals("Y")) {
            return(24);
        } else {
            return(Integer.parseInt(relevant));
        }
    }

}
