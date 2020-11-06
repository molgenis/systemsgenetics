package nl.systemsgenetics.downstreamer.summarystatistic;

import java.io.Serializable;

/**
 * The type Snp.
 */
public class SNP implements Serializable, OverlappableGenomicRange {

    /**
     * The primary variant id.
     */
    protected String primaryVariantId;
    /**
     * The Allele 1.
     */
    protected String allele1;
    /**
     * The Allele 2.
     */
    protected String allele2;
    /**
     * The Position.
     */
    protected int position;
    /**
     * The Sequence name.
     */
    protected String sequenceName;

    /**
     * Gets rs id.
     *
     * @return the rs id
     */
    public String getPrimaryVariantId() {
        return primaryVariantId;
    }

    /**
     * Sets rs id.
     *
     * @param primaryVariantId the rs id
     */
    public void setPrimaryVariantId(String primaryVariantId) {
        this.primaryVariantId = primaryVariantId;
    }

    /**
     * Gets allele 1.
     *
     * @return the allele 1
     */
    public String getAllele1() {
        return allele1;
    }

    /**
     * Sets allele 1.
     *
     * @param allele1 the allele 1
     */
    public void setAllele1(String allele1) {
        this.allele1 = allele1;
    }

    /**
     * Gets allele 2.
     *
     * @return the allele 2
     */
    public String getAllele2() {
        return allele2;
    }

    /**
     * Sets allele 2.
     *
     * @param allele2 the allele 2
     */
    public void setAllele2(String allele2) {
        this.allele2 = allele2;
    }

    /**
     * Gets position.
     *
     * @return the position
     */
    public int getPosition() {
        return position;
    }

    /**
     * Sets position.
     *
     * @param position the position
     */
    public void setPosition(int position) {
        this.position = position;
    }

    /**
     * Gets human chromosome.
     *
     * @return the human chromosome
     */
    public int getHumanChromosome() {

        if (sequenceName.toLowerCase().equals("x")) {
            return (23);
        } else if (sequenceName.toLowerCase().equals("y")) {
            return (24);
        } else if (sequenceName.toLowerCase().equals("mt")) {
            return (25);
        } else {
            String curSeq = sequenceName;
            return Integer.parseInt(curSeq.replace("chr", ""));
        }
    }

    /**
     * Sets sequence name.
     *
     * @param sequenceName the sequence name
     */
    public void setSequenceName(String sequenceName) {
        this.sequenceName = sequenceName;
    }

    @Override
    public int getStart() {
        return position;
    }

    @Override
    public int getEnd() {
        return position;
    }

    /**
     * Get sequence name string.
     *
     * @return the string
     */
    public String getSequenceName(){return sequenceName;}

    @Override
    public boolean isOverlapping(OverlappableGenomicRange other) {
        return LocusUtils.partialGenomicRangeOverlap(this, other);
    }

    @Override
    public boolean isOverlapping(OverlappableGenomicRange other, int window) {
        return LocusUtils.partialGenomicRangeOverlapWindow(this, other, window);
    }

    /**
     * Is transition boolean.
     *
     * @return the boolean
     */
    public boolean isTransition() {

        // Check if a variant is a transition
        if (allele1.equals("C") && allele2.equals("T") || allele2.equals("T") && allele1.equals("C")){
            return true;
        }

        if (allele1.equals("A") && allele2.equals("G") || allele2.equals("A") && allele1.equals("G")) {
            return true;
        }

        return false;

    }
}
