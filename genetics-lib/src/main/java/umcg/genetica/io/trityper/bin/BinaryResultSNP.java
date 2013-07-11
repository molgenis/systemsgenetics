/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper.bin;

/**
 *
 * @author harmjan
 */
public class BinaryResultSNP {
    private String name;
    private Byte chr;
    private Integer chrpos;
    private Integer numsamples;
    private Integer id;
    private Double hwe;
    private Double maf;
    private Double cr;
    private Byte minorAllele;
    private Long zScoreIndex;
    private byte[] alleles;
    private Byte assessedAllele;

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the chr
     */
    public Byte getChr() {
        return chr;
    }

    /**
     * @param chr the chr to set
     */
    public void setChr(Byte chr) {
        this.chr = chr;
    }

    /**
     * @return the chrpos
     */
    public Integer getChrpos() {
        return chrpos;
    }

    /**
     * @param chrpos the chrpos to set
     */
    public void setChrpos(Integer chrpos) {
        this.chrpos = chrpos;
    }

    /**
     * @return the numsamples
     */
    public Integer getNumsamples() {
        return numsamples;
    }

    /**
     * @param numsamples the numsamples to set
     */
    public void setNumsamples(Integer numsamples) {
        this.numsamples = numsamples;
    }

    /**
     * @return the id
     */
    public Integer getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(Integer id) {
        this.id = id;
    }

    /**
     * @return the hwe
     */
    public Double getHwe() {
        return hwe;
    }

    /**
     * @param hwe the hwe to set
     */
    public void setHwe(Double hwe) {
        this.hwe = hwe;
    }

    /**
     * @return the maf
     */
    public Double getMaf() {
        return maf;
    }

    /**
     * @param maf the maf to set
     */
    public void setMaf(Double maf) {
        this.maf = maf;
    }

    /**
     * @return the cr
     */
    public Double getCr() {
        return cr;
    }

    /**
     * @param cr the cr to set
     */
    public void setCr(Double cr) {
        this.cr = cr;
    }

    /**
     * @return the minorAllele
     */
    public Byte getMinorAllele() {
        return minorAllele;
    }

    /**
     * @param minorAllele the minorAllele to set
     */
    public void setMinorAllele(Byte minorAllele) {
        this.minorAllele = minorAllele;
    }

    /**
     * @return the zScoreIndex
     */
    public Long getzScoreIndex() {
        return zScoreIndex;
    }

    /**
     * @param zScoreIndex the zScoreIndex to set
     */
    public void setzScoreIndex(Long zScoreIndex) {
        this.zScoreIndex = zScoreIndex;
    }

    /**
     * @return the alleles
     */
    public byte[] getAlleles() {
        return alleles;
    }

    /**
     * @param alleles the alleles to set
     */
    public void setAlleles(byte[] alleles) {
        this.alleles = alleles;
    }

    /**
     * @return the assessedAllele
     */
    public Byte getAssessedAllele() {
        return assessedAllele;
    }

    /**
     * @param assessedAllele the assessedAllele to set
     */
    public void setAssessedAllele(Byte assessedAllele) {
        this.assessedAllele = assessedAllele;
    }

    public void clearMetaData(){
	name = null;
	chr = null;
	chrpos = null;
    }
}
