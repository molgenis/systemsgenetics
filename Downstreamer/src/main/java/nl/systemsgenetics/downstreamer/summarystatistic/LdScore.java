package nl.systemsgenetics.downstreamer.summarystatistic;

public class LdScore extends SNP {


    private double centimorgan;
    private double ldscore;

    public LdScore(String sequence, int pos, String variantId, double ldScore) {
        super(sequence, pos, pos);
        this.primaryVariantId = variantId;
        this.ldscore = ldScore;
    }


    public double getCentimorgan() {
        return centimorgan;
    }

    public void setCentimorgan(double centimorgan) {
        this.centimorgan = centimorgan;
    }

    public double getLdscore() {
        return ldscore;
    }

    public void setLdscore(double ldscore) {
        this.ldscore = ldscore;
    }
}
