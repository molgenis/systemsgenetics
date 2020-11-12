package nl.systemsgenetics.downstreamer.gene;

public class GenePathwayAssociationStatistic {
    private String pathway;
    private String gene;
    private double beta;
    private double standardError;
    private double tStatistic;
    private double zscore;
    private double analyticalPvalue;
    private double empiricalPvalue;


    public GenePathwayAssociationStatistic(double beta, double standardError, double tStatistic, double zscore, double analyticalPvalue) {
        this.beta = beta;
        this.standardError = standardError;
        this.tStatistic = tStatistic;
        this.analyticalPvalue = analyticalPvalue;
        this.zscore = zscore;
    }

    public GenePathwayAssociationStatistic(String pathway, String gene, double beta, double standardError, double tStatistic, double zscore, double analyticalPvalue) {
        this.pathway = pathway;
        this.gene = gene;
        this.beta = beta;
        this.standardError = standardError;
        this.tStatistic = tStatistic;
        this.zscore = zscore;
        this.analyticalPvalue = analyticalPvalue;
    }

/*    public GenePathwayAssociationStatistic(String pathway, String gene, double beta, double standardError, double tStatistic, double analyticalPvalue, double empiricalPvalue) {
        this.pathway = pathway;
        this.gene = gene;
        this.beta = beta;
        this.standardError = standardError;
        this.tStatistic = tStatistic;
        this.analyticalPvalue = analyticalPvalue;
        this.empiricalPvalue = empiricalPvalue;
    }*/

    public String getPathway() {
        return pathway;
    }

    public void setPathway(String pathway) {
        this.pathway = pathway;
    }

    public String getGene() {
        return gene;
    }

    public void setGene(String gene) {
        this.gene = gene;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

    public double getStandardError() {
        return standardError;
    }

    public void setStandardError(double standardError) {
        this.standardError = standardError;
    }

    public double gettStatistic() {
        return tStatistic;
    }

    public void settStatistic(double tStatistic) {
        this.tStatistic = tStatistic;
    }

    public double getZscore() {
        return zscore;
    }

    public void setZscore(double zscore) {
        this.zscore = zscore;
    }

    public double getAnalyticalPvalue() {
        return analyticalPvalue;
    }

    public void setAnalyticalPvalue(double analyticalPvalue) {
        this.analyticalPvalue = analyticalPvalue;
    }

    public double getEmpiricalPvalue() {
        return empiricalPvalue;
    }

    public void setEmpiricalPvalue(double empiricalPvalue) {
        this.empiricalPvalue = empiricalPvalue;
    }
}
