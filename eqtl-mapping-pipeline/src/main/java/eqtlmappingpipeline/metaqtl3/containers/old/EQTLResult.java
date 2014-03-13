///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package eqtlmappingpipeline.metaqtl3.containers.old;
//
///**
// *
// * @author harmjan
// */
//class EQTLResult implements Comparable<EQTLResult> {
//
////    private int m_p;
////    private Double m_pval;
////    private Double[] m_correlations;
////    private Double[] m_zscores;
////    private Double m_finalZScore;
////    private Double m_finalAbsoluteZScore;
////    private Integer[] m_numSamples;
////    private String m_alleles;
////    private String m_assessedAllele;
////
////    EQTLResult(int p, Double pval, Double[][] correlations, Double[][] zscores, Double aDouble, Double aDouble0, Integer[] numSamples, String alleles, String assessedAllele) {
////	m_p = p;
////	m_pval = pval;
////	int numDatsets = numSamples.length;
////	m_correlations = new Double[numDatsets];
////	m_zscores = new Double[numDatsets];
////
////	for (int i = 0; i < numDatsets; i++) {
////	    m_correlations[i] = correlations[i][p];
////	    correlations[i][p] = null;
////	    m_zscores[i] = zscores[i][p];
////	    zscores[i][p] = null;
////	}
////
////	m_numSamples = numSamples;
////	m_alleles = alleles;
////	m_assessedAllele = assessedAllele;
////    }
////
////    public int compareTo(EQTLResult o) {
////	if (this.m_pval.doubleValue() == o.getM_pval().doubleValue()) {
////	    if (Math.abs(this.m_finalZScore.doubleValue()) == Math.abs(o.getM_finalZScore().doubleValue())) {
////		return 0;
////	    } else if (Math.abs(this.m_finalZScore.doubleValue()) > Math.abs(o.getM_finalZScore().doubleValue())) {
////		return 1;
////	    } else {
////		return -1;
////	    }
////	} else if (this.m_pval.doubleValue() > o.getM_pval().doubleValue()) {
////	    return 1;
////	} else {
////	    return -1;
////	}
////
////    }
////
////    public boolean equals(EQTLResult o) {
////	if (this.m_pval.doubleValue() == o.getM_pval().doubleValue()) {
////	    if (Math.abs(this.m_finalZScore.doubleValue()) == Math.abs(o.getM_finalZScore().doubleValue())) {
////		return true;
////	    } else {
////		return false;
////	    }
////	} else {
////	    return false;
////	}
////    }
////
////    /**
////     * @return the m_p
////     */
////    public int getM_p() {
////	return m_p;
////    }
////
////    /**
////     * @param m_p the m_p to set
////     */
////    public void setM_p(int m_p) {
////	this.m_p = m_p;
////    }
////
////    /**
////     * @return the m_correlations
////     */
////    public Double[] getM_correlations() {
////	return m_correlations;
////    }
////
////    /**
////     * @param m_correlations the m_correlations to set
////     */
////    public void setM_correlations(Double[] m_correlations) {
////	this.m_correlations = m_correlations;
////    }
////
////    /**
////     * @return the m_zscores
////     */
////    public Double[] getM_zscores() {
////	return m_zscores;
////    }
////
////    /**
////     * @param m_zscores the m_zscores to set
////     */
////    public void setM_zscores(Double[] m_zscores) {
////	this.m_zscores = m_zscores;
////    }
////
////    /**
////     * @return the m_finalZScore
////     */
////    public Double getM_finalZScore() {
////	return m_finalZScore;
////    }
////
////    /**
////     * @param m_finalZScore the m_finalZScore to set
////     */
////    public void setM_finalZScore(Double m_finalZScore) {
////	this.m_finalZScore = m_finalZScore;
////    }
////
////    /**
////     * @return the m_finalAbsoluteZScore
////     */
////    public Double getM_finalAbsoluteZScore() {
////	return m_finalAbsoluteZScore;
////    }
////
////    /**
////     * @param m_finalAbsoluteZScore the m_finalAbsoluteZScore to set
////     */
////    public void setM_finalAbsoluteZScore(Double m_finalAbsoluteZScore) {
////	this.m_finalAbsoluteZScore = m_finalAbsoluteZScore;
////    }
////
////    /**
////     * @return the m_numSamples
////     */
////    public Integer[] getM_numSamples() {
////	return m_numSamples;
////    }
////
////    /**
////     * @param m_numSamples the m_numSamples to set
////     */
////    public void setM_numSamples(Integer[] m_numSamples) {
////	this.m_numSamples = m_numSamples;
////    }
////
////    /**
////     * @return the m_alleles
////     */
////    public String getM_alleles() {
////	return m_alleles;
////    }
////
////    /**
////     * @param m_alleles the m_alleles to set
////     */
////    public void setM_alleles(String m_alleles) {
////	this.m_alleles = m_alleles;
////    }
////
////    /**
////     * @return the m_assessedAllele
////     */
////    public String getM_assessedAllele() {
////	return m_assessedAllele;
////    }
////
////    /**
////     * @param m_assessedAllele the m_assessedAllele to set
////     */
////    public void setM_assessedAllele(String m_assessedAllele) {
////	this.m_assessedAllele = m_assessedAllele;
////    }
////
////    /**
////     * @return the m_pval
////     */
////    public Double getM_pval() {
////	return m_pval;
////    }
////
////    /**
////     * @param m_pval the m_pval to set
////     */
////    public void setM_pval(Double m_pval) {
////	this.m_pval = m_pval;
////    }
//}
