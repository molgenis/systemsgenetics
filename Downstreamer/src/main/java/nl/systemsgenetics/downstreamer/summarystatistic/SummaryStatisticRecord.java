package nl.systemsgenetics.downstreamer.summarystatistic;

import org.molgenis.genotype.variant.GeneticVariant;

import java.io.Serializable;

/**
 * Created by olivier on 06/12/2017.
 */
public class SummaryStatisticRecord extends SNP implements Comparable<SummaryStatisticRecord>, Serializable {

	//TODO: refactoring to more appropriate names
	private double beta;
	private double pvalue;
	private double standardError;
	private double falseDiscoveryRate;
	private String traitName;

	public SummaryStatisticRecord(SummaryStatisticRecord other) {
		this.setAllele1(other.getAllele1());
		this.setAllele2(other.getAllele2());
		this.setPosition(other.getPosition());
		this.setSequenceName(other.getSequenceName());
		this.setPrimaryVariantId(other.getPrimaryVariantId());
		this.beta = other.getBeta();
		this.pvalue = other.getPvalue();
		this.standardError = other.getStandardError();
		this.falseDiscoveryRate = other.getFalseDiscoveryRate();
		this.traitName = other.getTraitName();
	}

	/**
	 * Instantiates a new minimal Summary statistic record. based on gentic
	 * variant
	 */
	public SummaryStatisticRecord(GeneticVariant variant, double pvalue) {
		this.pvalue = pvalue;
		this.primaryVariantId = variant.getPrimaryVariantId();
		this.position = variant.getStartPos();
		this.sequenceName = variant.getSequenceName();
	}

	public SummaryStatisticRecord(String id, String chr, int pos, double pvalue) {
		this.pvalue = pvalue;
		this.primaryVariantId = id;
		this.position = pos;
		this.sequenceName = chr;
	}

	/**
	 * Instantiates a new Summary statistic record.
	 */
	public SummaryStatisticRecord() {
	}

	/**
	 * Sets effect allele.
	 *
	 * @param effectAllele the effect allele
	 */
	public void setEffectAllele(String effectAllele) {
		this.allele1 = effectAllele;
	}

	/**
	 * Gets effect allele.
	 *
	 * @return the effect allele
	 */
	public String getEffectAllele() {
		return allele1;
	}

	/**
	 * Sets other allele.
	 *
	 * @param otherAllele the other allele
	 */
	public void setOtherAllele(String otherAllele) {
		this.allele2 = otherAllele;
	}

	/**
	 * Gets other allele.
	 *
	 * @return the other allele
	 */
	public String getOtherAllele() {
		return allele2;
	}

	/**
	 * Gets beta.
	 *
	 * @return the beta
	 */
	public double getBeta() {
		return beta;
	}

	/**
	 * Sets beta.
	 *
	 * @param beta the beta
	 */
	public void setBeta(String beta) {
		this.beta = Double.parseDouble(beta);
	}

	/**
	 * Sets beta.
	 *
	 * @param beta the beta
	 */
	public void setBeta(double beta) {
		this.beta = beta;
	}

	/**
	 * Gets pvalue.
	 *
	 * @return the pvalue
	 */
	public double getPvalue() {
		return pvalue;
	}

	/**
	 * Sets pvalue.
	 *
	 * @param pvalue the pvalue
	 */
	public void setPvalue(String pvalue) {
		this.pvalue = Double.parseDouble(pvalue);
	}

	/**
	 * Sets pvalue.
	 *
	 * @param pvalue the pvalue
	 */
	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}

	/**
	 * Gets trait name.
	 *
	 * @return the trait name
	 */
	public String getTraitName() {
		return traitName;
	}

	/**
	 * Sets trait name.
	 *
	 * @param traitName the trait name
	 */
	public void setTraitName(String traitName) {
		this.traitName = traitName;
	}

	/**
	 * Gets standard error.
	 *
	 * @return the standard error
	 */
	public double getStandardError() {
		return standardError;
	}

	/**
	 * Sets standard error.
	 *
	 * @param standardError the standard error
	 */
	public void setStandardError(double standardError) {
		this.standardError = standardError;
	}

	/**
	 * Gets false discovery rate.
	 *
	 * @return the false discovery rate
	 */
	public double getFalseDiscoveryRate() {
		return falseDiscoveryRate;
	}

	/**
	 * Sets false discovery rate.
	 *
	 * @param falseDiscoveryRate the false discovery rate
	 */
	public void setFalseDiscoveryRate(double falseDiscoveryRate) {
		this.falseDiscoveryRate = falseDiscoveryRate;
	}

	/**
	 * To output string string.
	 *
	 * @param seperator the seperator
	 * @return the string
	 */
	public String toOutputString(String seperator) {
		return toOutputString(seperator, new boolean[]{true, true, true, true, true, true, true});
	}

	/**
	 * To output string string.
	 *
	 * @param separator the separator
	 * @param includeInOutput the include in output
	 * @return the string
	 */
	public String toOutputString(String separator, boolean[] includeInOutput) {
		StringBuilder output = new StringBuilder();

		if (includeInOutput[0]) {
			output.append(primaryVariantId)
					.append(separator);
		}
		if (includeInOutput[1]) {
			output.append(sequenceName)
					.append(separator);
		}
		if (includeInOutput[2]) {
			output.append(Long.toString(position))
					.append(separator);
		}
		if (includeInOutput[3]) {
			output.append(allele1)
					.append(separator);
		}
		if (includeInOutput[4]) {
			output.append(allele2)
					.append(separator);
		}
		if (includeInOutput[5]) {
			output.append(beta)
					.append(separator);
		}
		if (includeInOutput[6]) {
			output.append(pvalue)
					.append(separator);
		}

		return output.toString();
	}

	@Override
	public int compareTo(SummaryStatisticRecord o) {
		// In case on ties sort on the effect size
		if (this.pvalue == o.getPvalue() && this.beta != 0) {
			return Double.compare(Math.abs(this.beta), Math.abs(o.getBeta()));
		}
		return Double.compare(this.pvalue, o.getPvalue());
	}
}
