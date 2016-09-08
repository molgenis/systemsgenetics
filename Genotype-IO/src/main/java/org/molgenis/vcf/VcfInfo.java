package org.molgenis.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaInfo;

public class VcfInfo {

	private final VcfMeta vcfMeta;
	private String key;
	private String val;

	public VcfInfo(VcfMeta vcfMeta) {
		this(vcfMeta, null, null);
	}

	public VcfInfo(VcfMeta vcfMeta, String key, String val) {
		if (vcfMeta == null) {
			throw new IllegalArgumentException("vcfMeta is null");
		}
		this.vcfMeta = vcfMeta;
		this.key = key;
		this.val = val;
	}

	public String getKey() {
		return key;
	}

	public Object getVal() {
		VcfMetaInfo vcfMetaInfo = vcfMeta.getInfoMeta(key);
		if (vcfMetaInfo == null) {
			return val;
		}

		// A: one value per alternate allele
		// R: one value for each possible allele (including the reference)
		// G: one value for each possible genotype
		// .: number of possible values varies, is unknown, or is unbounded
		String number = vcfMetaInfo.getNumber();
		boolean isListValue;
		try {
			isListValue = number.equals("A") || number.equals("R") || number.equals("G") || number.equals(".") || Integer.valueOf(number) > 1;
		} catch (NumberFormatException ex) {
			throw new GenotypeDataException("Error parsing length of vcf info field. " + number + " is not a valid int or expected preset (A, R, G, .)", ex);
		}
		switch (vcfMetaInfo.getType()) {
			case CHARACTER:
				if (isListValue) {
					String[] valTokens = StringUtils.split(val, ',');
					if (valTokens.length == 0) {
						return Collections.<Character>emptyList();
					}
					List<Character> valList = new ArrayList<Character>(valTokens.length);
					for (String valToken : valTokens) {
						valList.add(Character.valueOf(valToken.charAt(0)));
					}
					return valList;
				} else {
					return Character.valueOf(val.charAt(0));
				}
			case FLAG:
				return null;
			case FLOAT:
				if (isListValue) {
					String[] valTokens = StringUtils.split(val, ',');
					if (valTokens.length == 0) {
						return Collections.<Float>emptyList();
					}
					List<Float> valList = new ArrayList<Float>(valTokens.length);
					for (String valToken : valTokens) {
						if(valToken.equals("."))
						{
							valList.add(null);
						}
						else if (valToken.equalsIgnoreCase("nan") || valToken.equalsIgnoreCase("Na")) {
							valList.add(Float.NaN);
						} else {
							try {
								valList.add(Float.valueOf(valToken));
							} catch (NumberFormatException ex) {
								throw new GenotypeDataException("Error parsing VCF info column value: " + valToken + " is not a float in field: " + key);
							}
						}

					}
					return valList;
				} else {
					if (val.equalsIgnoreCase("nan") || val.equalsIgnoreCase("Na")) {
						return Float.NaN;
					} else {
						try {
							return Float.valueOf(val);
						} catch (NumberFormatException ex) {
							throw new GenotypeDataException("Error parsing VCF info column value: " + val + " is not a float in field: " + key);
						}
					}
				}
			case INTEGER:
				if (isListValue) {
					String[] valTokens = StringUtils.split(val, ',');
					if (valTokens.length == 0) {
						return Collections.<Integer>emptyList();
					}
					List<Integer> valList = new ArrayList<Integer>(valTokens.length);
					for (String valToken : valTokens) {
						if (valToken.equals(".")) {
							valList.add(null);
						} else {
							try {
								valList.add(Integer.valueOf(valToken));
							} catch (NumberFormatException ex) {
								throw new GenotypeDataException("Error parsing VCF info column value: " + val + " is not a int in field: " + key);
							}
						}
					}
					return valList;
				} else {

					try {
						return Integer.valueOf(val);
					} catch (NumberFormatException ex) {
						throw new GenotypeDataException("Error parsing VCF info column value: " + val + " is not a int in field: " + key);
					}

				}
			case STRING:
				if (isListValue) {
					String[] valTokens = StringUtils.split(val, ',');
					if (valTokens.length == 0) {
						return Collections.<String>emptyList();
					}
					return Arrays.asList(valTokens);
				} else {
					return val;
				}
			default:
				throw new GenotypeDataException("unknown vcf info type [" + vcfMetaInfo.getType() + "] in field: " + key);
		}
	}

	public void reset(String key, String val) {
		this.key = key;
		this.val = val;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((key == null) ? 0 : key.hashCode());
		result = prime * result + ((val == null) ? 0 : val.hashCode());
		result = prime * result + ((vcfMeta == null) ? 0 : vcfMeta.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		VcfInfo other = (VcfInfo) obj;
		if (key == null) {
			if (other.key != null) {
				return false;
			}
		} else if (!key.equals(other.key)) {
			return false;
		}
		if (val == null) {
			if (other.val != null) {
				return false;
			}
		} else if (!val.equals(other.val)) {
			return false;
		}
		if (vcfMeta == null) {
			if (other.vcfMeta != null) {
				return false;
			}
		} else if (!vcfMeta.equals(other.vcfMeta)) {
			return false;
		}
		return true;
	}
}
