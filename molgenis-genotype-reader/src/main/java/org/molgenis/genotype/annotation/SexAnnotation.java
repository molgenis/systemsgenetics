package org.molgenis.genotype.annotation;

public enum SexAnnotation {

	MALE((byte) 1, "Male"), FEMALE((byte) 2, "Female"), UNKNOWN((byte) 0, "Unknown");
	private final byte plinkSex;
	private final String gender;

	private SexAnnotation(byte plinkSex, String gender) {
		this.plinkSex = plinkSex;
		this.gender = gender;
	}

	/**
	 * @return the plinkSex
	 */
	public byte getPlinkSex() {
		return plinkSex;
	}

	public static SexAnnotation getSexAnnotationForPlink(byte plinkSex) {
		switch (plinkSex) {
			case 1:
				return SexAnnotation.MALE;
			case 2:
				return SexAnnotation.FEMALE;
			default:
				return SexAnnotation.UNKNOWN;
		}
	}

	public static SexAnnotation getSexAnnotationForTriTyper(String ttSex) {
		if (ttSex == null) {
			return SexAnnotation.UNKNOWN;
		} else if (ttSex.toLowerCase().equals("female")) {
			return SexAnnotation.FEMALE;
		} else if (ttSex.toLowerCase().equals("male")) {
			return SexAnnotation.MALE;
		} else {
			return SexAnnotation.UNKNOWN;
		}
	}

	@Override
	public String toString() {
		return gender;
	}

	public String getGender() {
		return gender;
	}
}
