package org.molgenis.genotype.vcf;

import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaFormat;

public class VcfGeneticVariantMeta implements GeneticVariantMeta {

	private final VcfMeta vcfMeta;
	private final Iterable<String> vcfRecordFormat;

	public VcfGeneticVariantMeta(VcfMeta vcfMeta, Iterable<String> vcfRecordFormat) {
		if (vcfMeta == null) {
			throw new IllegalArgumentException("vcfMeta is null");
		}
		if (vcfRecordFormat == null) {
			throw new IllegalArgumentException("vcfRecord is null");
		}
		this.vcfMeta = vcfMeta;
		this.vcfRecordFormat = vcfRecordFormat;
	}

	@Override
	public Iterable<String> getRecordIds() {
		return vcfRecordFormat;
	}

	@Override
	public Type getRecordType(String recordId) {

		boolean found = false;
		for (String record : vcfRecordFormat) {
			found = record.equals(recordId);
			if (found) {
				break;
			}
		}
		if (!found) {
			return null;
		}

		VcfMetaFormat format = vcfMeta.getFormatMeta(recordId);
		if (format == null) {
			return null;
		}

		// A: one value per alternate allele
		// R: one value for each possible allele (including the reference)
		// G: one value for each possible genotype
		// .: number of possible values varies, is unknown, or is unbounded
		String number = format.getNumber();
		boolean isListValue = number.equals("A") || number.equals("R") || number.equals("G") || number.equals(".") || Integer.valueOf(number) > 1;

		switch (format.getType()) {
			case CHARACTER:
				return isListValue ? Type.CHAR_LIST : Type.CHAR;
			case FLOAT:
				return isListValue ? Type.FLOAT_LIST : Type.FLOAT;
			case INTEGER:
				return isListValue ? Type.INTEGER_LIST : Type.INTEGER;
			case STRING:
				return isListValue ? Type.STRING_LIST : Type.STRING;
			default:
				throw new IllegalArgumentException("invalid vcf format type [" + format.getType() + "]");

		}
	}
}
