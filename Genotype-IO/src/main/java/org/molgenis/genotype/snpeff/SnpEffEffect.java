package org.molgenis.genotype.snpeff;

/**
 * Note: very incomplete
 *
 *
 * @author Patrick Deelen
 */
public class SnpEffEffect {

	//Enum code copied from: ca.mcgill.mcb.pcingola.snpEffect
	public enum Coding {

		CODING, NON_CODING
	}

	public enum EffectImpact {

		HIGH, MODERATE, LOW, MODIFIER
	}

	public enum EffectType {

		NONE //
		, CHROMOSOME //
		, CHROMOSOME_LARGE_DELETION //
		, INTERGENIC //
		, UPSTREAM //
		, UTR_5_PRIME //
		, UTR_5_DELETED //
		, START_GAINED //
		, SPLICE_SITE_ACCEPTOR //
		, SPLICE_SITE_BRANCH //
		, SPLICE_SITE_BRANCH_U12 //
		, SPLICE_SITE_DONOR //
		, SPLICE_SITE_REGION //
		, START_LOST //
		, SYNONYMOUS_START //
		, NON_SYNONYMOUS_START //
		, CDS //
		, GENE //
		, GENOME //
		, TRANSCRIPT //
		, EXON //
		, EXON_DELETED //
		, NON_SYNONYMOUS_CODING //
		, SYNONYMOUS_CODING //
		, FRAME_SHIFT //
		, CODON_CHANGE //
		, CODON_INSERTION //
		, CODON_CHANGE_PLUS_CODON_INSERTION //
		, CODON_DELETION //
		, CODON_CHANGE_PLUS_CODON_DELETION //
		, RARE_AMINO_ACID //
		, STOP_GAINED //
		, SYNONYMOUS_STOP //
		, NON_SYNONYMOUS_STOP //
		, STOP_LOST //
		, INTRON //
		, UTR_3_PRIME //
		, UTR_3_DELETED //
		, DOWNSTREAM //
		, INTRON_CONSERVED //
		, INTERGENIC_CONSERVED //
		, INTRAGENIC //
		, REGULATION //
		, MOTIF //
		, MICRO_RNA //
		, CUSTOM //
		, NEXT_PROT //
		;
	}

	public enum FunctionalClass {

		NONE, SILENT, MISSENSE, NONSENSE
	}
	private final EffectType effectType;
	private final EffectImpact effectImpact;
	private final FunctionalClass functionalClass;

	public SnpEffEffect(EffectType effectType, EffectImpact effectImpact, FunctionalClass functionalClass) {
		this.effectType = effectType;
		this.effectImpact = effectImpact;
		this.functionalClass = functionalClass;
	}

	public EffectType getEffectType() {
		return effectType;
	}

	public EffectImpact getEffectImpact() {
		return effectImpact;
	}

	public FunctionalClass getFunctionalClass() {
		return functionalClass;
	}
}