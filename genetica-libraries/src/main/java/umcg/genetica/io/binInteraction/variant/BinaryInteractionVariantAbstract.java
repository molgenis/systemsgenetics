package umcg.genetica.io.binInteraction.variant;

import org.molgenis.genotype.Allele;

/**
 *
 * @author Patrick Deelen
 */
public abstract class BinaryInteractionVariantAbstract implements BinaryInteractionVariant{

	private final String name;
	private final String chr;
	private final int pos;
	private final Allele refAllele;
	private final Allele altAllele;

	public BinaryInteractionVariantAbstract(String name, String chr, int pos, Allele refAllele, Allele altAllele) {
		this.name = name;
		this.chr = chr;
		this.pos = pos;
		this.refAllele = refAllele;
		this.altAllele = altAllele;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public String getChr() {
		return chr;
	}

	@Override
	public int getPos() {
		return pos;
	}

	@Override
	public Allele getRefAllele() {
		return refAllele;
	}

	@Override
	public Allele getAltAllele() {
		return altAllele;
	}
	
}
