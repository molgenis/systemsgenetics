/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.variant;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.util.MafCalculator;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class ReadOnlyGeneticVariantTriTyper extends AbstractGeneticVariant {

	private final GeneticVariantId variantId;
	private final int startPos;
	private final String sequenceName;
	private final SampleVariantsProvider sampleVariantsProvider;
	private Alleles alleles;

	public ReadOnlyGeneticVariantTriTyper(String variantId, int startPos, String sequenceName, SampleVariantsProvider sampleVariantsProvider) {
		this.variantId = GeneticVariantId.createVariantId(variantId);
		this.startPos = startPos;
		this.sequenceName = sequenceName.intern();
		this.sampleVariantsProvider = sampleVariantsProvider;
	}

	@Override
	public String getPrimaryVariantId() {
		return variantId.getPrimairyId();
	}

	@Override
	public List<String> getAlternativeVariantIds() {
		return variantId.getAlternativeIds();
	}

	@Override
	public List<String> getAllIds() {
		return variantId.getVariantIds();
	}

	@Override
	public GeneticVariantId getVariantId() {
		return variantId;
	}

	@Override
	public int getStartPos() {
		return startPos;
	}

	@Override
	public String getSequenceName() {
		return sequenceName;
	}

	@Override
	public final Alleles getVariantAlleles() {
		if (alleles == null) {
			getSampleVariants();
		}
		return alleles;
	}

	@Override
	public int getAlleleCount() {
		return this.getVariantAlleles().getAlleleCount();
	}

	@Override
	public Allele getRefAllele() {
		return null;
	}

	@Override
	public List<Alleles> getSampleVariants() {
		List<Alleles> SampleVariantAlleles = Collections.unmodifiableList(sampleVariantsProvider.getSampleVariants(this));

		if (this.alleles == null) {
			//set alleles here

			HashSet<Allele> variantAlleles = new HashSet<Allele>(2);

			for (Alleles alleles2 : SampleVariantAlleles) {
				for (Allele allele : alleles2) {
					variantAlleles.add(allele);
				}
			}

			this.alleles = Alleles.createAlleles(new ArrayList<Allele>(variantAlleles));


		}
		return SampleVariantAlleles;
	}

	@Override
	public Map<String, ?> getAnnotationValues() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public double getMinorAlleleFrequency() {
		return MafCalculator.calculateMaf(this.getVariantAlleles(), this.getRefAllele(), getSampleVariants()).getFreq();
	}

	@Override
	public Allele getMinorAllele() {
		return MafCalculator.calculateMaf(this.getVariantAlleles(), this.getRefAllele(), getSampleVariants()).getMinorAllele();
	}

	@Override
	public float[] getSampleDosages() {
		return sampleVariantsProvider.getSampleDosage(this);
	}

	@Override
	public SampleVariantsProvider getSampleVariantsProvider() {
		return sampleVariantsProvider;
	}

	@Override
	public byte[] getSampleCalledDosages() {
		return sampleVariantsProvider.getSampleCalledDosage(this);
	}

	@Override
	public List<Boolean> getSamplePhasing() {
		return sampleVariantsProvider.getSamplePhasing(this);
	}

	@Override
	public int hashCode() {
		//This works good enough for trityper data.
		return variantId.getPrimairyId().hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (!(obj instanceof ReadOnlyGeneticVariantTriTyper)) {
			return false;
		}
		ReadOnlyGeneticVariantTriTyper other = (ReadOnlyGeneticVariantTriTyper) obj;
		if (getSequenceName() == null) {
			if (other.getSequenceName() != null) {
				return false;
			}
		} else if (!getSequenceName().equals(other.getSequenceName())) {
			return false;
		}
		if (getStartPos() != other.getStartPos()) {
			return false;
		}

		// If we get here pos and sequence are identical
		if (getSampleVariantsProvider() == null) {
			if (other.getSampleVariantsProvider() != null) {
				return false;
			}
		} else if (!getSampleVariantsProvider().equals(other.getSampleVariantsProvider())) {
			//For trityper id must be identical
			return this.variantId.getPrimairyId().equals(other.getPrimaryVariantId());
		}
		return true;
	}
}
