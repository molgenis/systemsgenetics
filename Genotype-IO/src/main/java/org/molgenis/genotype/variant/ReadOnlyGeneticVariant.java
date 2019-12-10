package org.molgenis.genotype.variant;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.MafCalculator;
import org.molgenis.genotype.util.MafResult;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class ReadOnlyGeneticVariant extends AbstractGeneticVariant {

	private final GeneticVariantId variantId;
	private final int startPos;
	private final String sequenceName;
	private final SampleVariantsProvider sampleVariantsProvider;
	private final Alleles alleles;
	private final Allele refAllele;
	private MafResult mafResult = null;
	protected final Map<String, ?> annotationValues;
	private final GeneticVariantMeta variantMeta;

	private ReadOnlyGeneticVariant(GeneticVariantMeta variantMeta, GeneticVariantId variantId, int startPos, String sequenceName,
			Map<String, ?> annotationValues, SampleVariantsProvider sampleVariantsProvider, Alleles alleles,
			Allele refAllele) {
		
		this.variantMeta = variantMeta;
		
		alleles = alleles.createCopyWithoutDuplicates();

		if (refAllele != null) {
			if (!alleles.contains(refAllele)) {
				throw new GenotypeDataException("Supplied ref allele (" + refAllele
						+ ") is not a found in supplied alleles " + alleles.getAllelesAsString()
						+ " for variant with ID: " + variantId.getPrimairyId() + " at: " + sequenceName + ":"
						+ startPos);
			}
			if (alleles.get(0) != refAllele) {
				// ref allele is not first in alleles. We need to change this
				ArrayList<Allele> allelesWithoutRef = new ArrayList<Allele>(alleles.getAlleles());
				allelesWithoutRef.remove(refAllele);
				allelesWithoutRef.add(0, refAllele);
				alleles = Alleles.createAlleles(allelesWithoutRef);
			}
		}

		this.variantId = variantId;
		this.startPos = startPos;
		this.sequenceName = sequenceName.intern();
		this.sampleVariantsProvider = sampleVariantsProvider;
		this.alleles = alleles;
		this.refAllele = refAllele;
		this.annotationValues = annotationValues != null ? annotationValues : Collections.<String, Object>emptyMap();

	}

	public static GeneticVariant createSnp(GeneticVariantMeta variantMeta, String snpId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(snpId), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), null);
	}

	public static GeneticVariant createSnp(GeneticVariantMeta variantMeta, String snpId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, char refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(snpId), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), Allele.create(refAllele));
	}

	public static GeneticVariant createSnp(GeneticVariantMeta variantMeta, List<String> snpIds, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(snpIds), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), null);
	}

	public static GeneticVariant createSnp(GeneticVariantMeta variantMeta, List<String> snpIds, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, char refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(snpIds), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), Allele.create(refAllele));
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, String variantId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), null);
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, String variantId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, Allele allele1, Allele allele2) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createAlleles(allele1, allele2), null);
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, GeneticVariantId variantId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, Allele allele1, Allele allele2) {
		return new ReadOnlyGeneticVariant(variantMeta, variantId, pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createAlleles(allele1, allele2), null);
	}
	
	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, String variantId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, String refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), Allele.create(refAllele));
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, List<String> variantIds, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), null);
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, List<String> variantIds, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, String refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), Allele.create(refAllele));
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, String variantId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, List<String> alleles) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(alleles), null);
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, String variantId, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, List<String> alleles, String refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(alleles), Allele.create(refAllele));
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, List<String> variantIds, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, List<String> alleles) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(alleles), null);
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, List<String> variantIds, int pos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, List<String> alleles, String refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
				sampleVariantsProvider, Alleles.createBasedOnString(alleles), Allele.create(refAllele));
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, List<String> variantIds, int pos, String sequenceName,
			Map<String, ?> annotationValues, SampleVariantsProvider sampleVariantsProvider, List<String> alleles,
			String refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantIds), pos, sequenceName,
				annotationValues, sampleVariantsProvider, Alleles.createBasedOnString(alleles),
				Allele.create(refAllele));
	}
	
	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, List<String> variantIds, int pos, String sequenceName,
			Map<String, ?> annotationValues, SampleVariantsProvider sampleVariantsProvider, Alleles alleles,
			Allele refAllele) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantIds), pos, sequenceName,
				annotationValues, sampleVariantsProvider, alleles, refAllele);
	}

	public static GeneticVariant createVariant(GeneticVariantMeta variantMeta, String variantId, int startPos, String sequenceName,
			SampleVariantsProvider sampleVariantsProvider, Alleles alleles) {
		return new ReadOnlyGeneticVariant(variantMeta, GeneticVariantId.createVariantId(variantId), startPos, sequenceName, null,
				sampleVariantsProvider, alleles, null);
	}
	
	@Override
	public GeneticVariantMeta getVariantMeta()
	{
		return variantMeta;
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
		return alleles;
	}

	@Override
	public int getAlleleCount() {
		return this.getVariantAlleles().getAlleleCount();
	}

	@Override
	public Allele getRefAllele() {
		return refAllele;
	}

	@Override
	public final List<Alleles> getSampleVariants() {
		return Collections.unmodifiableList(sampleVariantsProvider.getSampleVariants(this));
	}

	@Override
	public List<Boolean> getSamplePhasing() {
		return sampleVariantsProvider.getSamplePhasing(this);
	}

	@Override
	public float[][] getSampleGenotypeProbilities() {
		return sampleVariantsProvider.getSampleProbilities(this);
	}

	@Override
	public double[][] getSampleGenotypeProbabilitiesComplex() {
		return sampleVariantsProvider.getSampleProbabilitiesComplex(this);
	}

	@Override
	public Map<String, ?> getAnnotationValues() {
		return Collections.unmodifiableMap(annotationValues);
	}

	@Override
	public double getMinorAlleleFrequency() {
		if (mafResult == null) {
			try {
				mafResult = MafCalculator.calculateMaf(this.getVariantAlleles(), this.getRefAllele(), this.getSampleVariants());
			} catch (NullPointerException e) {
				throw new GenotypeDataException("NullPointerException in maf caculation. " + getVariantAlleles() + " ref: "
						+ getRefAllele(), e);
			}
		}

		return mafResult.getFreq();

	}

	@Override
	public Allele getMinorAllele() {
		if (mafResult == null) {
			mafResult = MafCalculator.calculateMaf(this.getVariantAlleles(), this.getRefAllele(), this.getSampleVariants());
		}
		return mafResult.getMinorAllele();
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
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords() {
		return sampleVariantsProvider.getSampleGenotypeRecords(this);
	}
	
	@Override
	public Alleles getAlternativeAlleles() {
		ArrayList<Allele> altAlleles = new ArrayList<>(this.getVariantAlleles().getAlleles());
		altAlleles.remove(this.getRefAllele());
		return Alleles.createAlleles(altAlleles);
	}

	@Override
	public double[][][] getSampleGenotypeProbabilitiesPhased() {
		return this.sampleVariantsProvider.getSampleProbabilitiesPhased(this);
	}
	
}
