package org.molgenis.genotype.variant;

import java.util.*;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.MafCalculator;
import org.molgenis.genotype.util.MafResult;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantProviderBgen;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class ReadOnlyGeneticVariantBgen extends AbstractGeneticVariant {

    private GeneticVariantId variantId;
    private final int startPos;
    private final String sequenceName;
    private final SampleVariantProviderBgen sampleVariantsProvider;
    private Alleles alleles;
    private final Allele refAllele;
    private MafResult mafResult = null;
    private final GeneticVariantMeta variantMeta = GeneticVariantMetaMap.getGeneticVariantMetaGp();
    private final long variantReadingPosition;
    private final int alleleCount;

    private ReadOnlyGeneticVariantBgen(GeneticVariantId variantId,
                                       int startPos, String sequenceName,
                                       SampleVariantProviderBgen sampleVariantsProvider, int alleleCount,
                                       Alleles alleles, Allele refAllele,
                                       Long variantReadingPosition, Long variantDataSizeInBytes) {

        if (alleles != null) {
            alleles = alleles.createCopyWithoutDuplicates();
        }

        if (refAllele != null) {
            if (alleles == null) {
                throw new GenotypeDataException("A ref allele was supplied while alleles are equal to null");
            }
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
        this.alleleCount = alleleCount;
        this.refAllele = refAllele;
        this.variantReadingPosition = variantReadingPosition;
    }

    public static ReadOnlyGeneticVariantBgen createVariant(List<String> variantId, int pos, String sequenceName,
                                                           SampleVariantProviderBgen sampleVariantsProvider, List<String> alleles,
                                                           long variantReadingPosition) {

        return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
                sampleVariantsProvider, alleles.size(), Alleles.createBasedOnString(alleles), null,
                variantReadingPosition, null);
    }

    public static ReadOnlyGeneticVariantBgen createVariant(String variantId, int pos, String sequenceName,
                                                           SampleVariantProviderBgen sampleVariantsProvider,
                                                           int numberOfAlleles, String allele1, String allele2,
                                                           long variantReadingPosition, long variantDataSizeInBytes) {
        Alleles alleles = null;
        if (numberOfAlleles < 2) {
            alleles = Alleles.createAlleles(Allele.create(allele1));
        } else if (numberOfAlleles == 2) {
            alleles = Alleles.createBasedOnString(allele1, allele2);
        }
        return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
                sampleVariantsProvider, numberOfAlleles, alleles, null,
                variantReadingPosition, variantDataSizeInBytes);
    }

    public void extendWithAdditionalVariantData() {
        ReadOnlyGeneticVariantBgen variant = sampleVariantsProvider.extendReadOnlyGeneticVariantBgen(this);
        updateWithAdditionalVariantData(variant);
    }

    private void updateWithAdditionalVariantData(ReadOnlyGeneticVariantBgen variant) {
        this.alleles = variant.alleles;
        this.variantId = variant.variantId;
    }

    @Override
    public GeneticVariantMeta getVariantMeta() {
        return variantMeta;
    }

    @Override
    public String getPrimaryVariantId() {
        return variantId.getPrimairyId();
    }

    @Override
    public List<String> getAlternativeVariantIds() {
        return getVariantId().getAlternativeIds();
    }

    @Override
    public List<String> getAllIds() {
        return getVariantId().getVariantIds();
    }

    @Override
    public GeneticVariantId getVariantId() {
        extendWithAdditionalVariantData();
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
            this.extendWithAdditionalVariantData();
        }
        return alleles;
    }

    @Override
    public int getAlleleCount() {
        return alleleCount;
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
        return this.sampleVariantsProvider.getSampleProbabilitiesComplex(this);
    }

    @Override
    public double[][][] getSampleGenotypeProbabilitiesPhased() {
        return this.sampleVariantsProvider.getSampleProbabilitiesPhased(this);
    }

    @Override
    public Map<String, ?> getAnnotationValues() {
        return Collections.emptyMap();
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

    public long getVariantReadingPosition() {
        return variantReadingPosition;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;
        ReadOnlyGeneticVariantBgen that = (ReadOnlyGeneticVariantBgen) o;
        return variantReadingPosition == that.variantReadingPosition &&
                sampleVariantsProvider.equals(that.sampleVariantsProvider);
    }

    @Override
    public int hashCode() {
        return (int)(variantReadingPosition ^ (variantReadingPosition >>> 32));
    }
}
