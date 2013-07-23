package org.molgenis.genotype.variant;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculator;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.util.MafCalculator;
import org.molgenis.genotype.util.MafResult;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class ReadOnlyGeneticVariant extends AbstractGeneticVariant {

    private final GeneticVariantId variantId;
    private final int startPos;
    private final String sequenceName;
    private Map<String, ?> annotationValues = null;
    private final SampleVariantsProvider sampleVariantsProvider;
    private Alleles alleles;
    private final Allele refAllele;
    private MafResult mafResult = null;

    private ReadOnlyGeneticVariant(GeneticVariantId variantId, int startPos, String sequenceName,
            Map<String, ?> annotationValues, SampleVariantsProvider sampleVariantsProvider, Alleles alleles,
            Allele refAllele) {
        super();

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
        this.sequenceName = sequenceName;
        this.sampleVariantsProvider = sampleVariantsProvider;
        this.alleles = alleles;
        this.refAllele = refAllele;
        this.annotationValues = annotationValues != null ? annotationValues : Collections.<String, Object>emptyMap();
    }

    public static GeneticVariant createSnp(String snpId, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(snpId), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), null);
    }

    public static GeneticVariant createSnp(String snpId, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(snpId), pos, sequenceName, null,
                sampleVariantsProvider, null, null);
    }

    public static GeneticVariant createSnp(String snpId, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, char refAllele) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(snpId), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), Allele.create(refAllele));
    }

    public static GeneticVariant createSnp(List<String> snpIds, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(snpIds), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), null);
    }

    public static GeneticVariant createSnp(List<String> snpIds, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, char refAllele) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(snpIds), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), Allele.create(refAllele));
    }

    public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), null);
    }

    public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, String refAllele) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), Allele.create(refAllele));
    }

    public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), null);
    }

    public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, String refAllele) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), Allele.create(refAllele));
    }

    public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, List<String> alleles) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(alleles), null);
    }

    public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, List<String> alleles, String refAllele) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantId), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(alleles), Allele.create(refAllele));
    }

    public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, List<String> alleles) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(alleles), null);
    }

    public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, List<String> alleles, String refAllele) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantIds), pos, sequenceName, null,
                sampleVariantsProvider, Alleles.createBasedOnString(alleles), Allele.create(refAllele));
    }

    public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
            Map<String, ?> annotationValues, SampleVariantsProvider sampleVariantsProvider, List<String> alleles,
            String refAllele) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantIds), pos, sequenceName,
                annotationValues, sampleVariantsProvider, Alleles.createBasedOnString(alleles),
                Allele.create(refAllele));
    }

    public static GeneticVariant createVariant(String variantId, int startPos, String sequenceName,
            SampleVariantsProvider sampleVariantsProvider, Alleles alleles) {
        return new ReadOnlyGeneticVariant(GeneticVariantId.createVariantId(variantId), startPos, sequenceName, null,
                sampleVariantsProvider, alleles, null);
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
    public Alleles getVariantAlleles() {
        return alleles;
    }

    @Override
    public int getAlleleCount() {
        if (alleles == null) {
            getSampleVariants();
        }
        return alleles.getAlleleCount();
    }

    @Override
    public Allele getRefAllele() {
        return refAllele;
    }

    @Override
    public List<Alleles> getSampleVariants() {

        List<Alleles> SampleVariantAlleles = Collections.unmodifiableList(sampleVariantsProvider.getSampleVariants(this));

        if (alleles == null) {
            //set alleles here

            HashSet<Allele> variantAlleles = new HashSet<Allele>(2);

            for (Alleles alleles2 : SampleVariantAlleles) {
                for (Allele allele : alleles2) {
                    variantAlleles.add(allele);
                }
            }

            this.alleles = Alleles.createAlleles(new ArrayList(variantAlleles));


        }
        return SampleVariantAlleles;
    }

    @Override
    public List<Boolean> getSamplePhasing() {
        return sampleVariantsProvider.getSamplePhasing(this);
    }

    @Override
    public Map<String, ?> getAnnotationValues() {
        return Collections.unmodifiableMap(annotationValues);
    }

    @Override
    public double getMinorAlleleFrequency() {
        if (mafResult == null) {
            try {
                mafResult = MafCalculator.calculateMaf(alleles, refAllele, getSampleVariants());
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
            mafResult = MafCalculator.calculateMaf(alleles, refAllele, getSampleVariants());
        }
        return mafResult.getMinorAllele();
    }

    @Override
    public boolean isSnp() {
        if (alleles == null) {
            getSampleVariants();
        }
        return alleles.isSnp();
    }

    @Override
    public boolean isAtOrGcSnp() {
        if (alleles == null) {
            getSampleVariants();
        }
        return alleles.isAtOrGcSnp();
    }

    @Override
    public Ld calculateLd(GeneticVariant other) throws LdCalculatorException {
        return LdCalculator.calculateLd(this, other);
    }

    @Override
    public boolean isBiallelic() {
        if (alleles == null) {
            getSampleVariants();
        }
        return alleles.getAlleleCount() == 2;
    }

    @Override
    public float[] getSampleDosages() {
        if (alleles == null) {
            getSampleVariants();
        }
        return sampleVariantsProvider.getSampleDosage(this);
    }

    @Override
    public SampleVariantsProvider getSampleVariantsProvider() {
        return sampleVariantsProvider;
    }

    @Override
    public byte[] getSampleCalledDosages() {
        if (alleles == null) {
            getSampleVariants();
        }
        return sampleVariantsProvider.getSampleCalledDosage(this);

    }

    /**
     * @param annotationValues the annotationValues to set
     */
    public void setAnnotationValues(Map<String, ?> annotationValues) {
        this.annotationValues = annotationValues;
    }
}
