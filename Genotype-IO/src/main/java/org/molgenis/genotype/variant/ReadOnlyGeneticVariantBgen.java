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
    private Long variantDataSizeInBytes;
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
        this.variantDataSizeInBytes = variantDataSizeInBytes;
    }

//	public static GeneticVariant createSnp(String snpId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(snpId), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createSnp(String snpId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, char refAllele, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(snpId), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), Allele.create(refAllele), variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createSnp(List<String> snpIds, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(snpIds), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createSnp(List<String> snpIds, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, char allele1, char allele2, char refAllele, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(snpIds), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnChars(allele1, allele2), Allele.create(refAllele), variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, Allele allele1, Allele allele2, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createAlleles(allele1, allele2), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(GeneticVariantId variantId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, Allele allele1, Allele allele2, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(variantId, pos, sequenceName,
//				sampleVariantsProvider, Alleles.createAlleles(allele1, allele2), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, String refAllele, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), Allele.create(refAllele), variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantIds), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, String allele1, String allele2, String refAllele, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantIds), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(allele1, allele2), Allele.create(refAllele), variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, List<String> alleles, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(alleles), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, List<String> alleles, String refAllele, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(alleles), Allele.create(refAllele), variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, List<String> alleles, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantIds), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(alleles), null, variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(List<String> variantIds, int pos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, List<String> alleles, String refAllele, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantIds), pos, sequenceName,
//				sampleVariantsProvider, Alleles.createBasedOnString(alleles), Allele.create(refAllele), variantReadingPosition, variantDataSizeInBytes);
//	}
//
//	public static GeneticVariant createVariant(String variantId, int startPos, String sequenceName,
//			SampleVariantsProvider sampleVariantsProvider, Alleles alleles, long variantReadingPosition, int variantDataSizeInBytes) {
//		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), startPos, sequenceName,
//				sampleVariantsProvider, alleles, null, variantReadingPosition, variantDataSizeInBytes);
//	}

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

    public long getVariantDataSizeInBytes() {
        if (variantDataSizeInBytes == null) {
            throw new IllegalStateException("variantDataSizeInBytes was not set. Use .setVariantDataSizeInBytes(...)");
        }
        return variantDataSizeInBytes;
    }

//    public long getVariantProbabilitiesStartPosition() {
//        if (this.variantProbabilitiesStartPosition == null) {
//            readAdditionalVariantData();
//        }
//        return this.variantProbabilitiesStartPosition;
//    }
//
//    public int getVariantGenotypeDataBlockLength() {
//        if (this.variantGenotypeDataBlockLength == null) {
//            readAdditionalVariantData();
//        }
//        return this.variantGenotypeDataBlockLength;
//    }
//
//    public int getVariantGenotypeDataDecompressedBlockLength() {
//        if (this.variantGenotypeDataDecompressedBlockLength == null) {
//            readAdditionalVariantData();
//        }
//        return this.variantGenotypeDataDecompressedBlockLength;
//    }
    /**
     * @param alleles The alleles for the variant to create combinations for.
     * @param ploidy  The ploidity of the sample to create combinations for.
     * @return the combinations of alleles that represent all possible haplotypes
     */
    public static List<List<Integer>> getAlleleCountsPerProbability(List<Integer> alleles, int ploidy) {
        // Construct nested lists
        List<List<Integer>> combinations = new ArrayList<>();

        // Set the maximum value of an allele, which is the number of alleles minus one, because we want to count
        // from 0 to n-1
        int maxAlleleValue = alleles.size() - 1;

        Map<Integer, Integer> initialCounter = new HashMap<>();
        for (int pe : alleles) {
            if (initialCounter.put(pe, 0) != null) {
                throw new IllegalStateException("Duplicate key");
            }
        }

        // Get the combinations
        getAlleleCountsPerProbabilityRecursively(combinations,
                initialCounter, alleles,
                maxAlleleValue, ploidy);

        return combinations;
    }

    /**
     * Method that recursively fills a list of allAlleleCounts.
     * Combinations are always ordered.
     *
     * @param allAlleleCounts   The list of allAlleleCounts to fill.
     * @param alleleCounts    The current combination that is being constructed.
     * @param alleles        The alleles to put into the allAlleleCounts
     * @param maxAlleleValue The number of different values to fit into a allAlleleCounts.
     * @param ploidy         The size of a combination.
     */
    private static void getAlleleCountsPerProbabilityRecursively(
            List<List<Integer>> allAlleleCounts, Map<Integer, Integer> alleleCounts, List<Integer> alleles, int maxAlleleValue, int ploidy) {
        // If the combination is complete, the size of the combination equals the required size.
        // Add the combination and return
        if (alleleCounts.values().stream().reduce(0, Integer::sum) == ploidy) {
            allAlleleCounts.add(0, new ArrayList<>(alleleCounts.values())); // Add in the beginning to maintain the correct order.
            return;
        }

        // Loop through the possible values from high to low.
        for (int newAlleleValue = maxAlleleValue; newAlleleValue >= 0; newAlleleValue--) {
            // Copy the preliminary combination
            LinkedHashMap<Integer, Integer> alleleCountsCopy = new LinkedHashMap<>(alleleCounts);
            // Add a new value to the combination
            alleleCountsCopy.merge(alleles.get(newAlleleValue), 1, Integer::sum);
            getAlleleCountsPerProbabilityRecursively(allAlleleCounts, alleleCountsCopy, alleles, newAlleleValue, ploidy);
        }
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

    public void setVariantDataSizeInBytes(long variantDataSizeInBytes) {
        this.variantDataSizeInBytes = variantDataSizeInBytes;
    }
}
