package org.molgenis.genotype.variant;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.bgen.BgenGenotypeData;
import org.molgenis.genotype.bgen.VariantGenotypeDataBlockInfo;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.MafCalculator;
import org.molgenis.genotype.util.MafResult;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

public class ReadOnlyGeneticVariantBgen extends AbstractGeneticVariant {

	private GeneticVariantId variantId;
	private final int startPos;
	private final String sequenceName;
	private final SampleVariantsProvider sampleVariantsProvider; // Change for sampleVariantsProviderBgen
	private Alleles alleles;
	private final Allele refAllele;
	private MafResult mafResult = null;
	private static final GeneticVariantMeta variantMeta = GeneticVariantMetaMap.getGeneticVariantMetaGp();
	private long variantReadingPosition;
	private int variantDataSizeInBytes;
	private Long variantProbabilitiesStartPosition;
	private Integer variantGenotypeDataBlockLength;
	private Integer variantGenotypeDataDecompressedBlockLength;
	private int alleleCount;

	private ReadOnlyGeneticVariantBgen(GeneticVariantId variantId,
									   int startPos, String sequenceName,
									   SampleVariantsProvider sampleVariantsProvider, int alleleCount,
									   Alleles alleles, Allele refAllele,
									   long variantReadingPosition, int variantDataSizeInBytes,
									   Long variantProbabilitiesStartPosition, Integer variantGenotypeDataBlockLength,
									   Integer variantGenotypeDataDecompressedBlockLength) {

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
		this.variantProbabilitiesStartPosition = variantProbabilitiesStartPosition;
		this.variantGenotypeDataBlockLength = variantGenotypeDataBlockLength;
		this.variantGenotypeDataDecompressedBlockLength = variantGenotypeDataDecompressedBlockLength;
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

	public static GeneticVariant createVariant(List<String> variantId, int pos, String sequenceName,
											   SampleVariantsProvider sampleVariantsProvider, List<String> alleles,
											   VariantGenotypeDataBlockInfo variantGenotypeDataBlockInfo) {


		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
				sampleVariantsProvider, alleles.size(), Alleles.createBasedOnString(alleles), null,
				variantGenotypeDataBlockInfo.getVariantReadingPosition(),
				variantGenotypeDataBlockInfo.getVariantDataSizeInBytes(),
				variantGenotypeDataBlockInfo.getVariantProbabilitiesStartPosition(),
				(int) variantGenotypeDataBlockInfo.getBlockLength(),
				(int) variantGenotypeDataBlockInfo.getDecompressedBlockLength());
	}

	public static GeneticVariant createVariant(String variantId, int pos, String sequenceName,
											   SampleVariantsProvider sampleVariantsProvider,
											   int numberOfAlleles, String allele1, String allele2,
											   long variantReadingPosition, int variantDataSizeInBytes) {
		Alleles alleles = Alleles.createBasedOnString(allele1, allele2);
		if (numberOfAlleles > 2) {
			alleles = null;
		}
		return new ReadOnlyGeneticVariantBgen(GeneticVariantId.createVariantId(variantId), pos, sequenceName,
				sampleVariantsProvider, numberOfAlleles, alleles, null,
				variantReadingPosition, variantDataSizeInBytes,
				null, null, null);
	}

	private void readFromSource() {
		try {
			ReadOnlyGeneticVariantBgen variant = ((BgenGenotypeData) sampleVariantsProvider).
					processVariantIdentifyingData(getVariantReadingPosition());

			this.alleles = variant.alleles;
			this.variantId = variant.variantId;
			this.variantProbabilitiesStartPosition = variant.variantProbabilitiesStartPosition;
			this.variantGenotypeDataBlockLength = variant.variantGenotypeDataBlockLength;
			this.variantGenotypeDataDecompressedBlockLength = variant.variantGenotypeDataDecompressedBlockLength;

		} catch (IOException e) {
			throw new GenotypeDataException(String.format(
					"Could not read variant data %s at position %d%n",
					getPrimaryVariantId(), getVariantReadingPosition()));
		}
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
		return getVariantId().getAlternativeIds();
	}

	@Override
	public List<String> getAllIds() {
		return getVariantId().getVariantIds();
	}

	@Override
	public GeneticVariantId getVariantId() {
		readFromSource();
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
			this.readFromSource();
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

	public double[][] getSampleGenotypeProbabilitiesBgen() {
		return ((BgenGenotypeData) this.sampleVariantsProvider).getSampleGenotypeProbabilitiesBgen(this);
	}

	public long getVariantProbabilitiesStartPosition() {
		if (this.variantProbabilitiesStartPosition == null) {
			readFromSource();
		}
		return this.variantProbabilitiesStartPosition;
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

	@Override
	public int hashCode() {
		return (int) (this.getVariantReadingPosition() ^ (this.getVariantReadingPosition() >>> 32));
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
		final ReadOnlyGeneticVariantBgen other = (ReadOnlyGeneticVariantBgen) obj;
		return this.getVariantReadingPosition() == other.getVariantReadingPosition();
	}

	public long getVariantReadingPosition() {
		return variantReadingPosition;
	}

	public int getVariantDataSizeInBytes() {
		return variantDataSizeInBytes;
	}

	public int getVariantGenotypeDataBlockLength() {
		if (this.variantGenotypeDataBlockLength == null) {
			readFromSource();
		}
		return this.variantGenotypeDataBlockLength;
	}

	public int getVariantGenotypeDataDecompressedBlockLength() {
		if (this.variantGenotypeDataDecompressedBlockLength == null) {
			readFromSource();
		}
		return this.variantGenotypeDataDecompressedBlockLength;
	}
}
