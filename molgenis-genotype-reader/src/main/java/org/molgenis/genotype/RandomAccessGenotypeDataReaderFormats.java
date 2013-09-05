package org.molgenis.genotype;

import java.io.File;
import java.io.IOException;

import org.molgenis.genotype.impute2.Impute2GenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.plink.BedBimFamGenotypeData;
import org.molgenis.genotype.plink.PedMapGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleFilterableGenotypeDataDecorator;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;
import org.molgenis.genotype.vcf.VcfGenotypeData;

public enum RandomAccessGenotypeDataReaderFormats {

	PED_MAP("PED / MAP files", "plink PED / MAP files"),
	VCF("VCF file", "gziped vcf with tabix index file"),
	VCF_FOLDER("VCF folder", "Matches all gziped vcf files + tabix index in a folder"),
	SHAPEIT2("Shapeit2 output", ".haps.gz, haps.gz.tbi and .samples with phased haplotypes as outputted by Shapeit2 converted to tab separated and bgziped with tabix index"),
	PLINK_BED("Plink BED / BIM / FAM files", "Plink BED / BIM / FAM files"),
	TRITYPER("TriTyper folder", "Folder with files in trityper format: GenotypeMatrix.dat, Individuals.txt, PhenotypeInformation.txt, SNPMappings.txt, SNPs.txt. Optionally: ImputedDosageMatrix.dat");
	private final String name;
	private final String description;

	RandomAccessGenotypeDataReaderFormats(String name, String description) {
		this.name = name;
		this.description = description;
	}

	public String getName() {
		return name;
	}

	public String getDescription() {
		return description;
	}

	public RandomAccessGenotypeData createGenotypeData(String path, int cacheSize) throws IOException,
			IncompatibleMultiPartGenotypeDataException {

		switch (this) {
			case PED_MAP:
				return new PedMapGenotypeData(new File(path + ".ped"), new File(path + ".map"));
			case VCF:
				return new VcfGenotypeData(new File(path + ".vcf.gz"), cacheSize);
			case VCF_FOLDER:
				return MultiPartGenotypeData.createFromVcfFolder(new File(path), cacheSize);
			case SHAPEIT2:
				return new Impute2GenotypeData(new File(path + ".haps"), new File(
						path + ".sample"), cacheSize);
			case PLINK_BED:
				return new BedBimFamGenotypeData(new File(path + ".bed"), new File(path + ".bim"), new File(path + ".fam"), cacheSize);
			case TRITYPER:
				return new TriTyperGenotypeData(path, cacheSize);
			default:
				throw new RuntimeException("This should not be reachable. Please contact the authors");
		}
	}

	/**
	 * Samples are filtered first then the variant filter is applied.
	 *
	 * @param path
	 * @param cacheSize
	 * @param variantFilter
	 * @param sampleFilter
	 * @return
	 * @throws IOException
	 */
	public RandomAccessGenotypeData createFilteredGenotypeData(String path, int cacheSize, VariantFilter variantFilter, SampleFilter sampleFilter) throws IOException {

		switch (this) {
			case TRITYPER:
				return new TriTyperGenotypeData(new File(path), cacheSize, variantFilter, sampleFilter);
			default:
				RandomAccessGenotypeData genotypeData = createGenotypeData(path, cacheSize);
				if (sampleFilter != null) {
					genotypeData = new SampleFilterableGenotypeDataDecorator(genotypeData, sampleFilter);
				}
				if (variantFilter != null) {
					genotypeData = new VariantFilterableGenotypeDataDecorator(genotypeData, variantFilter);
				}
				return genotypeData;
		}

	}
}
