package org.molgenis.genotype;

import java.io.File;
import java.io.IOException;
import org.apache.log4j.Logger;
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
	VCF_FOLDER("VCF folder", "Matches all gziped vcf files + tabix index in a folder. Each file must contain a sigle chromosome. Each chromsome should be unique to a single file."),
	SHAPEIT2("Shapeit2 output", ".haps.gz, haps.gz.tbi and .samples with phased haplotypes as outputted by Shapeit2 converted to tab separated and bgziped with tabix index"),
	PLINK_BED("Plink BED / BIM / FAM files", "Plink BED / BIM / FAM files"),
	TRITYPER("TriTyper folder", "Folder with files in trityper format: GenotypeMatrix.dat, Individuals.txt, PhenotypeInformation.txt, SNPMappings.txt, SNPs.txt. Optionally: ImputedDosageMatrix.dat");
	private final String name;
	private final String description;
	private static final Logger LOGGER = Logger.getLogger(RandomAccessGenotypeDataReaderFormats.class);

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

	public static RandomAccessGenotypeDataReaderFormats matchFormatToPath(String path) {

		File pathFile = new File(path);

		if (path.endsWith(".vcf.zg") && pathFile.exists()) {
			return VCF;
		}

		if (new File(path + ".vcf.gz").exists()) {
			return VCF;
		}

		if (new File(path + ".ped").exists() && new File(path + ".map").exists()) {
			return PED_MAP;
		}

		if (new File(path + ".bed").exists() && new File(path + ".bim").exists() && new File(path + ".fam").exists()) {
			return PLINK_BED;
		}

		if (new File(path + ".haps").exists() && new File(path + ".sample").exists()) {
			return SHAPEIT2;
		}

		if (pathFile.isDirectory()) {
			for (File file : pathFile.listFiles()) {
				if (file.getName().endsWith(".vcg.gz")) {
					return VCF_FOLDER;
				}
			}
			if (new File(pathFile, "GenotypeMatrix.dat").exists()
					&& (new File(pathFile, "SNPs.txt").exists() || new File(pathFile, "SNPs.txt.gz").exists())
					&& (new File(pathFile, "SNPMappings.txt").exists() || new File(pathFile, "SNPMappings.txt.gz").exists())
					&& (new File(pathFile, "Individuals.txt").exists() || new File(pathFile, "Individuals.txt.gz").exists())
					&& (new File(pathFile, "PhenotypeInformation.txt").exists() || new File(pathFile, "PhenotypeInformation.txt.gz").exists())) {
				
				return TRITYPER;				
				
			}
		}



		throw new GenotypeDataException("Cannot find any suitable genotype data format based on path");
	}

	public RandomAccessGenotypeData createGenotypeData(String path) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(path, 1000);
	}

	public RandomAccessGenotypeData createGenotypeData(String path, int cacheSize) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(path, cacheSize, null);
	}

	public RandomAccessGenotypeData createGenotypeData(String path, int cacheSize, String forcedSequence) throws IOException,
			IncompatibleMultiPartGenotypeDataException {

		switch (this) {
			case PED_MAP:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return new PedMapGenotypeData(new File(path + ".ped"), new File(path + ".map"));
			case VCF:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				File vcfFile;
				if (path.endsWith(".vcf.gz")) {
					vcfFile = new File(path);
				} else {
					vcfFile = new File(path + ".vcf.gz");
				}
				return new VcfGenotypeData(vcfFile, cacheSize);
			case VCF_FOLDER:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return MultiPartGenotypeData.createFromVcfFolder(new File(path), cacheSize);
			case SHAPEIT2:
				return new Impute2GenotypeData(new File(path + ".haps"), new File(
						path + ".sample"), cacheSize, forcedSequence);
			case PLINK_BED:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return new BedBimFamGenotypeData(new File(path + ".bed"), new File(path + ".bim"), new File(path + ".fam"), cacheSize);
			case TRITYPER:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
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
