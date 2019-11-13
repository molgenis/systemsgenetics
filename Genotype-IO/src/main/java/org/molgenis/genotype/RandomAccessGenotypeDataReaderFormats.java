package org.molgenis.genotype;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.EnumSet;

import org.molgenis.genotype.bgen.BgenGenotypeData;
import org.molgenis.genotype.oxford.HapsGenotypeData;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.oxford.GenGenotypeData;
import org.molgenis.genotype.plink.BedBimFamGenotypeData;
import org.molgenis.genotype.plink.PedMapGenotypeData;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleFilterableGenotypeDataDecorator;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import static org.molgenis.genotype.GenotypeFileType.*;

public enum RandomAccessGenotypeDataReaderFormats {

	PED_MAP("PED / MAP files", "plink PED / MAP files", EnumSet.of(PED, MAP)),
	VCF("VCF file", "gziped vcf with tabix index file", EnumSet.of(GenotypeFileType.VCF)),
	VCF_FOLDER("VCF folder", "Matches all gziped vcf files + tabix index in a folder. Each file must contains only a sigle chromosome. Each chromsome must be unique to a single file.", EnumSet.of(GenotypeFileType.VCF_FOLDER)),
	SHAPEIT2("Shapeit2 output", ".haps, and .samples with phased haplotypes as outputted by Shapeit2", EnumSet.of(HAPS, SAMPLE)),
	PLINK_BED("Plink BED / BIM / FAM files", "Plink BED / BIM / FAM files", EnumSet.of(BED, BIM, FAM)),
	TRITYPER("TriTyper folder", "Folder with files in trityper format: GenotypeMatrix.dat, Individuals.txt, PhenotypeInformation.txt, SNPMappings.txt, SNPs.txt and optionally: ImputedDosageMatrix.dat", EnumSet.of(TRITYPER_GENOTYPE, TRITYPER_IND, TRITYPER_MAPPINGS, TRITYPER_MAPPINGS, TRITYPER_PHENO, TRITYPER_SNPS)),
	GEN("Oxford GEN / SAMPLE files", "Oxford .gen and .sample", EnumSet.of(GenotypeFileType.GEN, SAMPLE)),
	GEN_FOLDER("Oxford GEN folder", "Folder with oxford .gen and .sample files", EnumSet.of(GenotypeFileType.GEN_FOLDER)),
	BGEN("Binary Oxford GEN", "Oxford .bgen file", EnumSet.of(GenotypeFileType.BGEN));
	private final String name;
	private final String description;
	private final EnumSet<GenotypeFileType> requiredFiles;

	private RandomAccessGenotypeDataReaderFormats(String name, String description, EnumSet<GenotypeFileType> requiredFiles) {
		this.name = name;
		this.description = description;
		this.requiredFiles = requiredFiles;
	}

	public String getName() {
		return name;
	}

	public String getDescription() {
		return description;
	}

	public EnumSet<GenotypeFileType> getRequiredFiles() {
		return requiredFiles;
	}
	

	public static RandomAccessGenotypeDataReaderFormats matchFormatToPath(String... paths) {

		if (paths.length == 0) {
			throw new GenotypeDataException("Cannot find any suitable genotype data format based on path, please make sure the files exist");
		} else if (paths.length == 1) {

			String path = paths[0];

			File pathFile = new File(path);

			if (GenotypeFileType.getTypeForPath(path) == GenotypeFileType.VCF && pathFile.exists()) {
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

			if (new File(path + ".gen").exists() && new File(path + ".sample").exists()) {
				return GEN;
			}

			if (pathFile.exists() && pathFile.isFile() && new File(path + ".sample").exists()) {
				return GEN;
			}

			if (new File(path + ".bgen").exists()) {
				return BGEN;
			}

			if (GenotypeFileType.getTypeForPath(path) == GenotypeFileType.BGEN && pathFile.exists()) {
				return BGEN;
			}

			if (pathFile.isDirectory()) {
				for (File file : pathFile.listFiles()) {
					if (file.getName().endsWith(".vcf.gz")) {
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

		}





		throw new GenotypeDataException("Cannot find any suitable genotype data format based on path, please make sure the files exist");
	}

	public RandomAccessGenotypeData createGenotypeData(String path) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(new String[]{path}, 1000);
	}

	public RandomAccessGenotypeData createGenotypeData(String path, int cacheSize) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(new String[]{path}, cacheSize, null);
	}

	public RandomAccessGenotypeData createGenotypeData(String path, int cacheSize, String forcedSequence) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(new String[]{path}, cacheSize, forcedSequence, 0.34f);
	}

	public RandomAccessGenotypeData createGenotypeData(String path, int cacheSize, String forcedSequence, double minimumPosteriorProbabilityToCall) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(new String[]{path}, cacheSize, forcedSequence, minimumPosteriorProbabilityToCall);
	}

	public RandomAccessGenotypeData createGenotypeData(String[] paths) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(paths, 1000);
	}

	public RandomAccessGenotypeData createGenotypeData(String[] paths, int cacheSize) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(paths, cacheSize, null);
	}

	public RandomAccessGenotypeData createGenotypeData(String[] paths, int cacheSize, String forcedSequence) throws IOException,
			IncompatibleMultiPartGenotypeDataException {
		return createGenotypeData(paths, cacheSize, forcedSequence, 0.34f);
	}

	public RandomAccessGenotypeData createGenotypeData(String[] paths, int cacheSize, String forcedSequence, double minimumPosteriorProbabilityToCall) throws IOException,
			IncompatibleMultiPartGenotypeDataException {

		if (paths == null || paths.length == 0) {
			throw new GenotypeDataException("Error no paths specified");
		}

		switch (this) {
			case PED_MAP:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return new PedMapGenotypeData(new File(paths[0] + ".ped"), new File(paths[0] + ".map"));
			case VCF:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				File vcfFile;
				if (paths[0].endsWith(".vcf.gz")) {
					vcfFile = new File(paths[0]);
				} else {
					vcfFile = new File(paths[0] + ".vcf.gz");
				}
				return new VcfGenotypeData(vcfFile, cacheSize, minimumPosteriorProbabilityToCall);
			case VCF_FOLDER:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return MultiPartGenotypeData.createFromVcfFolder(new File(paths[0]), cacheSize, minimumPosteriorProbabilityToCall);
			case SHAPEIT2:
				return new HapsGenotypeData(new File(paths[0] + ".haps"), new File(
						paths[0] + ".sample"), cacheSize, forcedSequence);
			case PLINK_BED:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return new BedBimFamGenotypeData(new File(paths[0] + ".bed"), new File(paths[0] + ".bim"), new File(paths[0] + ".fam"), cacheSize);
			case TRITYPER:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return new TriTyperGenotypeData(paths[0], cacheSize);

			case GEN:
				if (paths.length == 1) {
					if (new File(paths[0] + ".gen").exists()) {
						return new GenGenotypeData(new File(paths[0] + ".gen"), new File(
								paths[0] + ".sample"), cacheSize, forcedSequence, minimumPosteriorProbabilityToCall);
					} else if (new File(paths[0]).exists() && new File(paths[0] + ".sample").exists()) {
						return new GenGenotypeData(new File(paths[0]), new File(
								paths[0] + ".sample"), cacheSize, forcedSequence, minimumPosteriorProbabilityToCall);
					} else {
						throw new FileNotFoundException("Unable to load gen and sample file at: " + paths[0]);
					}
				} else if (paths.length == 2) {
					File genFile = null;
					File sampleFile = null;
					for (String path : paths) {
						if (GenotypeFileType.getTypeForPath(path) == GenotypeFileType.SAMPLE) {
							sampleFile = new File(path);
						} else {
							genFile = new File(path);
						}
					}
					if (sampleFile == null) {
						throw new GenotypeDataException("Path to .sample file not specified for oxford gen data");
					}
					return new GenGenotypeData(genFile, sampleFile, cacheSize, forcedSequence, minimumPosteriorProbabilityToCall);
				} else {
					throw new GenotypeDataException("Expected 2 files for oxford gen data but found: " + paths.length);
				}

			case GEN_FOLDER:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				return MultiPartGenotypeData.createFromGenFolder(new File(paths[0]), cacheSize, minimumPosteriorProbabilityToCall);

			case BGEN:
				if (forcedSequence != null) {
					throw new GenotypeDataException("Cannot force sequence for " + this.getName());
				}
				if (paths.length == 1) {
					if (new File(paths[0] + ".bgen").exists()) {
						return new BgenGenotypeData(new File(paths[0] + ".bgen"), new File(
								paths[0] + ".sample"), cacheSize, minimumPosteriorProbabilityToCall);
					} else if (new File(paths[0]).exists() && new File(paths[0] + ".sample").exists()) {
						return new BgenGenotypeData(new File(paths[0]), new File(
								paths[0] + ".sample"), cacheSize, minimumPosteriorProbabilityToCall);
					} else {
						throw new FileNotFoundException("Unable to load .bgen and .sample file at: " + paths[0]);
					}
				} else if (paths.length == 2) {
					File bgenFile = null;
					File sampleFile = null;
					for (String path : paths) {
						if (GenotypeFileType.getTypeForPath(path) == GenotypeFileType.SAMPLE) {
							sampleFile = new File(path);
						} else {
							bgenFile = new File(path);
						}
					}
					if (sampleFile == null) {
						throw new GenotypeDataException("Path to .sample file not specified for oxford BGEN data");
					}
					return new BgenGenotypeData(bgenFile, sampleFile, cacheSize, minimumPosteriorProbabilityToCall);
				} else {
					throw new GenotypeDataException("Expected 2 files for oxford BGEN data but found: " + paths.length);
				}

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

		return createFilteredGenotypeData(path, cacheSize, variantFilter, sampleFilter, null, 0.34f);
		
	}
	
	/**
	 * Samples are filtered first then the variant filter is applied.
	 *
	 * @param paths
	 * @param cacheSize
	 * @param variantFilter
	 * @param sampleFilter
	 * @param forcedSequence null if not used
	 * @param minimumPosteriorProbabilityToCall 
	 * @return
	 * @throws IOException
	 */
	public RandomAccessGenotypeData createFilteredGenotypeData(String paths[], int cacheSize, VariantFilter variantFilter, SampleFilter sampleFilter, String forcedSequence, double minimumPosteriorProbabilityToCall) throws IOException {

		switch (this) {
			case TRITYPER:
				return new TriTyperGenotypeData(new File(paths[0]), cacheSize, variantFilter, sampleFilter);
			default:
				RandomAccessGenotypeData genotypeData = createGenotypeData(paths, cacheSize, forcedSequence, minimumPosteriorProbabilityToCall);
				if (sampleFilter != null) {
					genotypeData = new SampleFilterableGenotypeDataDecorator(genotypeData, sampleFilter);
				}
				if (variantFilter != null) {
					genotypeData = new VariantFilterableGenotypeDataDecorator(genotypeData, variantFilter);
				}
				return genotypeData;
		}

	}
	
	/**
	 * Samples are filtered first then the variant filter is applied.
	 *
	 * @param path
	 * @param cacheSize
	 * @param variantFilter
	 * @param sampleFilter
	 * @param forcedSequence null if not used
	 * @param minimumPosteriorProbabilityToCall 
	 * @return
	 * @throws IOException
	 */
	public RandomAccessGenotypeData createFilteredGenotypeData(String path, int cacheSize, VariantFilter variantFilter, SampleFilter sampleFilter, String forcedSequence, double minimumPosteriorProbabilityToCall) throws IOException {

		return createFilteredGenotypeData(new String[]{path}, cacheSize, variantFilter, sampleFilter, forcedSequence, minimumPosteriorProbabilityToCall);

	}
	
	public static RandomAccessGenotypeDataReaderFormats valueOfSmart(String value){
		
		value = value.toUpperCase();
		if(value.equals("VCFFOLDER")){
			return VCF_FOLDER;
		} else if (value.equals("PLINKBED")){
			return PLINK_BED;
		} else if (value.equals("GENFOLDER")){
			return GEN_FOLDER;
		} else if (value.equals("BPLINK")){
			return PLINK_BED;
		} else if (value.equals("B_PLINK")){
			return PLINK_BED;
		} else if (value.equals("PLINKB")){
			return PLINK_BED;
		} else if (value.equals("PLINK_B")){
			return PLINK_BED;
		}
		return RandomAccessGenotypeDataReaderFormats.valueOf(value);
		
	}
	
}
