package org.molgenis.genotype;

import java.io.File;
import java.util.EnumSet;
import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public enum GenotypeFileType {

	BIM(Pattern.compile(".*\\.bim$", Pattern.CASE_INSENSITIVE), ".bim", "Plink bim"),
	FAM(Pattern.compile(".*\\.fam$", Pattern.CASE_INSENSITIVE), ".fam", "Plink fam"),
	BED(Pattern.compile(".*\\.bed$", Pattern.CASE_INSENSITIVE), ".bed", "Plink bed"),
	PED(Pattern.compile(".*\\.ped$", Pattern.CASE_INSENSITIVE), ".ped", "Plink ped"),
	MAP(Pattern.compile(".*\\.map$", Pattern.CASE_INSENSITIVE), ".map", "Plink map"),
	GEN(Pattern.compile(".*\\.gen$", Pattern.CASE_INSENSITIVE), ".gen", "Oxford gen"),
	HAPS(Pattern.compile(".*\\.haps$", Pattern.CASE_INSENSITIVE), ".haps", "Oxford haps"),
	SAMPLE(Pattern.compile(".*\\.sample$", Pattern.CASE_INSENSITIVE), ".sample", "Oxford sample"),
	VCF(Pattern.compile(".*\\.vcf\\.gz$", Pattern.CASE_INSENSITIVE), ".vcf.gz", "gzipped vcf"),
	BGEN(Pattern.compile(".*\\.bgen$", Pattern.CASE_INSENSITIVE), ".bgen", "Oxford Binary gen"),
	BGENIX(Pattern.compile(".*\\.bgen\\.bgi$", Pattern.CASE_INSENSITIVE), ".bgen.bgi", "Oxford Binary gen index file"),
	VCF_FOLDER(null, "", "folder with gzipped vcf"),
	GEN_FOLDER(null, "", "folder with oxford gen files"),
	TRITYPER_GENOTYPE(Pattern.compile(".*GenotypeMatrix.dat$", Pattern.CASE_INSENSITIVE), "GenotypeMatrix.dat", "Trityper GenotypeMatrix.dat"),//Do not change prefix, it is also used as default file name
	TRITYPER_DOSAGE(Pattern.compile(".*ImputedDosageMatrix.dat$", Pattern.CASE_INSENSITIVE), "ImputedDosageMatrix.dat", "Trityper ImputedDosageMatrix.dat"),//Do not change prefix, it is also used as default file name
	TRITYPER_SNPS(Pattern.compile(".*SNPs.txt(\\.gz)?$", Pattern.CASE_INSENSITIVE), "SNPs.txt", "Trityper SNPs.txt"),//Do not change prefix, it is also used as default file name
	TRITYPER_MAPPINGS(Pattern.compile(".*SNPMappings.txt(\\.gz)?$", Pattern.CASE_INSENSITIVE), "SNPMappings.txt", "Trityper SNPMappings.txt"),//Do not change prefix, it is also used as default file name
	TRITYPER_IND(Pattern.compile(".*Individuals.txt(\\.gz)?$", Pattern.CASE_INSENSITIVE), "Individuals.txt", "Trityper Individuals.txt"),//Do not change prefix, it is also used as default file name
	TRITYPER_PHENO(Pattern.compile(".*PhenotypeInformation.txt(\\.gz)?$", Pattern.CASE_INSENSITIVE), "PhenotypeInformation.txt", "Trityper PhenotypeInformation.txt"),//Do not change prefix, it is also used as default file name
	UNKNOWN(null, "", "unspecified type, potentially oxford gen without extention");

	private final Pattern namePattern;
	private final String suffix;
	private final String friendlyName;
	
	private static final EnumSet<GenotypeFileType> TRITYPER_TYPES = EnumSet.of(TRITYPER_DOSAGE, TRITYPER_GENOTYPE, TRITYPER_IND, TRITYPER_MAPPINGS, TRITYPER_PHENO, TRITYPER_SNPS);

	private GenotypeFileType(Pattern namePattern, String suffix, String friendlyName) {
		this.namePattern = namePattern;
		this.suffix = suffix;
		this.friendlyName = friendlyName;
	}

	public static GenotypeFileType getTypeForPath(String path) {

		for (GenotypeFileType type : values()) {
			if (type.namePattern != null && type.namePattern.matcher(path).matches()) {
				return type;
			}
		}

		return null;

	}

	public String getSuffix() {
		return suffix;
	}

	@Override
	public String toString() {
		return this.friendlyName;
	}

	public String getFriendlyName() {
		return friendlyName;
	}

	public static EnumSet<GenotypeFileType> getTrityperFileTypes() {
		return TRITYPER_TYPES;
	}

	public boolean matches(File file) {
		return namePattern.matcher(file.getAbsolutePath()).matches();
	}
}
