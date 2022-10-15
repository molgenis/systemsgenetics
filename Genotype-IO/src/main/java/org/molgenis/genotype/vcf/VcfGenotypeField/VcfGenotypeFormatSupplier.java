package org.molgenis.genotype.vcf.VcfGenotypeField;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.vcf.VcfRecord;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;

/**
 * Class that picks Vcf Genotype Formats and identifiers from a list of Vcf Genotype Formats,
 * possible preferring a preset, preferred genotype format, possibly with a synonymous identifier.
 */
public class VcfGenotypeFormatSupplier {
    private VcfGenotypeFormat preferredGenotypeFormat;
    private String preferredGenotypeFormatIdentifier;
    private boolean raiseExceptionIfUnavailable;

    public VcfGenotypeFormatSupplier(VcfGenotypeFormat preferredGenotypeFormat) {
        this(preferredGenotypeFormat, preferredGenotypeFormat.toString(), false);
    }

    public VcfGenotypeFormatSupplier(VcfGenotypeFormat preferredGenotypeFormat, String formatIdentifier) {
        this(preferredGenotypeFormat, formatIdentifier, false);
    }

    public VcfGenotypeFormatSupplier(VcfGenotypeFormat preferredGenotypeFormat, boolean raiseExceptionIfUnavailable) {
        this(preferredGenotypeFormat, preferredGenotypeFormat.toString(), false);
    }

    public VcfGenotypeFormatSupplier(VcfGenotypeFormat preferredGenotypeFormat, String formatIdentifier, boolean raiseExceptionIfUnavailable) {

        this.preferredGenotypeFormat = preferredGenotypeFormat;
        this.preferredGenotypeFormatIdentifier = formatIdentifier;
        this.raiseExceptionIfUnavailable = raiseExceptionIfUnavailable;
    }

    public VcfGenotypeFormatSupplier() {
        this(null, null, false);
    }

    /**
     * @param vcfRecord record, row, within a VCF file. corresponding to a particular variant.
     * @param genotypeDosageFieldPrecedence LinkedHashSet that lists all formats that can be read,
     *                                      in order of precedence (high precedence to low precedence).
     * @return The preferred genotype format if this is available from the vcf record and the list of
     * formats that can be read according to the genotype dosage field precedence hash set. If not, will return
     * the first format from the genotype dosage field precedence list that can be read from the vcf record.
     * If nothing matches these conditions, null is returned.
     */
    public VcfGenotypeFormat getVcfGenotypeFormat(
            VcfRecord vcfRecord,
            LinkedHashSet<VcfGenotypeFormat> genotypeDosageFieldPrecedence) {

        List<String> formatIdentifiers = Arrays.asList(vcfRecord.getFormat());

        if (preferredGenotypeFormat != null
                && genotypeDosageFieldPrecedence.contains(preferredGenotypeFormat)
                && formatIdentifiers.contains(this.getGenotypeFormatIdentifier(preferredGenotypeFormat))) {
            return preferredGenotypeFormat;
        }

        if (this.raiseExceptionIfUnavailable) {
            throw new GenotypeDataException(String.format(
                    "Preferred genotype format field (%s) is unavailable for vcf record: %n%s (%s:%s). " +
                            "Available format fields: %s",
                    preferredGenotypeFormatIdentifier,
                    String.join(", ", vcfRecord.getIdentifiers()),
                    vcfRecord.getChromosome(), vcfRecord.getPosition(),
                    String.join(", ", vcfRecord.getFormat())));
        }

        for (VcfGenotypeFormat genotypeFormat: genotypeDosageFieldPrecedence) {
            if (formatIdentifiers.contains(this.getGenotypeFormatIdentifier(genotypeFormat))) {
                return genotypeFormat;
            }
        }

        return null;
    }

    public VcfGenotypeFormat getPreferredGenotypeFormat() {
        return preferredGenotypeFormat;
    }

    public void setPreferredGenotypeFormat(VcfGenotypeFormat preferredGenotypeFormat) {
        this.preferredGenotypeFormat = preferredGenotypeFormat;
    }

    public String getPreferredGenotypeFormatIdentifier() {
        return preferredGenotypeFormatIdentifier;
    }

    public void setPreferredGenotypeFormatIdentifier(String preferredGenotypeFormatIdentifier) {
        this.preferredGenotypeFormatIdentifier = preferredGenotypeFormatIdentifier;
    }

    /**
     * Returns the format identifier corresponding to the genotype format.
     *
     * @param genotypeFormat The genotype format to return a string identifier for.
     * @return The string identifier that matches the genotype format.
     * returns the given synonymous identifier for the preferred genotype format if
     * this is set and given. Otherwise returns the string representation of the genotype format.
     */
    public String getGenotypeFormatIdentifier(VcfGenotypeFormat genotypeFormat) {
        String formatIdentifier;
        if (genotypeFormat.equals(this.preferredGenotypeFormat)) {
            formatIdentifier = this.preferredGenotypeFormatIdentifier;
        } else {
            formatIdentifier = genotypeFormat.toString();
        }

        return formatIdentifier;
    }
}
