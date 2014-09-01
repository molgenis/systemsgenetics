/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.plink;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeData;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeWriter;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.util.Utils;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.NotASnpException;

/**
 *
 * @author Patrick Deelen
 */
public class BedBimFamGenotypeWriter implements GenotypeWriter {

    private static final byte MAGIC_NUMBER_1 = 108;
    private static final byte MAGIC_NUMBER_2 = 27;
    private static final byte MODE = 1; //We only write snp major mode
    private static final int HOMOZYGOTE_SECOND_BITMASK = 192;
    private static final int HETEROZYGOTE_BITMASK = 128;
    private static final int MISSING_BIT_MASK = 64;
    private static final Charset FILE_ENCODING = Charset.forName("UTF-8");
    private static final char SEPARATOR = ' ';
    private static final DecimalFormat PHENO_FORMATTER = new DecimalFormat("0.#####");
    private final GenotypeData genotypeData;
    private int writtenSamplesCounter;
    private int writtenVariantsCounter;
    private int excludedVariantsCounter;
    private static final Logger LOGGER = Logger.getLogger(BedBimFamGenotypeWriter.class);

    public BedBimFamGenotypeWriter(GenotypeData genotypeData) {
        this.genotypeData = genotypeData;
    }

    @Override
    public void write(String path) throws IOException, NotASnpException {
        write(new File(path + ".bed"), new File(path + ".bim"), new File(path + ".fam"));
    }

    public void write(File bedFile, File bimFile, File famFile) throws IOException {

        if (bedFile == null) {
            throw new IllegalArgumentException("No bed file specified to write to");
        }
        if (bimFile == null) {
            throw new IllegalArgumentException("No bim file specified to write to");
        }
        if (famFile == null) {
            throw new IllegalArgumentException("No fam file specified to write to");
        }

        writtenSamplesCounter = 0;
        writtenVariantsCounter = 0;
        excludedVariantsCounter = 0;

        writeBimBedFile(bimFile, bedFile);
		writeFamFile(famFile);

        LOGGER.info("Binary plink data write completed.\n"
                + " - Number of samples: " + writtenSamplesCounter + "\n"
                + " - Number of SNPs: " + writtenVariantsCounter + "\n"
                + " - Excluded non biallelic SNPs: " + excludedVariantsCounter);


    }

    private void writeFamFile(File famFile) throws IOException {
        Utils.createEmptyFile(famFile, "fam");

        BufferedWriter famFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(famFile), FILE_ENCODING));

        for (Sample sample : genotypeData.getSamples()) {
            famFileWriter.append(sample.getFamilyId() != null ? sample.getFamilyId() : "0");
            famFileWriter.append(SEPARATOR);
            famFileWriter.append(sample.getId());
            famFileWriter.append(SEPARATOR);
            famFileWriter.append(sample.getFatherId());
            famFileWriter.append(SEPARATOR);
            famFileWriter.append(sample.getMotherId());
            famFileWriter.append(SEPARATOR);
            famFileWriter.append(Byte.toString(sample.getSex().getPlinkSex()));
            famFileWriter.append(SEPARATOR);
            famFileWriter.append(PHENO_FORMATTER.format(getPhenotype(sample)));
            famFileWriter.append('\n');

            ++writtenSamplesCounter;
        }

        famFileWriter.close();

    }

	private void writeBimBedFile(File bimFile, File bedFile) throws IOException {
		Utils.createEmptyFile(bimFile, "bim");
		Utils.createEmptyFile(bedFile, "bed");

        BufferedWriter bimFileWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(bimFile), FILE_ENCODING));
		
		BufferedOutputStream bedStreamWriter;
        FileOutputStream stream;
        try {
            stream = new FileOutputStream(bedFile);
            bedStreamWriter = new BufferedOutputStream(stream);
        } catch (FileNotFoundException ex) {
            throw new RuntimeException("This should never happen. File will be present at this time", ex);
        }

        bedStreamWriter.write(MAGIC_NUMBER_1);
        bedStreamWriter.write(MAGIC_NUMBER_2);
        bedStreamWriter.write(MODE);

        for (GeneticVariant variant : genotypeData) {

			Alleles variantAlleles = variant.getVariantAlleles();
			
            if (variantAlleles.getAlleleCount() > 2 || !variantAlleles.isSnp()) {
                LOGGER.warn("Skipping variant: " + variant.getPrimaryVariantId() + ", it is not a biallelic SNP.");
                ++excludedVariantsCounter;
                continue;
            }
            
            if (variantAlleles.getAlleleCount() == 0) {
                LOGGER.warn("Skipping variant: " + variant.getPrimaryVariantId() + ", this SNP has no alles.");
                ++excludedVariantsCounter;
                continue;
            }

			bimFileWriter.append(variant.getSequenceName());
            bimFileWriter.append(SEPARATOR);
            bimFileWriter.append(variant.getPrimaryVariantId() == null ? variant.getSequenceName() + ":" + variant.getStartPos() : variant.getPrimaryVariantId());
            bimFileWriter.append(SEPARATOR);
            bimFileWriter.append('0');
            bimFileWriter.append(SEPARATOR);
            bimFileWriter.append(String.valueOf(variant.getStartPos()));
            bimFileWriter.append(SEPARATOR);
            bimFileWriter.append(variantAlleles.getAlleleCount() == 0 ? Allele.ZERO.toString() : variantAlleles.get(0).toString());
            bimFileWriter.append(SEPARATOR);
            bimFileWriter.append(variantAlleles.getAlleleCount() <= 1 ? Allele.ZERO.toString() : variantAlleles.get(1).toString());
            bimFileWriter.append('\n');			

            Alleles homozygoteFirst = Alleles.createAlleles(variantAlleles.get(0), variantAlleles.get(0));
            Alleles homozygoteSecond = null;
            if (variantAlleles.getAlleleCount() == 2) {
                homozygoteSecond = Alleles.createAlleles(variantAlleles.get(1), variantAlleles.get(1));
            }

            int currentByte = 0; //Bit operations are on int level, but we only write the last byte
            byte counterCurrentByte = 0;

            for (Alleles alleles : variant.getSampleVariants()) {
                if (alleles == homozygoteFirst) {
                    //Do nothing, already 00
                } else if (variantAlleles.getAlleleCount() == 2 && alleles.sameAlleles(variantAlleles)) {
                    currentByte = currentByte | HETEROZYGOTE_BITMASK;
                } else if (variantAlleles.getAlleleCount() == 2 && alleles == homozygoteSecond) {
                    currentByte = currentByte | HOMOZYGOTE_SECOND_BITMASK;
                } else if (alleles.contains(Allele.ZERO)) {
                    currentByte = currentByte | MISSING_BIT_MASK;
                } else {
                    throw new GenotypeDataException("Trying to write alleles " + alleles.getAllelesAsString() + " for " + variantAlleles + " SNP");
                }
                ++counterCurrentByte;
                if (counterCurrentByte == 4) {
                    bedStreamWriter.write(currentByte);
                    currentByte = 0;
                    counterCurrentByte = 0;
                } else {
                    currentByte = currentByte >>> 2;
                }
            }

            if (counterCurrentByte != 0) {
                while (counterCurrentByte < 3) {
                    ++counterCurrentByte;
                    currentByte = currentByte >>> 2;
                }
                bedStreamWriter.write(currentByte);
            }
			
			++writtenVariantsCounter;

        }

		bimFileWriter.close();
        bedStreamWriter.close();
        stream.close();
		
	}

    private double getPhenotype(Sample sample) {

        Object value = sample.getAnnotationValues().get(GenotypeData.DOUBLE_PHENOTYPE_SAMPLE_ANNOTATION_NAME);
        if (value == null) {
            return -9;
        }

        if (value instanceof Double) {
            return (Double) value;
        }

        return Double.valueOf(value.toString());
    }
}
