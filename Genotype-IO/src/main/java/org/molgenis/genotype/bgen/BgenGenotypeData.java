/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import com.facebook.presto.orc.zstd.ZstdDecompressor;
import com.google.common.math.IntMath;
import org.apache.log4j.Logger;
import org.molgenis.genotype.*;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.util.RecordIteratorCreators;
import org.molgenis.genotype.variant.*;
import org.molgenis.genotype.variant.range.GeneticVariantRange;
import org.molgenis.genotype.variant.sampleProvider.*;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * @author Patrick Deelen
 */
public class BgenGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantProviderBgen {

    public enum BlockRepresentation {
        compression_0, compression_1, compression_2
    }

    public enum Layout {
        layOut_1, layOut_2
    }

    /**
     * Index is the number of bits used from the last byte. (these are the rightmost bits)
     */
    private static final int[] LAST_BYTE_MASK = {255, 1, 3, 7, 15, 31, 63, 127, 255};

    /**
     * Index is the number of bits used from the first byte. (these are the leftmost bits)
     */
    private static final int[] FIRST_BYTE_MASK = {0, 128, 192, 224, 240, 248, 252, 254, 255};

    private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL = 0.4f;
    private static final Logger LOGGER = Logger.getLogger(BgenGenotypeData.class);
    private static final Charset CHARSET = Charset.forName("UTF-8");

    private final RandomAccessFile bgenFile;
    private final byte[] byteArray4 = new byte[4]; //resuable 4 byte array
    private final byte[] byteArray2 = new byte[2]; //resuable 2 byte array
    private final List<Sample> samples = new ArrayList<>();
    private final Inflater gzipInflater = new Inflater();
    private final ZstdDecompressor zstdInflater = new ZstdDecompressor();
    private final LinkedHashSet<String> sequenceNames = new LinkedHashSet<String>();
    private final BlockRepresentation snpBlockRepresentation;
    private final Layout fileLayout;
    private boolean sampleIdentifiersPresent;
    ;
    private final SampleVariantProviderBgen sampleVariantProvider;
    private final BgenixReader bgenixReader;
    private final double minimumPosteriorProbabilityToCall;
    private final int sampleVariantProviderUniqueId;
    private final int sampleCount;

    public BgenGenotypeData(File bgenFile) throws IOException {
        this(bgenFile, 1000);
    }

    public BgenGenotypeData(File bgenFile, File bgenixFile) throws IOException {
        this(bgenFile, bgenixFile, 1000);
    }

    public BgenGenotypeData(File bgenFile, int cacheSize) throws IOException {
        this(bgenFile, cacheSize, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
    }

    public BgenGenotypeData(File bgenFile, File bgenixFile, int cacheSize) throws IOException {
        this(bgenFile, bgenixFile, cacheSize, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
    }

    public BgenGenotypeData(File bgenFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException {
        this(bgenFile, new File(bgenFile.getAbsolutePath() + ".bgi"),
                cacheSize, minimumPosteriorProbabilityToCall);
    }

    public BgenGenotypeData(File bgenFile, File bgenixFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException {
        this.minimumPosteriorProbabilityToCall = minimumPosteriorProbabilityToCall;
        this.bgenFile = new RandomAccessFile(bgenFile, "r");

        // Get offset of variants in the first four bytes.
        long snpOffset = readFourBytesAsUInt32(
                "SNP offset",
                "Error reading bgen file header. File is corrupt");

        // Get offset the size of the header in the following four bytes.
        long headerSize = readFourBytesAsUInt32(
                "Header size (Lh)",
                "Error reading bgen file header. File is corrupt");

        // Throw an exception if the size of the header is smaller than.
        if (headerSize > snpOffset) {
            throw new GenotypeDataException(
                    "Error reading bgen file header. Header information is bigger than expected offset.");
        }

        // Get the number of variants in the file.
        long variantCount = readFourBytesAsUInt32(
                "Number of SNPs",
                "Error reading bgen file header. File is corrupt");

//               Not needed when using the bgenix indexing file.
//		if (snpCount > Integer.MAX_VALUE) {
//			throw new GenotypeDataException("Found more than (2^31)-1 SNPs in bgen file. This is not supported");
//		}

        // Get the number of samples in the file.
        long sampleCountLong = readFourBytesAsUInt32(
                "Number of samples",
                "Error reading bgen file header. File is corrupt");

        // Throw an exception if the number of samples exceeds the maximum value of an integer.
        if (sampleCountLong > Integer.MAX_VALUE) {
            throw new GenotypeDataException("Found more than (2^31)-6 samples in bgen file. This is not supported");
        }

        this.sampleCount = (int) sampleCountLong;

        // Magic number is in the next four bytes but skipped over.

        // Skip over reserved and free area,
        // this seek works because there are 4 (interesting) bytes
        // in the beginning of the file and in the end of the header block.
        this.bgenFile.seek(headerSize);

        //Read flags
        if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
            throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
        }

        // Changed storing of version and data.
        // Read the SNP block representation (compression state).
        snpBlockRepresentation = readSnpBlockRepresentation();

        // Read the file layout.
        fileLayout = readFileLayout(byteArray4);

        // Throw an exception if the SNP block representation is 2 while layout 1 is used.
        if (snpBlockRepresentation == BlockRepresentation.compression_2 && fileLayout == Layout.layOut_1) {
            throw new GenotypeDataException("Invalid compression method for layout one observed. Trying to use ZSTD compression on layout one file, which is not supported.");
        }

        // Read the sample identifiers presence; set sampleIdentifiersPresent to true if present, false if absent.
        sampleIdentifiersPresent = readSampleIdentifiersPresence(byteArray4[3]);

        if (sampleIdentifiersPresent) {
            // Process sample identifier block.
            processSampleIdentifierBlock(snpOffset, headerSize);
        } else {
            // Set dummy samples.
            IntStream.range(0, sampleCount).boxed()
                    .map(i -> new Sample(i.toString(), null, null))
                    .collect(Collectors.toCollection(() -> samples));
        }

        // Get the start of the variant data block
        long pointerFirstSnp = snpOffset + 4;

        // Check if BGENIX file is present.
        // If it is non-existent, make a new one.
        if (!bgenixFile.exists()) {
            LOGGER.info("Creating bgenix file at: " + bgenixFile.getAbsolutePath());
            BgenixWriter bgenixWriter = new BgenixWriter(bgenixFile);
            createBgenixFile(
                    bgenFile,
                    bgenixWriter,
                    pointerFirstSnp);

            bgenixWriter.finalizeIndex();
        }

        // Not sure what bgenix reader should do in this object other than read an existing BGENIX file.
        bgenixReader = new BgenixReader(bgenixFile);
        readExistingBgenixFile(bgenFile);

        sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();
        // If the specified cache size is greater than 0, construct a new SampleVariantProvider
        if (cacheSize > 0) {
            sampleVariantProvider = new CachedSampleVariantProviderBgen(this, cacheSize);
        } else {
            sampleVariantProvider = this;
        }
    }

    /**
     * Method responsible for reading a BGENIX file. This BGENIX file stores indexes of variants
     * for quick random access to genotype data.
     *
     * @param bgenFile The BGENIX file to read.
     * @throws IOException if an I/O error occurs
     */
    private void readExistingBgenixFile(File bgenFile) throws IOException {
        BgenixMetadata metadata = bgenixReader.getMetadata();
        if (metadata != null) {
            // Check if the metadata of the read BGENIX file is equal to
            // the metadata of the BGEN file.
            if (!metadata.getFileName().equals(bgenFile.getName())) {
                throw new GenotypeDataException("Sample name between bgenix and bgen is not equal. Invalid Bgen and Bgenix combination.");
            }

            // Check if the BGEN file length corresponds to the expected file length in the BGENIX file.
            if (metadata.getFileSize() != bgenFile.length()) {
                throw new GenotypeDataException("Number of expected bytes is different to the numer of observed bytes. Invalid Bgen and Bgenix combination.");
            }

            // Read the first 1000 bytes to check if this is equal to that in the metadata in the
            // BGENIX file.
            this.bgenFile.seek(0);
            byte[] firstBytes = new byte[1000];
            this.bgenFile.read(firstBytes, 0, 1000);
            if (!Arrays.equals(metadata.getFirst1000bytes(), firstBytes)) {
                throw new GenotypeDataException("First 1000 bytes of meta data and actual data are not equal. Invalid Bgen and Bgenix combination.");
            }
        } else {
            // Also throw an exception whenever the BGENIX files metadata is null as returned by the getMetadata method
            // in the bgenixReader object.
            // Not sure if we're ok with null?
            throw new GenotypeDataException("Metadata could not be obtained from the BGENIX file.");
        }
    }

    /**
     * Method that processes the sample identifier block of the BGEN file.
     *
     * @param snpOffset  The number of bytes the variant data is offset from the start of the file header.
     * @param headerSize The size of the header within the BGEN file.
     * @throws IOException if an I/O error has occurred.
     */
    private void processSampleIdentifierBlock(long snpOffset, long headerSize) throws IOException {
        // Get the length of the sample identifier data block.
        long byteSizeSampleIds = readFourBytesAsUInt32(
                "Sample identifier block size",
                "Error in sample identifier block. File is corrupt.");

        // Throw an exception if the sample identifier block size plus the header size is more than the snp offset
        if ((byteSizeSampleIds + headerSize) > snpOffset) {
            throw new GenotypeDataException(
                    "Error reading bgen file header. Combination of header & sample id information is bigger than expected offset.");
        }

        // Read the number of samples represented in this file
        long sampleCountFromSampleIdBlock = readFourBytesAsUInt32(
                "Number of samples with sample-ids",
                "Error in sample identifier block. File is corrupt.");

        // Throw an exception if the number of samples given by the metadata in the sample identifier block and the
        // header are not equal
        if (sampleCountFromSampleIdBlock != sampleCount) {
            throw new GenotypeDataException("Number of samples in metadata and sample id data is not equal. File is corrupt.");
        }

        // Initialize an array of Strings representing the sample identifiers
        String[] sampleIds = new String[(int) sampleCountFromSampleIdBlock];

        // Read the sample identifiers sample by sample.
        for (int i = 0; i < sampleCountFromSampleIdBlock; i++) {
            // Read the sample id length within the next two bytes
            if (this.bgenFile.read(byteArray2, 0, 2) != 2) {
                throw new GenotypeDataException("Error in sample Id. File is corrupt.");
            }
            int sampleIdLength = getUInt16(byteArray2, 0);
            // Initialize the sample identifier.
            byte[] sampleName = new byte[sampleIdLength];
            // Read the sample identifier.
            this.bgenFile.read(sampleName, 0, sampleIdLength);
            // Append the sample identifier to the array of sample ids.
            sampleIds[i] = new String(sampleName, CHARSET);
            samples.add(new Sample(sampleIds[i], null, null));
        }
    }

    /**
     * Read the field that represents the presence of sample identifiers within the BGEN file.
     *
     * @param sampleIdentifiersField The field representing the presence of sample identifiers
     *                               within the BGEN file.
     * @return true if the smaple identifiers are present, false if not.
     */
    private boolean readSampleIdentifiersPresence(byte sampleIdentifiersField) {
        if ((sampleIdentifiersField & 128) == 128) {
            LOGGER.debug("SampleIdentifiers present in bgen file");
            return true;
        } else {
            LOGGER.debug("SampleIdentifiers not present in bgen file");
            return false;
        }
    }

    /**
     * Read the field that indicates which layout is used for this specific BGEN file.
     *
     * @param flagArray Byte array with flags of length 4 according to BGEN specifications.
     * @return the layout type.
     */
    private Layout readFileLayout(byte[] flagArray) {
        // Layout flag is located within the first byte.
        byte byte_tmp = flagArray[0];

        // Shift the layout flag, located on bit 2-5 relative to the first bit.
        byte_tmp = (byte) (byte_tmp >> 2); // relocate to bit 0-3.

        // Perform bitwise operation to mask everything outside bit 0-3.
        int layoutFlag = byte_tmp & 7; // mask using 7 (00000111) and convert to an integer.

        // Find the layout that corresponds to the layout flag value.
        switch (layoutFlag) {
            case 1:
                LOGGER.debug("Genotype data is in layout 1 (BGEN 1.1)");
                return Layout.layOut_1;
            case 2:
                LOGGER.debug("Genotype data is in layout 2 (BGEN 1.2 & 1.3)");
                return Layout.layOut_2;
            default:
                throw new GenotypeDataException("Invalid layout, observed: " + (byteArray4[0] & 3));
        }
    }

    /**
     * Read the flag that indicates how the SNP probability block is represented.
     * Data can be zlib compressed, zstd compressed (for layout 2), or uncompressed.
     *
     * @return the representation of the SNP probabilities block.
     */
    private BlockRepresentation readSnpBlockRepresentation() {

        // Perform bitwise operation to mask everything outside the 0, 1 bits.
        int blockRepresentationFlag = byteArray4[0] & 3; // mask with 7 (00000011) and convert to an integer.

        // Find the representation that corresponds to the flag.
        switch (blockRepresentationFlag) {
            case 0:
                LOGGER.debug("Genotype data is not compressed.");
                return BlockRepresentation.compression_0;
            case 1:
                LOGGER.debug("Genotype data is zlib compressed");
                return BlockRepresentation.compression_1;
            case 2:
                LOGGER.debug("Genotype data is zstd compressed");
                return BlockRepresentation.compression_2;
            default:
                throw new GenotypeDataException("Invalid compression method, observed: " + (blockRepresentationFlag));
        }
    }

    /**
     * @param fieldName        The name of the field to read.
     * @param exceptionMessage What to report if the length of bytes read is not equal to four.
     * @return the bytes as a long.
     * @throws IOException if an I/O error has occurred.
     */
    private long readFourBytesAsUInt32(String fieldName, String exceptionMessage) throws IOException {
        // Throw an exception if the read number of bytes is not equal to 4.
        if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
            throw new GenotypeDataException(exceptionMessage);
        }
        // Convert the read bytes to a long
        long value = getUInt32(byteArray4, 0);

        LOGGER.debug(fieldName + ": " + value);
        return value;
    }

    /**
     * Method responsible for creating a BGENIX file using the BGEN file as guide.
     *
     * @param bgen            The bgen file that is being read.
     * @param bgenixWriter    The bgenixWriter to use for creating this BGENIX file.
     * @param pointerFirstSnp The byte index of the file to start reading variants.
     * @throws IOException if an I/O error has occurred.
     */
    private void createBgenixFile(File bgen,
                                  BgenixWriter bgenixWriter,
                                  long pointerFirstSnp) throws IOException {
        // Go to the first byte...
        this.bgenFile.seek(0);

        // Read the first 1000 bytes of the bgen file.
        byte[] firstBytes = new byte[1000];
        this.bgenFile.read(firstBytes, 0, 1000);

        //Add current time in int.
        System.out.println((System.currentTimeMillis() / 1000L));
        // Create and write new metadata.
        BgenixMetadata m = new BgenixMetadata(
                bgen.getName(),
                (int) this.bgenFile.length(),
                (int) (bgen.lastModified() / 1000L),
                firstBytes,
                (System.currentTimeMillis() / 1000L));
        bgenixWriter.writeMetadata(m);

        //Loop through the start of the file
        long variantReadingPosition;
        variantReadingPosition = pointerFirstSnp;
        while ((variantReadingPosition) < bgenFile.length()) {
            //Loop through variants.

            // Read the variant data creating a variant.
            ReadOnlyGeneticVariantBgen variant = readGeneticVariantBgen(variantReadingPosition);

            // Add the read variant to the BGENIX file so that it can quickly be retrieved.
            bgenixWriter.addVariantToIndex(
                    variant,
                    variant.getVariantReadingPosition(),
                    variant.getVariantDataSizeInBytes(),
                    variant.getPrimaryVariantId());

            // Get the position in the file to start reading the next variant.
            variantReadingPosition = variant.getVariantReadingPosition() +
                    variant.getVariantDataSizeInBytes();
        }
    }

    private ReadOnlyGeneticVariantBgen readGeneticVariantBgen(long variantReadingPosition) throws IOException {
        ReadOnlyGeneticVariantBgen variant = processVariantIdentifyingData(variantReadingPosition);
        VariantGenotypeBlockInfo variantGenotypeBlockInfo = extractVariantGenotypeDataBlockInfo();
        variant.setVariantDataSizeInBytes(variantGenotypeBlockInfo.getVariantDataSizeInBytes(variantReadingPosition));
        return variant;
    }

    /**
     * Processes variant identifying data from the BGENIX file, using a the variant start
     * position to read from the right location.
     *
     * @param variantStartPosition The position to start reading the variant from.
     * @return a genetic variant.
     * @throws IOException if an I/O error has occurred.
     */
    private ReadOnlyGeneticVariantBgen processVariantIdentifyingData(long variantStartPosition) throws IOException {
        long filePointer = variantStartPosition;

        // If layout is equal to 1 then the variant identifying data starts with 4 bytes describing the
        // number of individuals within the row.
        if (fileLayout == Layout.layOut_1) {
            // We chose to ignore this apparently ...
            filePointer += 4;
        }

        // Go to the start of the variant to begin reading there.
        this.bgenFile.seek(filePointer);

        // Proposing not to do a buffer search here as the maximum number of possible bytes is very large
        // (16 + 4K + Lid + Lrsid + Lchr + the sum of the allele lengths (maximum of 2^32 for every allele))
        byte[] byteArray2 = new byte[2];
        byte[] byteArray4 = new byte[4];

//		if (snpInfoBufferSize < 20) {
//			throw new GenotypeDataException("Error reading bgen snp data. File is corrupt");
//		}
        // Need to check that it is correct with the block in front of the snp id.

        // Read the variant identifiers
        ArrayList<String> variantIds = new ArrayList<>();
        String snpId = readVariantInfo(byteArray2);
        String snpRsId = readVariantInfo(byteArray2);

        // add the variant identifiers in the variantIds list so that the RSID is the primary variantID
        variantIds.add(snpRsId);
        if (!snpId.isEmpty()) {
            variantIds.add(snpId);
        }

        // Read the sequence identifier
        String seqName = readVariantInfo(byteArray2);
        sequenceNames.add(seqName);

        // Get the position of the variant.
        bgenFile.read(byteArray4);
        int variantPosition = getVariantPosition(byteArray4);

        // Get the alleles for this variant.
        int numberOfAlleles = 2; // Default is two. (layout one)
        if (fileLayout.equals(Layout.layOut_2)) {
            bgenFile.read(byteArray2);
            numberOfAlleles = getUInt16(byteArray2, 0);
        }

        // Read the alleles
        List<String> alleles = new ArrayList<>();
        for (int i = 0; i < numberOfAlleles; i++) {
            bgenFile.read(byteArray4);
            readAllele(byteArray4, alleles);
        }

        // Log this variant
        LOGGER.debug(String.format("reading %s, %s at %d | seq:pos = %s:%d, %d alleles",
                snpRsId, snpId, variantStartPosition, seqName, variantPosition, numberOfAlleles));

        return ReadOnlyGeneticVariantBgen.createVariant(
                variantIds, variantPosition, seqName, sampleVariantProvider, alleles, variantStartPosition);
        // Not providing first allele as a reference allele
        // now in order to test against gen data.
    }

    private String readVariantInfo(byte[] byteArray2) throws IOException {
        int fieldLength;
        byte[] variableByteArray;
        bgenFile.read(byteArray2);
        fieldLength = getUInt16(byteArray2, 0);

        variableByteArray = new byte[fieldLength];
        bgenFile.read(variableByteArray);
        return new String(variableByteArray, CHARSET);
    }

    /**
     * Method for extracting the position of the current variant.
     *
     * @param snpInfoBuffer The byte array buffer starting from the start of the variant block.
     * @return the position of the variant.
     */
    private int getVariantPosition(byte[] snpInfoBuffer) {
        long variantPosLong = getUInt32(snpInfoBuffer, 0);
        // Throw an exception if the variant position exceeds the maximum value of an integer.
        if (variantPosLong > Integer.MAX_VALUE) {
            throw new GenotypeDataException("SNP pos larger than (2^31)-6 not supported");
        }
        return (int) variantPosLong;
    }

    /**
     * Read an allele and add this to the list of alleles
     *
     * @param snpInfoBuffer A byte array buffer starting from the start of the variant block.
     * @param alleles       A list of alleles.
     */
    private void readAllele(byte[] snpInfoBuffer, List<String> alleles) throws IOException {

        // Length of the allele
        long fieldLengthLong = getUInt32(snpInfoBuffer, 0);

        // Throw an exception if the allele length is longer
        if (fieldLengthLong > Integer.MAX_VALUE - 5) {
            throw new GenotypeDataException("SNP with allele longer than (2^31)-6 characters not supported");
        }

        // Create a new buffer with the correct size.
        byte[] alleleByteArray = new byte[(int) fieldLengthLong];
        bgenFile.read(alleleByteArray);
        // Get the allele from the buffer.
        String allele = new String(alleleByteArray, CHARSET);
        alleles.add(allele);
    }

    private double[][] readGenotypeDataFromVariant(ReadOnlyGeneticVariantBgen variant) throws IOException {
        double[][] probabilities = new double[sampleCount][];
        // Get the decompressed variant data of length D
        byte[] variantBlockData = getDecompressedBlockData(variant);

        // Read
        if (fileLayout == Layout.layOut_1) {
            for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
                double[] sampleProbabilities = new double[3];

//				System.out.println(getUInt16(variantBlockData, 0) / 32768f + " " + getUInt16(variantBlockData, 2) / 32768f + " " + getUInt16(variantBlockData, 4) / 32768f);
                sampleProbabilities[0] = getUInt16(variantBlockData, 0) / 32768f;
                sampleProbabilities[1] = getUInt16(variantBlockData, 2) / 32768f;
                sampleProbabilities[2] = getUInt16(variantBlockData, 4) / 32768f;
                probabilities[sampleIndex] = sampleProbabilities;
            }
        } else if (fileLayout.equals(Layout.layOut_2)) {

            int blockBufferOffset = 0;
            //must equal data before.
            if ((int) getUInt32(variantBlockData, blockBufferOffset) != sampleCount) {
                throw new GenotypeDataException(String.format(
                        "BGEN file format error. " +
                                "The variant's sample count (%d) does not match with the header (%d).",
                        sampleCount, (int) getUInt32(variantBlockData, blockBufferOffset)));
            }
            blockBufferOffset += 4;

            //must equal data before.
            int numberOfAlleles = getUInt16(variantBlockData, blockBufferOffset);
            blockBufferOffset += 2;

            // Read the minimum ploidy
            int minPloidy = getUInt8(variantBlockData, blockBufferOffset);
            blockBufferOffset += 1;

            // Read the maximum ploidy
            int maxPloidy = getUInt8(variantBlockData, blockBufferOffset);
            blockBufferOffset += 1;

            // Loop through every individual
            List<Boolean> isMissing = new ArrayList<>();
            List<Integer> ploidies = new ArrayList<>();
            ReadPloidiesAndMissingnessByte(variantBlockData, blockBufferOffset, isMissing, ploidies);
            // Add the number of individuals to the buffer.
            blockBufferOffset += sampleCount;

            boolean phased = isPhased(variantBlockData, blockBufferOffset);
            blockBufferOffset += 1;

            int probabilitiesLengthInBits = getUInt8(variantBlockData, blockBufferOffset);
            blockBufferOffset += 1;

            LOGGER.debug(String.format("%s | alleles: %d, ploidy: %d-%d, %s, %d bit representation",
                    variant.getPrimaryVariantId(), numberOfAlleles, minPloidy, maxPloidy,
                    phased ? "phased" : "unphased", probabilitiesLengthInBits));

            byte[] probabilitiesArray = Arrays.copyOfRange(
                    variantBlockData, blockBufferOffset,
                    variantBlockData.length);

            if (phased) {
                double[][][] haplotypeProbabilities = readHaplotypeProbabilities(
                        probabilitiesArray,
                        probabilitiesLengthInBits,
                        numberOfAlleles,
                        isMissing, ploidies);
                probabilities = convertHaplotypeProbabilitiesToGenotypeProbabilities(numberOfAlleles, ploidies, haplotypeProbabilities);
            } else {
                probabilities = readGenotypeProbabilities(
                        probabilitiesArray,
                        probabilitiesLengthInBits,
                        numberOfAlleles,
                        isMissing, ploidies);
            }
        }
        return probabilities;
    }

    private double[][][] readPhasedGenotypeDataFromVariant(ReadOnlyGeneticVariantBgen variant) throws IOException {
        if (fileLayout != Layout.layOut_2) {
            throw new GenotypeDataException("Phased data not available");
        }
        // Get the decompressed variant data of length D
        byte[] variantBlockData = getDecompressedBlockData(variant);

        int blockBufferOffset = 0;
        //must equal data before.
        if ((int) getUInt32(variantBlockData, blockBufferOffset) != sampleCount) {
            throw new GenotypeDataException(String.format(
                    "BGEN file format error. " +
                            "The variant's sample count (%d) does not match with the header (%d).",
                    sampleCount, (int) getUInt32(variantBlockData, blockBufferOffset)));
        }
        blockBufferOffset += 4;

        //must equal data before.
        int numberOfAlleles = getUInt16(variantBlockData, blockBufferOffset);
        blockBufferOffset += 2;

        // Read the minimum ploidy
        int minPloidy = getUInt8(variantBlockData, blockBufferOffset);
        blockBufferOffset += 1;

        // Read the maximum ploidy
        int maxPloidy = getUInt8(variantBlockData, blockBufferOffset);
        blockBufferOffset += 1;

        // Loop through every individual
        List<Boolean> isMissing = new ArrayList<>();
        List<Integer> ploidies = new ArrayList<>();
        ReadPloidiesAndMissingnessByte(variantBlockData, blockBufferOffset, isMissing, ploidies);
        // Add the number of individuals to the buffer.
        blockBufferOffset += sampleCount;

        boolean phased = isPhased(variantBlockData, blockBufferOffset);
        if (!phased) {
            throw new GenotypeDataException("Phased data not available");
        }
        blockBufferOffset += 1;

        int probabilitiesLengthInBits = getUInt8(variantBlockData, blockBufferOffset);
        blockBufferOffset += 1;

        LOGGER.debug(String.format("%s | alleles: %d, ploidy: %d-%d, %s, %d bit representation",
                variant.getPrimaryVariantId(), numberOfAlleles, minPloidy, maxPloidy,
                "phased", probabilitiesLengthInBits));

        byte[] probabilitiesArray = Arrays.copyOfRange(
                variantBlockData, blockBufferOffset,
                variantBlockData.length);

        return readHaplotypeProbabilities(
                probabilitiesArray,
                probabilitiesLengthInBits,
                numberOfAlleles,
                isMissing, ploidies);
    }

    private double[][] convertHaplotypeProbabilitiesToGenotypeProbabilities(int numberOfAlleles, List<Integer> ploidies, double[][][] haplotypeProbabilities) {
        double[][] probabilities;
        probabilities = new double[haplotypeProbabilities.length][];
        // Calculate the probability for homozygous genotype 'AA',
        // the probability for 'AB', and the probability for 'BB'
        for (int sampleIndex = 0; sampleIndex < haplotypeProbabilities.length; sampleIndex++) {
            // Convert the probabilities for this particular sample
            probabilities[sampleIndex] = phasedProbabilitiesToGenotypeProbabilities(
                    haplotypeProbabilities[sampleIndex], ploidies.get(sampleIndex), numberOfAlleles);
        }
        return probabilities;
    }

    private double[] computeApproximateProbabilities(byte[] probabilitiesArray, int probabilitiesLengthInBits, int bitOffset, int numberOfProbabilities) {
        // Keep track of the probabilities
        double[] probabilities = new double[numberOfProbabilities];
        // Keep track of the sum of all probabilities as the last value is calculated by
        // subtracting this sum from 1.
        double sumOfProbabilities = 0;
        for (int j = 0; j < numberOfProbabilities - 1; j++) {
            // Get the probability

            // Each probability is stored in B bits.
            // Values are interpreted by linear interpolation between 0 and 1;
            // value b corresponds to probability b / ((2^B)-1).

            // To interpret a stored value x as a probability:
            // Convert x to an integer in floating point representation (double because values are stored with the double precision)
            double probability = readProbabilityValue(probabilitiesArray, bitOffset, probabilitiesLengthInBits);

            // Divide by (2^B)-1
            double maxValue = Math.pow(2, probabilitiesLengthInBits) - 1;
            probability = probability / maxValue;
            probabilities[j] = probability;
            sumOfProbabilities += probability;
            // Update the offset of the current bit
            bitOffset += probabilitiesLengthInBits;
        }
        // Calculate probability of Kth allele (P_iK)
        probabilities[numberOfProbabilities - 1] = 1 - sumOfProbabilities;
        return probabilities;
    }

    private double[][] readGenotypeProbabilities(
            byte[] probabilitiesArray,
            int probabilitiesLengthInBits,
            int numberOfAlleles, List<Boolean> isMissing,
            List<Integer> ploidies) {

        // Get bit offset
        int bitOffset = 0;

        // Initialize an array of probabilities.
        double[][] probabilities = new double[getSamples().size()][];

        for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {

            // Get the list of combinations that should have probabilities in this probability-block.
            List<List<Integer>> combinations = getGenotypeCombinations(
                    numberOfAlleles,
                    ploidies.get(sampleIndex), true);

            // If the probabilities are missing for this sample, read zero and continue with the
            // next sample.
            if (isMissing.get(sampleIndex)) {
                // If this is missing, the probability is zero.
                bitOffset += probabilitiesLengthInBits * (combinations.size() - 1);
                probabilities[sampleIndex] = new double[combinations.size()];
                continue;
            }

            // Get all probabilities for this sample.
            double[] genotypeProbabilities = computeApproximateProbabilities(
                    probabilitiesArray, probabilitiesLengthInBits, bitOffset, combinations.size());
            probabilities[sampleIndex] = genotypeProbabilities;

            // Update the bit to read the next probabilities from.
            bitOffset += probabilitiesLengthInBits * (combinations.size() - 1);
        }
        return probabilities;
    }

    private double[][][] readHaplotypeProbabilities(
            byte[] probabilitiesArray,
            int probabilitiesLengthInBits,
            int numberOfAlleles, List<Boolean> isMissing,
            List<Integer> ploidies) {

        // Define an array consisting of an array of posterior probabilities for each genotype
        double[][][] haplotypeProbabilities = new double[getSamples().size()][][];

        // Get bit offset
        int bitOffset = 0;

        // Each probability is stored in B bits.
        // Values are interpreted by linear interpolation between 0 and 1;
        // value b corresponds to probability b / ((2^B)-1).

        for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
            Integer ploidy = ploidies.get(sampleIndex);
            // If the probabilities are missing for this sample, read zero and continue with the
            // next sample.
            if (isMissing.get(sampleIndex)) {
                // If this is missing, the probability is zero.
                bitOffset += probabilitiesLengthInBits * (ploidy * (numberOfAlleles - 1));
                haplotypeProbabilities[sampleIndex] = new double[ploidy][numberOfAlleles];
                continue;
            }

            double[][] phasedSampleProbabilities = new double[ploidy][numberOfAlleles];

            for (int i = 0; i < ploidy; i++) {
                // Get the probabilities for every allele in this haplotype.
                double[] alleleProbabilities = computeApproximateProbabilities(
                        probabilitiesArray, probabilitiesLengthInBits, bitOffset, numberOfAlleles);
                bitOffset += ((numberOfAlleles - 1) * probabilitiesLengthInBits);

                phasedSampleProbabilities[i] = alleleProbabilities;
            }
            haplotypeProbabilities[sampleIndex] = phasedSampleProbabilities;
        }
        return haplotypeProbabilities;
    }

    /**
     * Converts an array of arrays with probabilities per allele for particular haplotypes, to
     * regular genotype probability values per genotype.
     *
     * @param phasedProbabilities An array, with its size equal to the number haplotypes, with arrays
     *                            with a probability for every possible allele for the haplotype.
     * @param ploidy              The ploidy for this particular sample.
     * @param numberOfAlleles     The number of alleles.
     * @return The number of probabilities for every possible genotype.
     */
    private double[] phasedProbabilitiesToGenotypeProbabilities(double[][] phasedProbabilities,
                                                                Integer ploidy,
                                                                int numberOfAlleles) {
        // Get all possible combinations of alleles for the haplotypes (all permutations)
        List<List<Integer>> haplotypeCombinations = getGenotypeCombinations(
                numberOfAlleles, ploidy, false);

        // Initialize an array of probabilities.
        double[] genotypeProbabilities = new double[numberOfOrderedCombinationsWithRepetition(ploidy,
                numberOfAlleles - 1)]; // number of alleles is equal to n, n-1 = numberOfAlleles - 1

        // Loop through the combinations of haplotypes
        for (List<Integer> haplotypeCombination : haplotypeCombinations) {
            double probability = 1;
            // Multiply the probabilities that correspond to the indices within the haplotype combination
            for (int haplotypeIndex = 0; haplotypeIndex < haplotypeCombination.size(); haplotypeIndex++) {
                probability *= phasedProbabilities[haplotypeIndex][haplotypeCombination.get(haplotypeIndex)];
            }
            // Add the multiplied probability of this combination to the other combinations with the same genotype
            Collections.sort(haplotypeCombination);
            genotypeProbabilities[getIndexOfOrderedCombination(
                    Collections.unmodifiableList(haplotypeCombination))] += probability;
        }
        return genotypeProbabilities;
    }

    /**
     * Static method that returns the index of a particular genotype combination within
     * all the possible genotypes for a particular variant and individual (ploidy).
     *
     * @param combination The ordered combination sequence of allele indices. (with increasing value from left to right)
     * @return the index that corresponds to the given ordered combination.
     */
    private static int getIndexOfOrderedCombination(List<Integer> combination) {
        // Initiate the index as 0.
        int index = 0;
        // Loop through the values of the combinations, starting from the last value.
        for (int r = combination.size(); r > 0; r--) {
            // Get n - 1 and r in order to calculate the number of possible ordered combinations, with repetition.
            int nMinOne = combination.get(r - 1) - 1;
            // If n minus 1 is below zero, return the current index, as other values in the combinations left of this
            // are also zero.
            if (nMinOne < 0) {
                return index;
            }
            // Calculate the number of possible ordered combinations, with repetition, and add this to the index value.
            index += numberOfOrderedCombinationsWithRepetition(r, nMinOne);
        }
        return index;
    }

    private static int numberOfOrderedCombinationsWithRepetition(int r, int nMinOne) {
        return IntMath.factorial(r + nMinOne) /
                (IntMath.factorial(r) * IntMath.factorial(nMinOne));
    }

    /**
     * Read the ploidity and missingness.
     *
     * @param snpBlockData The byte buffer array to read the genotype data.
     * @param byteOffset   the byte number to start reading ploidies from within the input byte array.
     * @param isMissing    A list of booleans representing if the corresponding sample has missing probabilities
     * @param ploidity     A list representing the ploidy of corresponding samples.
     */
    private void ReadPloidiesAndMissingnessByte(
            byte[] snpBlockData, int byteOffset,
            List<Boolean> isMissing, List<Integer> ploidity) {

        for (int i = byteOffset; i < sampleCount + byteOffset; i++) {
            // Here we need to handle ploidity and missigness of prababilities.
            // Missingness is encoded by the most significant bit (MSB);
            // thus a value of 1 for the most significant bit
            // indicates that no probability data is stored for this sample.
            // Ploidity is encoded by the 6 least significant bits.

            // Get the MSB and check if the result equals 128.
            isMissing.add((snpBlockData[i] & 0x80) == 128);
            // Get the 6 least significant bits.
            ploidity.add(snpBlockData[i] & 0x3F);

//			String s1 = String.format("%8s", Integer.toBinaryString(snpBlockData[i] & 0xFF)).replace(' ', '0');
//			System.out.println("value = " + s1 + " | " + (snpBlockData[i] & 0x3F) + " | " + (snpBlockData[i] & 0x80));
        }
    }

    /**
     * Read if the variant data is phased or not.
     *
     * @param snpBlockData          The byte buffer array to read the genotype data.
     * @param isPhasedFieldPosition The position in the byte buffer array of the field representing the phasing
     * @return true if the data is phased, false if not.
     */
    private boolean isPhased(byte[] snpBlockData, int isPhasedFieldPosition) {
        boolean phased;
        // Get the field with the phasing flag.
        switch (getUInt8(snpBlockData, isPhasedFieldPosition)) {
            case 0:
                phased = false;
                break;
            case 1:
                phased = true;
                break;
            default:
                throw new GenotypeDataException(
                        "Bgen file format error. Unsupported value for the phased flag observed."
                );
        }
        return phased;
    }

    /**
     * Method for obtaining decompressed data for layout 2.
     */
    private byte[] getDecompressedBlockData(ReadOnlyGeneticVariantBgen variant) throws IOException {
        // First extend the genetic variant with ids, alleles, etc.
        variant.extendWithAdditionalVariantData();

        // Extract the variant genotype data block info starting from the current position of the
        // file pointer, which should be right after all alleles for a specific variant.
        VariantGenotypeBlockInfo variantGenotypeBlockInfo = extractVariantGenotypeDataBlockInfo();
        variant.setVariantDataSizeInBytes(variantGenotypeBlockInfo.getVariantDataSizeInBytes(
                variant.getVariantReadingPosition()));

        long decompressedVariantBlockLength = variantGenotypeBlockInfo.getDecompressedBlockLength();
        long variantBlockLength = variantGenotypeBlockInfo.getBlockLength();
        long variantProbabilitiesStartPosition = variantGenotypeBlockInfo.getVariantProbabilitiesStartPosition();

        if (decompressedVariantBlockLength > Integer.MAX_VALUE - 5) {
            throw new GenotypeDataException(String.format(
                    "Length of decompressed genotype data exceeds maximum supported value of (2^31)-6 (%d)",
                    decompressedVariantBlockLength));
        }

        if (variantBlockLength > Integer.MAX_VALUE - 5) {
            throw new GenotypeDataException(String.format(
                    "Length of compressed genotype data exceeds maximum supported value of (2^31)-6 (%d)", variantBlockLength));
        }

        // Initialize byte arrays.
        byte[] compressedBlockData = new byte[(int) variantBlockLength];
        byte[] decompressedBlockData = new byte[(int) decompressedVariantBlockLength];

        // Read the compressed / uncompressed data starting from the correct location.
        this.bgenFile.seek(variantProbabilitiesStartPosition);
        this.bgenFile.read(compressedBlockData, 0, (int) variantBlockLength);

        switch (snpBlockRepresentation) {

            case compression_1:
                decompressVariantBlockGzip(
                        compressedBlockData,
                        decompressedBlockData);
                break;

            case compression_2:
                zstdInflater.decompress(
                        compressedBlockData,
                        0,
                        (int) variantBlockLength,
                        decompressedBlockData,
                        0,
                        (int) decompressedVariantBlockLength);
                break;

            default:
                decompressedBlockData = compressedBlockData;
                break;
        }
        return decompressedBlockData;
    }

    /**
     * Method that decompresses data using the gzipInflater.
     *
     * @param compressedVariantDataBlock The input byte array to decompress.
     * @param outputVariantDataBlock     The decompressed output byte array.
     */
    private void decompressVariantBlockGzip(
            byte[] compressedVariantDataBlock, byte[] outputVariantDataBlock) {

        // Set the input for the gzip inflater.
        gzipInflater.setInput(compressedVariantDataBlock);

        // Try to decompress the data.
        try {
            gzipInflater.inflate(outputVariantDataBlock);
        } catch (DataFormatException e) {
            throw new GenotypeDataException("Error decompressing bgen data", e);
        }
        gzipInflater.reset();
    }

    /**
     * Method that extracts info about the variant genotype data block.
     * Method must be called with the file pointer at the start of a variant genotype block (the C field if present)
     *
     * @return an object with info of this block.
     * @throws IOException if an I/O error has occurred.
     */
    private VariantGenotypeBlockInfo extractVariantGenotypeDataBlockInfo() throws IOException {

        // Get the file pointer (bgenFile.getFilePointer();)
        long variantGenotypeStartPosition = this.bgenFile.getFilePointer();

        // Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
        byte[] snpInfoBuffer = new byte[8];
        this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
        int snpInfoBufferPos = 0;

        VariantGenotypeBlockInfo variantGenotypeBlockInfo = null;

        // Separate methods for layout 1 and 2.
        // Separation like this indicates that subclasses of this reader class should probably be made...
        if (fileLayout == Layout.layOut_1) {
            variantGenotypeBlockInfo = getVariantGenotypeDataBlockInfoForLayoutOne(
                    snpInfoBuffer, snpInfoBufferPos, variantGenotypeStartPosition);
        } else if (fileLayout == Layout.layOut_2) {
            variantGenotypeBlockInfo = determineVariantGenotypeDataBlockSizeForLayoutTwo(
                    snpInfoBuffer, snpInfoBufferPos, variantGenotypeStartPosition);
        }
        return variantGenotypeBlockInfo;
    }

    /**
     * Method for obtaining info of a genotype data block for layout two.
     *
     * @param snpInfoBuffer                A byte array buffer starting from the start of the variant block.
     * @param snpInfoBufferPos             The position of an genotype data block in the snp info buffer.
     * @param variantGenotypeStartPosition The position of the genotype data block in the BGEN file.
     * @return an object with info of this block.
     */
    private VariantGenotypeBlockInfo determineVariantGenotypeDataBlockSizeForLayoutTwo(
            byte[] snpInfoBuffer, int snpInfoBufferPos, long variantGenotypeStartPosition) {
        long variantBlockSize;
        long snpBlockSizeDecompressed;
        // If the file layout is 2,
        // The genotype data for this variant is equal to the next 4 bytes
        variantBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
        // Shift the buffer position 4 bytes
        snpInfoBufferPos += 4;

        // If the probabilities block is compressed read the size of the decompressed data within
        // the next four bytes.
        if (!(snpBlockRepresentation == BlockRepresentation.compression_0)) {
            snpBlockSizeDecompressed = getUInt32(snpInfoBuffer, snpInfoBufferPos);
            return new VariantGenotypeBlockInfo(
                    variantGenotypeStartPosition,
                    fileLayout,
                    variantBlockSize,
                    snpBlockSizeDecompressed, true); // The data is compressed.
        }
        return new VariantGenotypeBlockInfo(
                variantGenotypeStartPosition,
                fileLayout,
                variantBlockSize,
                false); // The data is not compressed.
    }

    /**
     * Method for obtaining info of a genotype data block for layout one.
     *
     * @param snpInfoBuffer                A byte array buffer starting from the start of the variant block.
     * @param snpInfoBufferPos             The position of an genotype data block in the snp info buffer.
     * @param variantGenotypeStartPosition The position of the genotype data block in the BGEN file.
     * @return an object with info of this block.
     */
    private VariantGenotypeBlockInfo getVariantGenotypeDataBlockInfoForLayoutOne(
            byte[] snpInfoBuffer, int snpInfoBufferPos, long variantGenotypeStartPosition) {
        long variantBlockSize;
        if (snpBlockRepresentation == BlockRepresentation.compression_1) {
            // If the file layout is 1 and the variants are zlib compressed
            // The snp block starts in the 4 bytes after the block size field
            variantBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
            return new VariantGenotypeBlockInfo(
                    variantGenotypeStartPosition,
                    fileLayout, variantBlockSize, true);
        } else {
            // If the file layout is 1 and the variants are not zlib compressed,
            // the genotype data for this variant is 6 times the number of samples
            variantBlockSize = 6 * sampleCount;
            return new VariantGenotypeBlockInfo(
                    variantGenotypeStartPosition, fileLayout, variantBlockSize, false);
        }
    }

    /**
     * Getter for the presence of sample identifiers within the BGEN file.
     *
     * @return true if sample identifiers are present in the BGEN file, false if otherwise.
     */
    public boolean areSampleIdentifiersPresent() {
        return sampleIdentifiersPresent;
    }

    /**
     * Getter for the compression type of the genotype data block.
     * Equivalent to the first two bits of the set of <i>flags</i> within the header of the BGEN file.
     *
     * @return the compression type of the genotype data block.
     */
    public BlockRepresentation getGenotypeDataBlockRepresentation() {
        return snpBlockRepresentation;
    }

    /**
     * Getter for the layout type of the BGEN file that corresponds to this object.
     * This is equivalent to the 2-5 bits of the set of <i>flags</i> within the header of the BGEN file.
     *
     * @return the layout type of the BGEN file.
     */
    public Layout getFileLayout() {
        return fileLayout;
    }

    @Override
    public List<Sample> getSamples() {
        return Collections.unmodifiableList(samples);
    }

    @Override
    public Map<String, Annotation> getVariantAnnotationsMap() {
        return Collections.emptyMap();
    }

    @Override
    public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
        return new LinkedHashMap<>();
    }

    @Override
    public boolean isOnlyContaingSaveProbabilityGenotypes() {
        return true;
    }

    @Override
    public void close() throws IOException {
        bgenFile.close();
    }

    @Override
    public List<String> getSeqNames() {
        return new ArrayList<>(sequenceNames);
    }

    @Override
    public Iterable<Sequence> getSequences() {
        List<Sequence> sequences = new ArrayList<>();
        for (String seqName : getSeqNames()) {
            sequences.add(new SimpleSequence(seqName, null, this));
        }
        return sequences;
    }

    @Override
    public Iterator<GeneticVariant> iterator() {
        return getGeneticVariants(bgenixReader.getVariants()).iterator();
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
        return getGeneticVariants(bgenixReader.getVariantsPostion(seqName, startPos));
    }

    @Override
    public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
        return getGeneticVariants(bgenixReader.getVariantsChromosome(seqName));
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
        return getGeneticVariants(bgenixReader.getVariantsRange(seqName, rangeStart, rangeEnd));
    }

    @Override
    public List<Alleles> getSampleVariants(GeneticVariant variant) {
        return ProbabilitiesConvertor.convertProbabilitiesToAlleles(
                variant.getSampleGenotypeProbilities(),
                variant.getVariantAlleles(),
                minimumPosteriorProbabilityToCall);
    }

    @Override
    public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
        return RecordIteratorCreators.createIteratorFromProbs(
                variant.getSampleGenotypeProbilities());
    }

    @Override
    public List<Boolean> getSamplePhasing(GeneticVariant variant) {
        ReadOnlyGeneticVariantBgen bgenVariant = getCastedBgenVariant(variant);
        // Check if the file layout is equal to layout 2
        if (fileLayout != Layout.layOut_2) {
            // If it is not the case, return a list of booleans with a false value.
            return Collections.unmodifiableList(Collections.nCopies(sampleCount, false));
        }
        try {
            // Get the decompressed variant data of length D
            byte[] variantBlockData = getDecompressedBlockData(bgenVariant);
            boolean phased = isPhased(variantBlockData, 8 + sampleCount);
            return Collections.unmodifiableList(Collections.nCopies(sampleCount, phased));
        } catch (IOException e) {
            throw new GenotypeDataException(String.format(
                    "Could not read variant data %s at position %d%n",
                    variant.getPrimaryVariantId(), bgenVariant.getVariantReadingPosition()));
        }
    }

    @Override
    public int cacheSize() {
        return 0;
    }

    @Override
    public int getSampleVariantProviderUniqueId() {
        return sampleVariantProviderUniqueId;
    }

    @Override
    public ReadOnlyGeneticVariantBgen extendReadOnlyGeneticVariantBgen(ReadOnlyGeneticVariantBgen variant) {
        try {
            return processVariantIdentifyingData(variant.getVariantReadingPosition());
        } catch (IOException e) {
            throw new GenotypeDataException(String.format(
                    "Could not read variant data %s at position %d%n",
                    variant.getPrimaryVariantId(), variant.getVariantReadingPosition()));
        }
    }

    @Override
    public byte[] getSampleCalledDosage(GeneticVariant variant) {
        return CalledDosageConvertor.convertCalledAllelesToCalledDosage(variant.getSampleVariants(), variant.getVariantAlleles(), null);
    }

    @Override
    public float[] getSampleDosage(GeneticVariant variant) {
        return ProbabilitiesConvertor.convertProbabilitiesToDosage(variant.getSampleGenotypeProbilities(), DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
    }

    @Override
    public float[][] getSampleProbilities(GeneticVariant variant) {
        // Make sure that probabilities for other than biallelic variants return missingness
        double[][] sampleGenotypeProbabilitiesBgen = getSampleGenotypeProbabilitiesBgen(variant);
        if (variant.isBiallelic()) {
            return ProbabilitiesConvertor.convertBiallelicBgenProbabilitiesToProbabilities(
                    sampleGenotypeProbabilitiesBgen);
        } else {
            return new float[sampleGenotypeProbabilitiesBgen.length][3];
        }
    }

    @Override
    public double[][] getSampleGenotypeProbabilitiesBgen(GeneticVariant variant) {
        ReadOnlyGeneticVariantBgen bgenVariant = getCastedBgenVariant(variant);
        try {
            return readGenotypeDataFromVariant(bgenVariant);
        } catch (IOException e) {
            throw new GenotypeDataException(String.format(
                    "Could not read variant data %s at position %d%n",
                    variant.getPrimaryVariantId(), bgenVariant.getVariantReadingPosition()));
        }
    }

    @Override
    public double[][][] getSampleGenotypeProbabilitiesBgenPhased(GeneticVariant variant) {
        ReadOnlyGeneticVariantBgen bgenVariant = getCastedBgenVariant(variant);
        try {
            return readPhasedGenotypeDataFromVariant(bgenVariant);
        } catch (IOException e) {
            throw new GenotypeDataException(String.format(
                    "Could not read variant data %s at position %d%n",
                    variant.getPrimaryVariantId(), bgenVariant.getVariantReadingPosition()));
        }
    }

    private ReadOnlyGeneticVariantBgen getCastedBgenVariant(GeneticVariant variant) {
        if (!(variant instanceof ReadOnlyGeneticVariantBgen)) {
            throw new GenotypeDataException("Variant is not of type 'ReadOnlyGeneticVariantBgen' and thus cannot" +
                    "be used within this sampleVariantProvider");
        }
        return (ReadOnlyGeneticVariantBgen) variant;
    }

    private Iterable<GeneticVariant> getGeneticVariants(BgenixVariantQueryResult variantQueryResult) {
        GeneticVariantRange.GeneticVariantRangeCreate variantRangeFactory = GeneticVariantRange.createRangeFactory();

        while (variantQueryResult.hasNext()) {
            BgenixVariantData variantData = variantQueryResult.next();
            // Create a read only genetic variant bgen.
            GeneticVariant variant = ReadOnlyGeneticVariantBgen.createVariant(
                    variantData.getRsid(),
                    variantData.getPosition(),
                    variantData.getChromosome(),
                    this,
                    variantData.getNumber_of_alleles(),
                    variantData.getAllele1(), variantData.getAllele2(),
                    variantData.getFile_start_position(),
                    variantData.getSize_in_bytes());
            variantRangeFactory.addVariant(variant);
        }
        return variantRangeFactory.createRange();
    }

    /**
     * Method that returns the list of all possible allele combinations, for the given ploidy.
     * These allele combinations, or genotypes represent the way in which genotypes are stored in a layout 2 BGEN file
     * for a specific sample, if onlyOrdered is true, otherwise also permutations of these combinations are returned.
     * <p>
     * For example, the method returns 4 combinations with 2 alleles and a ploidy of 2
     * in case of all combinations:
     * [1, 1]
     * [2, 1]
     * [1, 2]
     * [2, 2]
     * The method returns 3 combinations with 2 alleles and a ploidy of 2
     * in case only ordered combinations are requested:
     * [1, 1]
     * [1, 2]
     * [2, 2]
     *
     * @param numberOfAlleles The number of alleles for the variant to create combinations for.
     * @param ploidy          The ploidity of the sample to create combinations for.
     * @param onlyOrdered     Flag indicating if the combinations should only be the sorted combinations according to the
     *                        BGEN file specifications.
     * @return the combinations of alleles that represent all possible genotypes or haplotypes
     * for an unphased variant and sample, in this specified order.
     */
    private static List<List<Integer>> getGenotypeCombinations(int numberOfAlleles, int ploidy, boolean onlyOrdered) {
        // Construct nested lists
        List<List<Integer>> combinations = new ArrayList<>();

        // Set the maximum value of an allele, which is the number of alleles minus one, because we want to count
        // from 0 to n-1
        int maxAlleleValue = numberOfAlleles - 1;

        // Get the combinations
        if (onlyOrdered) {
            getCombinationsRecursively(combinations,
                    new ArrayList<>(),
                    maxAlleleValue,
                    ploidy);
        } else {
            getAllCombinationsRecursively(combinations,
                    new ArrayList<>(),
                    maxAlleleValue,
                    ploidy);
        }

        return combinations;
    }

    /**
     * Method that recursively fills a list of combinations.
     * Every permutation of a possible combination are returned.
     *
     * @param combinations   The list of combinations to fill.
     * @param combination    The current combination that is being constructed.
     * @param maxAlleleValue The number of different values to fit into a combinations.
     * @param ploidy         The size of a combination.
     */
    private static void getAllCombinationsRecursively(
            List<List<Integer>> combinations, List<Integer> combination, int maxAlleleValue, int ploidy) {
        // If the combination is complete, the size of the combination equals the required size.
        // Add the combination and return
        if (combination.size() == ploidy) {
            combinations.add(combination); // Add in the beginning to maintain the correct order.
            return;
        }

        // Loop through the possible values from high to low.
        for (int alleleValue = 0; alleleValue <= maxAlleleValue; alleleValue++) {
            // Copy the preliminary combination
            List<Integer> newCombination = new ArrayList<>(combination);
            // Add a new value to the combination
            newCombination.add(alleleValue); // Add in the beginning to maintain the correct order.
            getAllCombinationsRecursively(combinations, newCombination, maxAlleleValue, ploidy);
        }
    }

    /**
     * Method that recursively fills a list of combinations.
     * Combinations are always ordered.
     *
     * @param combinations   The list of combinations to fill.
     * @param combination    The current combination that is being constructed.
     * @param maxAlleleValue The number of different values to fit into a combinations.
     * @param ploidy         The size of a combination.
     */
    private static void getCombinationsRecursively(
            List<List<Integer>> combinations, List<Integer> combination, int maxAlleleValue, int ploidy) {
        // If the combination is complete, the size of the combination equals the required size.
        // Add the combination and return
        if (combination.size() == ploidy) {
            combinations.add(0, combination); // Add in the beginning to maintain the correct order.
            return;
        }

        // Loop through the possible values from high to low.
        for (int newAlleleValue = maxAlleleValue; newAlleleValue >= 0; newAlleleValue--) {
            // Copy the preliminary combination
            List<Integer> newCombination = new ArrayList<>(combination);
            // Add a new value to the combination
            newCombination.add(0, newAlleleValue); // Add in the beginning to maintain the correct order.
            getCombinationsRecursively(combinations, newCombination, newAlleleValue, ploidy);
        }
    }

    /**
     * Convert 4 bytes to unsigned 32 bit int from index. Returns long since
     * java does not have unsigned int
     * <p>
     * https://stackoverflow.com/questions/13203426/convert-4-bytes-to-an-unsigned-32-bit-integer-and-storing-it-in-a-long
     *
     * @return
     * @throws EOFException
     * @throws IOException
     */
    private long getUInt32(byte[] bytes, int startIndex) {
        long value = bytes[0 + startIndex] & 0xFF;
        value |= (bytes[1 + startIndex] << 8) & 0xFFFF;
        value |= (bytes[2 + startIndex] << 16) & 0xFFFFFF;
        value |= (bytes[3 + startIndex] << 24) & 0xFFFFFFFF;
        return value;
    }

    /**
     * Convert 2 bytes to unsigned 16 bit int from start index.
     * <p>
     * https://stackoverflow.com/questions/13203426/convert-4-bytes-to-an-unsigned-32-bit-integer-and-storing-it-in-a-long
     *
     * @return
     * @throws EOFException
     * @throws IOException
     */
    private int getUInt16(byte[] bytes, int startIndex) {
        int value = bytes[0 + startIndex] & 0xFF;
        value |= (bytes[1 + startIndex] << 8) & 0xFFFF;
        return value;
    }

    /**
     * Convert 1 byte to unsigned 8 bit int from start index.
     * <p>
     * https://stackoverflow.com/questions/13203426/convert-4-bytes-to-an-unsigned-32-bit-integer-and-storing-it-in-a-long
     *
     * @return
     * @throws EOFException
     * @throws IOException
     */
    private int getUInt8(byte[] bytes, int startIndex) {
        int value = bytes[0 + startIndex] & 0xFF;
        return value;
    }

    /**
     * Converts a specified range of bits, up to a total of 32, from a byte array into a 'long' value.
     *
     * @param bytes           The byte array to get the result value from.
     * @param bitOffset       The bit in the byte array to start reading from.
     * @param totalBitsToRead The total number of bits to read from the byte array.
     * @return The specified bits converted to a 'long' value.
     */
    private static long readProbabilityValue(byte[] bytes, int bitOffset, int totalBitsToRead) {

        // Get the byte to start reading from.
        int firstByteIndex = bitOffset / 8;
        // Get the bit within the byte to start reading from.
        int remainingBitOffset = bitOffset % 8;

        int totalBytesMin1 = (remainingBitOffset + totalBitsToRead - 1) / 8;

        int nBitsFromLastByte = totalBitsToRead - (totalBytesMin1 * 8) + remainingBitOffset;

        int bitShiftAfterFirstByte = 8 - remainingBitOffset;

        long value; // long because values are stored as unsigned int

        //Switch because this is very fast compared to loops or if statements
        switch (totalBytesMin1) {
            case 0:
                value
                        = bytes[firstByteIndex] >> (remainingBitOffset) & LAST_BYTE_MASK[totalBitsToRead];
                break;
            case 1:
                //Last byte parsing is different then longer encoding
                value
                        = (bytes[firstByteIndex + 1] & LAST_BYTE_MASK[nBitsFromLastByte]) << bitShiftAfterFirstByte
                        | (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
                break;
            case 2:
                value
                        = (bytes[firstByteIndex + 2] & LAST_BYTE_MASK[nBitsFromLastByte]) << (8 + bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex + 1] & 255) << (bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
                break;
            case 3:
                value
                        = ((long) bytes[firstByteIndex + 3] & LAST_BYTE_MASK[nBitsFromLastByte]) << (16 + bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex + 2] & 255) << (8 + bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex + 1] & 255) << (bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
                break;
            case 4:
                value
                        = ((long) bytes[firstByteIndex + 4] & LAST_BYTE_MASK[nBitsFromLastByte]) << (24 + bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex + 3] & 255) << (16 + bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex + 2] & 255) << (8 + bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex + 1] & 255) << (bitShiftAfterFirstByte)
                        | (bytes[firstByteIndex] & FIRST_BYTE_MASK[bitShiftAfterFirstByte]) >> remainingBitOffset;
                break;

            default:
                throw new GenotypeDataException("Error parsing bgen file. Debug info: totalBits=" + totalBitsToRead + "" +
                        " remainingBitOffset=" + remainingBitOffset + " totalBits=" + totalBitsToRead);
        }
        return value;
    }

    /**
     * Represents information from a genotype data block within a BGEN file format.
     *
     * @author Robert Warmerdam
     */
    public static class VariantGenotypeBlockInfo {
        private long variantGenotypeStartPosition;
        private final Layout bgenFileLayout;
        private final long blockLength;
        private long decompressedBlockLength;
        private final boolean isCompressed;

        /**
         * The number of bytes that comprise the field representing total length C
         * of the remaining part of the variant block (probability data).
         */
        private static final short BLOCK_LENGTH_FIELD_SIZE = 4;
        /**
         * The number of bytes that comprise the field representing total length D
         * of the probability data after decompression.
         */
        private static final short DECOMPRESSED_BLOCK_LENGTH_FIELD_SIZE = 4;

        /**
         * Default complete constructor.
         *
         * @param bgenFileLayout               The layout of the BGEN file
         * @param variantGenotypeStartPosition The position of the genotype data block in the BGEN file.
         * @param blockLength                  The total length C of the rest of the data for the variant
         * @param decompressedBlockLength      The total length of the probability data after decompression.
         * @param isCompressed                 Flag indicating if the probability data is compressed
         */
        VariantGenotypeBlockInfo(
                long variantGenotypeStartPosition,
                Layout bgenFileLayout,
                long blockLength,
                long decompressedBlockLength,
                boolean isCompressed) {
            this(variantGenotypeStartPosition, bgenFileLayout, blockLength, isCompressed);
            this.decompressedBlockLength = decompressedBlockLength;
        }

        /**
         * Constructor that should be used for layout 1 as a compressed state in layout 2
         * requires a decompressed length D.
         *
         * @param variantGenotypeStartPosition The position of the genotype data block in the BGEN file.
         * @param bgenFileLayout               The layout of the BGEN file
         * @param blockLength                  The total length C of the rest of the data for the variant
         * @param isCompressed                 Flag indicating if the probability data is compressed
         */
        VariantGenotypeBlockInfo(
                long variantGenotypeStartPosition,
                Layout bgenFileLayout,
                long blockLength,
                boolean isCompressed) {

            this.variantGenotypeStartPosition = variantGenotypeStartPosition;
            this.bgenFileLayout = bgenFileLayout;
            this.blockLength = blockLength;
            this.isCompressed = isCompressed;
        }

        /**
         * Getter for the length of the probability data.
         * Is equal to C or C-4 when the length of the decompressed probability data is given
         *
         * @return length of the probability data in number of bytes
         */
        long getBlockLength() {
            long actualBlockLength = blockLength;
            // If the decompressed block length is given, subtract the length of the field representing
            // this value.
            if (isDecompressedBlockLengthGiven()) {
                actualBlockLength -= DECOMPRESSED_BLOCK_LENGTH_FIELD_SIZE;
            }
            return actualBlockLength;
        }

        private boolean isDecompressedBlockLengthGiven() {
            return bgenFileLayout == Layout.layOut_2 && isCompressed;
        }

        /**
         * Getter for the length of the probability data in addition to the lengths of field C
         * (representing the variants block length) if it is present.
         *
         * @return the total length of the genotype data block
         */
        private long getBlockLengthHeaderInclusive() {
            long blockLengthHeaderInclusive = blockLength;
            // Account for the 4 extra bytes that are always present
            // at the start of the genotype data block in layout 2 and present if the data in layout 1 is compressed
            if (isCompressedBlockLengthGiven()) {
                blockLengthHeaderInclusive += BLOCK_LENGTH_FIELD_SIZE;
            }
            return blockLengthHeaderInclusive;
        }

        private boolean isCompressedBlockLengthGiven() {
            return bgenFileLayout == Layout.layOut_2 ||
                    (bgenFileLayout == Layout.layOut_1 && isCompressed);
        }

        /**
         * Getter for the length of the decompressed probability data, D, in bytes.
         * Returns 0 if not set.
         *
         * @return the length of the decompressed probability data.
         */
        long getDecompressedBlockLength() {
            return decompressedBlockLength != 0 ? decompressedBlockLength : blockLength;
        }

        /**
         * Getter for the offset of the probability data within the genotype data block.
         *
         * @return the offset of the probability data within the genotype data block.
         */
        private int getBlockOffset() {
            int blockOffset = 0;
            // If the compressed block length (C) is given, take this field length into account.
            if (isCompressedBlockLengthGiven()) {
                blockOffset += BLOCK_LENGTH_FIELD_SIZE;
                // If additionally the decompressed block length (D) is given, take this field
                // into account as well.
                if (isDecompressedBlockLengthGiven()) {
                    blockOffset += DECOMPRESSED_BLOCK_LENGTH_FIELD_SIZE;
                }
            }
            return blockOffset;
        }

        long getVariantProbabilitiesStartPosition() {
            return variantGenotypeStartPosition + getBlockOffset();
        }

        long getVariantDataSizeInBytes(long variantReadingPosition) {
            return variantGenotypeStartPosition - variantReadingPosition + getBlockLengthHeaderInclusive();
        }
    }
}
