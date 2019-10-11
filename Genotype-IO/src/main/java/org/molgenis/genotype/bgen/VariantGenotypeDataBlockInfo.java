package org.molgenis.genotype.bgen;

/**
 * Represents information from a genotype data block within a BGEN file format.
 * @author Robert Warmerdam
 */
class VariantGenotypeDataBlockInfo {
    private final BgenGenotypeData.Layout bgenFileLayout;
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
     * @param bgenFileLayout The layout of the BGEN file
     * @param blockLength The total length C of the rest of the data for the variant
     * @param decompressedBlockLength The total length of the probability data after decompression.
     * @param isCompressed Flag indicating if the probability data is compressed
     */
    VariantGenotypeDataBlockInfo(
            BgenGenotypeData.Layout bgenFileLayout,
            long blockLength,
            long decompressedBlockLength,
            boolean isCompressed) {
        this(bgenFileLayout, blockLength, isCompressed);
        this.decompressedBlockLength = decompressedBlockLength;
    }

    /**
     * Constructor that should be used for layout 1 as a compressed state in layout 2
     * requires a decompressed length D.
     *
     * @param bgenFileLayout The layout of the BGEN file
     * @param blockLength The total length C of the rest of the data for the variant
     * @param isCompressed Flag indicating if the probability data is compressed
     */
    VariantGenotypeDataBlockInfo(
            BgenGenotypeData.Layout bgenFileLayout,
            long blockLength,
            boolean isCompressed) {

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
        return bgenFileLayout == BgenGenotypeData.Layout.layOut_2 && isCompressed;
    }

    /**
     * Getter for the length of the probability data in addition to the lengths of field C
     * (representing the variants block length) is it is present.
     *
     * @return the total length of the genotype data block
     */
    long getBlockLengthHeaderInclusive() {
        long blockLengthHeaderInclusive = blockLength;
        // Account for the 4 extra bytes that are always present
        // at the start of the genotype data block in layout 2 and present if the data in layout 1 is compressed
        if (isCompressedBlockLengthGiven()) {
            blockLengthHeaderInclusive += BLOCK_LENGTH_FIELD_SIZE;
        }
        return blockLengthHeaderInclusive;
    }

    private boolean isCompressedBlockLengthGiven() {
        return bgenFileLayout == BgenGenotypeData.Layout.layOut_2 ||
                (bgenFileLayout == BgenGenotypeData.Layout.layOut_1 && isCompressed);
    }

    /**
     * Getter for the length of the decompressed probability data, D, in bytes.
     * Returns 0 if not set.
     *
     * @return the length of the decompressed probability data.
     */
    long getDecompressedBlockLength() {
        return decompressedBlockLength;
    }

    /**
     * Getter for the offset of the probability data within the genotype data block.
     *
     * @return the offset of the probability data within the genotype data block.
     */
    int getBlockOffset() {
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
}
