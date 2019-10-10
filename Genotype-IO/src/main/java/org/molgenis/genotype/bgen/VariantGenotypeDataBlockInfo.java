package org.molgenis.genotype.bgen;

class VariantGenotypeDataBlockInfo {
    private final BgenGenotypeData.Layout bgenFileLayout;
    private final long blockLength;
    private long decompressedBlockLength;
    private final boolean isCompressed;

    private static final short BLOCK_LENGTH_FIELD_SIZE = 4;
    private static final short DECOMPRESSED_BLOCK_LENGTH_FIELD_SIZE = 4;

    private VariantGenotypeDataBlockInfo(
            BgenGenotypeData.Layout bgenFileLayout,
            long blockLength,
            long decompressedBlockLength,
            boolean isCompressed) {
        this(bgenFileLayout, blockLength, isCompressed);
        this.decompressedBlockLength = decompressedBlockLength;
    }

    VariantGenotypeDataBlockInfo(
            BgenGenotypeData.Layout bgenFileLayout,
            long blockLength,
            boolean isCompressed) {

        this.bgenFileLayout = bgenFileLayout;
        this.blockLength = blockLength;
        this.isCompressed = isCompressed;
    }

    VariantGenotypeDataBlockInfo(
            BgenGenotypeData.Layout fileLayout,
            Long blockLength,
            long decompressedBlockLength) {
        this(fileLayout, blockLength, decompressedBlockLength,true);
    }

    long getBlockLength() {
        long actualBlockLength = blockLength;
        if (bgenFileLayout == BgenGenotypeData.Layout.layOut_2 && isCompressed) {
            actualBlockLength -= DECOMPRESSED_BLOCK_LENGTH_FIELD_SIZE;
        }
        return actualBlockLength;
    }

    long getBlockLengthHeaderInclusive() {
        long blockLengthHeaderInclusive = blockLength;
        // Account for the 4 extra bytes that are always present
        // at the start of the genotype data block in layout 2 and present if the data in layout 1 is compressed
        if (bgenFileLayout == BgenGenotypeData.Layout.layOut_2 ||
                (bgenFileLayout == BgenGenotypeData.Layout.layOut_1 && isCompressed)) {
            blockLengthHeaderInclusive += BLOCK_LENGTH_FIELD_SIZE;
        }
        return blockLengthHeaderInclusive;
    }

    long getDecompressedBlockLength() {
        return decompressedBlockLength;
    }
}
