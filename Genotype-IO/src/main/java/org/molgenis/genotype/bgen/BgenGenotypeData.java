/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.bgen;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import com.facebook.presto.orc.zstd.ZstdDecompressor;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GeneticVariantMetaMap;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class BgenGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

	public enum blockRepresentation {
		compression_0, compression_1, compression_2
	}

	public enum Layout {
		layOut_1, layOut_2
	}
	private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL = 0.4f;
	private final RandomAccessFile bgenFile;
	private final byte[] byteArray4 = new byte[4]; //resuable 4 byte array
	private final byte[] byteArray2 = new byte[2]; //resuable 2 byte array
	private static final Logger LOGGER = Logger.getLogger(BgenGenotypeData.class);
	private final boolean sampleIdentifiersPresent;
	private final Inflater gzipInflater = new Inflater();
	private final ZstdDecompressor zstdInflater = new ZstdDecompressor();
	private static final Charset CHARSET = Charset.forName("UTF-8");
	private final blockRepresentation snpBlockRepresentation;
	private final Layout fileLayout;
	private static final GeneticVariantMeta GP_VARIANT = GeneticVariantMetaMap.getGeneticVariantMetaGp();
	private final SampleVariantsProvider sampleVariantProvider;
	private BgenixReader bgenixReader; // Was previously a final field
	private final long sampleCount;

	public BgenGenotypeData(File bgenFile, File sampleFile) throws IOException {
		this(bgenFile, sampleFile, 1000);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize) throws IOException {
		this(bgenFile, sampleFile, cacheSize, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws IOException {

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
        sampleCount = readFourBytesAsUInt32(
                "Number of samples",
                "Error reading bgen file header. File is corrupt");

        // Throw an exception if the number of samples exceeds the maximum value of an integer.
		if (sampleCount > Integer.MAX_VALUE) {
			throw new GenotypeDataException("Found more than (2^31)-1 samples in bgen file. This is not supported");
		}

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
        fileLayout = readFileLayout();

        // Throw an exception if the SNP block representation is 2 while layout 1 is used.
        if (snpBlockRepresentation.equals(blockRepresentation.compression_2) & fileLayout.equals(Layout.layOut_1)) {
			throw new GenotypeDataException("Invalid compression method for layout one observed. Trying to use ZSTD compression on layout one file, which is not supported.");
		}

        // Read the sample identifiers presence; set sampleIdentifiersPresent to true if present, false if absent.
        sampleIdentifiersPresent = readSampleIdentifiersPresence();

        if (sampleIdentifiersPresent) {
			// Process sample identifier block.
            processSampleIdentifierBlock(snpOffset, headerSize);
        }

//		System.out.println("offset:"+ (snpOffset + 4));
		File bgenixFile = new File(bgenFile.getAbsolutePath() + ".bgi");
		System.out.println(bgenFile.getAbsolutePath() + ".bgi");

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

        readExistingBgenixFile(bgenFile, bgenixFile);

        // If the specified cache size is greater than 0, construct a new SampleVariantProvider
		if (cacheSize > 0) {
			sampleVariantProvider = new CachedSampleVariantProvider(this, cacheSize);
		} else {
			sampleVariantProvider = this;
		}
		//Read the first snp to get into genotype-io.
//		readCompleteGeneticVariant(bgenFile, pointerFirstSnp, (int) sampleCount, this.fileLayout, this.snpBlockRepresentation);
	}

    private void readExistingBgenixFile(File bgenFile, File bgenixFile) throws IOException {
	    // Not sure what bgenix reader should do in this object other than read an existing BGENIX file.
        bgenixReader = new BgenixReader(bgenixFile);
        BgenixMetadata metadata = bgenixReader.getMetadata();
        if (metadata != null) {
            // Check if the metadata of the read BGENIX file is equal to
            // the metadata of the BGEN file.
            if (!metadata.getFileName().equals(bgenFile.getName())) {
                throw new GenotypeDataException("Sample name between bgenix and bgen is not equal. Invalid Bgen and Bgenix combination.");
            }

            if (metadata.getFileSize() != bgenFile.length()) {
                throw new GenotypeDataException("Number of expected bytes is different to the numer of observed bytes. Invalid Bgen and Bgenix combination.");
            }

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
            throw new GenotypeDataException("Numer of samples in metadata and sample id data is not equal. File is corrupt.");
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
//				System.out.println("\t"+sampleIds[i]);
        }
    }

    private boolean readSampleIdentifiersPresence() {
        if ((byteArray4[3] & 128) == 128) {
            LOGGER.debug("SampleIdentifiers present in bgen file");
            return true;
        } else {
            LOGGER.debug("SampleIdentifiers not present in bgen file");
            return false;
        }
    }

    private Layout readFileLayout() {
        byte byte_tmp = byteArray4[0];
        byte_tmp = (byte) (byte_tmp >> 2);
        switch (byte_tmp & 7) {
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

    private blockRepresentation readSnpBlockRepresentation() {
        switch (byteArray4[0] & 3) {
            case 0:
                LOGGER.debug("Genotype data is not compressed.");
                return blockRepresentation.compression_0;
            case 1:
                LOGGER.debug("Genotype data is zlib compressed");
                return blockRepresentation.compression_1;
            case 2:
                LOGGER.debug("Genotype data is zstd compressed");
                return blockRepresentation.compression_2;
            default:
                throw new GenotypeDataException("Invalid compression method, observed: " + (byteArray4[0] & 3));
        }
    }

    private long readFourBytesAsUInt32(String fieldName, String exceptionMessage) throws IOException {
        if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
            throw new GenotypeDataException(exceptionMessage);
        }
        long value = getUInt32(byteArray4, 0);
        LOGGER.debug(fieldName + ": " + value);
        return value;
    }

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
//		int stepToNextVariant = 0;
//		if (fileLayout.equals(layout.layOut_1)) {
//			stepToNextVariant += 4;
//		}
		long startSize = pointerFirstSnp;
		while ((startSize) < bgenFile.length()) {
			//Loop through variants.
//			long currentStart = pointerFirstSnp + stepToNextVariant;

			GeneticVariant var = readVariantIdentifyingData(startSize);
			long currentPointer = this.bgenFile.getFilePointer();
			long variantGenotypeDataBlockSize = determineVariantGenotypeDataBlockSize(
					currentPointer);

			// If layout 2, the block size C may include 4 bytes representing D
			// I suppose the following tries to correct for this??
			// No, should be substracting 4 instead of adding...
			if (fileLayout.equals(Layout.layOut_2)) {
				variantGenotypeDataBlockSize += 4;
			}
//			this.bgenFile.seek(currentPointer);
			startSize = currentPointer + variantGenotypeDataBlockSize;

			bgenixWriter.addVariantToIndex(
					var,
					pointerFirstSnp,
					(int) variantGenotypeDataBlockSize,
					var.getVariantId().getPrimairyId());
		}
	}

	private GeneticVariant readVariantIdentifyingData(long filePointer) throws IOException {
		long lastSnpStart = filePointer;
		// If layout is equal to 1 then the variant identifying data starts with 4 bytes describing the
		// number of individuals within the row.
		if (fileLayout == Layout.layOut_1) {
			// We chose to ignore this...
			lastSnpStart += 4;
		}
//		System.out.println("Seeking to:" +lastSnpStart);
		this.bgenFile.seek(lastSnpStart);
		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

//		if (snpInfoBufferSize < 20) {
//			throw new GenotypeDataException("Error reading bgen snp data. File is corrupt");
//		}
		// Need to check that it is correct with the block infront of the snp id.
		ArrayList<String> snpIds = new ArrayList<>();
		int fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp ID length " + fieldLength);
		snpInfoBufferPos += 2;
		String snpId = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
		LOGGER.debug("SNP ID: " + snpId);
//		System.out.println(snpId);
		snpInfoBufferPos += fieldLength; // skip id length and snp id

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp RS length " + fieldLength);

		snpInfoBufferPos += 2;
		String snpRsId = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
		LOGGER.debug("SNP RSID: " + snpRsId);

//		System.out.println(snpRsId);
		snpInfoBufferPos += fieldLength;
		snpIds.add(snpRsId);
		snpIds.add(snpId);

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 2;
		String seqName = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, CHARSET);
//		System.out.println(seqName);
		snpInfoBufferPos += fieldLength;

		long snpPosLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 4;
		if (snpPosLong > Integer.MAX_VALUE) {
			throw new GenotypeDataException("SNP pos larger than (2^31)-1 not supported");
		}

		int snpPos = (int) snpPosLong;

//		System.out.println("SNP pos " + snpPos);
		int numberOfAlleles = 2;
		if (fileLayout.equals(Layout.layOut_2)) {
			numberOfAlleles = getUInt16(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 2;
//			System.out.println("SNP Alleles " + numberOfAlleles);
		}

		ArrayList<String> alleles = new ArrayList<>();
		for (int i = 0; i < numberOfAlleles; i++) {
			snpInfoBufferPos = readAllele(snpInfoBuffer, snpInfoBufferPos, alleles);
		}

//		System.out.println(alleles.toString());
//			System.out.println("Location where genotypes start: "+this.bgenFile.getFilePointer());
		this.bgenFile.seek(snpInfoBufferPos + lastSnpStart);

		GeneticVariant var = ReadOnlyGeneticVariant.createVariant(
				GP_VARIANT, snpIds, snpPos, seqName, sampleVariantProvider, alleles, alleles.get(0));

		return var;
	}

	private int readAllele(byte[] snpInfoBuffer, int snpInfoBufferPos, ArrayList<String> alleles) {

		// Length of the allele
		long fieldLengthLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 4;

		// Throw an exception if the allele length is longer
		if (fieldLengthLong > Integer.MAX_VALUE) {
			throw new GenotypeDataException("SNP with allele longer than (2^31)-1 characters not supported");
		}

		// Get the allele from the buffer.
		String allele = new String(snpInfoBuffer, snpInfoBufferPos, (int) fieldLengthLong, CHARSET);
		snpInfoBufferPos += ((int) fieldLengthLong);
		alleles.add(allele);
		return snpInfoBufferPos;
	}

	private void readGenotypesFromVariant(long filePointer) throws IOException {
		this.bgenFile.seek(filePointer);

		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		// Read
		if (fileLayout == Layout.layOut_1) {
			long snpBlockSize = determineVariantGenotypeDataBlockSizeForLayoutOne(
					snpInfoBuffer,
					snpInfoBufferPos);

			byte[] snpBlockData = new byte[6 * (int) sampleCount];
			if (snpBlockRepresentation.equals(blockRepresentation.compression_1)) {

				snpInfoBufferPos += 4; // Shift reading frame 4 bytes to just after the variant size field

				decompressVariantBlockGzip(snpInfoBuffer,
						(int) snpBlockSize - 4,
						snpBlockData,
						snpInfoBufferPos);

				//Here we can parse from the inflater block for all the samples.
				System.out.println(getUInt16(snpBlockData, 0) / 32768f + " " + getUInt16(snpBlockData, 2) / 32768f + " " + getUInt16(snpBlockData, 4) / 32768f);
			} else {
				//Here we can directly parse from the snpInfoBuffer.
			}

		} else if (fileLayout.equals(Layout.layOut_2)) {
			long snpBlockSize = determineVariantGenotypeDataBlockSizeForLayoutTwo(snpInfoBuffer, snpInfoBufferPos);
			long snpBlockSizeDecompressed = determineDecompressedVariantGenotypeDataBlockSizeForLayoutTwo(
					snpInfoBuffer,
					snpInfoBufferPos + 4, // Shift reading frame 4 bytes to decompressed size field.
					snpBlockSize);
			System.out.println("Snp block size: " + snpBlockSize);
			System.out.println("Snp block size decompressed: " + snpBlockSizeDecompressed);

			snpInfoBufferPos += 4; // Shift reading frame 4 bytes to just after the decompressed size field.

			byte[] snpBlockData = new byte[(int) snpBlockSizeDecompressed];
			switch (snpBlockRepresentation) {

				case compression_1:
					decompressVariantBlockGzip(snpInfoBuffer, (int) snpBlockSize - 4, snpBlockData, snpInfoBufferPos);

					//At genotype / haplotype data
				case compression_2:
					zstdInflater.decompress(snpInfoBuffer, snpInfoBufferPos,(int) snpBlockSize - 4,snpBlockData, 0, (int) snpBlockSizeDecompressed);
				//Is this enough?

				default:
					break;
//					snpInfoBufferPos = ;
				//Not compressed.
			}
			int blockBuffer = 0;
			//must equal data before.
			int numberOfIndividuals = (int) getUInt32(snpBlockData, blockBuffer);
			System.out.println("Number of individuals: " + numberOfIndividuals);
			blockBuffer += 4;
			//must equal data before.
			int numberOfAlleles = (int) getUInt16(snpBlockData, blockBuffer);
			System.out.println("Number of Alleles: " + numberOfAlleles);
			blockBuffer += 2;

			int minPloidy = getUInt8(snpBlockData, blockBuffer);
			System.out.println("Min ploidy: " + minPloidy);
			blockBuffer += 1;
			int maxPloidy = getUInt8(snpBlockData, blockBuffer);
			System.out.println("Max ploidy: " + maxPloidy);
			blockBuffer += 1;
			for (int i = 0; i < numberOfIndividuals; i++) {
				//Here we need to handle missing ploidity.
				//Missingness is encoded by the most significant bit; thus a value of 1 for the most significant bit indicates that no probability data is stored for this sample.
//				System.out.println("ploidity: " + getUInt8(snpBlockData, blockBuffer));
				blockBuffer += 1;
			}
			int phased = getUInt8(snpBlockData, blockBuffer);
			System.out.println("phased: " + phased);
			blockBuffer += 1;
			if (phased > 1) {
				throw new GenotypeDataException("Bgen file format error. Unsupported value for the phased flag observed.");
			}
			int bitProbabilityRepresentation = getUInt8(snpBlockData, blockBuffer);
			System.out.println("Bit representation of probability: " + bitProbabilityRepresentation);
			blockBuffer += 1;

			if (phased == 1) {

			} else {

			}

		}
	}

	private void decompressVariantBlockGzip(byte[] snpInfoBuffer, int snpBlockSize, byte[] snpBlockData, int i2) {
		gzipInflater.setInput(snpInfoBuffer, i2, snpBlockSize);
		try {
			gzipInflater.inflate(snpBlockData);
		} catch (DataFormatException ex) {
			throw new GenotypeDataException("Error decompressing bgen data", ex);
		}
		gzipInflater.reset();
	}

//	private void readCompleteGeneticVariant(File bgen, long lastSnpStart, int sampleCount, layout currentFileLayout, blockRepresentation currentBlockRepresentation) throws IOException {
//		//Binary index writer.
////		for (int snpI = 0; snpI < variantCount; ++snpI) {
//
//		long currentPointer= this.bgenFile.getFilePointer();
//		long stepSize=0L;
//		
//		while((currentPointer+stepSize)<bgen.length()){
//			GeneticVariant var = readSnpInfo((currentPointer+stepSize));
//			System.out.println(var.getAllIds().get(0));
//			currentPointer= this.bgenFile.getFilePointer();
//			
//			stepSize = DetermineStepSize(currentPointer, currentFileLayout, currentBlockRepresentation, sampleCount);
//			if(currentFileLayout.equals(layout.layOut_2)){
//				stepSize+=4;
//			}
//			this.bgenFile.seek(currentPointer);
//			//Here we can process actual genotype info.
//			
//			
////			System.out.println(bgen.length() - (currentPointer+stepSize));
//			System.out.println("Step to skip, to next variant: " + (int) stepSize);
//		}
//
//	}

	private Long determineVariantGenotypeDataBlockSize(
			long filePointer) throws IOException {

		// Seek towards the file pointer (bgenFile.getFilePointer();)
		this.bgenFile.seek(filePointer);

		// Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		Long snpBlockSize = null;

		if (fileLayout == Layout.layOut_1) {
			snpBlockSize = determineVariantGenotypeDataBlockSizeForLayoutOne(snpInfoBuffer, snpInfoBufferPos);
		} else if (fileLayout == Layout.layOut_2) {
			snpBlockSize = determineVariantGenotypeDataBlockSizeForLayoutTwo(snpInfoBuffer, snpInfoBufferPos);
		}
		return snpBlockSize;
	}

	private Long determineVariantGenotypeDataBlockSizeForLayoutTwo(byte[] snpInfoBuffer, int snpInfoBufferPos) {
		Long snpBlockSize;
//		long snpBlockSizeDecompressed;
		// If the file layout is 2,
		// The genotype data for this variant is equal to the next 4 bytes
		snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		// Shift the buffer position 4 bytes
		snpInfoBufferPos += 4;

//		if (snpBlockRepresentation.equals(blockRepresentation.compression_0)) {
//			snpBlockSizeDecompressed = snpBlockSize;
//		} else {
//			snpBlockSizeDecompressed = getUInt32(snpInfoBuffer, snpInfoBufferPos);
//			snpInfoBufferPos += 4; // Has no function
//		}
//		System.out.println("Snp block size: " + snpBlockSize);
//		System.out.println("Snp block size decompressed: " + snpBlockSizeDecompressed);
		return snpBlockSize;
	}

	private Long determineDecompressedVariantGenotypeDataBlockSizeForLayoutTwo(byte[] snpInfoBuffer,
																			   int snpInfoBufferPos,
																			   long variantBlockSize) {
		long snpBlockSizeDecompressed;
		if (snpBlockRepresentation.equals(blockRepresentation.compression_0)) {
			snpBlockSizeDecompressed = variantBlockSize;
		} else {
			snpBlockSizeDecompressed = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 4; // Has no function
		}
		return snpBlockSizeDecompressed;
	}

	private Long determineVariantGenotypeDataBlockSizeForLayoutOne(byte[] snpInfoBuffer, int snpInfoBufferPos) {
		long snpBlockSize;
		if (snpBlockRepresentation.equals(blockRepresentation.compression_1)) {
			// If the file layout is 1 and the variants are zlib compressed
			// get the snp block size by within the next four bytes
			snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			snpBlockSize += 4;
		} else {
			// If the file layout is 1 and the variants are not zlib compressed,
			// the genotype data for this variant is 6 times the number of samples
			snpBlockSize = 6 * sampleCount;
		}
		System.out.println("Snp block size: " + snpBlockSize);
		return snpBlockSize;
	}

	@Override
	public List<Sample> getSamples() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Map<String, Annotation> getVariantAnnotationsMap() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public boolean isOnlyContaingSaveProbabilityGenotypes() {
		return true;
	}

	@Override
	public void close() throws IOException {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public List<String> getSeqNames() {
	    // Corresponds to
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Iterable<Sequence> getSequences() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByPos(String seqName, int startPos) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
		BgenixVariantQueryResult bgenixVariants = bgenixReader.getVariantsChromosome(seqName);

		return null;
	}

	@Override
	public Iterable<GeneticVariant> getVariantsByRange(String seqName, int rangeStart, int rangeEnd) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public List<Alleles> getSampleVariants(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public List<Boolean> getSamplePhasing(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public int cacheSize() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public int getSampleVariantProviderUniqueId() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public byte[] getSampleCalledDosage(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public float[] getSampleDosage(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	@Override
	public float[][] getSampleProbilities(GeneticVariant variant) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	/**
	 * Convert 4 bytes to unsigned 32 bit int from index. Returns long since
	 * java does not have unsigned int
	 *
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
	 *
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
	 *
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
}
