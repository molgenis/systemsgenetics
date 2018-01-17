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
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.id.GeneticVariantId;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;

/**
 *
 * @author Patrick Deelen
 */
public class BgenGenotypeData {

	public enum blockRepresentation {
		compression_0, compression_1, compression_2
	}

	public enum layout {
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
	private static final Charset charset = Charset.forName("UTF-8");
	private final blockRepresentation snpBlockRepresentation;
	private final layout fileLayout;

	public BgenGenotypeData(File bgenFile, File sampleFile) throws IOException {
		this(bgenFile, sampleFile, 1000);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize) throws IOException {
		this(bgenFile, sampleFile, cacheSize, DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL);
	}

	public BgenGenotypeData(File bgenFile, File sampleFile, int cacheSize, double minimumPosteriorProbabilityToCall)
			throws IOException {

		this.bgenFile = new RandomAccessFile(bgenFile, "r");

		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long snpOffset = getUInt32(byteArray4, 0);
		LOGGER.debug("SNP offset: " + snpOffset);

		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long headerSize = getUInt32(byteArray4, 0);
		LOGGER.debug("Header size (Lh): " + headerSize);
		if (headerSize > snpOffset) {
			throw new GenotypeDataException("Error reading bgen file header. Header information is bigger than expected offset.");
		}
		//add checks for header size and header size must be smaller than offset.
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long variantCount = getUInt32(byteArray4, 0);
		LOGGER.debug("Number of SNPs: " + variantCount);

//               Not needed when using the bgenix indexing file.  
//		if (snpCount > Integer.MAX_VALUE) {
//			throw new GenotypeDataException("Found more than (2^31)-1 SNPs in bgen file. This is not supported");
//		}
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long sampleCount = getUInt32(byteArray4, 0);
		LOGGER.debug("Number of samples: " + sampleCount);

		if (sampleCount > Integer.MAX_VALUE) {
			throw new GenotypeDataException("Found more than (2^31)-1 samples in bgen file. This is not supported");
		}
		//Magic number is in the next four bytes but skipped over.

		//Skip over reserved and free area, 
		//this seek works because there are 4 (interesting) bytes 
		//in the beginning of the file and in the end of the header block.
		this.bgenFile.seek(headerSize);

		//Read flags
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}

		//Changed storing of version and data.
		switch (byteArray4[0] & 3) {
			case 0:
				LOGGER.debug("Genotype data is not compressed.");
				snpBlockRepresentation = blockRepresentation.compression_0;
				break;
			case 1:
				LOGGER.debug("Genotype data is zlib compressed");
				snpBlockRepresentation = blockRepresentation.compression_1;
				break;
			case 2:
				LOGGER.debug("Genotype data is zstd compressed");
				snpBlockRepresentation = blockRepresentation.compression_2;
				break;
			default:
				throw new GenotypeDataException("Invalid compression method, observed: " + (byteArray4[0] & 3));
		}
		byte byte_tmp = byteArray4[0];
		byte_tmp = (byte) (byte_tmp >> 2);
		switch (byte_tmp & 7) {
			case 1:
				LOGGER.debug("Genotype data is in layout 1 (BGEN 1.1)");
				fileLayout = layout.layOut_1;
				break;
			case 2:
				LOGGER.debug("Genotype data is in layout 2 (BGEN 1.2 & 1.3)");
				fileLayout = layout.layOut_2;
				break;
			default:
				throw new GenotypeDataException("Invalid layout, observed: " + (byteArray4[0] & 3));
		}

		if (snpBlockRepresentation.equals(blockRepresentation.compression_2) & fileLayout.equals(layout.layOut_1)) {
			throw new GenotypeDataException("Invalid compression method for layout one observed. Trying to use ZSTD compression on layout one file, which is not supported.");
		}

		if ((byteArray4[3] & 128) == 128) {
			LOGGER.debug("SampleIdentifiers present in bgen file");
			sampleIdentifiersPresent = true;
		} else {
			LOGGER.debug("SampleIdentifiers not present in bgen file");
			sampleIdentifiersPresent = false;
		}

		if (sampleIdentifiersPresent) {
			// Proccess sample indentifier block
			if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
				throw new GenotypeDataException("Error in sample identifier block. File is corrupt.");
			}
			long byteSizeSampleIds = getUInt32(byteArray4, 0);
			LOGGER.debug("LSI size: " + byteSizeSampleIds);
			if ((byteSizeSampleIds + headerSize) > snpOffset) {
				throw new GenotypeDataException("Error reading bgen file header. Combination of header & sample id information is bigger than expected offset.");
			}

			// Proccess sample indentifier block
			if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
				throw new GenotypeDataException("Error in sample identifier block. File is corrupt.");
			}
			long sampleCount2 = getUInt32(byteArray4, 0);
			LOGGER.debug("Number of samples with sample-ids: " + sampleCount2);
			if (sampleCount2 != sampleCount) {
				throw new GenotypeDataException("Numer of samples in metadata and sample id data is not equal. File is corrupt.");
			}
			String[] sampleIds = new String[(int) sampleCount2];
			for (int i = 0; i < sampleCount2; i++) {
				if (this.bgenFile.read(byteArray2, 0, 2) != 2) {
					throw new GenotypeDataException("Error in sample Id. File is corrupt.");
				}
				int sampleIdLength = getUInt16(byteArray2, 0);
				byte[] sampleName = new byte[sampleIdLength];
				this.bgenFile.read(sampleName, 0, sampleIdLength);
				sampleIds[i] = new String(sampleName, charset);
//				System.out.println("\t"+sampleIds[i]);
			}
		}
//		System.out.println("offset:"+ (snpOffset + 4));
		//Check if bgenix file is present.
		File bgenixFile = new File(bgenFile.getAbsolutePath() + ".bgi");
		System.out.println(bgenFile.getAbsolutePath() + ".bgi");

		long lastSnpStart = snpOffset + 4;

		if (bgenixFile.exists()) {
			BgenixReader indexInformation = new BgenixReader(bgenixFile);
			BgenixMetadata metadata = indexInformation.getMetadata();
			if (metadata != null) {
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
			}
		} else {
			BgenixWriter b = new BgenixWriter(bgenixFile);
			createBgenixFile(bgenFile, b, lastSnpStart, (int) sampleCount, this.fileLayout);
//			throw new GenotypeDataException("Currently only bgen genotype data indexed using bgenix is supported.");
		}

		//Read the first snp to get into genotype-io.
		readGeneticVariant(lastSnpStart, (int) sampleCount, this.fileLayout, this.snpBlockRepresentation);
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

	private void readGeneticVariant(long lastSnpStart, int sampleCount, layout currentFileLayout, blockRepresentation currentBlockRepresentation) throws IOException {
		//Binary index writer.
//		for (int snpI = 0; snpI < variantCount; ++snpI) {
		if (currentFileLayout.equals(layout.layOut_1)) {
			lastSnpStart += 4;
		}

		this.bgenFile.seek(lastSnpStart);
		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

//		if (snpInfoBufferSize < 20) {
//			throw new GenotypeDataException("Error reading bgen snp data. File is corrupt");
//		}
		int fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp ID length " + fieldLength);

		snpInfoBufferPos += 2 + fieldLength; // skip id length and snp id

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		LOGGER.debug("Snp RS length " + fieldLength);

		snpInfoBufferPos += 2;
		String snpId = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, charset);
		System.out.println(snpId);

		snpInfoBufferPos += fieldLength;

		fieldLength = getUInt16(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 2;
		String seqName = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, charset);
		System.out.println(seqName);
		snpInfoBufferPos += fieldLength;

		long snpPosLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
		snpInfoBufferPos += 4;
		if (snpPosLong > Integer.MAX_VALUE) {
			throw new GenotypeDataException("SNP pos larger than (2^31)-1 not supported");
		}

		int snpPos = (int) snpPosLong;

		System.out.println("SNP pos " + snpPos);

		int numberOfAlleles = 2;
		if (currentFileLayout.equals(layout.layOut_2)) {
			numberOfAlleles = getUInt16(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 2;
			System.out.println("SNP Alleles " + numberOfAlleles);
		}

		ArrayList<String> alleles = new ArrayList<>();
		for (int i = 0; i < numberOfAlleles; i++) {
			long fieldLengthLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 4;
			if (fieldLengthLong > Integer.MAX_VALUE) {
				throw new GenotypeDataException("SNP with allele longer than (2^31)-1 characters not supported");
			}
			String a = new String(snpInfoBuffer, snpInfoBufferPos, (int) fieldLengthLong, charset);
			snpInfoBufferPos += ((int) fieldLengthLong);
			alleles.add(a);
		}

		StringBuilder alleleBuffer = new StringBuilder();
		for (int i = 0; i < alleles.size(); i++) {
			alleleBuffer.append(alleles.get(i));
			alleleBuffer.append('/');
		}
		System.out.println(alleleBuffer.toString());

//			System.out.println("Location where genotypes start: "+this.bgenFile.getFilePointer());
		this.bgenFile.seek(snpInfoBufferPos + lastSnpStart);
		System.out.println("Location where genotypes start: " + this.bgenFile.getFilePointer());
		//readGenotypesFromVariant(this.bgenFile.getFilePointer(), currentFileLayout, currentBlockRepresentation, sampleCount);
		//move to next SNP
		//Determine skipping based on block and type.
		lastSnpStart = lastSnpStart + snpInfoBufferPos + snpBlockSize;
		this.bgenFile.seek(lastSnpStart);
//
//		}

	}

	private void createBgenixFile(File bgen, BgenixWriter b, long pointerFirstSnp, int nSamples, layout fileLayout) throws IOException {
		this.bgenFile.seek(0);

		byte[] firstBytes = new byte[1000];
		this.bgenFile.read(firstBytes, 0, 1000);

		//Add current time in int.
		BgenixMetadata m = new BgenixMetadata(bgen.getName(), (int) this.bgenFile.length(), (int) bgen.lastModified(), firstBytes, 100000);
		b.writeMetadata(m);

		//Loop through the start of the file
		int stepToNextVariant = 0;
		if (fileLayout.equals(layout.layOut_1)) {
			stepToNextVariant += 4;
		}
		while ((pointerFirstSnp + stepToNextVariant) < bgenFile.length()) {
			//Loop through variants.
			long currentStart = pointerFirstSnp + stepToNextVariant;
			ReadOnlyGeneticVariant var = ReadOnlyGeneticVariant.createSnp(variantMeta, snpIds, nSamples, sequenceName, sampleVariantsProvider, 0, 0);

			String variantId = null;

			b.addVariantToIndex(var, pointerFirstSnp, stepToNextVariant, variantId);
		}

		b.finalizeIndex();
	}

	private void readGenotypesFromVariant(long filePointer, layout currentFileLayout, blockRepresentation currentBlockRepresentation, int sampleCount) throws IOException {
		this.bgenFile.seek(filePointer);
		//Not sure if we want to do the buffer search here. Or we might be able to take a smaller set.
		byte[] snpInfoBuffer = new byte[8096];
//		int snpInfoBufferSize = 
		this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
		int snpInfoBufferPos = 0;

		if (currentFileLayout.equals(fileLayout.layOut_1)) {
			long snpBlockSize;
			if (currentBlockRepresentation.equals(blockRepresentation.compression_1)) {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpBlockSize += 4;
			} else {
				snpBlockSize = 6 * sampleCount;
			}
			System.out.println("Snp block size: " + snpBlockSize);
			byte[] snpBlockData = new byte[6 * (int) sampleCount];
			if (currentBlockRepresentation.equals(blockRepresentation.compression_1)) {
				gzipInflater.setInput(snpInfoBuffer, snpInfoBufferPos + 4, (int) snpBlockSize - 4);
				try {
					gzipInflater.inflate(snpBlockData);
				} catch (DataFormatException ex) {
					throw new GenotypeDataException("Error decompressing bgen data", ex);
				}
				gzipInflater.reset();

				//Here we can parse from the inflater block for all the samples.
				System.out.println(getUInt16(snpBlockData, 0) / 32768f + " " + getUInt16(snpBlockData, 2) / 32768f + " " + getUInt16(snpBlockData, 4) / 32768f);
			} else {
				//Here we can directly parse from the snpInfoBuffer.
			}

		} else if (currentFileLayout.equals(fileLayout.layOut_2)) {
			long snpBlockSize;
			long snpBlockSizeDecompressed;
			if (!currentBlockRepresentation.equals(blockRepresentation.compression_0)) {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
				snpBlockSizeDecompressed = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
			} else {
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos);
				snpInfoBufferPos += 4;
				snpBlockSizeDecompressed = snpBlockSize;
			}
			System.out.println("Snp block size: " + snpBlockSize);
			System.out.println("Snp block size decompressed: " + snpBlockSizeDecompressed);

			byte[] snpBlockData = new byte[(int) snpBlockSizeDecompressed];
			switch (currentBlockRepresentation) {

				case compression_1:
					gzipInflater.setInput(snpInfoBuffer, snpInfoBufferPos, (int) snpBlockSize - 4);
					try {
						gzipInflater.inflate(snpBlockData);
					} catch (DataFormatException ex) {
						throw new GenotypeDataException("Error decompressing bgen data", ex);
					}
					gzipInflater.reset();

				//At genotype / haplotype data
				case compression_2:
//					zstdInflater.decompress(snpInfoBuffer, snpInfoBufferPos,(int) snpBlockSize - 4,snpBlockData, 0, (int) snpBlockSizeDecompressed);
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
				System.out.println("ploidity: " + getUInt8(snpBlockData, blockBuffer));
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

}
