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
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.variant.range.GeneticVariantRange;

/**
 *
 * @author Patrick Deelen
 */
public class BgenGenotypeData {

	private static final double DEFAULT_MINIMUM_POSTERIOR_PROBABILITY_TO_CALL = 0.4f;
	private final RandomAccessFile bgenFile;
	private final byte[] byteArray4 = new byte[4]; //resuable 4 byte array
	private static final Logger LOGGER = Logger.getLogger(BgenGenotypeData.class);
	private final boolean compressedGenotypes;
	private final boolean v1_1;
	private final Inflater inflater = new Inflater();
	private static final Charset charset = Charset.forName("UTF-8") ;

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
		LOGGER.debug("Header size: " + headerSize);

		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long snpCount = getUInt32(byteArray4, 0);
		LOGGER.debug("Number of SNPs: " + snpCount);

		if (snpCount > Integer.MAX_VALUE) {
			throw new GenotypeDataException("Found more than (2^31)-1 SNPs in bgen file. This is not supported");
		}

		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		long sampleCount = getUInt32(byteArray4, 0);
		LOGGER.debug("Number of samples: " + sampleCount);

		//Skip over reserved and free area
		this.bgenFile.seek(headerSize);

		//Read flags
		if (this.bgenFile.read(byteArray4, 0, 4) != 4) {
			throw new GenotypeDataException("Error reading bgen file header. File is corrupt");
		}
		if ((byteArray4[0] & 1) == 1) {
			LOGGER.debug("Genotype data is compressed");
			compressedGenotypes = true;
		} else {
			LOGGER.debug("Genotype data is not compressed");
			compressedGenotypes = false;
		}

		if ((byteArray4[0] & 3) == 1) {
			LOGGER.debug("v1.1 gen file");
			v1_1 = true;
		} else {
			LOGGER.debug("v1.0 gen file");
			v1_1 = false;
			throw new GenotypeDataException("bgen version v1.0 is not supported at the moment. Feel free to contact the developers of Genotype IO");
		}

		GeneticVariantRange.createRangeFactory((int) snpCount);
		
		byte[] snpInfoBuffer = new byte[8096];
		this.bgenFile.seek(snpOffset + 4);
		long lastSnpStart = snpOffset + 4;

		for (int snpI = 0; snpI < snpCount; ++snpI) {
		
			int snpInfoBufferSize = this.bgenFile.read(snpInfoBuffer, 0, snpInfoBuffer.length);
			int snpInfoBufferPos = 0;
		
			if (snpInfoBufferSize < 20) {
				throw new GenotypeDataException("Error reading bgen snp data. File is corrupt");
			}
			
			long sampleCountSnp = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			System.out.println("test " + sampleCountSnp);
			snpInfoBufferPos +=4;
			
			if(sampleCountSnp != sampleCount){
				throw new GenotypeDataException("Error reading bgen variant data. General sample count is not the same as current variant sample count");
			}
			
			
			int fieldLength = getUInt16(snpInfoBuffer ,snpInfoBufferPos);
			LOGGER.debug("Snp ID length " + fieldLength);
			this.bgenFile.skipBytes(fieldLength);
			
			snpInfoBufferPos += 2 + fieldLength; // skip id length and snp id
			
			fieldLength = getUInt16(snpInfoBuffer ,snpInfoBufferPos);
			LOGGER.debug("Snp RS length " + fieldLength);
			
			snpInfoBufferPos += 2;
			String snpId = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, charset);
			System.out.println(snpId);
			
			snpInfoBufferPos += fieldLength;
			
			fieldLength = getUInt16(snpInfoBuffer ,snpInfoBufferPos);
			snpInfoBufferPos += 2;
			String seqName = new String(snpInfoBuffer, snpInfoBufferPos, fieldLength, charset);
			System.out.println(seqName);
			snpInfoBufferPos += fieldLength;
			
			long snpPosLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 4;
			if(snpPosLong > Integer.MAX_VALUE){
				throw new GenotypeDataException("SNP pos larger than (2^31)-1 not supported");
			}
			
			int snpPos = (int) snpPosLong;
			
			System.out.println("SNP pos " + snpPos);
			
			long fieldLengthLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 4;
			if(fieldLengthLong > Integer.MAX_VALUE){
				throw new GenotypeDataException("SNP with allele longer than (2^31)-1 characters not supported");
			}
			String a1 = new String(snpInfoBuffer, snpInfoBufferPos, (int)fieldLengthLong, charset);
			snpInfoBufferPos += ((int) fieldLengthLong);
			
			fieldLengthLong = getUInt32(snpInfoBuffer, snpInfoBufferPos);
			snpInfoBufferPos += 4;
			if(fieldLengthLong > Integer.MAX_VALUE){
				throw new GenotypeDataException("SNP with allele longer than (2^31)-1 characters not supported");
			}
			String a2 = new String(snpInfoBuffer, snpInfoBufferPos, (int)fieldLengthLong, charset);
			snpInfoBufferPos += ((int) fieldLengthLong);
			
			System.out.println(a1 + "/" + a2);
			
			long snpBlockSize;
			if(compressedGenotypes){
				snpBlockSize = getUInt32(snpInfoBuffer, snpInfoBufferPos) + 4; 
			} else {
				snpBlockSize = 6 * sampleCount;
			}
			
			byte[] snpBlockData = new byte[6 * (int) sampleCount];
			inflater.setInput(snpInfoBuffer, snpInfoBufferPos + 4, (int) snpBlockSize - 4);
			try {
				inflater.inflate(snpBlockData);
			} catch (DataFormatException ex) {
				throw new GenotypeDataException("Error decompressing bgen data", ex);
			}
			inflater.reset();
			System.out.println(getUInt16(snpBlockData, 0) / 32768f + " " + getUInt16(snpBlockData, 2) / 32768f + " " + getUInt16(snpBlockData, 4) / 32768f );
			
			//move to next SNP
			lastSnpStart = lastSnpStart + snpInfoBufferPos + snpBlockSize;
			this.bgenFile.seek(lastSnpStart);
				
		}

		

	}

	/**
	 * Convert 4 bytes to unsigned 32 bit int from index. Returns long since java
	 * does not have unsigned int
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
	
}
