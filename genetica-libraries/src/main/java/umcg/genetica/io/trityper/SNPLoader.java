/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.trityper;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;


/**
 * @author harm-jan
 */
public class SNPLoader {
	
	private RandomAccessFile m_genotypehandle;
	private RandomAccessFile m_dosagehandle;
	
	private int m_numIndividuals;
	
	
	private final Boolean[] m_isIncluded, m_isFemale;
	
	private MappedByteBuffer mappedGenotypeHandle = null;
	private MappedByteBuffer mappedDosageHandle = null;
	
	private long currentGtMapStart;
	private long currentGtMapEnd;
	private long currentDosageMapStart;
	private long currentDosageMapEnd;
	private byte[] bDs;
	private byte[] bGt;
	
	public SNPLoader(RandomAccessFile genotypehandle, Boolean[] indIsIncluded, Boolean[] isFemale) {
		this(genotypehandle, null, indIsIncluded, isFemale, 1000);
	}
	
	public SNPLoader(RandomAccessFile genotypehandle, RandomAccessFile dosagehandle, Boolean[] indIsIncluded, Boolean[] isFemale) {
		this(genotypehandle, dosagehandle, indIsIncluded, isFemale, 1000);
	}
	
	int gtmaplen;
	int dsmaplen;
	int numberOfVariantsInMemoryMap = 1000;
	
	public SNPLoader(RandomAccessFile genotypeHandle, RandomAccessFile dosageHandle, Boolean[] isIncluded, Boolean[] isFemale, int numberOfVariantsToBuffer) {
		
		m_genotypehandle = genotypeHandle;
		m_dosagehandle = dosageHandle;
		m_isIncluded = isIncluded;
		m_isFemale = isFemale;
		this.numberOfVariantsInMemoryMap = numberOfVariantsToBuffer;
	}
	
	public void loadGenotypes(SNP snp) throws IOException {
		byte[] allele1 = new byte[m_numIndividuals];
		byte[] allele2 = new byte[m_numIndividuals];
		
		int nrBytesToRead = m_numIndividuals * 2;
		long seekLoc = (long) snp.getId() * (long) nrBytesToRead;

//		byte[] alleles = new byte[nrBytesToRead];
		
		// initiate buffer if it doesn't exist, or if location we're looking for is beyond the current map
		long seekEnd = seekLoc + (m_numIndividuals * 2);
//		System.out.println("Seekloc: " + seekLoc);
//		System.out.println("SeekEnd: " + seekEnd);
		if (mappedGenotypeHandle == null || seekLoc < currentGtMapStart || seekLoc > currentGtMapEnd || seekEnd > currentGtMapEnd) {
			// 32 megabytes worth of variants; (32*1048576)/(m_numIndividuals * 2) bytes
			if (mappedGenotypeHandle == null) {
				int bytesPerVariant = (m_numIndividuals * 2);
				int nrBytesPerBuffer = bytesPerVariant * numberOfVariantsInMemoryMap; //(32 * 1048576);
//				int remainder = nrBytesPerBuffer % bytesPerVariant;
//				nrBytesPerBuffer += remainder;
				gtmaplen = nrBytesPerBuffer;
				dsmaplen = nrBytesPerBuffer / 2;
//				System.out.println("bytes in buffer1: " + gtmaplen);
//				System.out.println("bytes in buffer2: " + dsmaplen);
			}
			long maplentouse = gtmaplen;
			// prevent overflow
			if (seekLoc + maplentouse > m_genotypehandle.length()) {
				maplentouse = m_genotypehandle.length() - seekLoc;
			}
			
			
			mappedGenotypeHandle = m_genotypehandle.getChannel().map(FileChannel.MapMode.READ_ONLY, seekLoc, maplentouse);
			mappedGenotypeHandle.load();
			bGt = new byte[(int) maplentouse];
			mappedGenotypeHandle.get(bGt);
			currentGtMapStart = seekLoc;
			currentGtMapEnd = currentGtMapStart + maplentouse;
//			System.out.println("Reload buffer:\t" + seekLoc + "\tstart\t" + currentGtMapStart + "\tstop\t" + currentGtMapEnd + "\tlen\t" + maplentouse);
			
		}

//		if (m_genotypehandle.getFilePointer() != seekLoc) {
//			m_genotypehandle.seek(seekLoc);
//		}

//		m_genotypehandle.read(alleles, 0, bytesize);
//		mappedGenotypeHandle(alleles, 0, bytesize);
		
		// recalculate where we should be looking for this particular snp


//		mappedGenotypeHandle.get(alleles, offset, nrBytesToRead);
//		mappedGenotypeHandle.slice();
		
		int offset = (int) (seekLoc - currentGtMapStart);
//		System.out.println("Seekloc: " + seekLoc);
//		System.out.println("offset: " + offset);
//		System.out.println("btr: " + nrBytesToRead);
//		System.out.println("capacity: " + mappedGenotypeHandle.capacity());
		System.arraycopy(bGt, offset, allele1, 0, m_numIndividuals);
		System.arraycopy(bGt, offset + m_numIndividuals, allele2, 0, m_numIndividuals);
//		System.arraycopy(alleles, 0, allele1, 0, m_numIndividuals);
//		System.arraycopy(alleles, m_numIndividuals, allele2, 0, m_numIndividuals);

//		alleles = null;
		
		snp.setAlleles(allele1, allele2, m_isIncluded, m_isFemale);
		
	}
	
	public void loadDosage(SNP snp) throws IOException {
		if (m_dosagehandle != null) {
			byte[] dosageValues = new byte[m_numIndividuals];
			//if (loadedSNP.getGcScores()==null||loadedSNP.getThetaValues()==null||loadedSNP.getRValues()==null) {
			long seekLoc = (long) snp.getId() * (long) m_numIndividuals;
			
			// initiate buffer if it doesn't exist, or if location we're looking for is beyond the current map
			long seekEnd = seekLoc + (m_numIndividuals);
			if (mappedDosageHandle == null || seekLoc < currentDosageMapStart || seekLoc > currentDosageMapEnd || seekEnd > currentDosageMapEnd) {
				long maplentouse = dsmaplen;
				// prevent overflow
				if (seekLoc + maplentouse > m_dosagehandle.length()) {
					maplentouse = m_dosagehandle.length() - seekLoc;
				}
				
				mappedDosageHandle = m_dosagehandle.getChannel().map(FileChannel.MapMode.READ_ONLY, seekLoc, maplentouse);
				mappedDosageHandle.load();
				bDs = new byte[(int) maplentouse];
				mappedDosageHandle.get(bDs);
				
				currentDosageMapStart = seekLoc;
				currentDosageMapEnd = currentDosageMapStart + maplentouse;
			}

//			m_dosagehandle.seek(seekLoc);
//			m_dosagehandle.read(dosageValues, 0, m_numIndividuals);
//
			int offset = (int) (seekLoc - currentDosageMapStart);
			System.arraycopy(bDs, offset, dosageValues, 0, m_numIndividuals);
			
			
			byte[] genotypes = snp.getGenotypes();
			
			boolean takeComplement = false;
			for (int ind = 0; ind < dosageValues.length; ind++) {
				double dosagevalue = ((double) (-Byte.MIN_VALUE + dosageValues[ind])) / 100;
				if (genotypes[ind] == 0 && dosagevalue > 1) {
					takeComplement = true;
					break;
				}
				if (genotypes[ind] == 2 && dosagevalue < 1) {
					takeComplement = true;
					break;
				}
			}
			if (takeComplement) {
				for (int ind = 0; ind < dosageValues.length; ind++) {
					byte dosageValue = (byte) (200 - (-Byte.MIN_VALUE + dosageValues[ind]) + Byte.MIN_VALUE);
					dosageValues[ind] = dosageValue;
				}
			}
			
			snp.setDosage(dosageValues);
		}
	}
	
	/**
	 * @return the numIndividuals
	 */
	public int getNumIndividuals() {
		return m_numIndividuals;
	}
	
	/**
	 * @param numIndividuals the numIndividuals to set
	 */
	public void setNumIndividuals(int numIndividuals) {
		this.m_numIndividuals = numIndividuals;
		
	}
	
	public boolean hasDosageInformation() {
		return (m_dosagehandle != null);
	}
	
	public double getAverageSNPSize(int numSNPs) throws IOException {
		long size = 0;
		
		size += m_genotypehandle.length();
		if (m_dosagehandle != null) {
			size += m_dosagehandle.length();
		}
		
		
		double avgSNPSize = 0;
		if (size > 0) {
			avgSNPSize = (double) size / numSNPs;
		}
		
		return avgSNPSize;
	}
	
	public void close() throws IOException {
		if (m_dosagehandle != null) {
			m_dosagehandle.close();
		}
		m_genotypehandle.close();
	}
}
