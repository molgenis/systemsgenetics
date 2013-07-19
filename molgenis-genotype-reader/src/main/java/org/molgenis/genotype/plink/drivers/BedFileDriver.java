package org.molgenis.genotype.plink.drivers;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import org.molgenis.genotype.GenotypeDataException;

/**
 * Driver to query BED (binary Plink genotype) files. See:
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
 * 
 * PLEASE NOTE THAT: this driver at the moment works ONLY on SNP-major mode
 * files!
 * 
 */
public class BedFileDriver
{
	private int mode;
	private long nrOfElements;
	private File bedFile;

	/**
	 * Get the mode: mode 1 = SNP-major, mode 0 = individual-major
	 * 
	 * @return
	 */
	public int getMode()
	{
		return mode;
	}

	/**
	 * Get the number of retrievable genotype elements of this BED file. Does
	 * not account for trailing null elements because they are indistinguishable
	 * from null genotypes in this format alone.
	 * 
	 * @return
	 */
	public long getNrOfElements()
	{
		return nrOfElements;
	}

	/**
	 * Construct new convertGenoCoding on this file
	 * 
	 * @param bedFile
	 * @throws Exception
	 */
	public BedFileDriver(File bedFile) throws IOException
	{
		RandomAccessFile raf = new RandomAccessFile(bedFile, "r");
		try
		{
			this.bedFile = bedFile;

			byte mn1 = raf.readByte();
			byte mn2 = raf.readByte();

			if (mn1 == 108 && mn2 == 27) // tested, bit code 01101100 00011011
			{
				// System.out.println("Plink magic number valid");
			}
			else
			{
				throw new GenotypeDataException("Invalid Plink magic number");
			}

			byte bmode = raf.readByte();

			if (bmode == 1) // tested, bit code 00000001
			{
				// System.out.println("mode 1: SNP-major");
			}
			else if (bmode == 0) // assumed... bit code 00000000
			{
				// System.out.println("mode 0: individual-major");
				throw new GenotypeDataException("BED file individual-major mode not yet supported!");
			}
			else
			{
				throw new GenotypeDataException("Mode not recognized: " + bmode);
			}

			this.mode = bmode;
			this.nrOfElements = (raf.length() - 3) * 4;
		}
		finally
		{
			raf.close();
		}
	}

	/**
	 * Convert bit coding in custom genotype coding.
	 * 
	 * @param in
	 * @param hom1
	 * @param hom2
	 * @param het
	 * @param _null
	 * @return
	 * @throws Exception
	 */
	public String convertGenoCoding(String in, String hom1, String hom2, String het, String _null) throws Exception
	{
		if (in.equals("00"))
		{
			return hom1;
		}
		if (in.equals("01"))
		{
			return het;
		}
		if (in.equals("11"))
		{
			return hom2;
		}
		if (in.equals("10"))
		{
			return _null;
		}
		throw new Exception("Input '" + in + "' not recognized");
	}

	/**
	 * Convert bit coding in common genotype signs: A & B for homozygotes, H for
	 * heterozygote, N for null.
	 * 
	 * @param in
	 * @return
	 * @throws Exception
	 */
	public String genoCodingCommon(String in) throws Exception
	{
		return convertGenoCoding(in, "A", "B", "H", "N");
	}

	/**
	 * Get a single element from the BED file
	 * 
	 * @param index
	 * @return
	 * @throws Exception
	 */
	public String getElement(long index) throws Exception
	{
		// throw new Exception("fixme!");
		RandomAccessFile raf = new RandomAccessFile(bedFile, "r");
		raf.seek((index / 4) + 3);
		String byteString = reverse(bits(raf.readByte()));
		raf.close();
		int bitpair = (int) (index % 4) * 2;
		return byteString.substring(bitpair, bitpair + 2);
	}

	/**
	 * Get the SNP collection for this SNP index in the BIM file Returns all
	 * individuals for this set of SNPs You need to supply the numnber of
	 * individuals so that the reader knows 1) when to stop reading bytes and 2)
	 * how many padding bits to compensate for when reading any index > 0.
	 * 
	 * @param index
	 * @param nrOfIndividualsInFAMfile
	 * @return
	 * @throws Exception
	 */
	public String[] getSNPs(long index, int nrOfIndividualsInFAMfile) throws Exception
	{
		// calculate the number of individuals in the byte that is potentially
		// padded
		int nrOfIndividualsInPaddedByte = nrOfIndividualsInFAMfile % 4;

		// if nrOfIndividualsInPaddedByte is 0, there are no padding bit pairs
		// (paddingIndividuals = 0)
		// else, its 4 minus the amount of individuals in the padding byte (1, 2
		// or 3)
		int paddingIndividuals = nrOfIndividualsInPaddedByte == 0 ? 0 : 4 - nrOfIndividualsInPaddedByte;

		// check if we got the numbers right.. we want a multiplication of 4
		if ((nrOfIndividualsInFAMfile + paddingIndividuals) % 4 != 0)
		{
			throw new Exception("nrOfIndividuals + paddingIndividuals) % 4 must be 0");
		}

		int bytesPerIndividual = (nrOfIndividualsInFAMfile + paddingIndividuals) / 4;

		// inclusive: read this byte index
		// add 3 because of the reserved bytes in plink format
		long startByte = (index * bytesPerIndividual) + 3;

		long stopByte = startByte + bytesPerIndividual;

		byte[] res = new byte[(int) (stopByte - startByte)];

		RandomAccessFile raf = new RandomAccessFile(bedFile, "r");
		raf.seek(startByte);
		raf.read(res);
		raf.close();

		String[] result = new String[nrOfIndividualsInFAMfile];
		int res_index = 0;

		for (int i = 0; i < res.length; i++)
		{
			byte b = res[i];

			String byteString = reverse(bits(b));

			int toBit = 8; // normally we take the whole byte
			if (i == res.length - 1) // except at the end, when we correct for
										// padding 0's
			{
				// At the end, the string is padded with 0's -> check
				for (int j = paddingIndividuals * 2; j < 8; j++)
				{
					if (byteString.charAt(j) != '0')
					{
						throw new IOException("Fatal error: padding 0's not present where expected!");
					}
				}
				toBit -= (paddingIndividuals * 2);
			}
			for (int pair = 0; pair < toBit; pair += 2)
			{
				result[res_index++] = byteString.substring(pair, pair + 2);
			}
		}
		return result;
	}

	/**
	 * 
	 * UNIMPLEMENTED
	 * 
	 * Get a String[] of elements from the BED file. This function returns the
	 * elements in their intended order, with bits already reversed, and taking
	 * into account the padding bits at the end of each sequence of SNPs. (SNP
	 * major mode: lists ALL individuals for 1 SNP, then moves on to the next
	 * SNP)
	 * 
	 * From -> to is NOT whole SNPs (all individuals), but single SNP elements
	 * within the file Reason is that this function can be used for reading
	 * batches on the smallest elements to get a complete SNP set (all
	 * individuals for 1 snp) using this function, set your from / to arguments
	 * to match up with your number of individuals (in which case, the
	 * nrOfIndividuals argument is (to minus from)
	 * 
	 * @param from
	 *            : the starting bitpair, inclusive
	 * @param to
	 *            : last bitpair to read, exclusive (e.g. reading first element
	 *            is 0, 1)
	 * @param nrOfIndividualsInFAMfile
	 *            : need to correct in two ways: 1) when to adjust for padding
	 *            when reading SNPs from the bytes and 2) to know the number of
	 *            padding bytes: by taking modulo 4 (4 bitpairs in a byte, 0
	 *            padding when the amount of individuals is a multiplication of
	 *            4, etc)
	 * @return
	 * @throws Exception
	 */
	public String[] getElements(long from, long to, int nrOfIndividualsInFAMfile) throws Exception
	{
		if (1 == 1) throw new UnsupportedOperationException("Not yet implemented");

		// calculate the number of individuals in the byte that is potentially
		// padded
		int nrOfIndividualsInPaddedByte = nrOfIndividualsInFAMfile % 4;

		// if nrOfIndividualsInPaddedByte is 0, there are no padding bit pairs
		// (paddingIndividuals = 0)
		int paddingIndividuals = nrOfIndividualsInPaddedByte == 0 ? 0 : 4 - nrOfIndividualsInPaddedByte;

		// check if we got the numbers right.. we want a multiplication of 4
		if ((nrOfIndividualsInFAMfile + paddingIndividuals) % 4 != 0)
		{
			throw new Exception("nrOfIndividuals + paddingIndividuals) % 4 must be 0");
		}

		// calculate the amount of bytes that must be read to get the SNPs set
		// (1 snp for all individuals)
		int totalBytesToReadPerIndividual = (nrOfIndividualsInFAMfile + paddingIndividuals) / 4;

		// inclusive: read this byte index
		long startByte = 0;
		// exclusive: stop when reaching this index
		long stopByte = 0;

		byte[] res = new byte[(int) (stopByte - startByte)];

		String[] placeholder = new String[0];
		return placeholder;
	}

	/**
	 * Get a String[] of elements from the BED file.
	 * 
	 * NOTE: the 'pass' argument is vague.. it has to do with the multiplication
	 * of the padding, and is implicit on e.g. a for loop iteration retrieving
	 * all individuals, so for the next 'pass', more padding must be added to
	 * correct for the offset.
	 * 
	 * NOTE2: this function is difficult to use correctly and counter intuitive,
	 * please use public String[] getElements(long from, long to, int
	 * nrOfIndividuals) instead!
	 * 
	 * from = inclusive to = exclusive
	 * 
	 * @param from
	 * @param to
	 * @return
	 * @throws IOException
	 */
	public String[] getElements(long from, long to, int paddingBitpairs, int pass) throws IOException
	{
		double paddingFraction = paddingBitpairs / 4.0;
		// Start byte = byte position of start individual, corrected for padding
		// 0's that get added at every SNP:
		long start = (long) ((from / 4.0) - (pass * paddingFraction) + pass + 3);
		// Stop byte = byte position after last individual, corrected for
		// padding 0's that get added at every SNP:
		long stop = (long) ((to / 4.0) - ((pass + 1) * paddingFraction) + (pass + 1) + 3);
		System.out.println("GOING TO READ FROM BYTE " + start + " TO " + stop);
		byte[] res = new byte[(int) (stop - start)];
		int res_index = 0;
		String[] result = new String[(int) (to - from)]; // to - from = nr. of
															// individuals
		RandomAccessFile raf = new RandomAccessFile(bedFile, "r");
		raf.seek(start);
		raf.read(res);
		raf.close();

		for (int i = 0; i < res.length; i++)
		{
			byte b = res[i];
			String byteString = reverse(bits(b));

			int toBit = 8; // normally we take the whole byte
			if (i == res.length - 1) // except at the end, when we correct for
										// padding 0's
			{
				// At the end, the string is padded with 0's -> check
				for (int j = paddingBitpairs * 2; j < 8; j++)
				{
					if (byteString.charAt(j) != '0')
					{
						throw new IOException("Fatal error: padding 0's not present where expected!");
					}
				}
				toBit -= (paddingBitpairs * 2);
			}

			for (int pair = 0; pair < toBit; pair += 2)
			{
				result[res_index++] = byteString.substring(pair, pair + 2);
			}
		}
		return result;
	}

	/**
	 * Helper function to get the bit values
	 * 
	 * @param b
	 * @return
	 */
	private String bits(byte b)
	{
		StringBuilder bitsBuilder = new StringBuilder();
		for (int bit = 7; bit >= 0; --bit)
		{
			bitsBuilder.append(((b >>> bit) & 1));
		}
		return bitsBuilder.toString();
	}

	/**
	 * Helper function to reverse a string
	 * 
	 * @param string
	 * @return
	 */
	private String reverse(String string)
	{
		return new StringBuffer(string).reverse().toString();
	}

}
