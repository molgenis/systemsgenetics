package org.molgenis.genotype.plink.writers;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.plink.datatypes.TpedEntry;
import org.testng.Assert;
import org.testng.annotations.Test;

public class TpedFileWriterTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void TpedFileWriter() throws IOException
	{
		TpedFileWriter fileWriter = null;
		try
		{
			fileWriter = new TpedFileWriter(null);
		}
		finally
		{
			if (fileWriter != null) fileWriter.close();
		}
	}

	@Test
	public void writeTpedEntry() throws IOException
	{
		File file0 = File.createTempFile("TpedFileWriterTest_file0", null);
		try
		{
			TpedFileWriter fileWriter = null;
			try
			{
				fileWriter = new TpedFileWriter(file0);
				Alleles b1 = Alleles.createBasedOnChars('A', 'A');
				Alleles b2 = Alleles.createBasedOnChars('A', 'C');
				Alleles b3 = Alleles.createBasedOnChars('C', 'C');
				Alleles b4 = Alleles.createBasedOnChars('A', 'C');
				Alleles b5 = Alleles.createBasedOnChars('C', 'C');
				Alleles b6 = Alleles.createBasedOnChars('C', 'C');
				fileWriter.write(new TpedEntry("1", "snp1", 0.0, 5000650, Arrays.asList(b1, b2, b3, b4, b5, b6)));
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 snp1 0.0 5000650 A A A C C C A C C C C C\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}

	@Test
	public void writeIterableTpedEntry() throws IOException
	{
		List<TpedEntry> entryList = new ArrayList<TpedEntry>();
		Alleles b0_1 = Alleles.createBasedOnChars('A', 'A');
		Alleles b0_2 = Alleles.createBasedOnChars('A', 'C');
		Alleles b0_3 = Alleles.createBasedOnChars('C', 'C');
		Alleles b0_4 = Alleles.createBasedOnChars('A', 'C');
		Alleles b0_5 = Alleles.createBasedOnChars('C', 'C');
		Alleles b0_6 = Alleles.createBasedOnChars('C', 'C');
		entryList.add(new TpedEntry("1", "snp1", 0.0, 5000650, Arrays.asList(b0_1, b0_2, b0_3, b0_4, b0_5, b0_6)));
		Alleles b1_1 = Alleles.createBasedOnChars('G', 'T');
		Alleles b1_2 = Alleles.createBasedOnChars('G', 'T');
		Alleles b1_3 = Alleles.createBasedOnChars('G', 'G');
		Alleles b1_4 = Alleles.createBasedOnChars('T', 'T');
		Alleles b1_5 = Alleles.createBasedOnChars('G', 'T');
		Alleles b1_6 = Alleles.createBasedOnChars('T', 'T');
		entryList.add(new TpedEntry("1", "snp2", 0.0, 5000830, Arrays.asList(b1_1, b1_2, b1_3, b1_4, b1_5, b1_6)));

		File file0 = File.createTempFile("TpedFileWriterTest_file0", null);
		try
		{
			TpedFileWriter fileWriter = null;
			try
			{
				fileWriter = new TpedFileWriter(file0);
				fileWriter.write(entryList);
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 snp1 0.0 5000650 A A A C C C A C C C C C\n1 snp2 0.0 5000830 G T G T G G T T G T T T\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}
}
