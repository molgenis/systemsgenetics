package org.molgenis.genotype.plink.writers;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.plink.datatypes.BimEntry;
import org.testng.Assert;
import org.testng.annotations.Test;

public class BimFileWriterTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void BimFileWriter() throws IOException
	{
		BimFileWriter fileWriter = null;
		try
		{
			fileWriter = new BimFileWriter(null);
		}
		finally
		{
			if (fileWriter != null) fileWriter.close();
		}
	}

	@Test
	public void writeBimEntry() throws IOException
	{
		File file0 = File.createTempFile("BimFileWriterTest_file0", null);
		try
		{
			BimFileWriter fileWriter = null;
			try
			{
				fileWriter = new BimFileWriter(file0);
				fileWriter.write(new BimEntry("1", "snp1", 0.0, 1, Alleles.createBasedOnChars('A', 'C')));
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 snp1 0.0 1 A C\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}

	@Test
	public void writeIterableBimEntry() throws IOException
	{
		List<BimEntry> entryList = new ArrayList<BimEntry>();
		entryList.add(new BimEntry("1", "snp1", 0.0, 1, Alleles.createBasedOnChars('A', 'C')));
		entryList.add(new BimEntry("2", "snp2", 1.2, 1, Alleles.createBasedOnChars('C', 'A')));

		File file0 = File.createTempFile("BimFileWriterTest_file0", null);
		try
		{
			BimFileWriter fileWriter = null;
			try
			{
				fileWriter = new BimFileWriter(file0);
				fileWriter.write(entryList);
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 snp1 0.0 1 A C\n2 snp2 1.2 1 C A\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}
}
