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
import org.molgenis.genotype.plink.datatypes.PedEntry;
import org.testng.Assert;
import org.testng.annotations.Test;

public class PedFileWriterTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void PedFileWriter() throws IOException
	{
		PedFileWriter fileWriter = null;
		try
		{
			fileWriter = new PedFileWriter(null);
		}
		finally
		{
			if (fileWriter != null) fileWriter.close();
		}
	}

	@Test
	public void writePedEntry() throws IOException
	{
		File file0 = File.createTempFile("PedFileWriterTest_file0", null);
		try
		{
			PedFileWriter fileWriter = null;
			try
			{
				fileWriter = new PedFileWriter(file0);
				fileWriter.write(new PedEntry("1", "1", "0", "0", (byte) 1, 1.0, Arrays.asList(
						Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'T')).iterator()));
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 1 0 0 1 1 A A G T\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}

	@Test
	public void writeIterablePedEntry() throws IOException
	{
		List<PedEntry> entryList = new ArrayList<PedEntry>();
		entryList.add(new PedEntry("1", "1", "0", "0", (byte) 1, 1.0, Arrays.asList(
				Alleles.createBasedOnChars('A', 'A'), Alleles.createBasedOnChars('G', 'T')).iterator()));
		entryList.add(new PedEntry("2", "1", "0", "0", (byte) 1, 1.0, Arrays.asList(
				Alleles.createBasedOnChars('A', 'C'), Alleles.createBasedOnChars('T', 'G')).iterator()));

		File file0 = File.createTempFile("PedFileWriterTest_file0", null);
		try
		{
			PedFileWriter fileWriter = null;
			try
			{
				fileWriter = new PedFileWriter(file0);
				fileWriter.write(entryList);
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 1 0 0 1 1 A A G T\n2 1 0 0 1 1.0 A C T G\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}
}
