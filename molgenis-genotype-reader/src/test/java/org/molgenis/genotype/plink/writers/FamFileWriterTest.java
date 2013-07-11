package org.molgenis.genotype.plink.writers;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.plink.datatypes.FamEntry;
import org.testng.Assert;
import org.testng.annotations.Test;

public class FamFileWriterTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void FamFileWriter() throws IOException
	{
		FamFileWriter fileWriter = null;
		try
		{
			fileWriter = new FamFileWriter(null);
		}
		finally
		{
			if (fileWriter != null) fileWriter.close();
		}
	}

	@Test
	public void writeFamEntry() throws IOException
	{
		File file0 = File.createTempFile("FamFileWriterTest_file0", null);
		try
		{
			FamFileWriter fileWriter = null;
			try
			{
				fileWriter = new FamFileWriter(file0);
				fileWriter.write(new FamEntry("1", "Oleksandr", "0", "0", (byte) 1, 1.0));
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 Oleksandr 0 0 1 1.0\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}

	@Test
	public void writeIterableFamEntry() throws IOException
	{
		List<FamEntry> entryList = new ArrayList<FamEntry>();
		entryList.add(new FamEntry("1", "Oleksandr", "0", "0", (byte) 1, 1.0));
		entryList.add(new FamEntry("2", "Maksym", "0", "0", (byte) 1, 1.0));

		File file0 = File.createTempFile("FamFileWriterTest_file0", null);
		try
		{
			FamFileWriter fileWriter = null;
			try
			{
				fileWriter = new FamFileWriter(file0);
				fileWriter.write(entryList);
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 Oleksandr 0 0 1 1.0\n2 Maksym 0 0 1 1.0\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}
}
