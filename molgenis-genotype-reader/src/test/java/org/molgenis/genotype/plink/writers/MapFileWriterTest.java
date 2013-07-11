package org.molgenis.genotype.plink.writers;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.plink.datatypes.MapEntry;
import org.testng.Assert;
import org.testng.annotations.Test;

public class MapFileWriterTest
{

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void MapFileWriter() throws IOException
	{
		MapFileWriter fileWriter = null;
		try
		{
			fileWriter = new MapFileWriter(null);
		}
		finally
		{
			if (fileWriter != null) fileWriter.close();
		}
	}

	@Test
	public void writeMapEntry() throws IOException
	{
		File file0 = File.createTempFile("MapFileWriterTest_file0", null);
		try
		{
			MapFileWriter fileWriter = null;
			try
			{
				fileWriter = new MapFileWriter(file0);
				fileWriter.write(new MapEntry("1", "snp1", 0, 1l));
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 snp1 0 1\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}

	@Test
	public void writeIterableMapEntry() throws IOException
	{
		List<MapEntry> entryList = new ArrayList<MapEntry>();
		entryList.add(new MapEntry("1", "snp1", 0, 1l));
		entryList.add(new MapEntry("1", "snp2", 0, 2l));

		File file0 = File.createTempFile("MapFileWriterTest_file0", null);
		try
		{
			MapFileWriter fileWriter = null;
			try
			{
				fileWriter = new MapFileWriter(file0);
				fileWriter.write(entryList);
			}
			finally
			{
				IOUtils.closeQuietly(fileWriter);
			}

			String expected = "1 snp1 0 1\n1 snp2 0 2\n";
			Assert.assertEquals(FileUtils.readFileToString(file0, Charset.forName("UTF-8")), expected);
		}
		finally
		{
			file0.delete();
		}
	}
}
