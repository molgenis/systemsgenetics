package org.molgenis.vcf;

import static org.mockito.Mockito.mock;
import static org.testng.Assert.assertEquals;

import java.io.Reader;
import java.io.StringReader;
import java.util.List;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.molgenis.vcf.meta.VcfMeta;
import org.testng.annotations.Test;

import com.google.common.collect.Lists;

public class VcfRecordReaderTest
{
	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfRecordReaderReaderVcfMeta()
	{
		new VcfRecordReader(mock(Reader.class), null);
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void VcfRecordReaderBlockCompressedInputStreamVcfMeta()
	{
		new VcfRecordReader(mock(BlockCompressedInputStream.class), null);
	}

	@Test
	public void iterator()
	{
		String str = "1	565286	rs1578391	C	T	.	flt	NS=1;DP=5;AF=1.000;ANNOT=INT;GI=LOC100131754\n" + 
				 "2	575286	rs1578391	C	T	.	flt	NS=1;DP=5;AF=1.000;ANNOT=INT;GI=LOC100131754";
	
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfRecordReader vcfRecordReader = new VcfRecordReader(new StringReader(str), vcfMeta );
		List<VcfRecord> records = Lists.newArrayList(vcfRecordReader.iterator());
		assertEquals(records.size(), 2);
	}
	
	@Test
	public void iterator_noItems()
	{	
		VcfMeta vcfMeta = mock(VcfMeta.class);
		VcfRecordReader vcfRecordReader = new VcfRecordReader(new StringReader(""), vcfMeta );
		List<VcfRecord> records = Lists.newArrayList(vcfRecordReader.iterator());
		assertEquals(records.size(), 0);
	}
}
