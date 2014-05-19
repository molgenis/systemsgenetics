package org.molgenis.vcf;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertTrue;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Iterator;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.molgenis.vcf.meta.VcfMeta;
import org.testng.annotations.Test;

import com.google.common.collect.Lists;

public class VcfReaderTest
{
	private static String vcfStr = "##fileformat=VCFv4.1\n" +
			"##mycustomheader=mycustomvalue\n" +
			"##contig=<ID=1,length=249240621>\n" +
			"##contig=<ID=2,length=259240621>\n" +
			"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" + 
			"##FILTER=<ID=flt,Description=\"Failing one of the filters\">\n" + 
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" + 
			"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n" + 
			"1	565286	rs1578391	C	T	.	flt	NS=1;DP=5;AF=1.000;ANNOT=INT;GI=LOC100131754\n";
	
	@Test
	public void close_Reader() throws IOException
	{
		Reader reader = mock(Reader.class);
		new VcfReader(reader).close();
		verify(reader).close();
	}

	@Test
	public void close_BlockCompressedInputStream() throws IOException
	{
		BlockCompressedInputStream blockCompressedInputStream = mock(BlockCompressedInputStream.class);
		new VcfReader(blockCompressedInputStream).close();
		verify(blockCompressedInputStream).close();
	}
	
	@Test
	public void getVcfMeta() throws IOException
	{
		VcfReader vcfReader = new VcfReader(new StringReader(vcfStr));
		try {
			VcfMeta vcfMeta = vcfReader.getVcfMeta();
			assertEquals(vcfMeta.getFileFormat(), "VCFv4.1");
			assertEquals(vcfMeta.get("mycustomheader"), "mycustomvalue");
			assertEquals(Lists.newArrayList(vcfMeta.getContigMeta()).size(), 2);
			assertEquals(Lists.newArrayList(vcfMeta.getInfoMeta()).size(), 1);
			assertEquals(Lists.newArrayList(vcfMeta.getFilterMeta()).size(), 1);
			assertEquals(Lists.newArrayList(vcfMeta.getFormatMeta()).size(), 1);
			assertEquals(vcfMeta.getColNames(), new String[]{"#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"});
		} finally {
			vcfReader.close();
		}
	}

	@Test
	public void iterator() throws IOException
	{
		VcfReader vcfReader = new VcfReader(new StringReader(vcfStr));
		try {
			
			Iterator<VcfRecord> it = vcfReader.iterator();
			assertTrue(it.hasNext());
			assertTrue(it.hasNext());
			assertNotNull(it.next());
			assertFalse(it.hasNext());
			assertFalse(it.hasNext());
		} finally {
			vcfReader.close();
		}
	}
}

