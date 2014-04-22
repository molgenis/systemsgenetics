package org.molgenis.vcf;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.Reader;
import java.util.Iterator;

import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaParser;

import net.sf.samtools.util.BlockCompressedInputStream;

/**
 * A high-performance VCF reader based on the the Variant Call Format (VCF) Version 4.2 Specification
 * 
 * In order to achieve high performance VcfRecord and VcfSample objects are recycled,
 * use createClone() if you need to store these object in a iteration.  
 */
public class VcfReader implements Iterable<VcfRecord>, Closeable {
	private final BufferedReader reader;
	private final BlockCompressedInputStream blockCompressedInputStream;
	private VcfMeta vcfMeta;
	
	public VcfReader(Reader reader) throws IOException {
		if(reader == null) throw new IllegalArgumentException("reader is null");
		this.reader = reader instanceof BufferedReader ? (BufferedReader) reader : new BufferedReader(reader);
		this.blockCompressedInputStream = null;
	}

	public VcfReader(BlockCompressedInputStream blockCompressedInputStream) {
		if(blockCompressedInputStream == null) throw new IllegalArgumentException("blockCompressedInputStream is null");
		this.blockCompressedInputStream = blockCompressedInputStream;
		this.reader = null;
	}
	
	@Override
	public Iterator<VcfRecord> iterator()
	{
		if (vcfMeta == null) {
			try {
				vcfMeta = parseVcfMeta();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		VcfRecordReader vcfRecordReader = reader != null ? new VcfRecordReader(reader, vcfMeta) : new VcfRecordReader(blockCompressedInputStream, vcfMeta);  
		return vcfRecordReader.iterator();
	}
	
	public VcfMeta getVcfMeta() throws IOException {
		if (vcfMeta == null) {
			vcfMeta = parseVcfMeta();
		}
		return vcfMeta;
	}
	
	@Override
	public void close() throws IOException
	{
		if(reader != null) reader.close();
		else if (blockCompressedInputStream != null) blockCompressedInputStream.close();
	}

	private VcfMeta parseVcfMeta() throws IOException {
		if(reader != null) return new VcfMetaParser(reader).parse();
		else return new VcfMetaParser(blockCompressedInputStream).parse();
	}
}
