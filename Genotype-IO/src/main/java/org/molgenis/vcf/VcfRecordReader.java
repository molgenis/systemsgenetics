package org.molgenis.vcf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Iterator;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.vcf.meta.VcfMeta;

public class VcfRecordReader implements Iterable<VcfRecord>
{
	private final BufferedReader reader;
	private final BlockCompressedInputStream blockCompressedInputStream;
	private final VcfMeta vcfMeta;

	public VcfRecordReader(Reader reader, VcfMeta vcfMeta) {
		if(reader == null) throw new IllegalArgumentException("reader is null");
		if(vcfMeta == null) throw new IllegalArgumentException("vcfMeta is null");
		this.reader = reader instanceof BufferedReader ? (BufferedReader) reader : new BufferedReader(reader);
		this.vcfMeta = vcfMeta;
		this.blockCompressedInputStream = null;
	}
	
	public VcfRecordReader(BlockCompressedInputStream blockCompressedInputStream, VcfMeta vcfMeta) {
		if(blockCompressedInputStream == null) throw new IllegalArgumentException("blockCompressedInputStream is null");
		if(vcfMeta == null) throw new IllegalArgumentException("vcfMeta is null");
		this.blockCompressedInputStream = blockCompressedInputStream;
		this.vcfMeta = vcfMeta;
		this.reader = null;
	}
	
	@Override
	public Iterator<VcfRecord> iterator()
	{			
		return new Iterator<VcfRecord>() {
			private final VcfRecord recycableVcfRecord = new VcfRecord(vcfMeta);
			private String[] tokens;

			@Override
			public boolean hasNext() {
				if (tokens == null) {
					try {
						tokens = StringUtils.split(readLine(), '\t');
					}
					catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
				return tokens != null;
			}

			@Override
			public VcfRecord next() {
				recycableVcfRecord.reset(tokens);
				tokens = null;
				return recycableVcfRecord;
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	private String readLine() throws IOException {
		String line; 
		if(reader != null) line = reader.readLine();
		else line = blockCompressedInputStream.readLine();
		return line;
	}
}
