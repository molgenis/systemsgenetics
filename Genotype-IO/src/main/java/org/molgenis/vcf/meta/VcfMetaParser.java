package org.molgenis.vcf.meta;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

import net.sf.samtools.util.BlockCompressedInputStream;

public class VcfMetaParser
{
	private static final String PREFIX_ALT = "##ALT";
	private static final String PREFIX_CONTIG = "##contig";
	private static final String PREFIX_FILTER = "##FILTER";
	private static final String PREFIX_FORMAT = "##FORMAT";
	private static final String PREFIX_INFO = "##INFO";
	private static final String PREFIX_PEDIGREE = "##PEDIGREE";
	private static final String PREFIX_SAMPLE = "##SAMPLE";
	
	private final BufferedReader reader;
	private final BlockCompressedInputStream blockCompressedInputStream;

	public VcfMetaParser(Reader reader) throws IOException {
		if(reader == null) throw new IllegalArgumentException("reader is null");
		this.reader = reader instanceof BufferedReader ? (BufferedReader) reader : new BufferedReader(reader);
		this.blockCompressedInputStream = null;
	}

	public VcfMetaParser(BlockCompressedInputStream blockCompressedInputStream) {
		if(blockCompressedInputStream == null) throw new IllegalArgumentException("blockCompressedInputStream is null");
		this.blockCompressedInputStream = blockCompressedInputStream;
		this.reader = null;
	}
	
	public VcfMeta parse() throws IOException {
		VcfMeta vcfMeta = new VcfMeta();
		
		String line;
		while((line = readLine()) != null && line.startsWith("##")) {
			if(line.startsWith(PREFIX_ALT)) {
				vcfMeta.addAltMeta(new VcfMetaAlt(parseMetaLine(line)));
			} else if(line.startsWith(PREFIX_CONTIG)) {
				vcfMeta.addContigMeta(new VcfMetaContig(parseMetaLine(line)));
			} else if(line.startsWith(PREFIX_FILTER)) {
				vcfMeta.addFilterMeta(new VcfMetaFilter(parseMetaLine(line)));
			} else if(line.startsWith(PREFIX_FORMAT)) {
				vcfMeta.addFormatMeta(new VcfMetaFormat(parseMetaLine(line)));
			} else if(line.startsWith(PREFIX_INFO)) {
				vcfMeta.addInfoMeta(new VcfMetaInfo(parseMetaLine(line)));
			} else if(line.startsWith(PREFIX_PEDIGREE)) {
				vcfMeta.addPedigreeMeta(new VcfMetaPedigree(parseMetaLine(line)));
			} else if(line.startsWith(PREFIX_SAMPLE)) {
				vcfMeta.addSampleMeta(new VcfMetaSample(parseMetaLine(line)));
			} else {
				int idx = line.indexOf('=');
				vcfMeta.add(line.substring(2, idx), line.substring(idx + 1));
			}
		}
		if(line == null || !line.startsWith("#CHROM")) {
			throw new IOException("missing column headers");
		}
		vcfMeta.setColNames(StringUtils.split(line, '\t'));

		return vcfMeta;
	}

	private Map<String, String> parseMetaLine(String line) {
		Map<String, String> properties = new LinkedHashMap<String, String>();
		
		// header block starts with < and ends with >
		boolean inHeaderBlock = false;
		// values can be quoted to allow for = and ,
		boolean inQuotes = false;
		// header is divided in key and value using '='
		boolean inKey = true;
		// to store a key while parsing
		String key = "";
		// to store a value while parsing
		String value = "";

		final int nrChars = line.length();
		for (int i = 0; i < nrChars; ++i)
		{
			char c = line.charAt(i);
			if (!inHeaderBlock)
			{
				if ('<' == c)
				{
					inHeaderBlock = true;
				}
			}
			else
			{
				// check for seperator between key and value
				if (inKey && '=' == c)
				{
					inKey = false;
				}
				// parse quotes
				else if (!inKey && '"' == c)
				{
					if (inQuotes)
					{
						inQuotes = false;
						properties.put(key, value);
						key = "";
						value = "";
						inKey = true;
					}
					else
						inQuotes = true;
				}
				// close the key/value pair
				else if (!inQuotes && !inKey && (',' == c || '>' == c))
				{
					properties.put(key, value);
					key = "";
					value = "";
					inKey = true;
				}
				// otherwise just add key/value char
				else
				{
					if (inKey) key += c;
					else
						value += c;
				}
			}
		}

		return properties;
	}
	
	private String readLine() throws IOException {
		String line;
		if(reader != null) line = reader.readLine();
		else line = blockCompressedInputStream.readLine();
		return line;
	}
}
