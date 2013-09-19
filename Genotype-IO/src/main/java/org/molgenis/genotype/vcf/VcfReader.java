package org.molgenis.genotype.vcf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.molgenis.io.TupleReader;
import org.molgenis.io.csv.CsvReader;
import org.molgenis.io.processor.CellProcessor;
import org.molgenis.util.tuple.Tuple;

import com.google.common.collect.Lists;

public class VcfReader implements TupleReader
{
	private static final Charset CHARSET_UTF8 = Charset.forName("UTF-8");
	private static final List<String> NORMAL_COL_NAMES = Arrays.asList(new String[]
	{ "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "" });

	private List<String> headers = new ArrayList<String>();
	private final CsvReader csvReader;

	public VcfReader(InputStream vcfInputStream) throws IOException
	{
		BufferedReader br = new BufferedReader(new InputStreamReader(vcfInputStream, CHARSET_UTF8));

		// First read the headers (starts with '##')
		br.mark(vcfInputStream.available());
		headers = readHeaders(br);
		br.reset();

		// Create CsvReader, skip the headers starting with '##'
		for (int i = 0; i < headers.size(); i++)
		{
			br.readLine();
		}

		csvReader = new CsvReader(br, '	', true);
	}

	public List<String> getSampleNames() throws IOException
	{
		List<String> sampleNames = new ArrayList<String>();

		Iterator<String> it = colNamesIterator();
		while (it.hasNext())
		{
			String colName = it.next();
			if (!NORMAL_COL_NAMES.contains(colName))
			{
				sampleNames.add(colName);
			}
		}

		return Collections.unmodifiableList(sampleNames);
	}

	public List<String> getColNames() throws IOException
	{
		return Collections.unmodifiableList(Lists.newArrayList(colNamesIterator()));
	}

	public List<VcfInfo> getInfos()
	{
		List<VcfInfo> infos = new ArrayList<VcfInfo>();

		for (String header : headers)
		{
			if (header.startsWith("##INFO"))
			{
				Tuple tuple = new VcfHeaderParser(header).parse();
				VcfInfo info = new VcfInfo(tuple);
				infos.add(info);
			}
		}

		return infos;
	}

	public List<VcfFormat> getFormats()
	{
		List<VcfFormat> formats = new ArrayList<VcfFormat>();

		for (String header : headers)
		{
			if (header.startsWith("##FORMAT"))
			{
				Tuple tuple = new VcfHeaderParser(header).parse();
				VcfFormat format = new VcfFormat(tuple);
				formats.add(format);
			}
		}

		return formats;
	}

	public List<VcfSample> getSamples()
	{
		List<VcfSample> samples = new ArrayList<VcfSample>();

		for (String header : headers)
		{
			if (header.startsWith("##SAMPLE"))
			{
				Tuple tuple = new VcfHeaderParser(header).parse();
				samples.add(new VcfSample(tuple));
			}
		}

		return samples;
	}

	public List<VcfContig> getContigs()
	{
		List<VcfContig> contigs = new ArrayList<VcfContig>();

		for (String header : headers)
		{
			if (header.startsWith("##contig"))
			{
				Tuple tuple = new VcfHeaderParser(header).parse();
				VcfContig contig = new VcfContig(tuple);
				contigs.add(contig);
			}
		}

		return contigs;
	}

	public List<VcfAlt> getAlts()
	{
		List<VcfAlt> alts = new ArrayList<VcfAlt>();

		for (String header : headers)
		{
			if (header.startsWith("##ALT"))
			{
				Tuple tuple = new VcfHeaderParser(header).parse();
				alts.add(new VcfAlt(tuple));
			}
		}

		return alts;
	}

	public Iterator<VcfRecord> recordIterator()
	{
		return new VcfRecordIterator(iterator());
	}

	@Override
	public void close() throws IOException
	{
		csvReader.close();
	}

	@Override
	public Iterator<Tuple> iterator()
	{
		return csvReader.iterator();
	}

	@Override
	public boolean hasColNames()
	{
		return true;
	}

	@Override
	public Iterator<String> colNamesIterator() throws IOException
	{
		return csvReader.colNamesIterator();
	}

	@Override
	public void addCellProcessor(CellProcessor cellProcessor)
	{
		csvReader.addCellProcessor(cellProcessor);
	}

	private List<String> readHeaders(BufferedReader br) throws IOException
	{
		List<String> headers = new ArrayList<String>();

		String header = br.readLine();
		while ((header != null) && header.startsWith("##"))
		{
			headers.add(header);
			header = br.readLine();
		}

		return headers;
	}

	private static class VcfRecordIterator implements Iterator<VcfRecord>
	{
		private final Iterator<Tuple> tupleIterator;

		public VcfRecordIterator(Iterator<Tuple> tupleIterator)
		{
			this.tupleIterator = tupleIterator;
		}

		@Override
		public boolean hasNext()
		{
			return tupleIterator.hasNext();
		}

		@Override
		public VcfRecord next()
		{
			Tuple tuple = tupleIterator.next();
			return new VcfRecord(tuple);
		}

		@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}
	}

}
