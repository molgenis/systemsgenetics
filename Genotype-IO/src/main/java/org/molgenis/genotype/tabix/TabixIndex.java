package org.molgenis.genotype.tabix;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.molgenis.genotype.GenotypeDataIndex;
import org.molgenis.genotype.RawLineQuery;
import org.molgenis.genotype.VariantQuery;
import org.molgenis.genotype.variant.VariantLineMapper;

/**
 * Tabix implementation of the GenotypeDataIndex The constructor takes a gz.tbi index file and a bzipped file.
 * 
 * Code is copied from the TabixReader wich comes with Tabix when you download is see
 * http://sourceforge.net/projects/samtools/
 * 
 * Use creatQuery() to create a query to get a subset of the data
 * 
 * @author erwin
 * 
 */
@edu.umd.cs.findbugs.annotations.SuppressWarnings("RR_NOT_CHECKED")
public class TabixIndex implements GenotypeDataIndex
{
	private static final Charset CHARSET_UTF8 = Charset.forName("UTF-8");
	private static int MAX_BIN = 37450;
	private static int TAD_LIDX_SHIFT = 14;

	private String[] seqNames;// distinct #CHROM column

	private int mPreset;
	private int mSc;
	private int mBc;
	private int mEc;
	private int mMeta;
	private TIndex[] mIndex;
	private HashMap<String, Integer> mChr2tid;

	private final File bzipFile;
	private final VariantLineMapper variantLineMapper;

	public TabixIndex(File tabixIndexFile, File bzipFile, VariantLineMapper variantLineMapper) throws IOException
	{
		this.bzipFile = bzipFile;
		this.variantLineMapper = variantLineMapper;

		readIndexFile(tabixIndexFile);
	}

	@Override
	public List<String> getSeqNames()
	{
		return Collections.unmodifiableList(Arrays.asList(seqNames));
	}

	public TIndex[] getIndex()
	{
		return mIndex.clone();
	}

	@Override
	public VariantQuery createQuery()
	{
		return new TabixQuery(bzipFile, this, variantLineMapper);
	}

	@Override
	public RawLineQuery createRawLineQuery()
	{
		return new TabixRawLineQuery(bzipFile, this);
	}

	private int chr2tid(final String chr)
	{
		if (mChr2tid.containsKey(chr)) return mChr2tid.get(chr);
		else return -1;
	}

	private void readIndexFile(File tabixIndexFile) throws IOException
	{
		BlockCompressedInputStream bciStream = new BlockCompressedInputStream(tabixIndexFile);

		try
		{
			byte[] buf = new byte[4];
			bciStream.read(buf, 0, 4); // read "TBI\1"
			seqNames = new String[readInt(bciStream)]; // # sequences
			mChr2tid = new HashMap<String, Integer>();
			mPreset = readInt(bciStream);
			mSc = readInt(bciStream);
			mBc = readInt(bciStream);
			mEc = readInt(bciStream);
			mMeta = readInt(bciStream);
			readInt(bciStream);
			// read sequence dictionary
			int i, j, k, l = readInt(bciStream);
			buf = new byte[l];
			bciStream.read(buf);
			for (i = j = k = 0; i < buf.length; ++i)
			{
				if (buf[i] == 0)
				{
					byte[] b = new byte[i - j];
					System.arraycopy(buf, j, b, 0, b.length);
					String s = new String(b, CHARSET_UTF8);
					mChr2tid.put(s, k);
					seqNames[k++] = s;
					j = i + 1;
				}
			}
			// read the index
			mIndex = new TIndex[seqNames.length];
			for (i = 0; i < seqNames.length; ++i)
			{
				// the binning index
				int n_bin = readInt(bciStream);
				mIndex[i] = new TIndex();
				mIndex[i].b = new HashMap<Integer, TPair64[]>();
				for (j = 0; j < n_bin; ++j)
				{
					int bin = readInt(bciStream);
					TPair64[] chunks = new TPair64[readInt(bciStream)];
					for (k = 0; k < chunks.length; ++k)
					{
						long u = readLong(bciStream);
						long v = readLong(bciStream);
						chunks[k] = new TPair64(u, v); // in C, this is
														// inefficient
					}
					mIndex[i].b.put(bin, chunks);
				}
				// the linear index
				mIndex[i].l = new long[readInt(bciStream)];
				for (k = 0; k < mIndex[i].l.length; ++k)
					mIndex[i].l[k] = readLong(bciStream);
			}
		}
		finally
		{
			IOUtils.closeQuietly(bciStream);
		}
	}

	private int readInt(InputStream in) throws IOException
	{
		byte[] buf = new byte[4];
		in.read(buf);

		return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getInt();
	}

	private long readLong(InputStream in) throws IOException
	{
		byte[] buf = new byte[8];
		in.read(buf);

		return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getLong();
	}

	private static class TPair64 implements Comparable<TPair64>
	{
		long u, v;

		public TPair64(final long _u, final long _v)
		{
			u = _u;
			v = _v;
		}

		public TPair64(final TPair64 p)
		{
			u = p.u;
			v = p.v;
		}

		@Override
		public int compareTo(final TPair64 p)
		{
			return u == p.u ? 0 : ((u < p.u) ^ (u < 0) ^ (p.u < 0)) ? -1 : 1; // unsigned
																				// 64-bit
																				// comparison
		}

		@Override
		public boolean equals(Object obj)
		{
			if (!(obj instanceof TPair64))
			{
				return false;
			}

			return compareTo((TPair64) obj) == 0;
		}

		@Override
		public int hashCode()
		{
			return new HashCodeBuilder(3, 17).append(u).append(v).toHashCode();
		}

	};

	private static boolean less64(final long u, final long v)
	{ // unsigned 64-bit comparison
		return (u < v) ^ (u < 0) ^ (v < 0);
	}

	private static class TIndex
	{
		HashMap<Integer, TPair64[]> b; // binning index
		long[] l; // linear index
	}

	private int reg2bins(final int beg, final int _end, final int[] list)
	{
		int i = 0, k, end = _end;
		if (beg >= end) return 0;
		if (end >= 1 << 29) end = 1 << 29;
		--end;
		list[i++] = 0;
		for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k)
			list[i++] = k;
		for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k)
			list[i++] = k;
		for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k)
			list[i++] = k;
		for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k)
			list[i++] = k;
		for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k)
			list[i++] = k;
		return i;
	}

	public TabixIterator queryTabixIndex(String sequence, final int beg, final int end,
			BlockCompressedInputStream bzipInputStream) throws IOException
	{
		TPair64[] off, chunks;
		long min_off;
		int tid = chr2tid(sequence);

		if (tid == -1)
		{
			return null;
		}

		TIndex idx = mIndex[tid];
		int[] bins = new int[MAX_BIN];
		int i, l, n_off, n_bins = reg2bins(beg, end, bins);
		if (idx.l.length > 0) min_off = (beg >> TAD_LIDX_SHIFT >= idx.l.length) ? idx.l[idx.l.length - 1] : idx.l[beg >> TAD_LIDX_SHIFT];
		else min_off = 0;
		for (i = n_off = 0; i < n_bins; ++i)
		{
			if ((chunks = idx.b.get(bins[i])) != null) n_off += chunks.length;
		}
		if (n_off == 0) return null;
		off = new TPair64[n_off];
		for (i = n_off = 0; i < n_bins; ++i)
			if ((chunks = idx.b.get(bins[i])) != null) for (int j = 0; j < chunks.length; ++j)
				if (TabixIndex.less64(min_off, chunks[j].v)) off[n_off++] = new TPair64(chunks[j]);
		if (n_off == 0) return null;
		Arrays.sort(off, 0, n_off);
		// resolve completely contained adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i)
		{
			if (TabixIndex.less64(off[l].v, off[i].v))
			{
				++l;
				off[l].u = off[i].u;
				off[l].v = off[i].v;
			}
		}
		n_off = l + 1;
		// resolve overlaps between adjacent blocks; this may happen due to the
		// merge in indexing
		for (i = 1; i < n_off; ++i)
			if (!TabixIndex.less64(off[i - 1].v, off[i].u)) off[i - 1].v = off[i].u;
		// merge adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i)
		{
			if (off[l].v >> 16 == off[i].u >> 16) off[l].v = off[i].v;
			else
			{
				++l;
				off[l].u = off[i].u;
				off[l].v = off[i].v;
			}
		}
		n_off = l + 1;
		// return
		TPair64[] ret = new TPair64[n_off];
		for (i = 0; i < n_off; ++i)
			ret[i] = new TPair64(off[i].u, off[i].v); // in C, this is
														// inefficient
		return new TabixIterator(tid, beg, end, ret, bzipInputStream);
	}

	public class TabixIterator
	{
		private int i;
		private final int tid, beg, end;
		private final TPair64[] off;
		private long curr_off;
		private boolean iseof;
		private final BlockCompressedInputStream inputStream;

		public TabixIterator(final int _tid, final int _beg, final int _end, final TPair64[] _off,
				final BlockCompressedInputStream inputStream)
		{
			i = -1;
			curr_off = 0;
			iseof = false;
			off = _off.clone();
			tid = _tid;
			beg = _beg;
			end = _end;
			this.inputStream = inputStream;
		}

		public String next() throws IOException
		{
			if (iseof) return null;
			for (;;)
			{
				if (curr_off == 0 || !less64(curr_off, off[i].v))
				{ // then jump to the next chunk
					if (i == off.length - 1) break; // no more chunks
					if (i >= 0) assert (curr_off == off[i].v); // otherwise bug
					if (i < 0 || off[i].v != off[i + 1].u)
					{ // not adjacent chunks; then seek
						inputStream.seek(off[i + 1].u);
						curr_off = inputStream.getFilePointer();
					}
					++i;
				}
				String s;
				if ((s = inputStream.readLine()) != null)
				{
					TIntv intv;
					char[] str = s.toCharArray();
					curr_off = inputStream.getFilePointer();
					if (str.length == 0 || str[0] == mMeta) continue;
					intv = getIntv(s);
					if (intv.tid != tid || intv.beg >= end) break; // no need to
																	// proceed
					else if (intv.end > beg && intv.beg < end) return s; // overlap;
																			// return
				}
				else break; // end of file
			}
			iseof = true;
			return null;
		}

		private TIntv getIntv(final String s)
		{
			TIntv intv = new TIntv();
			int col = 0, end = 0, beg = 0;
			while ((end = s.indexOf('\t', beg)) >= 0 || end == -1)
			{
				++col;
				if (col == mSc)
				{
					intv.tid = chr2tid(s.substring(beg, end));
				}
				else if (col == mBc)
				{
					intv.beg = intv.end = Integer.parseInt(s.substring(beg, end == -1 ? s.length() : end));
					if ((mPreset & 0x10000) != 0) ++intv.end;
					else --intv.beg;
					if (intv.beg < 0) intv.beg = 0;
					if (intv.end < 1) intv.end = 1;
				}
				else
				{ // FIXME: SAM supports are not tested yet
					if ((mPreset & 0xffff) == 0)
					{ // generic
						if (col == mEc) intv.end = Integer.parseInt(s.substring(beg, end));
					}
					else if ((mPreset & 0xffff) == 1)
					{ // SAM
						if (col == 6)
						{ // CIGAR
							int l = 0, i, j;
							String cigar = s.substring(beg, end);
							for (i = j = 0; i < cigar.length(); ++i)
							{
								if (cigar.charAt(i) > '9')
								{
									int op = cigar.charAt(i);
									if (op == 'M' || op == 'D' || op == 'N') l += Integer.parseInt(cigar
											.substring(j, i));
								}
							}
							intv.end = intv.beg + l;
						}
					}
					else if ((mPreset & 0xffff) == 2)
					{ // VCF
						String alt;
						alt = end >= 0 ? s.substring(beg, end) : s.substring(beg);
						if (col == 4)
						{ // REF
							if (alt.length() > 0) intv.end = intv.beg + alt.length();
						}
						else if (col == 8)
						{ // INFO
							int e_off = -1, i = alt.indexOf("END=");
							if (i == 0) e_off = 4;
							else if (i > 0)
							{
								i = alt.indexOf(";END=");
								if (i >= 0) e_off = i + 5;
							}
							if (e_off > 0)
							{
								i = alt.indexOf(";", e_off);
								intv.end = Integer.parseInt(i > e_off ? alt.substring(e_off, i) : alt.substring(e_off));
							}
						}
					}
				}
				if (end == -1) break;
				beg = end + 1;
			}
			return intv;
		}
	};

	private static class TIntv
	{
		int tid, beg, end;
	}

}
