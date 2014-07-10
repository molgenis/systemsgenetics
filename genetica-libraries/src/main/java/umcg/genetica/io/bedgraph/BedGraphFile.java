package umcg.genetica.io.bedgraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;
import umcg.genetica.io.gtf.GffElement;

/**
 *
 * @author Patrick Deelen
 */
public class BedGraphFile implements Iterable<BedGraphEntry> {

	private final File bedGraphFile;
	private final boolean omitChr;
	private final boolean makeOneBased;
	
	private static final Pattern CHR_PATTERN = Pattern.compile("^chr(.*)$", Pattern.CASE_INSENSITIVE);

	public BedGraphFile(String bedGraphFilePath) throws FileNotFoundException, IOException {
		this(new File(bedGraphFilePath), false, false);
	}

	public BedGraphFile(String bedGraphFilePath, boolean omitChr, boolean makeOneBased) throws FileNotFoundException, IOException {
		this(new File(bedGraphFilePath), omitChr, makeOneBased);
	}

	public BedGraphFile(File bedGraphFile) throws FileNotFoundException, IOException {
		this(bedGraphFile, false, false);
	}

	public BedGraphFile(File bedGraphFile, boolean omitChr, boolean makeOneBased) throws FileNotFoundException, IOException {
		this.bedGraphFile = bedGraphFile;
		this.omitChr = omitChr;
		this.makeOneBased = makeOneBased;


		if (!this.bedGraphFile.exists()) {
			throw new FileNotFoundException("BedGraph file not found at: " + bedGraphFile.getAbsolutePath());
		} else if (!this.bedGraphFile.isFile()) {
			throw new IOException("Error reading BedGraph file at: " + bedGraphFile.getAbsolutePath());
		} else if (!this.bedGraphFile.canRead()) {
			throw new IOException("Error reading BedGraph file at: " + bedGraphFile.getAbsolutePath());
		}
	}

	private static BedGraphEntry parseLine(String line, boolean omitChr, boolean makeOneBased) throws IOException {
		String[] lineElements = StringUtils.split(line);

		if (lineElements.length != 4) {
			throw new IOException("Error parsing BedGraph, did not find 4 fields on line: " + line);
		}
		
		final String chr;
		if(omitChr){
			chr = removeChr(lineElements[0]).intern();
		} else {
			chr = lineElements[0].intern();
		}
		

		int start;
		try {
			start = Integer.parseInt(lineElements[1]);
		} catch (NumberFormatException ex) {
			throw new IOException("Error parsing BedGraph, Start is not an int on line: " + line);
		}
		
		int stop;
		try {
			stop = Integer.parseInt(lineElements[2]);
		} catch (NumberFormatException ex) {
			throw new IOException("Error parsing BedGraph, Stop is not an int on line: " + line);
		}
		
		final double value;
		try {
			value = Double.parseDouble(lineElements[3]);
		} catch (NumberFormatException ex) {
			throw new IOException("Error parsing BedGraph, Value is not a double on line: " + line);
		}
		
		if (makeOneBased) {
			start = start + 1;
			stop = stop + 1;
		}
		
		
		
		return new BedGraphEntry(chr, start, stop, value);
	}
	
		/**
	 * This method removes a chromosome?
	 *
	 * @param chromosome
	 * @return
	 */
	private static String removeChr(String chromosome) {

		Matcher chrMatcher = CHR_PATTERN.matcher(chromosome);
		if (chrMatcher.find()) {
			return chrMatcher.group(1);
		} else {
			return chromosome;
		}

	}

	@Override
	public Iterator<BedGraphEntry> iterator() {
		try {

			return new Iterator<BedGraphEntry>() {
				private final BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(bedGraphFile), "UTF-8"));
				private BedGraphEntry next;
				private boolean atNext;
				private boolean atEnd = false;

				@Override
				public boolean hasNext() {

					if (atEnd) {
						return false;
					}

					if (atNext) {
						return true;
					}

					try {
						String line;
						while ((line = reader.readLine()) != null) {
							if (line.startsWith("browser") || line.startsWith("track") || line.charAt(0) == '#') {
								continue;
							} else {
								break;
							}
						}
						if (line == null) {
							atEnd = true;
							try {
								reader.close();
							} catch (IOException ex) {
							}
							return false;
						} else {
							next = parseLine(line, omitChr, makeOneBased);
							atNext = true;
							return true;
						}
					} catch (IOException ex) {
						throw new RuntimeException(ex);
					}

				}

				@Override
				public BedGraphEntry next() {
					if (!hasNext()) {
						throw new NoSuchElementException();
					}
					atNext = false;
					return next;
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException("Not supported yet.");
				}
			};
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
	}
	
	public PerChrIntervalTree<BedGraphEntry> createIntervalTree() throws Exception{
		
		return PerChrIntervalTree.createFromChrGroupedIterable(this, BedGraphEntry.class);
		
	}
}
