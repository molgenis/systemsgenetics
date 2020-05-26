package umcg.genetica.io.bed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;

/**
 *
 * @author Patrick Deelen
 */
public class BedFile implements Iterable<BedEntry> {

	private final File bedFile;
	private final boolean omitChr;
	private final boolean makeOneBased;
	private final Map<String, String> trackInfo;
	private static final Pattern CHR_PATTERN = Pattern.compile("^chr(.*)$", Pattern.CASE_INSENSITIVE);
	private static final Pattern TRACK_INFO_PATTERN = Pattern.compile("(\\w+)=([^\"]\\S+|\".+?\")\\s*");

	public BedFile(String bedFilePath) throws FileNotFoundException, IOException {
		this(new File(bedFilePath), false, false);
	}

	public BedFile(String bedFilePath, boolean omitChr, boolean makeOneBased) throws FileNotFoundException, IOException {
		this(new File(bedFilePath), omitChr, makeOneBased);
	}

	public BedFile(File bedFile) throws FileNotFoundException, IOException {
		this(bedFile, false, false);
	}

	public BedFile(File bedFile, boolean omitChr, boolean makeOneBased) throws FileNotFoundException, IOException {
		this.bedFile = bedFile;
		this.omitChr = omitChr;
		this.makeOneBased = makeOneBased;


		if (!this.bedFile.exists()) {
			throw new FileNotFoundException("Bed file not found at: " + bedFile.getAbsolutePath());
		} else if (!this.bedFile.isFile()) {
			throw new IOException("Error reading bed file at: " + bedFile.getAbsolutePath());
		} else if (!this.bedFile.canRead()) {
			throw new IOException("Error reading bed file at: " + bedFile.getAbsolutePath());
		}

		Map<String, String> trackInfo = Collections.EMPTY_MAP;
		final BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile), "UTF-8"));
		String line;
		while ((line = reader.readLine()) != null) {
			if (line.startsWith("track")) {
				trackInfo = parseTrackLine(line);
				break;
			} else if (StringUtils.split(line, '\t').length > 3) {
				break;
			}
		}
		reader.close();

		this.trackInfo = trackInfo;

	}

	private static BedEntry parseLine(String line, boolean omitChr, boolean makeOneBased) throws IOException {
		String[] lineElements = StringUtils.split(line, '\t');

		if (lineElements.length < 3) {
			throw new IOException("Error parsing bed, did not find 3 fields on line: " + line);
		}

		final String chr;
		if (omitChr) {
			chr = removeChr(lineElements[0]).intern();
		} else {
			chr = lineElements[0].intern();
		}


		int start;
		try {
			start = Integer.parseInt(lineElements[1]);
		} catch (NumberFormatException ex) {
			throw new IOException("Error parsing bed file, Start is not an int on line: " + line);
		}

		int stop;
		try {
			stop = Integer.parseInt(lineElements[2]);
		} catch (NumberFormatException ex) {
			throw new IOException("Error parsing bed file, Stop is not an int on line: " + line);
		}


		if (makeOneBased) {
			start = start + 1;
			stop = stop + 1;
		}

		final String name;
		if (lineElements.length >= 4) {
			name = lineElements[3];
		} else {
			name = null;
		}

		final double score;
		if (lineElements.length >= 5) {
			try {
				score = Double.parseDouble(lineElements[4]);
			} catch (NumberFormatException ex) {
				throw new IOException("Error parsing bed file, score is not a double: " + line);
			}
		} else {
			score = Double.NaN;
		}

		return new BedEntry(chr, start, stop, name, score);
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
	public Iterator<BedEntry> iterator() {
		try {

			return new Iterator<BedEntry>() {
				private final BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile), "UTF-8"));
				private BedEntry next;
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
				public BedEntry next() {
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

	public PerChrIntervalTree<BedEntry> createIntervalTree() throws Exception {

		return PerChrIntervalTree.createFromChrGroupedIterable(this, BedEntry.class);

	}

	public Map<String, String> getTrackInfo() {
		return trackInfo;
	}
	
	

	private static Map<String, String> parseTrackLine(String line) {

		HashMap<String, String> infoMap = new HashMap<String, String>();
		
		Matcher m = TRACK_INFO_PATTERN.matcher(line);
		while (m.find()) {
				String infoName = m.group(1);
				String info = m.group(2).replace("\"", "");
				infoMap.put(infoName, info);				
		}
		
		return Collections.unmodifiableMap(infoMap);
		
	}
}
