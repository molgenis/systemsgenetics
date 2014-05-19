package umcg.genetica.io.pileup;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.Iterator;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import org.molgenis.genotype.Allele;

/**
 *
 * @author Patrick Deelen
 */
public class PileupFile implements Iterable<PileupEntry> {

	File pileupFile;
	private static final Pattern TAB_PATTERN = Pattern.compile("\t");
	private final int minimumBaseQuality;

	public PileupFile(String pileupFilePath, int minimumBaseQuality) throws FileNotFoundException, IOException {
		this(new File(pileupFilePath), minimumBaseQuality);
	}
	
	public PileupFile(File pileupFile, int minimumBaseQuality) throws FileNotFoundException, IOException {
		this.pileupFile = pileupFile;
		this.minimumBaseQuality = minimumBaseQuality;

		if (!pileupFile.isFile()) {
			throw new FileNotFoundException("Cannot find pileup file at: " + pileupFile.getAbsolutePath());
		}
		if (!pileupFile.canRead()) {
			throw new IOException("Cannot read pileup file at: " + pileupFile.getAbsolutePath());
		}

	}
	
	public static PileupEntry parsePileupLine(String line, int minimumBaseQuality) throws PileupParseException{
		String[] elements = TAB_PATTERN.split(line);
		
		final int pos;
		try {
			pos = Integer.parseInt(elements[1]);
		} catch (NumberFormatException ex){
			throw new PileupParseException("Error parsing position: " + elements[1]);
		}
		final int depth;
		try {
			depth = Integer.parseInt(elements[3]);
		} catch (NumberFormatException ex){
			throw new PileupParseException("Error parsing position: " + elements[3]);
		}
		
		Allele allele = Allele.create(elements[2]);
		if(elements.length == 5){
			return new PileupEntry(elements[0], pos, allele, depth, elements[4], minimumBaseQuality);
		} else {
			return new PileupEntry(elements[0], pos, allele, depth, elements[4], elements[5], minimumBaseQuality);
		}
		
	}

	@Override
	public Iterator<PileupEntry> iterator() {

		BufferedReader reader = null;
		try {
			if(pileupFile.getAbsolutePath().endsWith(".gz")){
				reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pileupFile)), "UTF-8"));
			} else {
				reader = new BufferedReader(new InputStreamReader(new FileInputStream(pileupFile), "UTF-8"));
			}
			
			return new pileupFileIterator(reader);
		} catch (FileNotFoundException ex) {
			throw new RuntimeException(ex);
		} catch (UnsupportedEncodingException ex) {
			throw new RuntimeException(ex);
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		} 
	}
	
	public String getAbsolutePath(){
			return pileupFile.getAbsolutePath();
		}

	private class pileupFileIterator implements Iterator<PileupEntry> {

		private final BufferedReader reader;
		private PileupEntry next;

		public pileupFileIterator(BufferedReader reader) {
			this.reader = reader;
			try {
				String line = reader.readLine();
				next = line == null ? null : parsePileupLine(line, minimumBaseQuality);
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			} catch (PileupParseException ex) {
				throw new RuntimeException(ex);
			}
		}
				
		@Override
		public boolean hasNext() {
			return next != null;
		}

		@Override
		public PileupEntry next() {
			PileupEntry current = next;
			try {
				String line = reader.readLine();
				next = line == null ? null : parsePileupLine(line, minimumBaseQuality);
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			} catch (PileupParseException ex) {
				throw new RuntimeException(ex);
			}
			return current;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("Not supported.");
		}
	}
}
