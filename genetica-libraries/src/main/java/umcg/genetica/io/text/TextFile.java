/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.text;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.text.Strings;

/**
 * @author harmjan
 */
public class TextFile implements Iterable<String> {
// check 1,2,3
	
	public static final int DEFAULT_BUFFER_SIZE = 1024 * 64;
	public static final Pattern tab = Pattern.compile("\\t");//Using the \t from string.tab is technically not valid. I would not want to depend on this
	public static final Pattern space = Strings.space;
	public static final Pattern colon = Strings.colon;
	public static final Pattern semicolon = Strings.semicolon;
	public static final Pattern comma = Strings.comma;
	protected BufferedReader in;
	protected File file;
	public static final boolean W = true;
	public static final boolean R = false;
	protected BufferedWriter out;
	protected boolean writeable;
	protected static final String ENCODING = "ISO-8859-1";
	private boolean gzipped;
	private int buffersize;
	
	public TextFile(String file, boolean mode) throws IOException {
		this(new File(file), mode, DEFAULT_BUFFER_SIZE);
	}
	
	public TextFile(File file, boolean mode) throws IOException {
		this(file, mode, DEFAULT_BUFFER_SIZE);
	}
	
	public TextFile(File file, boolean mode, int buffersize) throws IOException {
		
		this.buffersize = buffersize;
		this.file = file;
		
		String loc = file.getAbsolutePath();
		if (loc.trim().length() == 0) {
			throw new IOException("Could not find file: no file specified");
		}
		this.writeable = mode;
		
		if (loc.endsWith(".gz")) {
			gzipped = true;
		}
		open();
	}
	
	public TextFile(String file, boolean mode, int buffersize) throws IOException {
		this(new File(file), mode, buffersize);
	}
	
	public final void open() throws IOException {
		
		if (!file.exists() && !writeable) {
			throw new IOException("Could not find file: " + file);
		} else {
			if (writeable) {
				if (gzipped) {
//					ParallelGZIPOutputStream gzipOutputStream = new ParallelGZIPOutputStream(new FileOutputStream(file), buffersize));
					GZIPOutputStream gzipOutputStream = new GZIPOutputStream(new FileOutputStream(file), buffersize){{def.setLevel(Deflater.BEST_SPEED);}};
					out = new BufferedWriter(new OutputStreamWriter(gzipOutputStream), buffersize);
					
				} else {
					out = new BufferedWriter(new FileWriter(file), buffersize);
				}
			} else {
				if (gzipped) {
					GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(file), buffersize);
					in = new BufferedReader(new InputStreamReader(gzipInputStream, "US-ASCII"));
				} else {
//                System.out.println("Opening file: "+file);
					in = new BufferedReader(new InputStreamReader(new FileInputStream(file), ENCODING), buffersize);
				}
			}
		}
	}
	
	public String readLine() throws IOException {
		return in.readLine();
	}
	
	public void write(String line) throws IOException {
		out.write(line);
	}
	
	public void close() throws IOException {
//        System.out.println("Closing "+file);
		if (writeable) {
			out.close();
		} else {
			in.close();
		}
	}
	
	/**
	 * This method is a wrapper for readLineElemsReturnReference: this method
	 * returns default substrings delimited by Pattern p.
	 *
	 * @param p The Pattern object to split with (e.g. TextFile.tab or
	 *          Strings.comma)
	 * @return New String objects for each substring delimited by Pattern p
	 * @throws IOException
	 */
	public String[] readLineElems(Pattern p) throws IOException {
		return readLineElemsReturnReference(p);
	}
	
	private Iterator<String[]> readLineElemsIterator(Pattern p) throws IOException {
		return new TextFileIteratorElements(this, p);
	}
	
	public Iterable<String[]> readLineElemsIterable(Pattern p) {
		return new TextFileIterableElements(this, p);
	}
	
	/**
	 * This method returns default substrings delimited by Pattern p. As such,
	 * this method may be more memory-efficient in some situations (for example
	 * when only a multiple columns should be loaded and stored in memory).
	 *
	 * @param p The Pattern object to split with (e.g. TextFile.tab or
	 *          Strings.comma)
	 * @return New String objects for each substring delimited by Pattern p
	 * @throws IOException
	 */
	public String[] readLineElemsReturnReference(Pattern p) throws IOException {
		if (in != null) {
			String ln = readLine();
			if (ln != null) {
				String[] elems = p.split(ln);
				ln = null;
				return elems;
			} else {
				return null;
			}
		} else {
			return null;
		}
	}
	
	/**
	 * This method returns a new object for each of the splitted elements,
	 * instead of the default action, in which each substring is backed by the
	 * original full-length String. As such, this method may be more
	 * memory-efficient in some situations (for example when only a single
	 * column should be loaded and stored in memory).
	 *
	 * @param p The Pattern object to split with (e.g. TextFile.tab or
	 *          Strings.comma)
	 * @return New String objects for each substring delimited by Pattern p
	 * @throws IOException
	 */
	public String[] readLineElemsReturnObjects(Pattern p) throws IOException {
		if (in != null) {
			String ln = readLine();
			if (ln != null) {
				String[] origelems = p.split(ln);
				String[] returnelems = new String[origelems.length];
				for (int i = 0; i < origelems.length; i++) {
					returnelems[i] = new String(origelems[i]);
				}
				ln = null;
				return returnelems;
			} else {
				return null;
			}
		} else {
			return null;
		}
		
		
	}
	
	public int countLines() throws IOException {
		String ln = readLine();
		int ct = 0;
		while (ln != null) {
			if (ln.trim().length() > 0) {
				ct++;
			}
			ln = readLine();
		}
		close();
		open();
		return ct;
	}
	
	public int countCols(Pattern p) throws IOException {
		String ln = readLine();
		int ct = 0;
		if (ln != null) {
			String[] elems = p.split(ln);
			ct = elems.length;
		}
		close();
		open();
		return ct;
	}
	
	public String[] readAsArray() throws IOException {
		int numLines = countLines();
		String ln = readLine();
		if (ln == null) {
			return new String[0];
		}
		String[] data = new String[numLines];
		int i = 0;
		while (ln != null) {
			if (ln.trim().length() > 0) {
				data[i] = ln;
				i++;
			}
			ln = in.readLine();
		}
		return data;
	}
	
	public String[] readAsArray(int col, Pattern p) throws IOException {
		int numLines = countLines();
		String[] data = new String[numLines];
		int i = 0;
		String[] elems = readLineElems(p);
		while (elems != null) {
			if (elems.length > col) {
				data[i] = elems[col];
			}
			i++;
			elems = readLineElems(p);
		}
		return data;
	}
	
	public ArrayList<String> readAsArrayList() throws IOException {
		ArrayList<String> data = new ArrayList<String>();
		
		String ln = readLine();
		while (ln != null) {
			if (ln.trim().length() > 0) {
				data.add(ln);
			}
			ln = in.readLine();
		}
		return data;
	}
	
	public ArrayList<String> readAsArrayList(int col, Pattern p) throws IOException {
		ArrayList<String> data = new ArrayList<String>();
		
		String[] elems = readLineElems(p);
		while (elems != null) {
			if (elems.length > col) {
				data.add(elems[col]);
			}
			elems = readLineElems(p);
		}
		return data;
	}
	
	public void writeln(CharSequence csq) throws IOException {
		out.append(csq);
		out.append('\n');
	}
	
	public void writeln(String line) throws IOException {
		out.append(line);
		out.append('\n');
	}
	
	public void writeln() throws IOException {
		out.newLine();
	}
	
	public void append(char c) throws IOException {
		out.append(c);
	}
	
	public void append(CharSequence csq) throws IOException {
		out.append(csq);
	}
	
	public synchronized void writelnsynced(String str) throws IOException {
		this.writeln(str);
	}
	
	public void writelnTabDelimited(Object[] vals) throws IOException {
		String delim = "";
		for (Object val : vals) {
			out.write(delim);
			out.write(val.toString());
			delim = "\t";
		}
		writeln();
	}
	
	public void writelnDelimited(Object[] vals, Pattern p) throws IOException {
		String delim = "";
		for (Object val : vals) {
			out.write(delim);
			out.write(val.toString());
			delim = p.pattern();
		}
		writeln();
	}
	
	public Map<String, String> readAsHashMap(int col1, int col2) throws IOException {
		Map<String, String> output = new HashMap<String, String>();
		String[] elems = readLineElems(tab);
		while (elems != null) {
			if (elems.length > col1 && elems.length > col2) {
				output.put(elems[col1], elems[col2]);
			}
			elems = readLineElems(tab);
		}
		return output;
	}
	
	public Map<String, String> readAsHashMap(int col1, int[] col2, String seperator) throws IOException {
		Map<String, String> output = new HashMap<String, String>();
		
		int maxInt = col1;
		for (int i : col2) {
			if (i > maxInt) {
				maxInt = i;
			}
		}
		
		String[] elems = readLineElems(tab);
		while (elems != null) {
			if (elems.length > maxInt) {
				StringBuilder value = new StringBuilder(elems[col2[0]]);
				for (int i = 1; i < col2.length; ++i) {
					value.append(seperator);
					value.append(elems[col2[i]]);
				}
				output.put(elems[col1], value.toString());
			}
			elems = readLineElems(tab);
		}
		return output;
	}
	
	public Map<String, String> readAsHashMap(int col1, int col2, Pattern p) throws IOException {
		Map<String, String> output = new HashMap<String, String>();
		String[] elems = readLineElems(p);
		while (elems != null) {
			if (elems.length > col1 && elems.length > col2) {
				output.put(elems[col1], elems[col2]);
			}
			elems = readLineElems(p);
		}
		return output;
	}
	
	public Set<String> readAsSet(int col, Pattern p) throws IOException {
		Set<String> output = new HashSet<String>();
		String[] elems = readLineElems(p);
		while (elems != null) {
			if (elems.length > col) {
				output.add(elems[col]);
			}
			elems = readLineElems(p);
		}
		return output;
	}
	
	public void writeList(List l) throws IOException {
		for (Object e : l) {
			this.writeln(e.toString());
		}
	}
	
	public String getFileName() {
		return file.getAbsolutePath();
	}
	
	public HashSet<Pair<String, String>> readAsPairs(int A, int B) throws IOException {
		
		HashSet<Pair<String, String>> output = new HashSet<Pair<String, String>>();
		String[] elems = this.readLineElemsReturnObjects(tab);
		
		while (elems != null) {
			if (elems.length > A && elems.length > B) {
				output.add(new Pair(elems[A], elems[B]));
			}
			elems = this.readLineElemsReturnObjects(tab);
		}
		return output;
	}
	
	public void flush() throws IOException {
		out.flush();
	}
	
	public HashSet<Triple<String, String, String>> readAsTriple(int A, int B, int C) throws IOException {
		HashSet<Triple<String, String, String>> output = new HashSet<Triple<String, String, String>>();
		String[] elems = this.readLineElemsReturnObjects(tab);
		
		while (elems != null) {
			if (elems.length > A && elems.length > B && elems.length > C) {
				output.add(new Triple(elems[A], elems[B], elems[C]));
			}
			elems = this.readLineElemsReturnObjects(tab);
		}
		return output;
	}
	
	@Override
	public Iterator<String> iterator() {
		try {
			return new TextFileIterator(this);
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}
	
	private static class TextFileIterator implements Iterator<String> {
		
		private final TextFile textFile;
		String next;
		
		public TextFileIterator(TextFile textFile) throws IOException {
			this.textFile = textFile;
			next = textFile.readLine();
		}
		
		@Override
		public boolean hasNext() {
			return next != null;
		}
		
		@Override
		public String next() {
			String current = next;
			try {
				next = textFile.readLine();
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
			return current;
		}
		
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Not supported yet.");
		}
	}
	
	private static class TextFileIteratorElements implements Iterator<String[]> {
		
		private final TextFile textFile;
		private final Pattern pattern;
		String[] next;
		
		public TextFileIteratorElements(TextFile textFile, Pattern p) throws IOException {
			this.textFile = textFile;
			this.pattern = p;
			next = textFile.readLineElemsReturnObjects(this.pattern);
		}
		
		@Override
		public boolean hasNext() {
			return next != null;
		}
		
		@Override
		public String[] next() {
			String[] current = next;
			try {
				next = textFile.readLineElemsReturnObjects(this.pattern);
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
			return current;
		}
		
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Not supported yet.");
		}
	}
	
	private static class TextFileIterableElements implements Iterable<String[]> {
		
		private final TextFile textFile;
		private final Pattern pattern;
		
		public TextFileIterableElements(TextFile textFile, Pattern pattern) {
			this.textFile = textFile;
			this.pattern = pattern;
		}
		
		@Override
		public Iterator<String[]> iterator() {
			try {
				return textFile.readLineElemsIterator(pattern);
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
		}
		
	}
}
