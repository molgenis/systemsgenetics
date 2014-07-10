/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.gtf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;
import org.apache.log4j.Logger;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;

/**
 *
 * 
 * 
 * @author Patrick Deelen
 */
public class GtfReader implements GffReaderInterface{

	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final Pattern SPACE_PATTERN = Pattern.compile(" ");
	private static final Pattern SEMICOLON_PATTERN = Pattern.compile(";");
	private static final Logger LOGGER = Logger.getLogger(GtfReader.class);
	
	private final BufferedReader gtfFileBufferedReader;

	public GtfReader(File gftFile) throws FileNotFoundException{
		this(new FileReader(gftFile));
	}
	
	public GtfReader(Reader gftReader) {
		
		if(gftReader == null){
			throw new IllegalArgumentException("Supplied gft reader = null");
		}
		
		gtfFileBufferedReader = new BufferedReader(gftReader);
		
	}
	
	@Override
	public GffElement nextElement() throws IOException{
		
		
		String line;
		do {
			
			line = gtfFileBufferedReader.readLine();
		
			if(line == null){
				return null;
			}
			
		} while (line.charAt(0) == '#');
		

		
		String[] lineElements = TAB_PATTERN.split(line);
		
		if(lineElements.length != 9){
			throw new IOException("Line does not contain 9 tab sepparated columns: " + line);
		}
		
		String seqname = lineElements[0];
		String source = lineElements[1];
		String feature = lineElements[2]; 
		int start = Integer.parseInt(lineElements[3]);
		int end = Integer.parseInt(lineElements[4]);
		float score = lineElements[5].charAt(0) == '.' ? Float.NaN : Float.parseFloat(lineElements[5]);
		char strand = lineElements[6].charAt(0);
		int frame = lineElements[7].charAt(0) == '.' ? -1 : Integer.parseInt(lineElements[7]);
		String attributesInfo = lineElements[8];
		
		LinkedHashMap<String, String> attributes = new LinkedHashMap<String, String>();
		
		for(String attributeInfo : SEMICOLON_PATTERN.split(attributesInfo)){
			
			attributeInfo = attributeInfo.trim();
			String[] attributeInfoElements = SPACE_PATTERN.split(attributeInfo);
			if(attributeInfoElements.length != 2){
				throw new IOException("Line attributes format error part: " + attributeInfo + " of line: " + line);
			}
			
			String attributeName = attributeInfoElements[0];
			String attributeValue = attributeInfoElements[1];
			
			if(attributeValue.charAt(0) == '\"' && attributeValue.charAt(attributeValue.length()-1) == '\"'){
				attributeValue = attributeValue.substring(1, attributeValue.length()-1);
			}
			
			attributes.put(attributeName, attributeValue);
			
		}
		
		return new GffElement(seqname, source, feature, start, end, score, strand, frame, attributes);
		
	}
	
	@Override
	public void close() throws IOException{
		gtfFileBufferedReader.close();
	}

	@Override
	public Iterator<GffElement> iterator() {
		return new GffElementIterator(this);
	}
	
	public PerChrIntervalTree<GffElement> createIntervalTree() throws Exception{
		
		return PerChrIntervalTree.createFromChrGroupedIterable(this, GffElement.class);
				
	}
	
}
