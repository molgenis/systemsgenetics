package umcg.genetica.genomicboundaries;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author PatrickDeelen
 */
public class GenomicBoundariesSimple {

		/**
	 * May <b>not</b> contain boundaries that falls completely within an other.
	 */
    private HashMap<String, ArrayList<GenomicBoundary>> genomicsBoundaries;
    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final Pattern CHR_PATTERN = Pattern.compile("^chr(.*)$", Pattern.CASE_INSENSITIVE);
    private int margin = 0;

    public GenomicBoundariesSimple(String genomicsBoundarysFilePath) throws FileNotFoundException, IOException {

        genomicsBoundaries = new HashMap<String, ArrayList<GenomicBoundary>>();
        BufferedReader genomicsBoundarysReader = new BufferedReader(new FileReader(genomicsBoundarysFilePath));
        String line;


        while( ( line = genomicsBoundarysReader.readLine() ) != null ){
            String[] cells = TAB_PATTERN.split(line);
			String chromosme = removeChr(cells[0]);
			Integer beginPoint = Integer.valueOf(cells[1]);
			int endPoint = Integer.parseInt(cells[2]);

			ArrayList<GenomicBoundary> chromosomeBoundaries;
			if(!genomicsBoundaries.containsKey(chromosme)){
				chromosomeBoundaries = new ArrayList<GenomicBoundary>();
				genomicsBoundaries.put(chromosme, chromosomeBoundaries);
			} else {
				chromosomeBoundaries = genomicsBoundaries.get(chromosme);
			}

            GenomicBoundary genomicBoundary = new GenomicBoundary(chromosme, beginPoint, endPoint);

			chromosomeBoundaries.add(genomicBoundary);


        }

    }

    public GenomicBoundariesSimple(String genomicsBoundarysFilePath, int margin) throws FileNotFoundException, IOException {
        this(genomicsBoundarysFilePath);
        this.margin = margin;
    }

    public boolean isChromosomeInBoundary(String chromosome){
        return genomicsBoundaries.containsKey(chromosome);
    }

    public Set<String> getChromosomes(){
        return genomicsBoundaries.keySet();
    }

    public boolean isInBoundary(String chromosome, int position){

        return getBoundary(chromosome, position) != null;

    }

	public GenomicBoundary getBoundary(String chromosome, int position){

        if(isChromosomeInBoundary(chromosome)){

			for(GenomicBoundary boundary : genomicsBoundaries.get(chromosome)){
				if(boundary.isInBoundarie(position, margin)){
					return boundary;
				}
			}

			return null;

        }

        return null;

    }

	private static String removeChr(String chromosome){

		Matcher chrMatcher = CHR_PATTERN.matcher(chromosome);
		if(chrMatcher.find()){
			return chrMatcher.group(1);
		} else {
			return chromosome;
		}

	}

}
