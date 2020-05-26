
package umcg.genetica.genomicboundaries;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author PatrickDeelen
 */
public class GenomicBoundaries<V> implements Iterable<GenomicBoundary<V>>{

	/**
	 * May <b>not</b> contain boundaries that falls completely within an other.
	 */
    private HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> genomicsBoundaries;
    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final Pattern COMMENT_PATTERN = Pattern.compile("^#");
	private static final Pattern CHR_PATTERN = Pattern.compile("^chr(.*)$", Pattern.CASE_INSENSITIVE);
    private int margin = 0;
	private boolean removedSubBoundaries = false;
	private String feature;
	
	/**
	 * Constructor for building Genomic Boundaries from a file.
	 * @param genomicsBoundarysFilePath
	 * @throws FileNotFoundException
	 * @throws IOException 
	 */
    public GenomicBoundaries(String genomicsBoundarysFilePath) throws FileNotFoundException, IOException {
       
        genomicsBoundaries = new HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>>();
        BufferedReader genomicsBoundarysReader = new BufferedReader(new FileReader(genomicsBoundarysFilePath));
        String line;


        while( ( line = genomicsBoundarysReader.readLine() ) != null ){
			
			if(COMMENT_PATTERN.matcher(line).find()){
				continue;
			}
			
            String[] cells = TAB_PATTERN.split(line);
			String chromosome = removeChr(cells[0]);
			Integer beginPoint = Integer.valueOf(cells[1]);
			int endPoint = Integer.parseInt(cells[2]);

			TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chromosomeBoundaries;
			if(!genomicsBoundaries.containsKey(chromosome)){
				chromosomeBoundaries = new TreeMap<Integer, ArrayList<GenomicBoundary<V>>>();
				genomicsBoundaries.put(chromosome, chromosomeBoundaries);
			} else {
				chromosomeBoundaries = genomicsBoundaries.get(chromosome);
			}

            GenomicBoundary genomicBoundary;
            

			genomicBoundary = new GenomicBoundary<V>(chromosome, beginPoint, endPoint);
			
			int n = chromosomeBoundaries.size();
			if(!chromosomeBoundaries.containsKey(beginPoint) || chromosomeBoundaries.get(beginPoint).get(n-1).getStop() < endPoint){
				
				ArrayList<GenomicBoundary<V>> boundaries;
				if(chromosomeBoundaries.containsKey(beginPoint)){
					boundaries = chromosomeBoundaries.get(beginPoint);
				} else {
					boundaries = new ArrayList<GenomicBoundary<V>>();
					chromosomeBoundaries.put(beginPoint, boundaries);
				}
				boundaries.add(genomicBoundary);
				
				//chromosomeBoundaries.put(beginPoint, genomicBoundary);
            }

			
        }
		
		
		removeSubBoundaries();
		
		
    }
	
	/**
	 * Empty constructor that only creates an empty structure for saving Genomic Boundaries.
	 */
	public GenomicBoundaries() {

		genomicsBoundaries = new HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>>();
		
	}
	
	/**
	 * Method to add a Genomic Boundary to the Genmoic Boundaries set.
	 * @param chromosome
	 * @param beginPoint
	 * @param endPoint 
	 */
	public final void addBoundary(String chromosome, Integer beginPoint, int endPoint){
		
		chromosome = removeChr(chromosome);
		
		TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chromosomeBoundaries;
		if(!genomicsBoundaries.containsKey(chromosome)){
			chromosomeBoundaries = new TreeMap<Integer, ArrayList<GenomicBoundary<V>>>();
			genomicsBoundaries.put(chromosome, chromosomeBoundaries);
		} else {
			chromosomeBoundaries = genomicsBoundaries.get(chromosome);
		}

		GenomicBoundary genomicBoundary;

		genomicBoundary = new GenomicBoundary(chromosome, beginPoint, endPoint);
		
		//Test if start point already contains boundary. If not so save new boundary. If start site contains 

		ArrayList<GenomicBoundary<V>> boundaries;
		if(chromosomeBoundaries.containsKey(beginPoint)){
			boundaries = chromosomeBoundaries.get(beginPoint);
		} else {
			boundaries = new ArrayList<GenomicBoundary<V>>();
			chromosomeBoundaries.put(beginPoint, boundaries);
		}
		boundaries.add(genomicBoundary);

		//chromosomeBoundaries.put(beginPoint, genomicBoundary);
		
		
		removedSubBoundaries = false;
		
	}
	
	
	/**
	 * Method to remove boundaries that are sub boundaries of other boundaries.
	 */
	public final void removeSubBoundaries(){
		// remove all the boundaries that fall within a boundary
		GenomicBoundary previousGenomicBoundary = null;
		GenomicBoundary genomicBoundary;
		
		Iterator<GenomicBoundary<V>> gboit = this.iterator();
		while(gboit.hasNext()){
			genomicBoundary = gboit.next();
			
			//System.out.println("cur: " + genomicBoundary.getStart() + "-" + genomicBoundary.getStop());
			//System.out.println("prev:" + previousGenomicBoundary.getStart() + "-" + previousGenomicBoundary.getStop());
			
			if(genomicBoundary.isPartOfBoundary(previousGenomicBoundary)){
				//System.out.println("REMOVING");
				gboit.remove();
				//genomicBoundaryList.remove(genomicBoundary);
			}

			previousGenomicBoundary = genomicBoundary;
		}
		
		removedSubBoundaries = true;
	}
	
	
	/**
	 * Method that merges overlapping Boundaries.
	 */
	public final void mergeOverlappingBoundaries(){
	
		removeSubBoundaries();
		
		HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> genomicsBoundariesMergedOverlappingBoundaries = new HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>>();
		
		for( Entry<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> chromosomeBoundariesEntry : genomicsBoundaries.entrySet()){
			
			TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chromosomeBoundaries = chromosomeBoundariesEntry.getValue();
			
			TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chromosomeBoundariesMergedOverlappingBoundaries = new TreeMap<Integer, ArrayList<GenomicBoundary<V>>>();
			
			GenomicBoundary<V> currentGenomicBoundary = null;
			
			chrBoundaries:
			for(ArrayList<GenomicBoundary<V>> genomicBoundaryList : chromosomeBoundaries.values()){
				
				for(GenomicBoundary<V> genomicBoundary : genomicBoundaryList){
					if(currentGenomicBoundary == null){
						currentGenomicBoundary = genomicBoundary;
					}

					if(!currentGenomicBoundary.isOverlaping(genomicBoundary)){
						//Check of al een ArrayList bestaat op die positie? Nee: maak een nieuwe lijst aan, Ja: voeg enkel toe.
						
						if(!chromosomeBoundariesMergedOverlappingBoundaries.containsKey(currentGenomicBoundary.getStart())){
							ArrayList<GenomicBoundary<V>> boundaries = new ArrayList<GenomicBoundary<V>>();
							boundaries.add(currentGenomicBoundary);
							chromosomeBoundariesMergedOverlappingBoundaries.put(currentGenomicBoundary.getStart(), boundaries);
						}
						
						else{
							ArrayList<GenomicBoundary<V>> boundaries = chromosomeBoundariesMergedOverlappingBoundaries.get(currentGenomicBoundary.getStart());
							boundaries.add(currentGenomicBoundary);
						}
						
						
						//chromosomeBoundariesMergedOverlappingBoundaries.put(currentGenomicBoundary.getStart(), currentGenomicBoundary);
						currentGenomicBoundary = genomicBoundary;
					} else {
						currentGenomicBoundary = new GenomicBoundary<V>(currentGenomicBoundary.getChromosome(), Math.min(currentGenomicBoundary.getStart(), genomicBoundary.getStart()), Math.max(currentGenomicBoundary.getStop(), genomicBoundary.getStop()));
					}
				}
				
			}
			
			
			if(!chromosomeBoundariesMergedOverlappingBoundaries.containsKey(currentGenomicBoundary.getStart())){
				ArrayList<GenomicBoundary<V>> boundaries = new ArrayList<GenomicBoundary<V>>();
				boundaries.add(currentGenomicBoundary);
				chromosomeBoundariesMergedOverlappingBoundaries.put(currentGenomicBoundary.getStart(), boundaries);
			}

			else{
				ArrayList<GenomicBoundary<V>> boundaries = chromosomeBoundariesMergedOverlappingBoundaries.get(currentGenomicBoundary.getStart());
				boundaries.add(currentGenomicBoundary);
			}
			
			
			//chromosomeBoundariesMergedOverlappingBoundaries.put(currentGenomicBoundary.getStart(), currentGenomicBoundary);
			
			genomicsBoundariesMergedOverlappingBoundaries.put(chromosomeBoundariesEntry.getKey(), chromosomeBoundariesMergedOverlappingBoundaries);
			
		}
		
		genomicsBoundaries = genomicsBoundariesMergedOverlappingBoundaries;
		
	}
	
	
    public GenomicBoundaries(String genomicsBoundarysFilePath, int margin) throws FileNotFoundException, IOException {
        this(genomicsBoundarysFilePath);
        this.margin = margin;
    }
	
	/**
	 * Method to check if a chromosome is contained in the Genomic Boundaries.
	 * @param chromosome
	 * @return 
	 */
    public boolean isChromosomeInBoundary(String chromosome){
        return genomicsBoundaries.containsKey(chromosome);
    }
	
	/**
	 * 
	 * @return 
	 */
    public Set<String> getChromosomes(){
        return genomicsBoundaries.keySet();
    }
	
	
    public boolean isInBoundary(String chromosome, int position, int listIndex){

        return getBoundary(chromosome, position, listIndex) != null;
        
    }
	
	/**
	 * This method returns a Boundary object from a chromosome at a specific position.
	 * 
	 * @param chromosome
	 * @param position
	 * @return 
	 */
	public GenomicBoundary getBoundary(String chromosome, int position, int listIndex){

		/*
		if(!removedSubBoundaries){
			throw new IllegalStateException("Can not get boundary is sub boundaries not removed.");
		}
		*/
		
        if(isChromosomeInBoundary(chromosome)){

			TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chromosomeBoudaries = genomicsBoundaries.get(chromosome);
			//if(position >= (chromosomeBoudaries.firstKey() - margin) && position <= (chromosomeBoudaries.lastKey() + margin)){
			Entry<Integer, ArrayList<GenomicBoundary<V>>> boundaryEntry = chromosomeBoudaries.floorEntry(position + margin);

			if(boundaryEntry != null){
				ArrayList<GenomicBoundary<V>> boundary = boundaryEntry.getValue();
				if (boundary.get(listIndex).isInBoundarie(position,margin)){
					return boundary.get(listIndex);
				}
			}
			

			

        }

        return null;

    }
	
	/**
	 * This method removes a chromosome?
	 * @param chromosome
	 * @return 
	 */
	private static String removeChr(String chromosome){

		Matcher chrMatcher = CHR_PATTERN.matcher(chromosome);
		if(chrMatcher.find()){
			return chrMatcher.group(1);
		} else {
			return chromosome;
		}

	}

//	public Iterator<GenomicBoundary> iterator() {
//		return new GenomicBoundariesIterator(genomicsBoundaries);
//	}

	public Iterator<GenomicBoundary<V>> iterator() {
		return new GenomicBoundariesIterator<V>(genomicsBoundaries);
	}

	
	public int getBoundaryCountChromosome(String chromosome){
		if(isChromosomeInBoundary(chromosome)){
			return genomicsBoundaries.get(chromosome).size();
		} else  {
			return 0;
		}
	}
	
	/**
	 * This method returns the count of Boundary objects.
	 * @return 
	 */
	public int getBoudaryCount(){
		
		int count = 0;
		
		for(TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chrBoundaries : genomicsBoundaries.values()){
			count += chrBoundaries.size();
		}
		
		return count;
	}

	
	/**
	 * Write the data of the Boundary objects to a BED file, specified by @param bedFilePath.
	 * @param bedFilePath
	 * @param trackName
	 * @throws IOException 
	 */
	public void writeBoundariesToBedFile(String bedFilePath, String trackName) throws IOException{
		BufferedWriter bedBufferWriter = new BufferedWriter(new FileWriter(bedFilePath));
		bedBufferWriter.append("track name=\"" + trackName + "\"\n");
		
		for(GenomicBoundary<V> genomicBoundary : this){
				bedBufferWriter.append("chr" + genomicBoundary.getChromosome() + "\t" + genomicBoundary.getStart() + "\t" + genomicBoundary.getStop() + "\n");
		}
		bedBufferWriter.close();
	}
	
	
	/**
	 * Returns a TreeMap<Integer, GenomicBoundary<V>> if the provided chromosome is present with genomicsBoundaries.
	 * Otherwise an exception is thrown.
	 * 
	 * @param chromosome
	 * @return
	 * @throws Exception 
	 */
	public TreeMap<Integer, ArrayList<GenomicBoundary<V>>> getGenomicBoundariesMap(String chromosome) throws Exception{
		
		if(!this.genomicsBoundaries.containsKey(chromosome)){
			throw new Exception("Chromosome: " + chromosome + " doesn't exist.");
		}
		return this.genomicsBoundaries.get(chromosome);
	}
	
	
	/**
	 * Returns the GenomicBoundaries
	 * @param chromosome
	 * @return
	 * @throws Exception 
	 */
	public Collection<ArrayList<GenomicBoundary<V>>> getGenomicBoundaries(String chromosome) throws Exception{
		
		return this.getGenomicBoundariesMap(chromosome).values();
		
	}
	
	
	
	/**
	 * Merges GenomicBoundaries by adding 
	 * @param otherGenomicBoundariesSet 
	 */
	public void mergeGenomicBoundaries(GenomicBoundaries otherGenomicBoundariesSet){
		
		HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> otherGenomicBoundaries = otherGenomicBoundariesSet.getGenomicBoundaries();
		
		
		for(Entry<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> chromosomeBoundariesEntry : otherGenomicBoundaries.entrySet()){
			TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chromosomeBoundaries = chromosomeBoundariesEntry.getValue();
			
			for(ArrayList<GenomicBoundary<V>> genomicBoundaryList : chromosomeBoundaries.values()){
				
				for(GenomicBoundary<V> genomicBoundary : genomicBoundaryList){
				
				//Check hier of het chromosoom nummer al bestaat. Zo niet, maak dan een nieuwe entry aan en voeg de Boundary toe. Bestaat het chromosoom al, voeg dan enkel de Boundary toe.
				if(!this.genomicsBoundaries.containsKey(genomicBoundary.getChromosome())){
					ArrayList<GenomicBoundary<V>> aap = new ArrayList<GenomicBoundary<V>>();
					aap.add(genomicBoundary);
					this.genomicsBoundaries.get(genomicBoundary.getChromosome()).put(genomicBoundary.getStart(), aap);
					//this.genomicsBoundaries.put(genomicBoundary.getChromosome(), new TreeMap<Integer, ArrayList<GenomicBoundary<V>>>());
				}
				
				else{
					this.genomicsBoundaries.get(genomicBoundary.getChromosome()).get(genomicBoundary.getStart()).add(genomicBoundary);
				}
				//this.genomicsBoundaries.get(genomicBoundary.getChromosome()).put(genomicBoundary.getStart(), genomicBoundary);
				}
			}
		}
		
	}
	
	
	/**
	 * Returns the GenomicBoundaries
	 * @return 
	 */
	public HashMap<String, TreeMap<Integer, ArrayList<GenomicBoundary<V>>>> getGenomicBoundaries(){
		return this.genomicsBoundaries;
	}
	
	
	
	/**
	 * Method to add a Genomic Boundary to the Genmoic Boundaries set.
	 * @param chromosome
	 * @param beginPoint
	 * @param endPoint 
	 */
	public final void addBoundary(String chromosome, Integer beginPoint, int endPoint, V annotation){
		
		chromosome = removeChr(chromosome);
		
		TreeMap<Integer, ArrayList<GenomicBoundary<V>>> chromosomeBoundaries;
		if(!genomicsBoundaries.containsKey(chromosome)){
			chromosomeBoundaries = new TreeMap<Integer, ArrayList<GenomicBoundary<V>>>();
			genomicsBoundaries.put(chromosome, chromosomeBoundaries);
		} else {
			chromosomeBoundaries = genomicsBoundaries.get(chromosome);
		}

		GenomicBoundary genomicBoundary;

		genomicBoundary = new GenomicBoundary(chromosome, beginPoint, endPoint, annotation);
		
		ArrayList<GenomicBoundary<V>> boundaries;
		if(chromosomeBoundaries.containsKey(beginPoint)){
			boundaries = chromosomeBoundaries.get(beginPoint);
		} else {
			boundaries = new ArrayList<GenomicBoundary<V>>();
			chromosomeBoundaries.put(beginPoint, boundaries);
		}
		boundaries.add(genomicBoundary);
		
		removedSubBoundaries = false;
		
	}
	
}
