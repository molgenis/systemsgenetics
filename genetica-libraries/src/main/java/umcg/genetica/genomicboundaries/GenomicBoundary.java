package umcg.genetica.genomicboundaries;

/**
 *
 * @author PatrickDeelen
 */
public class GenomicBoundary<V> implements Comparable<GenomicBoundary>{

    private String chromosome;
    private Integer start;
    private int stop;
	private V annotation;

    public GenomicBoundary(String chromosome, Integer start, int stop) {
        this.chromosome = chromosome;
        this.start = start;
        this.stop = stop;
    }
	
	public GenomicBoundary(String chromosome, Integer start, int stop, V annotation) {
        this.chromosome = chromosome;
        this.start = start;
        this.stop = stop;
		this.annotation = annotation;
    }

    public String getChromosome() {
        return chromosome;
    }

    public Integer getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    @Override
    public int compareTo(GenomicBoundary other) {
        if(this.chromosome.equals(other.chromosome)){
			if(this.start == other.start){
				return other.stop - this.stop;
			} else {
				return this.start - other.start;
			}
        } else {
            return this.chromosome.compareTo(other.chromosome);
        }
    }

    /**
     * Assumes that chromosome is correct
     *
     * @param position
     * @return
     */
    public boolean isInBoundarie(int position){
        return isInBoundarie(position, 0);
    }

    /**
     * Assumes that chromosome is correct
     *
     * @param position
	 * @param margin 
     * @return
     */
    public boolean isInBoundarie(int position, int margin){
        return (position >= (start - margin) && position <= (stop + margin));
    }

	/**
	 * Is this boundary a sub boundary of the other boundary.
	 * Assumes that chromosome is correct
	 *
	 * @param other
	 * @return
	 */
	public boolean isPartOfBoundary(GenomicBoundary other){

		if(other == null){
			return false;
		}

		return other.start <= this.start && other.stop >= this.stop ? true : false;
		
	}
	
	public boolean isOverlaping(GenomicBoundary other){
		
		if(other == null){
			return false;
		}
		
		if( ! this.chromosome.equals(other.chromosome)){
			return false;
		}

		
		if(this.start <= other.start && this.stop >= other.start){
			
			return true;
			
		} else if(other.start <= this.start && other.stop >= this.start){
			
			return true;
			
		} else {
			
			return false;
			
		}
		
	}

	public V getAnnotation() {
		return annotation;
	}
	
	
	public int getLength(boolean isInclusive){
		
		if(isInclusive){
			return ((this.getStop() + 1) - this.getStart());
		}
		
		else{
			return (this.getStop() - this.getStart());
		}
		
	}

}
