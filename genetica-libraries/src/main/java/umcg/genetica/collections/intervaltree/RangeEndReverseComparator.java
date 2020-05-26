package umcg.genetica.collections.intervaltree;

import java.util.Comparator;

/**
 *
 * @author Patrick Deelen
 */
public class RangeEndReverseComparator implements Comparator<Range>{

	@Override
	/**
	 * Note: we consider large end values to be smaller for reverse sort
	 */
	public int compare(Range o1, Range o2) {
		return o2.getEnd() - o1.getEnd();
	}

}
