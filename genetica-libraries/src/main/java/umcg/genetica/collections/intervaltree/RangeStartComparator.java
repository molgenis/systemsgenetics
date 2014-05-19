package umcg.genetica.collections.intervaltree;

import java.util.Comparator;

/**
 *
 * @author Patrick Deelen
 */
public class RangeStartComparator implements Comparator<Range>{

	@Override
	public int compare(Range o1, Range o2) {
		return o1.getStart() - o2.getStart();
	}

}
